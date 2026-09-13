#!/usr/bin/env bash
# Build the audited locus-level and sample-call Parquets for parameter validation.
# Real extraction is intentionally restricted to a Bunya SLURM allocation.
# Agreed compromise: evaluate reported GT correctness regardless of GP confidence.
# Retain GP and maxGP for auditing, omit HD; no GP-argmax correctness or calibration.
# Final mask: all 25 targets must have evaluable GT and truth in all six conditions.

#SBATCH --job-name=quilt2_parameter_extract
#SBATCH --account=a_qaafi_cas
#SBATCH --partition=general
#SBATCH --constraint=epyc4
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=72:00:00

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
R_HELPER="${SCRIPT_DIR}/build_quilt2_parameter_parquet.R"
MASK_POSITION_HELPER="${SCRIPT_DIR}/extract_array_evaluation_positions.R"

# Hard-coded Bunya configuration: no command-line arguments are required.
# The bundled manifest and sample map are resolved relative to this script.
# ARRAY_TRUTH and its index are listed in scratch_structure.txt; verify on Bunya.
# OUTPUT_DIR is relative to where bash is launched, not SCRIPT_DIR.
# Existing command-line options can override these declarations for a pilot run.
declare RUN_MANIFEST="${ROOT_DIR}/analysis/quilt2_parameter_validation/run_manifest.example.tsv"
declare ARRAY_TRUTH="/scratch/project_mnt/S0218/downsampling/truth_array.vcf.gz"
declare WGS_TRUTH_DIR="/QRISdata/Q8367/WGS_Reference_Panel/NCBI_truth_set/7.Consolidated_VCF"
declare REFERENCE_FASTA="/QRISdata/Q8367/Reference_Genome/GDDH13_1-1_formatted.fasta"
declare SAMPLE_MAP="${ROOT_DIR}/analysis/quilt2_parameter_validation/sample_map.tsv"
declare OUTPUT_DIR="./quilt2_parameter_validation"
declare CHR_ARG=""  # Empty means all chromosomes, Chr01-Chr17.
declare MIN_GQ="60"
declare MIN_DP="10"
declare GP_SUM_TOLERANCE="0.001"
declare DRY_RUN=false
declare FORCE=false
declare CONDA_ENV="${CONDA_ENV:-myenv_py310}"
declare MINIFORGE_MODULE="${MINIFORGE_MODULE:-miniforge/26.1.0-0}"
declare BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"

# Completed, standardised Array Group I/II 2x evaluations, keyed by treatment|panel.
# The Filtered Liao run uses the legacy dot-separated Parquet filename.
declare -A ARRAY_MASK_PARQUETS=(
    ["Filtered|Liao"]="/scratch/project_mnt/S0218/downsampling/Liao_sample_downsampling_imputed_w_Liao_reference/eval/dosage_eval.concordance.parquet"
    ["Filtered|NCBI"]="/scratch/project_mnt/S0218/downsampling/Liao_sample_downsampling_imputed_w_NCBI_reference/eval/dosage_eval/concordance.parquet"
    ["Filtered|Combined"]="/scratch/project_mnt/S0218/downsampling/Liao_sample_downsampling_imputed_w_Combined_reference/eval/dosage_eval/concordance.parquet"
    ["No_filter|Liao"]="/scratch/project_mnt/S0218/downsampling/No_filter/Liao_sample_downsampling_imputed_w_Liao_reference/eval/dosage_eval/concordance.parquet"
    ["No_filter|NCBI"]="/scratch/project_mnt/S0218/downsampling/No_filter/Liao_sample_downsampling_imputed_w_NCBI_reference/eval/dosage_eval/concordance.parquet"
    ["No_filter|Combined"]="/scratch/project_mnt/S0218/downsampling/No_filter/Liao_sample_downsampling_imputed_w_Combined_reference/eval/dosage_eval/concordance.parquet"
)

# WGS evaluator datasets: CHROM=ChrNN/part-000.parquet beneath each directory.
# All six holdout paths are listed in scratch_structure.txt. WGS no_filter runs
# are under NCBI_downsampling_imputed/no_filter, not the Array No_filter directory.
declare -A WGS_MASK_PARQUETS=(
    ["Filtered|Liao"]="/scratch/project_mnt/S0218/downsampling/NCBI_downsampling_imputed/filter/NCBI_WGS_Liao_holdout/eval/dosage_eval_wgs/metrics/per_variant_metrics"
    ["Filtered|NCBI"]="/scratch/project_mnt/S0218/downsampling/NCBI_downsampling_imputed/filter/NCBI_WGS_NCBI_holdout/eval/dosage_eval_wgs/metrics/per_variant_metrics"
    ["Filtered|Combined"]="/scratch/project_mnt/S0218/downsampling/NCBI_downsampling_imputed/filter/NCBI_WGS_Combined_holdout/eval/dosage_eval_wgs/metrics/per_variant_metrics"
    ["No_filter|Liao"]="/scratch/project_mnt/S0218/downsampling/NCBI_downsampling_imputed/no_filter/NCBI_WGS_no_filter_Liao_holdout/eval/dosage_eval_wgs/metrics/per_variant_metrics"
    ["No_filter|NCBI"]="/scratch/project_mnt/S0218/downsampling/NCBI_downsampling_imputed/no_filter/NCBI_WGS_no_filter_NCBI_holdout/eval/dosage_eval_wgs/metrics/per_variant_metrics"
    ["No_filter|Combined"]="/scratch/project_mnt/S0218/downsampling/NCBI_downsampling_imputed/no_filter/NCBI_WGS_no_filter_Combined_holdout/eval/dosage_eval_wgs/metrics/per_variant_metrics"
)

log_info() { printf '[INFO] %s\n' "$*" >&2; }
log_warn() { printf '[WARN] %s\n' "$*" >&2; }
die() { printf '[ERROR] %s\n' "$*" >&2; exit 1; }

usage() {
    cat <<'EOF'
Usage: extract_quilt2_parameter_validation.sh [options]

No arguments are required. Edit the hard-coded configuration near the top of
the script, or supply optional overrides below. Real extraction requires an
existing Bunya SLURM allocation; this script does not submit itself.

Options:
  --run-manifest FILE    Twelve-row manifest: Array and WGS physical inputs for
                         each Filtered/No_filter x Liao/NCBI/Combined condition.
                         Default: bundled run_manifest.example.tsv.
  --array-truth FILE     Array truth VCF/BCF containing Group I/II samples.
                         Default: /scratch/project_mnt/S0218/downsampling/
                                  truth_array.vcf.gz
  --output-dir DIR       Directory for both Parquets and audit sidecars.
                         Default: ./quilt2_parameter_validation
                         Created after preflight passes.
  --wgs-truth-dir DIR    Directory of Chr*_consolidated.vcf.gz truth files.
                         Default: /QRISdata/Q8367/WGS_Reference_Panel/
                                  NCBI_truth_set/7.Consolidated_VCF
  --reference-fasta FILE Reference FASTA used to verify/normalise SNP alleles.
                         Default: /QRISdata/Q8367/Reference_Genome/
                                  GDDH13_1-1_formatted.fasta
  --sample-map FILE      25-sample metadata table. Default: repository copy.
  --chr LIST             Comma/space-separated ChrNN names. Default: Chr01-Chr17.
  --min-gq NUMBER        WGS truth minimum GQ. Default: 60.
  --min-dp NUMBER        WGS truth minimum DP. Default: 10.
  --gp-sum-tolerance X   GP field-validity audit only; never filters GT accuracy.
                         Allowed deviation of GP sum from 1. Default: 0.001.
  --dry-run              Validate inputs, source headers, samples, and dependencies;
                         do not stage VCFs or write outputs. May run outside SLURM.
  --force                Replace final outputs only after a complete new extraction.
  --help                 Show this help.

Run-manifest columns (tab-separated, exactly in this order):
  source_id logical_run_id treatment panel truth_source input_type input_path

Input modes:
  vcf              Read one existing indexed VCF/BCF directly, such as
                   imputed.all_chroms.vcf.gz.
  chromosome_vcfs  Read existing indexed chromosome VCFs directly from
                   INPUT_PATH/ChrNN/imputed.ChrNN.vcf.gz.

The supplied manifest uses vcf for Array and chromosome_vcfs for WGS.
Every selected chromosome VCF and its index must already exist. Source VCFs
are read-only; raw chunks are never rebuilt, even with --force. Only the small
masked chromosome subsets are combined afterward for the Parquet builder.
For an older manifest, replace WGS input_type "chunks" with "chromosome_vcfs";
keep its input_path pointing to the run's chunks/imputed directory.
Paths may contain spaces but not tabs/newlines.

The six Array concordance Parquet paths and six WGS per_variant_metrics dataset
paths are fixed in ARRAY_MASK_PARQUETS and WGS_MASK_PARQUETS near the top of this
script. Their twelve-way CHROM/POS intersection is applied to every Array and
WGS input, including truth. Inputs already use standardised ChrNN names.
This starting mask depends only on position presence. WGS CHROM is read from
its Hive partition directories. After staging, the R helper retains only unique
allele-compatible SNPs with valid GT and aligned truth for all 25 targets in
all six treatment-panel conditions (150/150 comparisons per locus). WGS truth
GQ/DP thresholds and unambiguous Array A/B alignment still apply. Incorrect GT
calls remain eligible; GP, DS, INFO_SCORE, and HWE do not filter the final mask.
Approximately 8,000 final loci is an expectation, not an enforced count.

Analysis compromise: assess reported GT accuracy regardless of GP confidence.
There are no GP bins, GP-argmax comparisons, or probability-calibration plots.
Missing or malformed GP does not exclude an otherwise evaluable GT call.
Truth validity requirements still apply to GT comparisons.

The command writes the audited two-file schema in OUTPUT_DIR:
  quilt2_locus_parameters.parquet   One row per locus x physical run.
  quilt2_sample_calls.parquet       One row per locus x run x sample.
  Each Parquet has an adjacent .sha256 audit sidecar.
  quilt2_parameter_validation.common_loci.tsv.gz
  quilt2_parameter_validation.summary.tsv

INFO fields have no info__ prefix. FORMAT fields keep fmt__, except HD is
omitted. Calls retain the original GP vector, gp_valid, maxGP, and the same-run
info_score, hwe, eaf, erc, eac, paf, info_score_valid and hwe_valid. These site
annotations keep their original payload types and are shared within each locus x
treatment x panel x Array/WGS cohort, not necessarily between the two cohorts.
Provenance, fixed denominators, and constant VCF fields are metadata, not repeated
columns. Existing single-file v5 outputs are not modified or deleted.
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-manifest|--array-truth|--wgs-truth-dir|--reference-fasta|--sample-map|--output-dir|--chr|--min-gq|--min-dp|--gp-sum-tolerance)
            [[ $# -ge 2 ]] || die "Missing value for $1"
            case "$1" in
                --run-manifest) RUN_MANIFEST="$2" ;;
                --array-truth) ARRAY_TRUTH="$2" ;;
                --wgs-truth-dir) WGS_TRUTH_DIR="$2" ;;
                --reference-fasta) REFERENCE_FASTA="$2" ;;
                --sample-map) SAMPLE_MAP="$2" ;;
                --output-dir) OUTPUT_DIR="$2" ;;
                --chr) CHR_ARG="$2" ;;
                --min-gq) MIN_GQ="$2" ;;
                --min-dp) MIN_DP="$2" ;;
                --gp-sum-tolerance) GP_SUM_TOLERANCE="$2" ;;
            esac
            shift 2
            ;;
        --dry-run) DRY_RUN=true; shift ;;
        --force) FORCE=true; shift ;;
        --output) die "The two-file schema uses --output-dir DIR instead of --output FILE" ;;
        --help|-h) usage; exit 0 ;;
        *) die "Unknown option: $1" ;;
    esac
done

[[ -n "${OUTPUT_DIR}" ]] || die "--output-dir must not be empty"
[[ -f "${RUN_MANIFEST}" ]] || die "Run manifest not found: ${RUN_MANIFEST}"
[[ -f "${SAMPLE_MAP}" ]] || die "Sample map not found: ${SAMPLE_MAP}"
[[ -f "${ARRAY_TRUTH}" ]] || die "Array truth not found: ${ARRAY_TRUTH}"
[[ -d "${WGS_TRUTH_DIR}" ]] || die "WGS truth directory not found: ${WGS_TRUTH_DIR}"
[[ -f "${REFERENCE_FASTA}" ]] || die "Reference FASTA not found: ${REFERENCE_FASTA}"
[[ -f "${REFERENCE_FASTA}.fai" ]] || die "Reference FASTA index not found: ${REFERENCE_FASTA}.fai"
[[ -f "${R_HELPER}" ]] || die "Missing R helper: ${R_HELPER}"
[[ -f "${MASK_POSITION_HELPER}" ]] || die "Missing mask-position helper: ${MASK_POSITION_HELPER}"
awk -v value="${MIN_GQ}" 'BEGIN {exit !(value ~ /^[0-9]+([.][0-9]+)?$/ && value + 0 >= 0)}' \
    || die "--min-gq must be non-negative"
awk -v value="${MIN_DP}" 'BEGIN {exit !(value ~ /^[0-9]+([.][0-9]+)?$/ && value + 0 >= 0)}' \
    || die "--min-dp must be non-negative"
awk -v value="${GP_SUM_TOLERANCE}" 'BEGIN {exit !(value ~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ && value + 0 > 0)}' \
    || die "--gp-sum-tolerance must be positive"

if [[ "${DRY_RUN}" != "true" ]]; then
    [[ "$(uname -s)" == "Linux" ]] || die "Real extraction is Bunya/Linux-only; use --dry-run for local preflight"
    [[ -n "${SLURM_JOB_ID:-}" ]] || die "Real extraction requires a Bunya SLURM allocation (SLURM_JOB_ID is unset)"
fi

setup_tools() {
    if command -v module >/dev/null 2>&1; then
        module load "${MINIFORGE_MODULE}" >/dev/null 2>&1 || die "Could not load ${MINIFORGE_MODULE}"
        if [[ -n "${ROOTMINIFORGE:-}" && -f "${ROOTMINIFORGE}/etc/profile.d/conda.sh" ]]; then
            # shellcheck source=/dev/null
            source "${ROOTMINIFORGE}/etc/profile.d/conda.sh"
            conda activate "${CONDA_ENV}" >/dev/null 2>&1 || die "Could not activate ${CONDA_ENV}"
        else
            die "Conda initialisation script was not found under ROOTMINIFORGE"
        fi
        command -v bcftools >/dev/null 2>&1 || module load "${BCFTOOLS_MODULE}" >/dev/null 2>&1 || true
    fi
    command -v bcftools >/dev/null 2>&1 || die "bcftools is required"
    command -v Rscript >/dev/null 2>&1 || die "Rscript is required"
    command -v sha256sum >/dev/null 2>&1 || die "sha256sum is required"
    Rscript -e 'packages <- c("arrow", "data.table", "dplyr", "digest", "jsonlite"); quit(status=as.integer(any(!vapply(packages, requireNamespace, logical(1), quietly=TRUE))))' \
        || die "R packages arrow, data.table, dplyr, digest, and jsonlite are required in ${CONDA_ENV}"
}

setup_tools

declare -a CHROMS=()
if [[ -n "${CHR_ARG}" ]]; then
    read -r -a RAW_CHROMS <<< "${CHR_ARG//,/ }"
else
    RAW_CHROMS=()
    for number in $(seq 1 17); do
        printf -v chr 'Chr%02d' "${number}"
        RAW_CHROMS+=("${chr}")
    done
fi
declare -A CHR_SEEN=()
for raw_chr in "${RAW_CHROMS[@]}"; do
    chr="${raw_chr}"
    [[ "${chr}" =~ ^Chr(0[1-9]|1[0-7])$ ]] || die "--chr requires Chr01-Chr17 names: ${chr}"
    [[ -z "${CHR_SEEN[${chr}]:-}" ]] || die "Duplicate --chr argument: ${chr}"
    CHR_SEEN[${chr}]=1
    CHROMS+=("${chr}")
done
(( ${#CHROMS[@]} > 0 )) || die "No chromosomes selected"
mapfile -t CHROMS < <(printf '%s\n' "${CHROMS[@]}" | LC_ALL=C sort)
CHR_CSV="$(IFS=,; printf '%s' "${CHROMS[*]}")"

for chr in "${CHROMS[@]}"; do
    awk -v wanted="${chr}" '$1 == wanted {found=1} END {exit !found}' "${REFERENCE_FASTA}.fai" \
        || die "${chr} is absent from ${REFERENCE_FASTA}.fai"
done

EXPECTED_SAMPLE_HEADER=$'sample_id\ttruth_source\tarray_group'
[[ "$(head -n 1 "${SAMPLE_MAP}" | tr -d '\r')" == "${EXPECTED_SAMPLE_HEADER}" ]] \
    || die "Unexpected sample-map header"
awk -F'\t' '
    NR == 1 {next}
    NF != 3 || $1 == "" || $2 == "" || $3 == "" {exit 20}
    seen[$1]++ {exit 21}
    $2 == "array" && ($3 == "Group I" || $3 == "Group II") {array++ ; next}
    $2 == "wgs" && $3 == "WGS" {wgs++ ; next}
    {exit 22}
    END {if (NR != 26 || array != 18 || wgs != 7) exit 23}
' "${SAMPLE_MAP}" || die "Sample map must contain unique IDs for 18 Group I/II Array and 7 WGS samples"

EXPECTED_RUN_HEADER=$'source_id\tlogical_run_id\ttreatment\tpanel\ttruth_source\tinput_type\tinput_path'
[[ "$(head -n 1 "${RUN_MANIFEST}" | tr -d '\r')" == "${EXPECTED_RUN_HEADER}" ]] \
    || die "Unexpected run-manifest header"
mapfile -t RUN_ROWS < <(tail -n +2 "${RUN_MANIFEST}" | sed 's/\r$//')
(( ${#RUN_ROWS[@]} == 12 )) || die "Run manifest must contain exactly 12 physical sources"

declare -A CONDITION_SEEN=()
declare -A LOGICAL_RUN_FOR_CONDITION=()
declare -A SOURCE_ID_SEEN=()
declare -A LOGICAL_ID_SEEN=()
for row in "${RUN_ROWS[@]}"; do
    IFS=$'\t' read -r source_id logical_run_id treatment panel truth_source input_type input_path extra <<< "${row}"
    [[ -z "${extra:-}" ]] || die "Too many columns in run-manifest row: ${source_id:-<unknown>}"
    for value in "${source_id}" "${logical_run_id}" "${treatment}" "${panel}" "${truth_source}" "${input_type}" "${input_path}"; do
        [[ -n "${value}" ]] || die "Blank run-manifest field"
    done
    [[ "${treatment}" == "Filtered" || "${treatment}" == "No_filter" ]] || die "Invalid treatment: ${treatment}"
    [[ "${panel}" == "Liao" || "${panel}" == "NCBI" || "${panel}" == "Combined" ]] || die "Invalid panel: ${panel}"
    [[ "${truth_source}" == "array" || "${truth_source}" == "wgs" ]] || die "Invalid truth_source: ${truth_source}"
    [[ "${input_type}" == "vcf" || "${input_type}" == "chromosome_vcfs" ]] \
        || die "Invalid input_type: ${input_type}; use vcf or chromosome_vcfs (raw-chunk rebuilding is not supported)"
    [[ -z "${SOURCE_ID_SEEN[${source_id}]:-}" ]] || die "Duplicate source_id: ${source_id}"
    SOURCE_ID_SEEN[${source_id}]=1
    condition_key="${treatment}|${panel}|${truth_source}"
    [[ -z "${CONDITION_SEEN[${condition_key}]:-}" ]] || die "Duplicate condition: ${condition_key}"
    CONDITION_SEEN[${condition_key}]=1
    logical_key="${treatment}|${panel}"
    if [[ -n "${LOGICAL_RUN_FOR_CONDITION[${logical_key}]:-}" && "${LOGICAL_RUN_FOR_CONDITION[${logical_key}]}" != "${logical_run_id}" ]]; then
        die "Array/WGS rows disagree on logical_run_id for ${logical_key}"
    fi
    LOGICAL_RUN_FOR_CONDITION[${logical_key}]="${logical_run_id}"
    if [[ -n "${LOGICAL_ID_SEEN[${logical_run_id}]:-}" && "${LOGICAL_ID_SEEN[${logical_run_id}]}" != "${logical_key}" ]]; then
        die "logical_run_id ${logical_run_id} is assigned to more than one treatment-panel condition"
    fi
    LOGICAL_ID_SEEN[${logical_run_id}]="${logical_key}"
    if [[ "${input_type}" == "vcf" ]]; then
        [[ -f "${input_path}" ]] || die "Source VCF not found: ${input_path}"
        [[ -f "${input_path}.csi" || -f "${input_path}.tbi" ]] || die "Source VCF/BCF is not indexed: ${input_path}"
    else
        [[ -d "${input_path}" ]] || die "Chromosome-VCF directory not found: ${input_path}"
        for chr in "${CHROMS[@]}"; do
            chromosome_vcf="${input_path}/${chr}/imputed.${chr}.vcf.gz"
            [[ -f "${chromosome_vcf}" && -s "${chromosome_vcf}" && -r "${chromosome_vcf}" ]] \
                || die "Existing chromosome VCF is missing, empty, or unreadable: ${chromosome_vcf}; no rebuild will be attempted"
            [[ -f "${chromosome_vcf}.csi" || -f "${chromosome_vcf}.tbi" ]] \
                || die "Existing chromosome VCF is not indexed: ${chromosome_vcf}; no rebuild will be attempted"
        done
    fi
done
for treatment in Filtered No_filter; do
    for panel in Liao NCBI Combined; do
        for truth_source in array wgs; do
            key="${treatment}|${panel}|${truth_source}"
            [[ -n "${CONDITION_SEEN[${key}]:-}" ]] || die "Missing condition: ${key}"
        done
    done
done

declare -a MASK_ROWS=()
for treatment in Filtered No_filter; do
    for panel in Liao NCBI Combined; do
        key="${treatment}|${panel}"
        for truth_source in array wgs; do
            if [[ "${truth_source}" == "array" ]]; then
                mask_input_path="${ARRAY_MASK_PARQUETS[${key}]}"
                mask_input_type="concordance_parquet"
                [[ -f "${mask_input_path}" ]] || die "Array evaluator Parquet not found: ${mask_input_path}"
            else
                mask_input_path="${WGS_MASK_PARQUETS[${key}]:-}"
                [[ -n "${mask_input_path}" ]] || die "WGS mask path is not configured for ${key}"
                mask_input_type="per_variant_metrics"
                [[ -d "${mask_input_path}" ]] || die "WGS evaluator dataset not found: ${mask_input_path}"
            fi
            printf -v row '%s\t%s\t%s\t%s\t%s\t%s' \
                "${LOGICAL_RUN_FOR_CONDITION[${key}]}" "${treatment}" "${panel}" \
                "${truth_source}" "${mask_input_type}" "${mask_input_path}"
            MASK_ROWS+=("${row}")
        done
    done
done

ARRAY_SAMPLES="$(mktemp "${TMPDIR:-/tmp}/quilt2_array_samples.XXXXXX")"
WGS_SAMPLES="$(mktemp "${TMPDIR:-/tmp}/quilt2_wgs_samples.XXXXXX")"
trap 'rm -f "${ARRAY_SAMPLES}" "${WGS_SAMPLES}"' EXIT
awk -F'\t' 'NR > 1 && $2 == "array" {print $1}' "${SAMPLE_MAP}" > "${ARRAY_SAMPLES}"
awk -F'\t' 'NR > 1 && $2 == "wgs" {print $1}' "${SAMPLE_MAP}" > "${WGS_SAMPLES}"

declare -A WGS_TRUTH_FILES=()
for chr in "${CHROMS[@]}"; do
    truth_vcf="${WGS_TRUTH_DIR}/${chr}_consolidated.vcf.gz"
    [[ -f "${truth_vcf}" ]] || die "WGS truth VCF not found: ${truth_vcf}"
    [[ -f "${truth_vcf}.csi" || -f "${truth_vcf}.tbi" ]] || die "WGS truth VCF is not indexed: ${truth_vcf}"
    WGS_TRUTH_FILES[${chr}]="${truth_vcf}"
done

check_samples_present() {
    local input="$1" expected="$2" label="$3" observed sorted_expected missing
    observed="$(mktemp "${TMPDIR:-/tmp}/quilt2_observed_samples.XXXXXX")"
    sorted_expected="$(mktemp "${TMPDIR:-/tmp}/quilt2_expected_samples.XXXXXX")"
    bcftools query -l "${input}" | LC_ALL=C sort -u > "${observed}"
    LC_ALL=C sort -u "${expected}" > "${sorted_expected}"
    missing="$(comm -23 "${sorted_expected}" "${observed}" | paste -sd, -)"
    rm -f "${observed}" "${sorted_expected}"
    [[ -z "${missing}" ]] || die "${label} is missing target sample(s): ${missing}"
}

check_samples_present "${ARRAY_TRUTH}" "${ARRAY_SAMPLES}" "Array truth"
for chr in "${CHROMS[@]}"; do
    check_samples_present "${WGS_TRUTH_FILES[${chr}]}" "${WGS_SAMPLES}" "WGS truth ${chr}"
done
for row in "${RUN_ROWS[@]}"; do
    IFS=$'\t' read -r source_id logical_run_id treatment panel truth_source input_type input_path <<< "${row}"
    [[ "${truth_source}" == "array" ]] && expected_samples="${ARRAY_SAMPLES}" || expected_samples="${WGS_SAMPLES}"
    if [[ "${input_type}" == "vcf" ]]; then
        check_samples_present "${input_path}" "${expected_samples}" "Source ${source_id}"
    else
        for chr in "${CHROMS[@]}"; do
            chromosome_vcf="${input_path}/${chr}/imputed.${chr}.vcf.gz"
            bcftools index -n "${chromosome_vcf}" >/dev/null \
                || die "Cannot read the chromosome VCF index: ${chromosome_vcf}"
            check_samples_present "${chromosome_vcf}" "${expected_samples}" "Source ${source_id} ${chr}"
        done
    fi
done

if [[ "${DRY_RUN}" == "true" ]]; then
    for row in "${MASK_ROWS[@]}"; do
        IFS=$'\t' read -r logical_run_id treatment panel truth_source mask_input_type mask_input_path <<< "${row}"
        Rscript "${MASK_POSITION_HELPER}" \
            --input "${mask_input_path}" \
            --chr "${CHR_CSV}" \
            --validate-only >/dev/null
    done
    log_info "Preflight passed for 12 physical sources, 12 evaluator masks (6 Array + 6 WGS), 6 logical conditions, 25 samples, and ${#CHROMS[@]} chromosome(s)."
    log_info "No VCFs were staged and no outputs were written."
    exit 0
fi

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR="$(cd "${OUTPUT_DIR}" && pwd)"
LOCUS_OUTPUT="${OUTPUT_DIR}/quilt2_locus_parameters.parquet"
CALLS_OUTPUT="${OUTPUT_DIR}/quilt2_sample_calls.parquet"
MASK_OUTPUT="${OUTPUT_DIR}/quilt2_parameter_validation.common_loci.tsv.gz"
SUMMARY_OUTPUT="${OUTPUT_DIR}/quilt2_parameter_validation.summary.tsv"
LOCUS_CHECKSUM="${LOCUS_OUTPUT}.sha256"
CALLS_CHECKSUM="${CALLS_OUTPUT}.sha256"
if [[ "${FORCE}" != "true" ]]; then
    for existing in "${LOCUS_OUTPUT}" "${CALLS_OUTPUT}" "${MASK_OUTPUT}" "${SUMMARY_OUTPUT}" "${LOCUS_CHECKSUM}" "${CALLS_CHECKSUM}"; do
        [[ ! -e "${existing}" ]] || die "Output exists (use --force): ${existing}"
    done
fi

STAGE_DIR="$(mktemp -d "${TMPDIR:-/tmp}/quilt2_parameter_extract.XXXXXX")"
PARTIAL_LOCUS="${LOCUS_OUTPUT}.partial.${SLURM_JOB_ID}"
PARTIAL_CALLS="${CALLS_OUTPUT}.partial.${SLURM_JOB_ID}"
PARTIAL_MASK="${MASK_OUTPUT}.partial.${SLURM_JOB_ID}"
PARTIAL_SUMMARY="${SUMMARY_OUTPUT}.partial.${SLURM_JOB_ID}"
cleanup() {
    status=$?
    [[ -z "${STAGE_DIR}" || ! -d "${STAGE_DIR}" ]] || rm -rf "${STAGE_DIR}"
    rm -f "${PARTIAL_LOCUS}" "${PARTIAL_CALLS}" "${PARTIAL_MASK}" "${PARTIAL_SUMMARY}" "${ARRAY_SAMPLES}" "${WGS_SAMPLES}"
    exit "${status}"
}
trap cleanup EXIT

sha256_file() {
    sha256sum "$1" | awk '{print $1}'
}

hash_file_set() {
    local list_file="$1" hash_rows="$2" path
    : > "${hash_rows}"
    while IFS= read -r path; do
        [[ -f "${path}" ]] || die "Cannot hash missing source file: ${path}"
        printf '%s\t%s\n' "$(sha256_file "${path}")" "${path}" >> "${hash_rows}"
    done < "${list_file}"
    [[ -s "${hash_rows}" ]] || die "No source files were available for hashing"
    sha256_file "${hash_rows}"
}

sample_signature() {
    local input="$1" output="$2"
    bcftools query -l "${input}" | LC_ALL=C sort -u > "${output}"
    [[ -s "${output}" ]] || die "No samples found in source VCF: ${input}"
    SAMPLE_COUNT="$(wc -l < "${output}" | tr -d ' ')"
    SAMPLE_SHA256="$(sha256_file "${output}")"
}

combine_vcfs() {
    local list_file="$1" output="$2" count first
    count="$(wc -l < "${list_file}" | tr -d ' ')"
    (( count > 0 )) || die "Cannot combine an empty VCF list"
    if (( count == 1 )); then
        first="$(head -n 1 "${list_file}")"
        cp "${first}" "${output}"
    else
        bcftools concat -Oz -f "${list_file}" -o "${output}"
    fi
    bcftools index -f -c "${output}"
}

candidate_hashes() {
    local candidate="$1"
    CANDIDATE_SHA256="$(sha256_file "${candidate}")"
    HEADER_SHA256="$(bcftools view -h "${candidate}" | sha256sum | awk '{print $1}')"
}

log_info "Building the starting presence mask from six Array and six WGS evaluations"
EVAL_MASK_DIR="${STAGE_DIR}/evaluation_mask"
mkdir -p "${EVAL_MASK_DIR}"
STAGED_EVAL_MASK_MANIFEST="${EVAL_MASK_DIR}/staged_evaluation_masks.tsv"
printf 'logical_run_id\ttreatment\tpanel\ttruth_source\tmask_input_type\tmask_input_path\tsource_sha256\tposition_count\tpositions_sha256\n' \
    > "${STAGED_EVAL_MASK_MANIFEST}"
declare -a EVAL_MASK_POSITION_FILES=()
mask_index=0
for row in "${MASK_ROWS[@]}"; do
    mask_index=$((mask_index + 1))
    IFS=$'\t' read -r logical_run_id treatment panel truth_source mask_input_type mask_input_path <<< "${row}"
    mask_stem="${EVAL_MASK_DIR}/mask_$(printf '%02d' "${mask_index}")"
    mask_positions="${mask_stem}.positions.tsv"
    Rscript "${MASK_POSITION_HELPER}" \
        --input "${mask_input_path}" \
        --chr "${CHR_CSV}" \
        --output "${mask_positions}" >/dev/null
    mask_position_count="$(wc -l < "${mask_positions}" | tr -d ' ')"
    if [[ "${truth_source}" == "array" ]]; then
        mask_source_sha256="$(sha256_file "${mask_input_path}")"
    else
        find "${mask_input_path}" -type f -name '*.parquet' -print | LC_ALL=C sort > "${mask_stem}.source_files.txt"
        mask_source_sha256="$(hash_file_set "${mask_stem}.source_files.txt" "${mask_stem}.source_hashes.tsv")"
    fi
    mask_positions_sha256="$(sha256_file "${mask_positions}")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${logical_run_id}" "${treatment}" "${panel}" "${truth_source}" "${mask_input_type}" \
        "${mask_input_path}" "${mask_source_sha256}" "${mask_position_count}" \
        "${mask_positions_sha256}" >> "${STAGED_EVAL_MASK_MANIFEST}"
    EVAL_MASK_POSITION_FILES+=("${mask_positions}")
done

POSITIONS_UNSORTED="${EVAL_MASK_DIR}/common_positions.unsorted.tsv"
POSITIONS_FILE="${EVAL_MASK_DIR}/common_positions.tsv"
# Count presence once per input file. Start with an Array list so the in-memory
# key set stays small even when WGS datasets contain millions of positions.
awk -F'\t' -v OFS='\t' '
    FILENAME == ARGV[1] {count[$1 OFS $2] = 1; next}
    {
        position = $1 OFS $2
        if (position in count && last_file[position] != FILENAME) {
            count[position]++
            last_file[position] = FILENAME
        }
    }
    END {for (position in count) if (count[position] == ARGC - 1) print position}
' "${EVAL_MASK_POSITION_FILES[@]}" > "${POSITIONS_UNSORTED}"
LC_ALL=C sort -k1,1 -k2,2n "${POSITIONS_UNSORTED}" > "${POSITIONS_FILE}"
EVAL_MASK_POSITION_COUNT="$(wc -l < "${POSITIONS_FILE}" | tr -d ' ')"
log_info "Starting presence mask written: ${POSITIONS_FILE}"
log_info "Positions present in all 12 runs (6 Array + 6 WGS), before allele and sample-completeness checks: ${EVAL_MASK_POSITION_COUNT}"
[[ -s "${POSITIONS_FILE}" ]] || die "No positions are shared by all six Array and six WGS evaluator artifacts"
EVAL_MASK_SHA256="$(sha256_file "${POSITIONS_FILE}")"
log_info "Position-mask SHA-256: ${EVAL_MASK_SHA256}"
log_info "Applying this same position mask to Array truth, WGS truth, and all 12 QUILT2 sources"
log_info "Candidate staging follows: bcftools 'Lines total/split/realigned/skipped' counts are per chromosome/source, not final mask counts"
log_info "The R helper will report the final mask after requiring 150/150 valid GT-truth comparisons per locus"

log_info "Staging Array truth at the twelve-run evaluator positions"
ARRAY_TRUTH_DIR="${STAGE_DIR}/truth_array"
mkdir -p "${ARRAY_TRUTH_DIR}"
ARRAY_TRUTH_CANDIDATE="${ARRAY_TRUTH_DIR}/array_truth.candidate.vcf.gz"
bcftools view --no-update -S "${ARRAY_SAMPLES}" -m2 -M2 -v snps \
    -T "${POSITIONS_FILE}" -Oz -o "${ARRAY_TRUTH_CANDIDATE}" "${ARRAY_TRUTH}"
bcftools index -f -c "${ARRAY_TRUTH_CANDIDATE}"
# The R builder checks truth/source coverage after unique-SNP selection.

ARRAY_SOURCE_LIST="${ARRAY_TRUTH_DIR}/source_files.txt"
printf '%s\n' "${ARRAY_TRUTH}" > "${ARRAY_SOURCE_LIST}"
ARRAY_SOURCE_SHA256="$(hash_file_set "${ARRAY_SOURCE_LIST}" "${ARRAY_TRUTH_DIR}/source_hashes.tsv")"
sample_signature "${ARRAY_TRUTH}" "${ARRAY_TRUTH_DIR}/source_samples.txt"
ARRAY_SOURCE_SAMPLE_COUNT="${SAMPLE_COUNT}"
ARRAY_SOURCE_SAMPLES_SHA256="${SAMPLE_SHA256}"
candidate_hashes "${ARRAY_TRUTH_CANDIDATE}"
ARRAY_CANDIDATE_SHA256="${CANDIDATE_SHA256}"
ARRAY_HEADER_SHA256="${HEADER_SHA256}"

log_info "Staging WGS truth at the twelve-run evaluator positions"
WGS_TRUTH_STAGE="${STAGE_DIR}/truth_wgs"
mkdir -p "${WGS_TRUTH_STAGE}"
WGS_TRUTH_LIST="${WGS_TRUTH_STAGE}/candidate_chromosomes.txt"
WGS_SOURCE_LIST="${WGS_TRUTH_STAGE}/source_files.txt"
: > "${WGS_TRUTH_LIST}"
: > "${WGS_SOURCE_LIST}"
WGS_SOURCE_SAMPLE_COUNT=""
WGS_SOURCE_SAMPLES_SHA256=""
for chr in "${CHROMS[@]}"; do
    input="${WGS_TRUTH_FILES[${chr}]}"
    output_chr="${WGS_TRUTH_STAGE}/${chr}.candidate.vcf.gz"
    log_info "WGS truth ${chr}: normalising candidate records (not the final completeness mask)"
    bcftools view --no-update -r "${chr}" -S "${WGS_SAMPLES}" -m2 -M2 -v snps \
        -T "${POSITIONS_FILE}" -Ou "${input}" \
        | bcftools norm -f "${REFERENCE_FASTA}" -c e -Oz -o "${output_chr}"
    bcftools index -f -c "${output_chr}"
    sample_signature "${input}" "${WGS_TRUTH_STAGE}/${chr}.source_samples.txt"
    if [[ -z "${WGS_SOURCE_SAMPLE_COUNT}" ]]; then
        WGS_SOURCE_SAMPLE_COUNT="${SAMPLE_COUNT}"
        WGS_SOURCE_SAMPLES_SHA256="${SAMPLE_SHA256}"
    elif [[ "${WGS_SOURCE_SAMPLE_COUNT}" != "${SAMPLE_COUNT}" || "${WGS_SOURCE_SAMPLES_SHA256}" != "${SAMPLE_SHA256}" ]]; then
        die "WGS truth sample set differs between chromosome files"
    fi
    printf '%s\n' "${output_chr}" >> "${WGS_TRUTH_LIST}"
    printf '%s\n' "${input}" >> "${WGS_SOURCE_LIST}"
done
WGS_TRUTH_CANDIDATE="${WGS_TRUTH_STAGE}/wgs_truth.candidate.vcf.gz"
combine_vcfs "${WGS_TRUTH_LIST}" "${WGS_TRUTH_CANDIDATE}"
WGS_SOURCE_SHA256="$(hash_file_set "${WGS_SOURCE_LIST}" "${WGS_TRUTH_STAGE}/source_hashes.tsv")"
candidate_hashes "${WGS_TRUTH_CANDIDATE}"
WGS_CANDIDATE_SHA256="${CANDIDATE_SHA256}"
WGS_HEADER_SHA256="${HEADER_SHA256}"

STAGED_MANIFEST="${STAGE_DIR}/staged_runs.tsv"
TRUTH_MANIFEST="${STAGE_DIR}/staged_truth.tsv"
printf 'source_id\tlogical_run_id\ttreatment\tpanel\ttruth_source\tinput_type\tvcf\tsource_path\tsource_sample_count\tselected_sample_count\tsource_samples_sha256\tsource_signature_sha256\tcandidate_vcf_sha256\tcandidate_header_sha256\n' > "${STAGED_MANIFEST}"
printf 'truth_source\tvcf\tsource_path\tsource_sample_count\tselected_sample_count\tsource_samples_sha256\tsource_signature_sha256\tcandidate_vcf_sha256\tcandidate_header_sha256\n' > "${TRUTH_MANIFEST}"
printf 'array\t%s\t%s\t%s\t18\t%s\t%s\t%s\t%s\n' \
    "${ARRAY_TRUTH_CANDIDATE}" "${ARRAY_TRUTH}" "${ARRAY_SOURCE_SAMPLE_COUNT}" "${ARRAY_SOURCE_SAMPLES_SHA256}" \
    "${ARRAY_SOURCE_SHA256}" "${ARRAY_CANDIDATE_SHA256}" "${ARRAY_HEADER_SHA256}" >> "${TRUTH_MANIFEST}"
printf 'wgs\t%s\t%s\t%s\t7\t%s\t%s\t%s\t%s\n' \
    "${WGS_TRUTH_CANDIDATE}" "${WGS_TRUTH_DIR}" "${WGS_SOURCE_SAMPLE_COUNT}" "${WGS_SOURCE_SAMPLES_SHA256}" \
    "${WGS_SOURCE_SHA256}" "${WGS_CANDIDATE_SHA256}" "${WGS_HEADER_SHA256}" >> "${TRUTH_MANIFEST}"

log_info "Extracting common-mask records from 12 QUILT2 source runs"
source_index=0
for row in "${RUN_ROWS[@]}"; do
    source_index=$((source_index + 1))
    IFS=$'\t' read -r source_id logical_run_id treatment panel truth_source input_type input_path <<< "${row}"
    source_dir="${STAGE_DIR}/source_$(printf '%02d' "${source_index}")"
    mkdir -p "${source_dir}"
    chromosome_candidates="${source_dir}/candidate_chromosomes.txt"
    source_files="${source_dir}/source_files.txt"
    : > "${chromosome_candidates}"
    : > "${source_files}"
    [[ "${truth_source}" == "array" ]] && expected_samples="${ARRAY_SAMPLES}" || expected_samples="${WGS_SAMPLES}"

    if [[ "${input_type}" == "vcf" ]]; then
        printf '%s\n' "${input_path}" > "${source_files}"
    else
        for chr in "${CHROMS[@]}"; do
            printf '%s\n' "${input_path}/${chr}/imputed.${chr}.vcf.gz" >> "${source_files}"
        done
    fi

    source_signature="$(hash_file_set "${source_files}" "${source_dir}/source_hashes.tsv")"
    source_probe="$(head -n 1 "${source_files}")"
    sample_signature "${source_probe}" "${source_dir}/source_samples.txt"
    source_sample_count="${SAMPLE_COUNT}"
    source_samples_sha256="${SAMPLE_SHA256}"
    for chr in "${CHROMS[@]}"; do
        full_chr_vcf=""
        if [[ "${input_type}" == "chromosome_vcfs" ]]; then
            full_chr_vcf="${input_path}/${chr}/imputed.${chr}.vcf.gz"
            log_info "${source_id} ${chr}: reusing existing indexed chromosome VCF: ${full_chr_vcf}"
            sample_signature "${full_chr_vcf}" "${source_dir}/${chr}.source_samples.txt"
            if [[ "${SAMPLE_COUNT}" != "${source_sample_count}" || "${SAMPLE_SHA256}" != "${source_samples_sha256}" ]]; then
                die "Source ${source_id} sample set differs on ${chr}"
            fi
        else
            full_chr_vcf="${input_path}"
        fi
        output_chr="${source_dir}/${chr}.candidate.vcf.gz"
        log_info "${source_id} ${chr}: normalising candidate records (not the final completeness mask)"
        bcftools view --no-update -r "${chr}" -S "${expected_samples}" -m2 -M2 -v snps \
            -T "${POSITIONS_FILE}" -Ou "${full_chr_vcf}" \
            | bcftools norm -f "${REFERENCE_FASTA}" -c e -Oz -o "${output_chr}"
        bcftools index -f -c "${output_chr}"
        printf '%s\n' "${output_chr}" >> "${chromosome_candidates}"
    done

    staged_vcf="${source_dir}/${source_id}.candidate.vcf.gz"
    log_info "${source_id}: combining only the masked chromosome subsets for the Parquet builder"
    combine_vcfs "${chromosome_candidates}" "${staged_vcf}"
    check_samples_present "${staged_vcf}" "${expected_samples}" "Staged source ${source_id}"
    staged_sample_count="$(bcftools query -l "${staged_vcf}" | wc -l | tr -d ' ')"
    [[ "${truth_source}" == "array" ]] && expected_count=18 || expected_count=7
    [[ "${staged_sample_count}" == "${expected_count}" ]] || die "Staged source ${source_id} has ${staged_sample_count} samples; expected ${expected_count}"
    candidate_hashes "${staged_vcf}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${source_id}" "${logical_run_id}" "${treatment}" "${panel}" "${truth_source}" \
        "${input_type}" "${staged_vcf}" "${input_path}" "${source_sample_count}" "${expected_count}" \
        "${source_samples_sha256}" "${source_signature}" "${CANDIDATE_SHA256}" "${HEADER_SHA256}" >> "${STAGED_MANIFEST}"
done

log_info "Building the final all-samples/all-runs mask, locus/call tables, and correctness metrics"
log_info "Analysis compromise: reported GT correctness only, regardless of GP; probability calibration is not assessed"
Rscript "${R_HELPER}" \
    --staged-manifest "${STAGED_MANIFEST}" \
    --truth-manifest "${TRUTH_MANIFEST}" \
    --eval-mask-manifest "${STAGED_EVAL_MASK_MANIFEST}" \
    --eval-mask "${POSITIONS_FILE}" \
    --sample-map "${SAMPLE_MAP}" \
    --array-truth "${ARRAY_TRUTH_CANDIDATE}" \
    --wgs-truth "${WGS_TRUTH_CANDIDATE}" \
    --locus-output "${PARTIAL_LOCUS}" \
    --calls-output "${PARTIAL_CALLS}" \
    --mask-out "${PARTIAL_MASK}" \
    --summary-out "${PARTIAL_SUMMARY}" \
    --min-gq "${MIN_GQ}" \
    --min-dp "${MIN_DP}" \
    --gp-sum-tolerance "${GP_SUM_TOLERANCE}"

[[ -s "${PARTIAL_LOCUS}" ]] || die "R helper did not create the locus Parquet"
[[ -s "${PARTIAL_CALLS}" ]] || die "R helper did not create the calls Parquet"
[[ -s "${PARTIAL_MASK}" ]] || die "R helper did not create the mask audit table"
[[ -s "${PARTIAL_SUMMARY}" ]] || die "R helper did not create the summary audit table"
FINAL_MASK_LOCUS_COUNT="$(awk -F'\t' '$1 == "common_loci" {print $2}' "${PARTIAL_SUMMARY}")"
[[ "${FINAL_MASK_LOCUS_COUNT}" =~ ^[1-9][0-9]*$ ]] || die "R helper did not report a positive final-mask locus count"
log_info "FINAL MASK: ${FINAL_MASK_LOCUS_COUNT} loci evaluable in all 25 samples across all six conditions (150/150 comparisons per locus)"
mv -f "${PARTIAL_LOCUS}" "${LOCUS_OUTPUT}"
mv -f "${PARTIAL_CALLS}" "${CALLS_OUTPUT}"
mv -f "${PARTIAL_MASK}" "${MASK_OUTPUT}"
mv -f "${PARTIAL_SUMMARY}" "${SUMMARY_OUTPUT}"
sha256sum "${LOCUS_OUTPUT}" > "${LOCUS_CHECKSUM}"
sha256sum "${CALLS_OUTPUT}" > "${CALLS_CHECKSUM}"

log_info "Created locus table: ${LOCUS_OUTPUT}"
log_info "Created sample-call table: ${CALLS_OUTPUT}"
log_info "Audit mask: ${MASK_OUTPUT}"
log_info "Audit summary: ${SUMMARY_OUTPUT}"
log_info "Parquet checksums: ${LOCUS_CHECKSUM}, ${CALLS_CHECKSUM}"
