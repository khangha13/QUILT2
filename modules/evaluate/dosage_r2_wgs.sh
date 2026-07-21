#!/bin/bash
# WGS-truth dosage evaluation: QUILT2 FORMAT/DS versus filtered GATK FORMAT/GT.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ENV_FILE="${ROOT_DIR}/config/environment.sh"
ENV_TEMPLATE="${ROOT_DIR}/config/environment.template.sh"
R_HELPER="${SCRIPT_DIR}/dosage_r2_wgs.R"

usage() {
    cat <<'EOF'
Usage: dosage_r2_wgs.sh --imputed VCF --truth-dataset-dir DIR \
       --out-prefix DIR [options]

Required:
  --imputed PATH             Indexed QUILT2 VCF/BCF containing FORMAT/DS
  --truth-dataset-dir PATH   GATK 7.Consolidated_VCF directory containing
                             indexed Chr*_consolidated.vcf.gz files
  --out-prefix PATH          WGS evaluation output directory

Options:
  --chr LIST                 Comma/space-separated chromosome list
  --region STR               Canonical chromosome region, e.g. Chr01:1-1000000
  --samples FILE             One sample ID per line (default: sample intersection)
  --reference-fasta PATH     Override QUILT2_REFERENCE_FASTA from environment.sh
  --force                    Recompute an existing WGS-mode output
  --keep-temp                Retain temporary normalized VCFs and extraction tables
  -h, --help                 Show this help

WGS filter thresholds are read only from config/environment.sh. Array-specific
options are intentionally not accepted by this evaluator.
EOF
}

IMPUTED=""
TRUTH_DATASET_DIR=""
REFERENCE_FASTA=""
OUT_PREFIX=""
CHR_ARG=""
REGION=""
SAMPLES_FILE=""
FORCE=false
KEEP_TEMP=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --imputed) IMPUTED="$2"; shift 2 ;;
        --truth-dataset-dir) TRUTH_DATASET_DIR="$2"; shift 2 ;;
        --reference-fasta) REFERENCE_FASTA="$2"; shift 2 ;;
        --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
        --chr) CHR_ARG="$2"; shift 2 ;;
        --region) REGION="$2"; shift 2 ;;
        --samples) SAMPLES_FILE="$2"; shift 2 ;;
        --force) FORCE=true; shift ;;
        --keep-temp) KEEP_TEMP=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "[ERROR] Unknown WGS evaluation option: $1" >&2; usage >&2; exit 1 ;;
    esac
done

if [[ -f "${ENV_FILE}" ]]; then
    # shellcheck source=/dev/null
    source "${ENV_FILE}"
elif [[ -f "${ENV_TEMPLATE}" ]]; then
    # shellcheck source=/dev/null
    source "${ENV_TEMPLATE}"
fi

REFERENCE_FASTA="${REFERENCE_FASTA:-${QUILT2_REFERENCE_FASTA:-}}"

# Defaults are repeated here so older user-created environment.sh files remain usable.
QUILT2_WGS_TRUTH_FILTER_ENABLED="${QUILT2_WGS_TRUTH_FILTER_ENABLED:-true}"
QUILT2_WGS_TRUTH_MIN_QUAL="${QUILT2_WGS_TRUTH_MIN_QUAL:-30}"
QUILT2_WGS_TRUTH_MIN_QD="${QUILT2_WGS_TRUTH_MIN_QD:-2.0}"
QUILT2_WGS_TRUTH_MAX_SOR="${QUILT2_WGS_TRUTH_MAX_SOR:-3.0}"
QUILT2_WGS_TRUTH_MAX_FS="${QUILT2_WGS_TRUTH_MAX_FS:-60.0}"
QUILT2_WGS_TRUTH_MIN_MQ="${QUILT2_WGS_TRUTH_MIN_MQ:-40.0}"
QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM="${QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM:--12.5}"
QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM="${QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM:--8.0}"
QUILT2_WGS_TRUTH_MIN_GQ="${QUILT2_WGS_TRUTH_MIN_GQ:-60}"
QUILT2_WGS_TRUTH_MIN_DP="${QUILT2_WGS_TRUTH_MIN_DP:-10}"

is_numeric() {
    [[ "$1" =~ ^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$ ]]
}

validate_config() {
    if [[ "${QUILT2_WGS_TRUTH_FILTER_ENABLED}" != "true" && "${QUILT2_WGS_TRUTH_FILTER_ENABLED}" != "false" ]]; then
        echo "[ERROR] QUILT2_WGS_TRUTH_FILTER_ENABLED must be true or false." >&2
        return 1
    fi
    local name value
    for name in \
        QUILT2_WGS_TRUTH_MIN_QUAL QUILT2_WGS_TRUTH_MIN_QD QUILT2_WGS_TRUTH_MAX_SOR \
        QUILT2_WGS_TRUTH_MAX_FS QUILT2_WGS_TRUTH_MIN_MQ QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM \
        QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM; do
        value="${!name}"
        if ! is_numeric "${value}"; then
            echo "[ERROR] ${name} must be numeric (found '${value}')." >&2
            return 1
        fi
    done
    if [[ ! "${QUILT2_WGS_TRUTH_MIN_GQ}" =~ ^[0-9]+$ ]] || \
       (( QUILT2_WGS_TRUTH_MIN_GQ < 0 || QUILT2_WGS_TRUTH_MIN_GQ > 99 )); then
        echo "[ERROR] QUILT2_WGS_TRUTH_MIN_GQ must be an integer from 0 to 99." >&2
        return 1
    fi
    if [[ ! "${QUILT2_WGS_TRUTH_MIN_DP}" =~ ^[0-9]+$ ]]; then
        echo "[ERROR] QUILT2_WGS_TRUTH_MIN_DP must be a nonnegative integer." >&2
        return 1
    fi
}

validate_config

for required in IMPUTED TRUTH_DATASET_DIR OUT_PREFIX; do
    if [[ -z "${!required}" ]]; then
        echo "[ERROR] Missing required WGS argument: ${required}" >&2
        usage >&2
        exit 1
    fi
done
[[ -n "${REFERENCE_FASTA}" ]] || {
    echo "[ERROR] Set QUILT2_REFERENCE_FASTA in config/environment.sh or pass --reference-fasta." >&2
    exit 1
}
[[ -f "${IMPUTED}" ]] || { echo "[ERROR] Imputed VCF not found: ${IMPUTED}" >&2; exit 1; }
[[ -d "${TRUTH_DATASET_DIR}" ]] || { echo "[ERROR] Truth dataset directory not found: ${TRUTH_DATASET_DIR}" >&2; exit 1; }
[[ -f "${REFERENCE_FASTA}" ]] || { echo "[ERROR] Reference FASTA not found: ${REFERENCE_FASTA}" >&2; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "[ERROR] Reference FASTA index not found: ${REFERENCE_FASTA}.fai" >&2; exit 1; }
[[ -z "${SAMPLES_FILE}" || -f "${SAMPLES_FILE}" ]] || { echo "[ERROR] Sample file not found: ${SAMPLES_FILE}" >&2; exit 1; }
[[ -f "${R_HELPER}" ]] || { echo "[ERROR] Missing R helper: ${R_HELPER}" >&2; exit 1; }

# Match the existing array evaluator's Bunya environment setup.
CONDA_ENV="${CONDA_ENV:-myenv_py310}"
MINIFORGE_MODULE="${MINIFORGE_MODULE:-miniforge/25.3.0-3}"
BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"
if command -v module >/dev/null 2>&1; then
    if module load "${MINIFORGE_MODULE}" >/dev/null 2>&1; then
        echo "[INFO] Loaded ${MINIFORGE_MODULE} module" >&2
    else
        echo "[ERROR] Failed to load module: ${MINIFORGE_MODULE}" >&2
        exit 1
    fi
    if [[ -z "${ROOTMINIFORGE:-}" || ! -f "${ROOTMINIFORGE}/etc/profile.d/conda.sh" ]]; then
        echo "[ERROR] conda.sh was not found under ROOTMINIFORGE." >&2
        exit 1
    fi
    # shellcheck source=/dev/null
    source "${ROOTMINIFORGE}/etc/profile.d/conda.sh"
    if conda activate "${CONDA_ENV}" >/dev/null 2>&1; then
        echo "[INFO] Activated conda environment: ${CONDA_ENV}" >&2
    else
        echo "[ERROR] Failed to activate conda environment: ${CONDA_ENV}" >&2
        exit 1
    fi
fi
if ! command -v bcftools >/dev/null 2>&1; then
    if command -v module >/dev/null 2>&1; then
        module load "${BCFTOOLS_MODULE}"
    fi
fi
command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools is required." >&2; exit 1; }
command -v Rscript >/dev/null 2>&1 || { echo "[ERROR] Rscript is required." >&2; exit 1; }

has_index() {
    [[ -f "$1.csi" || -f "$1.tbi" ]]
}
has_index "${IMPUTED}" || { echo "[ERROR] Imputed VCF/BCF is not indexed: ${IMPUTED}" >&2; exit 1; }

canonical_chr() {
    local value="$1" number
    value="${value#chr}"
    value="${value#Chr}"
    [[ "${value}" =~ ^[0-9]+$ ]] || return 1
    number=$((10#${value}))
    (( number >= 0 && number <= 99 )) || return 1
    printf 'Chr%02d\n' "${number}"
}

declare -a TRUTH_CHROMS=()
declare -a TRUTH_FILES=()
truth_file_for_chr() {
    local wanted="$1" i
    for ((i=0; i<${#TRUTH_CHROMS[@]}; i++)); do
        if [[ "${TRUTH_CHROMS[$i]}" == "${wanted}" ]]; then
            printf '%s\n' "${TRUTH_FILES[$i]}"
            return 0
        fi
    done
    return 1
}
while IFS= read -r truth_vcf; do
    base="$(basename "${truth_vcf}")"
    if [[ "${base}" =~ ^(Chr[0-9]+)_consolidated[.]vcf[.]gz$ ]]; then
        chr="$(canonical_chr "${BASH_REMATCH[1]}")"
        has_index "${truth_vcf}" || { echo "[ERROR] Truth VCF is not indexed: ${truth_vcf}" >&2; exit 1; }
        if truth_file_for_chr "${chr}" >/dev/null 2>&1; then
            echo "[ERROR] Multiple truth VCFs resolve to ${chr}." >&2
            exit 1
        fi
        TRUTH_CHROMS+=("${chr}")
        TRUTH_FILES+=("${truth_vcf}")
    fi
done < <(find "${TRUTH_DATASET_DIR}" -maxdepth 1 -type f -name 'Chr*_consolidated.vcf.gz' -print)
(( ${#TRUTH_FILES[@]} > 0 )) || { echo "[ERROR] No indexed Chr*_consolidated.vcf.gz truth files found." >&2; exit 1; }

declare -a IMPUTED_CHROMS=()
while IFS=$'\t' read -r source_chr _; do
    if chr="$(canonical_chr "${source_chr}" 2>/dev/null)"; then
        IMPUTED_CHROMS+=("${chr}")
    fi
done < <(bcftools index -s "${IMPUTED}")

array_contains() {
    local wanted="$1" value
    shift
    for value in "$@"; do
        [[ "${value}" == "${wanted}" ]] && return 0
    done
    return 1
}

requested_chr00=false
declare -a REQUESTED_CHROMS=()
if [[ -n "${CHR_ARG}" ]]; then
    read -r -a raw_chroms <<< "${CHR_ARG//,/ }"
    for raw_chr in "${raw_chroms[@]}"; do
        chr="$(canonical_chr "${raw_chr}")" || { echo "[ERROR] Invalid chromosome: ${raw_chr}" >&2; exit 1; }
        [[ "${chr}" == "Chr00" ]] && requested_chr00=true
        REQUESTED_CHROMS+=("${chr}")
    done
fi

declare -a SELECTED_CHROMS=()
for number in $(seq 0 99); do
    printf -v chr 'Chr%02d' "${number}"
    truth_file_for_chr "${chr}" >/dev/null 2>&1 || continue
    array_contains "${chr}" "${IMPUTED_CHROMS[@]}" || continue
    if [[ ${#REQUESTED_CHROMS[@]} -gt 0 ]]; then
        include=false
        for requested in "${REQUESTED_CHROMS[@]}"; do
            [[ "${requested}" == "${chr}" ]] && include=true
        done
        [[ "${include}" == "true" ]] || continue
    elif [[ "${chr}" == "Chr00" ]]; then
        continue
    fi
    SELECTED_CHROMS+=("${chr}")
done
(( ${#SELECTED_CHROMS[@]} > 0 )) || { echo "[ERROR] No requested chromosomes are shared by imputed and truth inputs." >&2; exit 1; }
if [[ "${requested_chr00}" == "true" && ! " ${SELECTED_CHROMS[*]} " =~ " Chr00 " ]]; then
    echo "[WARN] Chr00 was requested but is not available on both input sides; it will not be evaluated." >&2
fi

SELECTED_CSV="$(IFS=,; echo "${SELECTED_CHROMS[*]}")"
NORMALIZED_REGION=""
if [[ -n "${REGION}" ]]; then
    region_chr="${REGION%%:*}"
    region_tail="${REGION#${region_chr}}"
    region_chr="$(canonical_chr "${region_chr}")" || { echo "[ERROR] Invalid region chromosome: ${REGION}" >&2; exit 1; }
    [[ " ${SELECTED_CHROMS[*]} " =~ " ${region_chr} " ]] || { echo "[ERROR] Region chromosome is not selected: ${region_chr}" >&2; exit 1; }
    NORMALIZED_REGION="${region_chr}${region_tail}"
fi

echo "[INFO] WGS truth filters: enabled=${QUILT2_WGS_TRUTH_FILTER_ENABLED}, QUAL>=${QUILT2_WGS_TRUTH_MIN_QUAL}, QD>=${QUILT2_WGS_TRUTH_MIN_QD}, SOR<=${QUILT2_WGS_TRUTH_MAX_SOR}, FS<=${QUILT2_WGS_TRUTH_MAX_FS}, MQ>=${QUILT2_WGS_TRUTH_MIN_MQ}, MQRankSum>=${QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM}, ReadPosRankSum>=${QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM}, GQ>=${QUILT2_WGS_TRUTH_MIN_GQ}, DP>=${QUILT2_WGS_TRUTH_MIN_DP}" >&2
echo "[INFO] WGS chromosomes: ${SELECTED_CHROMS[*]}" >&2

OUT_PREFIX="${OUT_PREFIX%/}"
MODE_FILE="${OUT_PREFIX}/.evaluation_mode"
COMPLETE_FILE="${OUT_PREFIX}/.complete"
if [[ -e "${MODE_FILE}" && "$(<"${MODE_FILE}")" != "wgs" ]]; then
    echo "[ERROR] Output belongs to a different evaluation mode: ${OUT_PREFIX}" >&2
    exit 1
fi
if [[ ! -e "${MODE_FILE}" && -d "${OUT_PREFIX}" && -n "$(find "${OUT_PREFIX}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
    echo "[ERROR] Non-empty output has no WGS mode marker; refusing cross-mode cache reuse: ${OUT_PREFIX}" >&2
    exit 1
fi
mkdir -p "${OUT_PREFIX}/intermediate"
printf 'wgs\n' > "${MODE_FILE}"

required_outputs=(
    per_variant_metrics.tsv per_sample_metrics.tsv filter_summary.tsv
    genotype_masking_summary.tsv matched_variants.tsv imputed_only_variants.tsv
    truth_only_variants.tsv allele_mismatches.tsv intermediate/imputed_ds.tsv
    intermediate/truth_gt_dosage.tsv run_manifest.tsv
)
if [[ "${FORCE}" != "true" && -f "${COMPLETE_FILE}" ]]; then
    complete=true
    for output in "${required_outputs[@]}"; do
        [[ -f "${OUT_PREFIX}/${output}" ]] || complete=false
    done
    if [[ "${complete}" == "true" ]]; then
        echo "[INFO] Complete WGS evaluation already exists; use --force to recompute: ${OUT_PREFIX}" >&2
        exit 0
    fi
fi
rm -f "${COMPLETE_FILE}"

TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/dosage_r2_wgs.XXXXXX")"
cleanup() {
    if [[ "${KEEP_TEMP}" == "true" ]]; then
        echo "[INFO] Temporary files retained: ${TMP_DIR}" >&2
    else
        rm -rf "${TMP_DIR}"
    fi
}
trap cleanup EXIT

make_rename_map() {
    local input="$1" map="$2" source_chr target_chr
    : > "${map}"
    while IFS=$'\t' read -r source_chr _; do
        if target_chr="$(canonical_chr "${source_chr}" 2>/dev/null)"; then
            printf '%s\t%s\n' "${source_chr}" "${target_chr}" >> "${map}"
        fi
    done < <(bcftools index -s "${input}")
}

normalize_vcf() {
    local input="$1" output="$2" prefix="$3"
    local rename_map="${TMP_DIR}/${prefix}.rename.tsv"
    make_rename_map "${input}" "${rename_map}"
    bcftools annotate --rename-chrs "${rename_map}" -Ou "${input}" \
        | bcftools view -t "${SELECTED_CSV}" -Ou \
        | bcftools view -m2 -M2 -v snps -Ou \
        | bcftools norm -f "${REFERENCE_FASTA}" -Ou \
        | bcftools norm -d exact -Oz -o "${output}"
    bcftools index -f "${output}"
    if [[ -n "${NORMALIZED_REGION}" ]]; then
        bcftools view -r "${NORMALIZED_REGION}" -Oz -o "${output}.region.vcf.gz" "${output}"
        mv "${output}.region.vcf.gz" "${output}"
        bcftools index -f "${output}"
    fi
}

NORMALIZED_IMPUTED="${TMP_DIR}/imputed.normalized.vcf.gz"
normalize_vcf "${IMPUTED}" "${NORMALIZED_IMPUTED}" imputed
if ! bcftools view -h "${NORMALIZED_IMPUTED}" | grep -q '^##FORMAT=<ID=DS,'; then
    echo "[ERROR] QUILT2 imputed VCF does not define FORMAT/DS." >&2
    exit 1
fi

declare -a NORMALIZED_TRUTH=()
for chr in "${SELECTED_CHROMS[@]}"; do
    truth_out="${TMP_DIR}/truth.${chr}.normalized.vcf.gz"
    truth_input="$(truth_file_for_chr "${chr}")"
    normalize_vcf "${truth_input}" "${truth_out}" "truth.${chr}"
    NORMALIZED_TRUTH+=("${truth_out}")
done

IMPUTED_SAMPLES="${TMP_DIR}/imputed.samples.txt"
TRUTH_SAMPLES="${TMP_DIR}/truth.samples.txt"
COMMON_SAMPLES="${TMP_DIR}/common.samples.txt"
bcftools query -l "${NORMALIZED_IMPUTED}" > "${IMPUTED_SAMPLES}"
bcftools query -l "${NORMALIZED_TRUTH[0]}" > "${TRUTH_SAMPLES}"
if [[ -n "${SAMPLES_FILE}" ]]; then
    awk 'NF {print $1}' "${SAMPLES_FILE}" > "${COMMON_SAMPLES}"
else
    while IFS= read -r sample; do
        grep -Fqx "${sample}" "${TRUTH_SAMPLES}" && printf '%s\n' "${sample}"
    done < "${IMPUTED_SAMPLES}" > "${COMMON_SAMPLES}"
fi
[[ -s "${COMMON_SAMPLES}" ]] || { echo "[ERROR] No common evaluation samples found." >&2; exit 1; }
if [[ "$(sort -u "${COMMON_SAMPLES}" | wc -l | tr -d ' ')" != "$(wc -l < "${COMMON_SAMPLES}" | tr -d ' ')" ]]; then
    echo "[ERROR] Sample list contains duplicates." >&2
    exit 1
fi
for sample_set in "${IMPUTED_SAMPLES}" "${TRUTH_SAMPLES}"; do
    while IFS= read -r sample; do
        grep -Fqx "${sample}" "${sample_set}" || { echo "[ERROR] Requested sample is absent from an input: ${sample}" >&2; exit 1; }
    done < "${COMMON_SAMPLES}"
done
for truth_vcf in "${NORMALIZED_TRUTH[@]}"; do
    current_samples="${TMP_DIR}/$(basename "${truth_vcf}").samples"
    bcftools query -l "${truth_vcf}" > "${current_samples}"
    while IFS= read -r sample; do
        grep -Fqx "${sample}" "${current_samples}" || { echo "[ERROR] Sample ${sample} is absent from ${truth_vcf}." >&2; exit 1; }
    done < "${COMMON_SAMPLES}"
done
SAMPLE_CSV="$(paste -sd, "${COMMON_SAMPLES}")"

IMPUTED_RAW="${TMP_DIR}/imputed.raw.tsv"
TRUTH_RAW="${TMP_DIR}/truth.raw.tsv"
{
    printf 'CHROM\tPOS\tREF\tALT\tID'
    while IFS= read -r sample; do printf '\t%s' "${sample}"; done < "${COMMON_SAMPLES}"
    printf '\n'
} > "${IMPUTED_RAW}"
bcftools query -s "${SAMPLE_CSV}" -f '%CHROM\t%POS\t%REF\t%ALT\t%ID[\t%DS]\n' "${NORMALIZED_IMPUTED}" >> "${IMPUTED_RAW}"

{
    printf 'CHROM\tPOS\tREF\tALT\tID\tQUAL\tQD\tSOR\tFS\tMQ\tMQRankSum\tReadPosRankSum'
    while IFS= read -r sample; do
        printf '\t%s.GT' "${sample}"
        if [[ "${QUILT2_WGS_TRUTH_FILTER_ENABLED}" == "true" ]]; then
            printf '\t%s.GQ\t%s.DP' "${sample}" "${sample}"
        fi
    done < "${COMMON_SAMPLES}"
    printf '\n'
} > "${TRUTH_RAW}"
if [[ "${QUILT2_WGS_TRUTH_FILTER_ENABLED}" == "true" ]]; then
    truth_format='%CHROM\t%POS\t%REF\t%ALT\t%ID\t%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum[\t%GT\t%GQ\t%DP]\n'
else
    truth_format='%CHROM\t%POS\t%REF\t%ALT\t%ID\t%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum[\t%GT]\n'
fi
for truth_vcf in "${NORMALIZED_TRUTH[@]}"; do
    bcftools query -u -s "${SAMPLE_CSV}" -f "${truth_format}" "${truth_vcf}" >> "${TRUTH_RAW}"
done

MANIFEST_TMP="${TMP_DIR}/run_manifest.tsv"
{
    printf 'key\tvalue\n'
    printf 'truth_mode\twgs\n'
    printf 'imputed\t%s\n' "${IMPUTED}"
    printf 'truth_dataset_dir\t%s\n' "${TRUTH_DATASET_DIR}"
    printf 'reference_fasta\t%s\n' "${REFERENCE_FASTA}"
    printf 'chromosomes\t%s\n' "${SELECTED_CSV}"
    printf 'region\t%s\n' "${NORMALIZED_REGION:-ALL}"
    printf 'samples\t%s\n' "${SAMPLE_CSV}"
    printf 'QUILT2_WGS_TRUTH_FILTER_ENABLED\t%s\n' "${QUILT2_WGS_TRUTH_FILTER_ENABLED}"
    printf 'QUILT2_WGS_TRUTH_MIN_QUAL\t%s\n' "${QUILT2_WGS_TRUTH_MIN_QUAL}"
    printf 'QUILT2_WGS_TRUTH_MIN_QD\t%s\n' "${QUILT2_WGS_TRUTH_MIN_QD}"
    printf 'QUILT2_WGS_TRUTH_MAX_SOR\t%s\n' "${QUILT2_WGS_TRUTH_MAX_SOR}"
    printf 'QUILT2_WGS_TRUTH_MAX_FS\t%s\n' "${QUILT2_WGS_TRUTH_MAX_FS}"
    printf 'QUILT2_WGS_TRUTH_MIN_MQ\t%s\n' "${QUILT2_WGS_TRUTH_MIN_MQ}"
    printf 'QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM\t%s\n' "${QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM}"
    printf 'QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM\t%s\n' "${QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM}"
    printf 'QUILT2_WGS_TRUTH_MIN_GQ\t%s\n' "${QUILT2_WGS_TRUTH_MIN_GQ}"
    printf 'QUILT2_WGS_TRUTH_MIN_DP\t%s\n' "${QUILT2_WGS_TRUTH_MIN_DP}"
} > "${MANIFEST_TMP}"

Rscript "${R_HELPER}" \
    --imputed-raw "${IMPUTED_RAW}" \
    --truth-raw "${TRUTH_RAW}" \
    --samples "${COMMON_SAMPLES}" \
    --out-dir "${OUT_PREFIX}" \
    --filter-enabled "${QUILT2_WGS_TRUTH_FILTER_ENABLED}" \
    --min-qual "${QUILT2_WGS_TRUTH_MIN_QUAL}" \
    --min-qd "${QUILT2_WGS_TRUTH_MIN_QD}" \
    --max-sor "${QUILT2_WGS_TRUTH_MAX_SOR}" \
    --max-fs "${QUILT2_WGS_TRUTH_MAX_FS}" \
    --min-mq "${QUILT2_WGS_TRUTH_MIN_MQ}" \
    --min-mq-rank-sum "${QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM}" \
    --min-read-pos-rank-sum "${QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM}" \
    --min-gq "${QUILT2_WGS_TRUTH_MIN_GQ}" \
    --min-dp "${QUILT2_WGS_TRUTH_MIN_DP}"

mv "${MANIFEST_TMP}" "${OUT_PREFIX}/run_manifest.tsv"
printf 'complete\n' > "${COMPLETE_FILE}"
echo "[INFO] WGS dosage evaluation complete: ${OUT_PREFIX}" >&2
