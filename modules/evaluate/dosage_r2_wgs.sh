#!/bin/bash
# WGS-truth evaluation: chromosome-scoped exact-allele QUILT2 GT versus filtered GATK GT.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
ENV_FILE="${ROOT_DIR}/config/environment.sh"
ENV_TEMPLATE="${ROOT_DIR}/config/environment.template.sh"
R_HELPER="${SCRIPT_DIR}/dosage_r2_wgs.R"
CONCAT_SCRIPT="${SCRIPT_DIR}/concat_imputed.sh"
OUTPUT_SCHEMA_VERSION="wgs-gt-isec-v4"

usage() {
    cat <<'EOF'
Usage: dosage_r2_wgs.sh (--imputed VCF | --chunks-dir DIR) --truth-dataset-dir DIR \
       --out-prefix DIR [options]

Required:
  --imputed PATH             Indexed QUILT2 VCF/BCF containing FORMAT/GT
  --chunks-dir PATH          Per-chunk QUILT2 output; mutually exclusive with --imputed
  --truth-dataset-dir PATH   GATK 7.Consolidated_VCF directory
  --out-prefix PATH          WGS evaluation output directory

Options:
  --chr LIST                 Comma/space-separated chromosome list
  --region STR               Canonical chromosome region, e.g. Chr01:1-1000000
  --samples FILE             One sample ID per line (default: sample intersection)
  --reference-fasta PATH     Override QUILT2_REFERENCE_FASTA from environment.sh
  --concat-force             Rebuild chromosome VCFs in --chunks-dir mode
  --force                    Recompute an existing WGS-mode output
  --keep-temp                Retain task-local normalization and extraction files
  -h, --help                 Show this help

Large WGS tables are chromosome-partitioned Parquet datasets. The array-compatible
per_sample_metrics.tsv and run_manifest.tsv are written at the evaluation root.
WGS filter thresholds are read only from config/environment.sh.
EOF
}

IMPUTED=""
CHUNKS_DIR=""
TRUTH_DATASET_DIR=""
REFERENCE_FASTA=""
OUT_PREFIX=""
CHR_ARG=""
REGION=""
SAMPLES_FILE=""
FORCE=false
CONCAT_FORCE=false
KEEP_TEMP=false

# Internal orchestration flags used by bin/dosage_r2_sbatch.sh.
PREPARE_ONLY=false
FINALIZE_ONLY=false
WORKER_TASK=""
WORKER_CHR=""
TASK_MANIFEST=""
ARRAY_JOB_ID=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --imputed) IMPUTED="$2"; shift 2 ;;
        --chunks-dir) CHUNKS_DIR="$2"; shift 2 ;;
        --truth-dataset-dir) TRUTH_DATASET_DIR="$2"; shift 2 ;;
        --reference-fasta) REFERENCE_FASTA="$2"; shift 2 ;;
        --out-prefix) OUT_PREFIX="$2"; shift 2 ;;
        --chr) CHR_ARG="$2"; shift 2 ;;
        --region) REGION="$2"; shift 2 ;;
        --samples) SAMPLES_FILE="$2"; shift 2 ;;
        --concat-force) CONCAT_FORCE=true; shift ;;
        --force) FORCE=true; shift ;;
        --keep-temp) KEEP_TEMP=true; shift ;;
        --prepare-only) PREPARE_ONLY=true; shift ;;
        --finalize-only) FINALIZE_ONLY=true; shift ;;
        --worker-task) WORKER_TASK="$2"; shift 2 ;;
        --worker-chr) WORKER_CHR="$2"; shift 2 ;;
        --task-manifest) TASK_MANIFEST="$2"; shift 2 ;;
        --array-job-id) ARRAY_JOB_ID="$2"; shift 2 ;;
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
QUILT2_WGS_TRUTH_FILTER_ENABLED="${QUILT2_WGS_TRUTH_FILTER_ENABLED:-true}"
if [[ -n "${QUILT2_WGS_KEEP_DOSAGE_MATRICES+x}" ]]; then
    echo "[WARN] QUILT2_WGS_KEEP_DOSAGE_MATRICES is obsolete and will be ignored; WGS GT evaluation does not write dosage matrices." >&2
    unset QUILT2_WGS_KEEP_DOSAGE_MATRICES
fi
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
        is_numeric "${value}" || { echo "[ERROR] ${name} must be numeric (found '${value}')." >&2; return 1; }
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
echo "[INFO] WGS truth filters: enabled=${QUILT2_WGS_TRUTH_FILTER_ENABLED}; QUAL>=${QUILT2_WGS_TRUTH_MIN_QUAL}; QD>=${QUILT2_WGS_TRUTH_MIN_QD}; SOR<=${QUILT2_WGS_TRUTH_MAX_SOR}; FS<=${QUILT2_WGS_TRUTH_MAX_FS}; MQ>=${QUILT2_WGS_TRUTH_MIN_MQ}; MQRankSum>=${QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM}; ReadPosRankSum>=${QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM}; GQ>=${QUILT2_WGS_TRUTH_MIN_GQ}; DP>=${QUILT2_WGS_TRUTH_MIN_DP}" >&2

if [[ -n "${IMPUTED}" && -n "${CHUNKS_DIR}" ]] || [[ -z "${IMPUTED}" && -z "${CHUNKS_DIR}" ]]; then
    echo "[ERROR] Specify exactly one of --imputed or --chunks-dir." >&2
    exit 1
fi
for required in TRUTH_DATASET_DIR OUT_PREFIX; do
    [[ -n "${!required}" ]] || { echo "[ERROR] Missing required WGS argument: ${required}" >&2; usage >&2; exit 1; }
done
[[ -n "${REFERENCE_FASTA}" ]] || { echo "[ERROR] Set QUILT2_REFERENCE_FASTA or pass --reference-fasta." >&2; exit 1; }
[[ -z "${IMPUTED}" || -f "${IMPUTED}" ]] || { echo "[ERROR] Imputed VCF not found: ${IMPUTED}" >&2; exit 1; }
[[ -z "${CHUNKS_DIR}" || -d "${CHUNKS_DIR}" ]] || { echo "[ERROR] Chunks directory not found: ${CHUNKS_DIR}" >&2; exit 1; }
[[ -d "${TRUTH_DATASET_DIR}" ]] || { echo "[ERROR] Truth dataset directory not found: ${TRUTH_DATASET_DIR}" >&2; exit 1; }
[[ -f "${REFERENCE_FASTA}" ]] || { echo "[ERROR] Reference FASTA not found: ${REFERENCE_FASTA}" >&2; exit 1; }
[[ -f "${REFERENCE_FASTA}.fai" ]] || { echo "[ERROR] Reference FASTA index not found: ${REFERENCE_FASTA}.fai" >&2; exit 1; }
[[ -z "${SAMPLES_FILE}" || -f "${SAMPLES_FILE}" ]] || { echo "[ERROR] Sample file not found: ${SAMPLES_FILE}" >&2; exit 1; }
[[ -f "${R_HELPER}" ]] || { echo "[ERROR] Missing R helper: ${R_HELPER}" >&2; exit 1; }
[[ -z "${CHUNKS_DIR}" || -f "${CONCAT_SCRIPT}" ]] || { echo "[ERROR] Missing concat helper: ${CONCAT_SCRIPT}" >&2; exit 1; }

OUT_PREFIX="${OUT_PREFIX%/}"
TRUTH_DATASET_DIR="$(cd "${TRUTH_DATASET_DIR}" && pwd)"
[[ -z "${CHUNKS_DIR}" ]] || CHUNKS_DIR="$(cd "${CHUNKS_DIR}" && pwd)"
MODE_FILE="${OUT_PREFIX}/.evaluation_mode"
SIGNATURE_FILE="${OUT_PREFIX}/.run_signature"
COMPLETE_FILE="${OUT_PREFIX}/.complete"
DEFAULT_TASK_MANIFEST="${OUT_PREFIX}/intermediate/chromosome_tasks.tsv"
TASK_MANIFEST="${TASK_MANIFEST:-${DEFAULT_TASK_MANIFEST}}"
COMMON_SAMPLES="${OUT_PREFIX}/intermediate/common_samples.txt"
PENDING_MANIFEST="${OUT_PREFIX}/intermediate/run_manifest.pending.tsv"

CONDA_ENV="${CONDA_ENV:-myenv_py310}"
MINIFORGE_MODULE="${MINIFORGE_MODULE:-miniforge/25.3.0-3}"
BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"

setup_tools() {
    if command -v module >/dev/null 2>&1; then
        if module load "${MINIFORGE_MODULE}" >/dev/null 2>&1; then
            echo "[INFO] Loaded ${MINIFORGE_MODULE} module" >&2
        else
            echo "[ERROR] Failed to load module: ${MINIFORGE_MODULE}" >&2
            exit 1
        fi
        [[ -n "${ROOTMINIFORGE:-}" && -f "${ROOTMINIFORGE}/etc/profile.d/conda.sh" ]] || {
            echo "[ERROR] conda.sh was not found under ROOTMINIFORGE." >&2; exit 1;
        }
        # shellcheck source=/dev/null
        source "${ROOTMINIFORGE}/etc/profile.d/conda.sh"
        conda activate "${CONDA_ENV}" >/dev/null 2>&1 || {
            echo "[ERROR] Failed to activate conda environment: ${CONDA_ENV}" >&2; exit 1;
        }
        echo "[INFO] Activated conda environment: ${CONDA_ENV}" >&2
    fi
    if ! command -v bcftools >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
        module load "${BCFTOOLS_MODULE}"
    fi
    command -v bcftools >/dev/null 2>&1 || { echo "[ERROR] bcftools is required." >&2; exit 1; }
    command -v Rscript >/dev/null 2>&1 || { echo "[ERROR] Rscript is required." >&2; exit 1; }
    Rscript -e 'quit(status=if (requireNamespace("data.table", quietly=TRUE) && requireNamespace("arrow", quietly=TRUE)) 0 else 1)' || {
        echo "[ERROR] R packages data.table and arrow are required in ${CONDA_ENV}." >&2
        exit 1
    }
}

setup_tools

has_index() { [[ -f "$1.csi" || -f "$1.tbi" ]]; }
[[ -z "${IMPUTED}" ]] || has_index "${IMPUTED}" || { echo "[ERROR] Imputed VCF/BCF is not indexed: ${IMPUTED}" >&2; exit 1; }

canonical_chr() {
    local value="$1" number
    value="${value#chr}"
    value="${value#Chr}"
    [[ "${value}" =~ ^[0-9]+$ ]] || return 1
    number=$((10#${value}))
    (( number >= 0 && number <= 99 )) || return 1
    printf 'Chr%02d\n' "${number}"
}

declare -a TRUTH_CHROMS=() TRUTH_FILES=()
truth_file_for_chr() {
    local wanted="$1" i
    for ((i=0; i<${#TRUTH_CHROMS[@]}; i++)); do
        [[ "${TRUTH_CHROMS[$i]}" == "${wanted}" ]] && { printf '%s\n' "${TRUTH_FILES[$i]}"; return 0; }
    done
    return 1
}

while IFS= read -r truth_vcf; do
    base="$(basename "${truth_vcf}")"
    if [[ "${base}" =~ ^(Chr[0-9]+)_consolidated[.]vcf[.]gz$ ]]; then
        chr="$(canonical_chr "${BASH_REMATCH[1]}")"
        has_index "${truth_vcf}" || { echo "[ERROR] Truth VCF is not indexed: ${truth_vcf}" >&2; exit 1; }
        truth_file_for_chr "${chr}" >/dev/null 2>&1 && { echo "[ERROR] Multiple truth VCFs resolve to ${chr}." >&2; exit 1; }
        TRUTH_CHROMS+=("${chr}")
        TRUTH_FILES+=("${truth_vcf}")
    fi
done < <(find "${TRUTH_DATASET_DIR}" -maxdepth 1 -type f -name 'Chr*_consolidated.vcf.gz' -print | sort)
(( ${#TRUTH_FILES[@]} > 0 )) || { echo "[ERROR] No indexed Chr*_consolidated.vcf.gz truth files found." >&2; exit 1; }

declare -a IMPUTED_CHROMS=() IMPUTED_CONTIGS=()
imputed_contig_for_chr() {
    local wanted="$1" i
    for ((i=0; i<${#IMPUTED_CHROMS[@]}; i++)); do
        [[ "${IMPUTED_CHROMS[$i]}" == "${wanted}" ]] && { printf '%s\n' "${IMPUTED_CONTIGS[$i]}"; return 0; }
    done
    return 1
}

CHUNK_MANIFEST=""
resolve_chunk_manifest() {
    local output_guess run_manifest detected
    output_guess="$(cd "${CHUNKS_DIR}/../.." 2>/dev/null && pwd || true)"
    run_manifest="${output_guess}/run_manifest.tsv"
    if [[ -f "${run_manifest}" ]]; then
        detected="$(awk -F'\t' '$1=="chunk_manifest"{print $2}' "${run_manifest}")"
        [[ -z "${detected}" || ! -f "${detected}" ]] || CHUNK_MANIFEST="${detected}"
    fi
}

if [[ -n "${IMPUTED}" ]]; then
    while IFS=$'\t' read -r source_chr _; do
        if chr="$(canonical_chr "${source_chr}" 2>/dev/null)"; then
            imputed_contig_for_chr "${chr}" >/dev/null 2>&1 && {
                echo "[ERROR] Multiple imputed contigs resolve to ${chr}." >&2; exit 1;
            }
            IMPUTED_CHROMS+=("${chr}")
            IMPUTED_CONTIGS+=("${source_chr}")
        fi
    done < <(bcftools index -s "${IMPUTED}")
else
    resolve_chunk_manifest
    if [[ -n "${CHUNK_MANIFEST}" ]]; then
        while IFS= read -r source_chr; do
            chr="$(canonical_chr "${source_chr}")" || continue
            if ! imputed_contig_for_chr "${chr}" >/dev/null 2>&1; then
                IMPUTED_CHROMS+=("${chr}")
                IMPUTED_CONTIGS+=("${source_chr}")
            fi
        done < <(awk -F'|' '{print $2}' "${CHUNK_MANIFEST}" | awk '!seen[$0]++')
    else
        while IFS= read -r chr_dir; do
            source_chr="$(basename "${chr_dir}")"
            chr="$(canonical_chr "${source_chr}")" || continue
            IMPUTED_CHROMS+=("${chr}")
            IMPUTED_CONTIGS+=("${source_chr}")
        done < <(find "${CHUNKS_DIR}" -mindepth 1 -maxdepth 1 -type d -print | sort)
    fi
fi
(( ${#IMPUTED_CHROMS[@]} > 0 )) || { echo "[ERROR] No imputed chromosomes were discovered." >&2; exit 1; }

declare -a REQUESTED_CHROMS=()
requested_chr00=false
if [[ -n "${CHR_ARG}" ]]; then
    read -r -a raw_chroms <<< "${CHR_ARG//,/ }"
    for raw_chr in "${raw_chroms[@]}"; do
        chr="$(canonical_chr "${raw_chr}")" || { echo "[ERROR] Invalid chromosome: ${raw_chr}" >&2; exit 1; }
        [[ "${chr}" == "Chr00" ]] && requested_chr00=true
        REQUESTED_CHROMS+=("${chr}")
    done
fi

NORMALIZED_REGION=""
REGION_CHR=""
if [[ -n "${REGION}" ]]; then
    region_source="${REGION%%:*}"
    region_tail="${REGION#${region_source}}"
    REGION_CHR="$(canonical_chr "${region_source}")" || { echo "[ERROR] Invalid region chromosome: ${REGION}" >&2; exit 1; }
    NORMALIZED_REGION="${REGION_CHR}${region_tail}"
    [[ "${REGION_CHR}" != "Chr00" ]] || requested_chr00=true
    if (( ${#REQUESTED_CHROMS[@]} > 0 )); then
        found=false
        for chr in "${REQUESTED_CHROMS[@]}"; do [[ "${chr}" == "${REGION_CHR}" ]] && found=true; done
        [[ "${found}" == "true" ]] || { echo "[ERROR] Region chromosome is not selected: ${REGION_CHR}" >&2; exit 1; }
    else
        REQUESTED_CHROMS=("${REGION_CHR}")
    fi
fi

is_requested() {
    local wanted="$1" value
    (( ${#REQUESTED_CHROMS[@]} == 0 )) && return 0
    for value in "${REQUESTED_CHROMS[@]}"; do [[ "${value}" == "${wanted}" ]] && return 0; done
    return 1
}

declare -a SELECTED_CHROMS=()
for number in $(seq 0 99); do
    printf -v chr 'Chr%02d' "${number}"
    truth_file_for_chr "${chr}" >/dev/null 2>&1 || continue
    imputed_contig_for_chr "${chr}" >/dev/null 2>&1 || continue
    is_requested "${chr}" || continue
    [[ "${chr}" != "Chr00" || "${requested_chr00}" == "true" ]] || continue
    SELECTED_CHROMS+=("${chr}")
done
(( ${#SELECTED_CHROMS[@]} > 0 )) || { echo "[ERROR] No requested chromosomes are shared by imputed and truth inputs." >&2; exit 1; }
if [[ "${requested_chr00}" == "true" && ! " ${SELECTED_CHROMS[*]} " =~ " Chr00 " ]]; then
    echo "[WARN] Chr00 was requested but is unavailable on both input sides." >&2
fi

TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/dosage_r2_wgs.XXXXXX")"
cleanup() {
    if [[ "${KEEP_TEMP}" == "true" ]]; then
        echo "[INFO] Temporary files retained: ${TMP_DIR}" >&2
    else
        rm -rf "${TMP_DIR}"
    fi
}
trap cleanup EXIT

check_output_mode() {
    if [[ -e "${MODE_FILE}" && "$(<"${MODE_FILE}")" != "wgs" ]]; then
        echo "[ERROR] Output belongs to a different evaluation mode: ${OUT_PREFIX}" >&2
        exit 1
    fi
    if [[ ! -e "${MODE_FILE}" && -d "${OUT_PREFIX}" && -n "$(find "${OUT_PREFIX}" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
        echo "[ERROR] Non-empty output has no WGS mode marker; refusing cross-mode cache reuse: ${OUT_PREFIX}" >&2
        exit 1
    fi
}

write_task_manifest() {
    local path="$1" task=0 chromosome truth_vcf imputed_contig
    printf 'task_id\tchromosome\ttruth_vcf\timputed_contig\n' > "${path}"
    for chromosome in "${SELECTED_CHROMS[@]}"; do
        task=$((task + 1))
        truth_vcf="$(truth_file_for_chr "${chromosome}")"
        imputed_contig="$(imputed_contig_for_chr "${chromosome}")"
        printf '%d\t%s\t%s\t%s\n' "${task}" "${chromosome}" "${truth_vcf}" "${imputed_contig}" >> "${path}"
    done
}

representative_imputed() {
    if [[ -n "${IMPUTED}" ]]; then
        printf '%s\n' "${IMPUTED}"
        return
    fi
    local first_chr source_chr
    first_chr="${SELECTED_CHROMS[0]}"
    source_chr="$(imputed_contig_for_chr "${first_chr}")"
    find "${CHUNKS_DIR}/${source_chr}" -maxdepth 1 -type f -name "quilt2.diploid.${source_chr}.*.vcf.gz" -print | sort -V | sed -n '1p'
}

file_metadata() {
    local path="$1" meta
    if meta="$(stat -c '%s:%Y' "${path}" 2>/dev/null)"; then
        printf '%s:%s\n' "${path}" "${meta}"
    else
        meta="$(stat -f '%z:%m' "${path}")"
        printf '%s:%s\n' "${path}" "${meta}"
    fi
}

hash_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    elif command -v shasum >/dev/null 2>&1; then
        shasum -a 256 "$1" | awk '{print $1}'
    else
        cksum "$1" | awk '{print $1":"$2}'
    fi
}

prepare_run() {
    check_output_mode
    local tasks_tmp="${TMP_DIR}/chromosome_tasks.tsv"
    local samples_tmp="${TMP_DIR}/common_samples.txt"
    local imputed_samples="${TMP_DIR}/imputed.samples.txt"
    local truth_samples="${TMP_DIR}/truth.samples.txt"
    local signature_payload="${TMP_DIR}/signature.txt"
    local representative truth_first selected_csv sample_csv new_signature existing_signature=""
    local cache_complete=true chromosome checkpoint path

    write_task_manifest "${tasks_tmp}"
    representative="$(representative_imputed)"
    [[ -n "${representative}" && -f "${representative}" ]] || { echo "[ERROR] No representative imputed VCF was found for sample discovery." >&2; exit 1; }
    truth_first="$(truth_file_for_chr "${SELECTED_CHROMS[0]}")"
    bcftools query -l "${representative}" > "${imputed_samples}"
    bcftools query -l "${truth_first}" > "${truth_samples}"
    if [[ -n "${SAMPLES_FILE}" ]]; then
        awk 'NF {print $1}' "${SAMPLES_FILE}" > "${samples_tmp}"
    else
        while IFS= read -r sample; do
            grep -Fqx "${sample}" "${truth_samples}" && printf '%s\n' "${sample}"
        done < "${imputed_samples}" > "${samples_tmp}"
    fi
    [[ -s "${samples_tmp}" ]] || { echo "[ERROR] No common evaluation samples found." >&2; exit 1; }
    [[ "$(sort -u "${samples_tmp}" | wc -l | tr -d ' ')" == "$(wc -l < "${samples_tmp}" | tr -d ' ')" ]] || {
        echo "[ERROR] Sample list contains duplicates." >&2; exit 1;
    }
    for sample_set in "${imputed_samples}" "${truth_samples}"; do
        while IFS= read -r sample; do
            grep -Fqx "${sample}" "${sample_set}" || { echo "[ERROR] Requested sample is absent from an input: ${sample}" >&2; exit 1; }
        done < "${samples_tmp}"
    done

    selected_csv="$(IFS=,; echo "${SELECTED_CHROMS[*]}")"
    sample_csv="$(paste -sd, "${samples_tmp}")"
    {
        printf 'schema=%s\n' "${OUTPUT_SCHEMA_VERSION}"
        printf 'input_mode=%s\n' "$([[ -n "${IMPUTED}" ]] && echo imputed || echo chunks)"
        [[ -z "${IMPUTED}" ]] || file_metadata "${IMPUTED}"
        if [[ -n "${CHUNKS_DIR}" ]]; then
            printf 'chunks_dir=%s\n' "${CHUNKS_DIR}"
            [[ -z "${CHUNK_MANIFEST}" ]] || file_metadata "${CHUNK_MANIFEST}"
            while IFS= read -r chunk; do file_metadata "${chunk}"; done < <(find "${CHUNKS_DIR}" -type f -name 'quilt2.diploid.*.vcf.gz' -print | sort)
        fi
        while IFS=$'\t' read -r _ _ truth_vcf _; do [[ "${truth_vcf}" == "truth_vcf" ]] || file_metadata "${truth_vcf}"; done < "${tasks_tmp}"
        file_metadata "${REFERENCE_FASTA}"
        file_metadata "${REFERENCE_FASTA}.fai"
        printf 'chromosomes=%s\nregion=%s\nsamples=%s\n' "${selected_csv}" "${NORMALIZED_REGION:-ALL}" "${sample_csv}"
        printf 'filter=%s|%s|%s|%s|%s|%s|%s|%s|%s|%s\n' \
            "${QUILT2_WGS_TRUTH_FILTER_ENABLED}" "${QUILT2_WGS_TRUTH_MIN_QUAL}" "${QUILT2_WGS_TRUTH_MIN_QD}" \
            "${QUILT2_WGS_TRUTH_MAX_SOR}" "${QUILT2_WGS_TRUTH_MAX_FS}" "${QUILT2_WGS_TRUTH_MIN_MQ}" \
            "${QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM}" "${QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM}" \
            "${QUILT2_WGS_TRUTH_MIN_GQ}" "${QUILT2_WGS_TRUTH_MIN_DP}"
        cat "${tasks_tmp}" "${samples_tmp}"
    } > "${signature_payload}"
    new_signature="$(hash_file "${signature_payload}")"
    [[ ! -f "${SIGNATURE_FILE}" ]] || existing_signature="$(<"${SIGNATURE_FILE}")"
    if [[ "${FORCE}" != "true" && -z "${existing_signature}" && \
          ( -f "${COMPLETE_FILE}" || -e "${OUT_PREFIX}/per_sample_metrics.tsv" || -e "${OUT_PREFIX}/per_variant_metrics.tsv" ) ]]; then
        echo "[ERROR] Existing WGS output predates run signatures; use --force to migrate it: ${OUT_PREFIX}" >&2
        exit 1
    fi
    if [[ "${FORCE}" != "true" && -n "${existing_signature}" && "${existing_signature}" != "${new_signature}" ]]; then
        echo "[ERROR] Existing WGS output has a different run signature; use --force: ${OUT_PREFIX}" >&2
        exit 1
    fi

    if [[ "${FORCE}" == "true" ]]; then
        rm -rf "${OUT_PREFIX}/metrics" "${OUT_PREFIX}/qc" "${OUT_PREFIX}/intermediate"
        rm -f \
            "${OUT_PREFIX}/per_sample_metrics.tsv" "${OUT_PREFIX}/run_manifest.tsv" "${COMPLETE_FILE}" \
            "${OUT_PREFIX}/per_variant_metrics.tsv" "${OUT_PREFIX}/filter_summary.tsv" \
            "${OUT_PREFIX}/genotype_masking_summary.tsv" "${OUT_PREFIX}/matched_variants.tsv" \
            "${OUT_PREFIX}/site_filtered_variants.tsv" \
            "${OUT_PREFIX}/imputed_only_variants.tsv" "${OUT_PREFIX}/truth_only_variants.tsv" \
            "${OUT_PREFIX}/allele_mismatches.tsv"
    fi
    mkdir -p "${OUT_PREFIX}/metrics" "${OUT_PREFIX}/qc/chromosomes" \
        "${OUT_PREFIX}/intermediate/chromosome_stats" "${OUT_PREFIX}/intermediate/checkpoints" \
        "${OUT_PREFIX}/slurm"
    printf 'wgs\n' > "${MODE_FILE}"
    cp "${tasks_tmp}" "${TASK_MANIFEST}.tmp"
    mv "${TASK_MANIFEST}.tmp" "${TASK_MANIFEST}"
    cp "${samples_tmp}" "${COMMON_SAMPLES}.tmp"
    mv "${COMMON_SAMPLES}.tmp" "${COMMON_SAMPLES}"
    printf '%s\n' "${new_signature}" > "${SIGNATURE_FILE}.tmp"
    mv "${SIGNATURE_FILE}.tmp" "${SIGNATURE_FILE}"

    {
        printf 'key\tvalue\n'
        printf 'truth_mode\twgs\n'
        printf 'output_schema\t%s\n' "${OUTPUT_SCHEMA_VERSION}"
        printf 'output_format\tpartitioned_parquet_snappy\n'
        printf 'comparison_field\tGT\n'
        printf 'genotype_encoding\tALT_COUNT_0_1_2\n'
        printf 'intersection_key\tCHROM:POS:REF:ALT\n'
        printf 'intersection_tool\tbcftools_isec\n'
        printf 'imputed\t%s\n' "${IMPUTED:-NA}"
        printf 'chunks_dir\t%s\n' "${CHUNKS_DIR:-NA}"
        printf 'truth_dataset_dir\t%s\n' "${TRUTH_DATASET_DIR}"
        printf 'reference_fasta\t%s\n' "${REFERENCE_FASTA}"
        printf 'chromosomes\t%s\n' "${selected_csv}"
        printf 'region\t%s\n' "${NORMALIZED_REGION:-ALL}"
        printf 'samples\t%s\n' "${sample_csv}"
        printf 'run_signature\t%s\n' "${new_signature}"
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
    } > "${PENDING_MANIFEST}.tmp"
    mv "${PENDING_MANIFEST}.tmp" "${PENDING_MANIFEST}"

    if [[ "${FORCE}" != "true" && -f "${COMPLETE_FILE}" && -s "${OUT_PREFIX}/per_sample_metrics.tsv" && -s "${OUT_PREFIX}/run_manifest.tsv" ]]; then
        for chromosome in "${SELECTED_CHROMS[@]}"; do
            checkpoint="${OUT_PREFIX}/intermediate/checkpoints/${chromosome}.done"
            [[ -f "${checkpoint}" && "$(<"${checkpoint}")" == "${new_signature}" ]] || cache_complete=false
            while IFS= read -r path; do [[ -s "${path}" ]] || cache_complete=false; done < <(required_partition_paths "${chromosome}")
        done
        if [[ "${cache_complete}" == "true" ]]; then
            echo "[INFO] Complete WGS evaluation already exists; use --force to recompute: ${OUT_PREFIX}" >&2
            return 0
        fi
    fi
    rm -f "${COMPLETE_FILE}"
    echo "[INFO] Prepared WGS evaluation with ${#SELECTED_CHROMS[@]} chromosome task(s): ${TASK_MANIFEST}" >&2
}

find_source_contig() {
    local input="$1" wanted="$2" source _ canonical
    while IFS=$'\t' read -r source _; do
        canonical="$(canonical_chr "${source}" 2>/dev/null || true)"
        [[ "${canonical}" == "${wanted}" ]] && { printf '%s\n' "${source}"; return 0; }
    done < <(bcftools index -s "${input}")
    return 1
}

normalize_chromosome() {
    local input="$1" chromosome="$2" output="$3" prefix="$4"
    local source_chr source_region canonical threads rename_map="${TMP_DIR}/${prefix}.rename.tsv" region_tail=""
    source_chr="$(find_source_contig "${input}" "${chromosome}")" || {
        echo "[ERROR] ${chromosome} was not found in ${input}." >&2; exit 1;
    }
    [[ -z "${NORMALIZED_REGION}" ]] || region_tail="${NORMALIZED_REGION#${chromosome}}"
    source_region="${source_chr}${region_tail}"
    : > "${rename_map}"
    while IFS=$'\t' read -r source _; do
        canonical="$(canonical_chr "${source}" 2>/dev/null || true)"
        [[ -z "${canonical}" ]] || printf '%s\t%s\n' "${source}" "${canonical}" >> "${rename_map}"
    done < <(bcftools index -s "${input}")
    threads="${SLURM_CPUS_PER_TASK:-1}"
    bcftools view -r "${source_region}" -Ou "${input}" \
        | bcftools annotate --rename-chrs "${rename_map}" -Ou \
        | bcftools view -m2 -M2 -v snps -Ou \
        | bcftools norm -f "${REFERENCE_FASTA}" -Ou \
        | bcftools norm -d exact --threads "${threads}" -Oz -o "${output}"
    bcftools index -f -c "${output}"
}

read_task() {
    local row
    if [[ -n "${WORKER_TASK}" ]]; then
        row="$(awk -F'\t' -v task="${WORKER_TASK}" 'NR>1 && $1==task {print; exit}' "${TASK_MANIFEST}")"
    else
        row="$(awk -F'\t' -v chr="${WORKER_CHR}" 'NR>1 && $2==chr {print; exit}' "${TASK_MANIFEST}")"
    fi
    [[ -n "${row}" ]] || { echo "[ERROR] Worker task was not found in ${TASK_MANIFEST}." >&2; exit 1; }
    IFS=$'\t' read -r TASK_ID TASK_CHR TASK_TRUTH_VCF TASK_IMPUTED_CONTIG <<< "${row}"
}

required_partition_paths() {
    local chromosome="$1"
    printf '%s\n' \
        "${OUT_PREFIX}/metrics/per_variant_metrics/CHROM=${chromosome}/part-000.parquet" \
        "${OUT_PREFIX}/metrics/site_filtered_variants/CHROM=${chromosome}/part-000.parquet" \
        "${OUT_PREFIX}/metrics/imputed_only_variants/CHROM=${chromosome}/part-000.parquet" \
        "${OUT_PREFIX}/metrics/truth_only_variants/CHROM=${chromosome}/part-000.parquet" \
        "${OUT_PREFIX}/metrics/allele_mismatches/CHROM=${chromosome}/part-000.parquet"
    printf '%s\n' \
        "${OUT_PREFIX}/intermediate/chromosome_stats/${chromosome}.sample_stats.tsv" \
        "${OUT_PREFIX}/qc/chromosomes/${chromosome}.filter_summary.tsv" \
        "${OUT_PREFIX}/qc/chromosomes/${chromosome}.genotype_masking_summary.tsv"
}

run_worker() {
    check_output_mode
    [[ -f "${TASK_MANIFEST}" && -f "${SIGNATURE_FILE}" && -f "${COMMON_SAMPLES}" ]] || {
        echo "[ERROR] WGS run is not prepared: ${OUT_PREFIX}" >&2; exit 1;
    }
    read_task
    local run_signature checkpoint="${OUT_PREFIX}/intermediate/checkpoints/${TASK_CHR}.done" complete=true path
    run_signature="$(<"${SIGNATURE_FILE}")"
    if [[ "${FORCE}" != "true" && -f "${checkpoint}" && "$(<"${checkpoint}")" == "${run_signature}" ]]; then
        while IFS= read -r path; do [[ -s "${path}" ]] || complete=false; done < <(required_partition_paths "${TASK_CHR}")
        if [[ "${complete}" == "true" ]]; then
            echo "[INFO] Complete chromosome checkpoint exists; skipping ${TASK_CHR}." >&2
            return 0
        fi
    fi
    rm -f "${checkpoint}"

    local worker_imputed="${IMPUTED}" concat_args normalized_imputed normalized_truth
    if [[ -n "${CHUNKS_DIR}" ]]; then
        concat_args=(--chunks-dir "${CHUNKS_DIR}" --chr "${TASK_IMPUTED_CONTIG}")
        [[ "${CONCAT_FORCE}" == "true" ]] && concat_args+=(--force)
        worker_imputed="$(bash "${CONCAT_SCRIPT}" "${concat_args[@]}")"
    fi
    has_index "${worker_imputed}" || { echo "[ERROR] Imputed VCF is not indexed: ${worker_imputed}" >&2; exit 1; }
    normalized_imputed="${TMP_DIR}/imputed.${TASK_CHR}.normalized.vcf.gz"
    normalized_truth="${TMP_DIR}/truth.${TASK_CHR}.normalized.vcf.gz"
    normalize_chromosome "${worker_imputed}" "${TASK_CHR}" "${normalized_imputed}" imputed
    normalize_chromosome "${TASK_TRUTH_VCF}" "${TASK_CHR}" "${normalized_truth}" truth
    for input in "${normalized_imputed}" "${normalized_truth}"; do
        bcftools view -h "${input}" > "${TMP_DIR}/$(basename "${input}").header.txt"
        grep -q '^##FORMAT=<ID=GT,' "${TMP_DIR}/$(basename "${input}").header.txt" || {
            echo "[ERROR] Input VCF does not define FORMAT/GT: ${input}" >&2
            exit 1
        }
    done

    local current_samples="${TMP_DIR}/current.samples.txt" sample_csv
    local isec_dir="${TMP_DIR}/isec.${TASK_CHR}" imputed_common_raw truth_common_raw
    local imputed_only_raw truth_only_raw truth_format threads input isec_part
    for input in "${normalized_imputed}" "${normalized_truth}"; do
        bcftools query -l "${input}" > "${current_samples}"
        while IFS= read -r sample; do
            grep -Fqx "${sample}" "${current_samples}" || { echo "[ERROR] Sample ${sample} is absent from ${input}." >&2; exit 1; }
        done < "${COMMON_SAMPLES}"
    done
    sample_csv="$(paste -sd, "${COMMON_SAMPLES}")"
    threads="${SLURM_CPUS_PER_TASK:-1}"
    bcftools isec --threads "${threads}" -c none -Oz -p "${isec_dir}" \
        "${normalized_imputed}" "${normalized_truth}"
    for isec_part in 0000 0001 0002 0003; do
        [[ -f "${isec_dir}/${isec_part}.vcf.gz" ]] || {
            echo "[ERROR] bcftools isec did not create ${isec_part}.vcf.gz for ${TASK_CHR}." >&2
            exit 1
        }
    done

    imputed_only_raw="${TMP_DIR}/imputed.only.tsv"
    truth_only_raw="${TMP_DIR}/truth.only.tsv"
    imputed_common_raw="${TMP_DIR}/imputed.common.tsv"
    truth_common_raw="${TMP_DIR}/truth.common.tsv"
    {
        printf 'CHROM\tPOS\tREF\tALT\tID\n'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' "${isec_dir}/0000.vcf.gz"
    } > "${imputed_only_raw}"
    {
        printf 'CHROM\tPOS\tREF\tALT\tID\n'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' "${isec_dir}/0001.vcf.gz"
    } > "${truth_only_raw}"
    {
        printf 'CHROM\tPOS\tREF\tALT\tID'
        while IFS= read -r sample; do printf '\t%s' "${sample}"; done < "${COMMON_SAMPLES}"
        printf '\n'
        bcftools query -s "${sample_csv}" -f '%CHROM\t%POS\t%REF\t%ALT\t%ID[\t%GT]\n' "${isec_dir}/0002.vcf.gz"
    } > "${imputed_common_raw}"
    {
        printf 'CHROM\tPOS\tREF\tALT\tID\tQUAL\tQD\tSOR\tFS\tMQ\tMQRankSum\tReadPosRankSum'
        while IFS= read -r sample; do
            printf '\t%s.GT' "${sample}"
            [[ "${QUILT2_WGS_TRUTH_FILTER_ENABLED}" != "true" ]] || printf '\t%s.GQ\t%s.DP' "${sample}" "${sample}"
        done < "${COMMON_SAMPLES}"
        printf '\n'
        if [[ "${QUILT2_WGS_TRUTH_FILTER_ENABLED}" == "true" ]]; then
            truth_format='%CHROM\t%POS\t%REF\t%ALT\t%ID\t%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum[\t%GT\t%GQ\t%DP]\n'
        else
            truth_format='%CHROM\t%POS\t%REF\t%ALT\t%ID\t%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum[\t%GT]\n'
        fi
        bcftools query -u -s "${sample_csv}" -f "${truth_format}" "${isec_dir}/0003.vcf.gz"
    } > "${truth_common_raw}"

    Rscript "${R_HELPER}" \
        --mode chromosome --chromosome "${TASK_CHR}" \
        --imputed-common-raw "${imputed_common_raw}" --truth-common-raw "${truth_common_raw}" \
        --imputed-only-raw "${imputed_only_raw}" --truth-only-raw "${truth_only_raw}" \
        --samples "${COMMON_SAMPLES}" --out-dir "${OUT_PREFIX}" \
        --filter-enabled "${QUILT2_WGS_TRUTH_FILTER_ENABLED}" \
        --min-qual "${QUILT2_WGS_TRUTH_MIN_QUAL}" --min-qd "${QUILT2_WGS_TRUTH_MIN_QD}" \
        --max-sor "${QUILT2_WGS_TRUTH_MAX_SOR}" --max-fs "${QUILT2_WGS_TRUTH_MAX_FS}" \
        --min-mq "${QUILT2_WGS_TRUTH_MIN_MQ}" --min-mq-rank-sum "${QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM}" \
        --min-read-pos-rank-sum "${QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM}" \
        --min-gq "${QUILT2_WGS_TRUTH_MIN_GQ}" --min-dp "${QUILT2_WGS_TRUTH_MIN_DP}"

    while IFS= read -r path; do [[ -s "${path}" ]] || { echo "[ERROR] Missing chromosome output: ${path}" >&2; exit 1; }; done < <(required_partition_paths "${TASK_CHR}")
    printf '%s\n' "${run_signature}" > "${checkpoint}.tmp.$$"
    mv "${checkpoint}.tmp.$$" "${checkpoint}"
    echo "[INFO] Completed WGS chromosome task ${TASK_ID}: ${TASK_CHR}" >&2
}

run_finalize() {
    check_output_mode
    [[ -f "${TASK_MANIFEST}" && -f "${SIGNATURE_FILE}" && -f "${PENDING_MANIFEST}" ]] || {
        echo "[ERROR] WGS run is not prepared for finalization: ${OUT_PREFIX}" >&2; exit 1;
    }
    local run_signature chromosome checkpoint path
    run_signature="$(<"${SIGNATURE_FILE}")"
    while IFS=$'\t' read -r _ chromosome _ _; do
        [[ "${chromosome}" == "chromosome" ]] && continue
        checkpoint="${OUT_PREFIX}/intermediate/checkpoints/${chromosome}.done"
        [[ -f "${checkpoint}" && "$(<"${checkpoint}")" == "${run_signature}" ]] || {
            echo "[ERROR] Missing or stale chromosome checkpoint: ${chromosome}" >&2; exit 1;
        }
        while IFS= read -r path; do [[ -s "${path}" ]] || { echo "[ERROR] Missing chromosome output: ${path}" >&2; exit 1; }; done < <(required_partition_paths "${chromosome}")
    done < "${TASK_MANIFEST}"

    Rscript "${R_HELPER}" \
        --mode finalize --out-dir "${OUT_PREFIX}" --chromosome-manifest "${TASK_MANIFEST}" \
        --filter-enabled "${QUILT2_WGS_TRUTH_FILTER_ENABLED}" \
        --min-qual "${QUILT2_WGS_TRUTH_MIN_QUAL}" --min-qd "${QUILT2_WGS_TRUTH_MIN_QD}" \
        --max-sor "${QUILT2_WGS_TRUTH_MAX_SOR}" --max-fs "${QUILT2_WGS_TRUTH_MAX_FS}" \
        --min-mq "${QUILT2_WGS_TRUTH_MIN_MQ}" --min-mq-rank-sum "${QUILT2_WGS_TRUTH_MIN_MQ_RANK_SUM}" \
        --min-read-pos-rank-sum "${QUILT2_WGS_TRUTH_MIN_READ_POS_RANK_SUM}" \
        --min-gq "${QUILT2_WGS_TRUTH_MIN_GQ}" --min-dp "${QUILT2_WGS_TRUTH_MIN_DP}"

    cp "${PENDING_MANIFEST}" "${OUT_PREFIX}/run_manifest.tsv.tmp"
    {
        printf 'array_job_id\t%s\n' "${ARRAY_JOB_ID:-NA}"
        printf 'finalizer_job_id\t%s\n' "${SLURM_JOB_ID:-NA}"
        printf 'per_sample_metrics\t%s\n' "${OUT_PREFIX}/per_sample_metrics.tsv"
        printf 'per_variant_metrics\t%s\n' "${OUT_PREFIX}/metrics/per_variant_metrics"
        printf 'site_filtered_variants\t%s\n' "${OUT_PREFIX}/metrics/site_filtered_variants"
        printf 'imputed_only_variants\t%s\n' "${OUT_PREFIX}/metrics/imputed_only_variants"
        printf 'truth_only_variants\t%s\n' "${OUT_PREFIX}/metrics/truth_only_variants"
        printf 'allele_mismatches\t%s\n' "${OUT_PREFIX}/metrics/allele_mismatches"
    } >> "${OUT_PREFIX}/run_manifest.tsv.tmp"
    mv "${OUT_PREFIX}/run_manifest.tsv.tmp" "${OUT_PREFIX}/run_manifest.tsv"
    printf 'complete\n' > "${COMPLETE_FILE}.tmp.$$"
    mv "${COMPLETE_FILE}.tmp.$$" "${COMPLETE_FILE}"
    echo "[INFO] WGS GT evaluation complete: ${OUT_PREFIX}" >&2
}

if [[ "${PREPARE_ONLY}" == "true" ]]; then
    prepare_run
    exit 0
fi
if [[ -n "${WORKER_TASK}" || -n "${WORKER_CHR}" ]]; then
    run_worker
    exit 0
fi
if [[ "${FINALIZE_ONLY}" == "true" ]]; then
    run_finalize
    exit 0
fi

# Direct/in-allocation compatibility path: prepare, process chromosomes sequentially, finalize.
prepare_run
if [[ -f "${COMPLETE_FILE}" && "${FORCE}" != "true" ]]; then
    exit 0
fi
while IFS=$'\t' read -r task_id chromosome _ _; do
    [[ "${task_id}" == "task_id" ]] && continue
    WORKER_TASK="${task_id}"
    WORKER_CHR="${chromosome}"
    run_worker
done < "${TASK_MANIFEST}"
run_finalize
