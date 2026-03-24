#!/bin/bash
# Submit helper: wraps modules/evaluate/dosage_r2.sh in sbatch with logs placed beside the outputs.
# Usage:
#   bash bin/dosage_r2_sbatch.sh --imputed /path/imputed.vcf.gz --truth /path/truth.vcf.gz [--out-prefix prefix] [-- extra flags]
# If run under SLURM (SLURM_JOB_ID set), it will execute dosage_r2.sh directly (no re-submit).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DOSAGE_SCRIPT="${ROOT_DIR}/modules/evaluate/dosage_r2.sh"
ENV_FILE="${ROOT_DIR}/config/environment.sh"
ENV_TEMPLATE="${ROOT_DIR}/config/environment.template.sh"

usage() {
    cat <<'EOF'
Usage:
  bash bin/dosage_r2_sbatch.sh --imputed IMPT_VCF --truth TRUTH_VCF [--out-prefix OUT_PREFIX] [-- extra args]

Notes:
  - OUT_PREFIX defaults to the imputed VCF basename (without .vcf.gz).
  - Extra args after -- are forwarded to modules/evaluate/dosage_r2.sh (e.g., --samples FILE --region chr1:1-1e6 --no-parquet --use-vcfpp).
  - Resources/logs honor QUILT2_* env if set (ACCOUNT, PARTITION, QOS, CPUS_PER_TASK, MEMORY, TIME_LIMIT).
  - Logs are written to <out_dir>/slurm/dosage_r2_%j.(out|err).
  - Requires miniforge module and conda env myenv_py310 (override MINIFORGE_MODULE/CONDA_ENV),
    plus bcftools and Rscript (data.table/arrow, vcfppR optional).
EOF
}

IMPUTED=""
TRUTH=""
OUT_PREFIX=""
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --imputed)
            IMPUTED="$2"; shift 2 ;;
        --truth)
            TRUTH="$2"; shift 2 ;;
        --out-prefix)
            OUT_PREFIX="$2"; shift 2 ;;
        --)
            shift
            EXTRA_ARGS+=("$@")
            break ;;
        *)
            EXTRA_ARGS+=("$1"); shift ;;
    esac
done

if [[ -z "${IMPUTED}" || -z "${TRUTH}" ]]; then
    usage
    exit 1
fi

if [[ ! -f "${IMPUTED}" ]]; then
    echo "[ERROR] Imputed VCF not found: ${IMPUTED}" >&2
    exit 1
fi
if [[ ! -f "${TRUTH}" ]]; then
    echo "[ERROR] Truth VCF not found: ${TRUTH}" >&2
    exit 1
fi

if [[ -z "${OUT_PREFIX}" ]]; then
    base="$(basename "${IMPUTED}")"
    dir="$(cd "$(dirname "${IMPUTED}")" && pwd)"
    case "${base}" in
        *.vcf.gz) OUT_PREFIX="${dir}/${base%.vcf.gz}" ;;
        *.vcf) OUT_PREFIX="${dir}/${base%.vcf}" ;;
        *.bcf) OUT_PREFIX="${dir}/${base%.bcf}" ;;
        *) OUT_PREFIX="${dir}/${base}" ;;
    esac
fi

# Load environment defaults if present.
if [[ -f "${ENV_FILE}" ]]; then
    # shellcheck source=/dev/null
    source "${ENV_FILE}"
elif [[ -f "${ENV_TEMPLATE}" ]]; then
    # shellcheck source=/dev/null
    source "${ENV_TEMPLATE}"
fi

if [[ ! -f "${DOSAGE_SCRIPT}" ]]; then
    echo "[ERROR] Missing helper: ${DOSAGE_SCRIPT}" >&2
    exit 1
fi

OUT_DIR="$(cd "$(dirname "${OUT_PREFIX}")" && pwd)"
LOG_DIR="${OUT_DIR}/slurm"
mkdir -p "${LOG_DIR}"

# If already running under SLURM, just execute directly (no re-submit).
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    cmd=(bash "${DOSAGE_SCRIPT}" --imputed "${IMPUTED}" --truth "${TRUTH}" --out-prefix "${OUT_PREFIX}")
    if [[ "${#EXTRA_ARGS[@]}" -gt 0 ]]; then
        cmd+=("${EXTRA_ARGS[@]}")
    fi
    echo "[INFO] Running inside SLURM job ${SLURM_JOB_ID}: ${cmd[*]}"
    exec "${cmd[@]}"
fi

# Build sbatch arguments from environment defaults if present.
sbatch_args=( --job-name=dosage_r2 --output="${LOG_DIR}/dosage_r2_%j.out" --error="${LOG_DIR}/dosage_r2_%j.err" )
[[ -n "${QUILT2_ACCOUNT:-}" ]] && sbatch_args+=( --account="${QUILT2_ACCOUNT}" )
[[ -n "${QUILT2_PARTITION:-}" ]] && sbatch_args+=( --partition="${QUILT2_PARTITION}" )
[[ -n "${QUILT2_QOS:-}" ]] && sbatch_args+=( --qos="${QUILT2_QOS}" )
[[ -n "${QUILT2_CPUS_PER_TASK:-}" ]] && sbatch_args+=( --cpus-per-task="${QUILT2_CPUS_PER_TASK}" )
[[ -n "${QUILT2_MEMORY:-}" ]] && sbatch_args+=( --mem="${QUILT2_MEMORY}" )
[[ -n "${QUILT2_TIME_LIMIT:-}" ]] && sbatch_args+=( --time="${QUILT2_TIME_LIMIT}" )

cmd=(bash "${DOSAGE_SCRIPT}" --imputed "${IMPUTED}" --truth "${TRUTH}" --out-prefix "${OUT_PREFIX}")
if [[ "${#EXTRA_ARGS[@]}" -gt 0 ]]; then
    cmd+=("${EXTRA_ARGS[@]}")
fi
cmd_str="$(printf " %q" "${cmd[@]}")"

echo "[INFO] Submitting: sbatch ${sbatch_args[*]} --wrap=\"${cmd_str}\""
sbatch "${sbatch_args[@]}" --wrap="${cmd_str}"
