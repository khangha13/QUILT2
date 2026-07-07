#!/bin/bash
# Submit helper: wraps modules/evaluate/dosage_r2.sh in sbatch with logs placed beside the outputs.
# Usage:
#   bash bin/dosage_r2_sbatch.sh --imputed /path/imputed.vcf.gz --truth /path/truth.vcf.gz [--out-prefix prefix] [-- extra flags]
# If run under SLURM (SLURM_JOB_ID set), it will execute dosage_r2.sh directly (no re-submit).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
DOSAGE_SCRIPT="${ROOT_DIR}/modules/evaluate/dosage_r2.sh"
CONCAT_SCRIPT="${ROOT_DIR}/modules/evaluate/concat_imputed.sh"
ENV_FILE="${ROOT_DIR}/config/environment.sh"
ENV_TEMPLATE="${ROOT_DIR}/config/environment.template.sh"

usage() {
    cat <<'EOF'
Usage:
  bash bin/dosage_r2_sbatch.sh --imputed IMPT_VCF --truth TRUTH_VCF [--out-prefix OUT_PREFIX] [-- extra args]
  bash bin/dosage_r2_sbatch.sh --chunks-dir CHUNKS_DIR --truth TRUTH_VCF [--chr LIST] [--concat-force] [--out-prefix OUT_PREFIX] [-- extra args]

Input (exactly one of):
  --imputed PATH        Already-concatenated imputed VCF/BCF
  --chunks-dir PATH     OUTPUT_DIR/chunks/imputed layout of per-chunk VCFs; concatenated via
                        modules/evaluate/concat_imputed.sh before evaluation

Options (--chunks-dir mode only):
  --chr LIST            Comma/space list of chromosomes to include (default: every chromosome
                        found via the run_quilt2.sh chunk manifest, or every --chunks-dir subdir)
  --concat-force        Re-concatenate even if concat_imputed.sh outputs already exist

Notes:
  - OUT_PREFIX defaults to the imputed VCF basename (without .vcf.gz) in --imputed mode, or
    OUTPUT_DIR/eval/dosage_eval in --chunks-dir mode (falls back to <chunks-dir>/dosage_eval
    if OUTPUT_DIR can't be resolved from a run_manifest.tsv).
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
CHUNKS_DIR=""
CHR_ARG=""
CONCAT_FORCE=false
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --imputed)
            IMPUTED="$2"; shift 2 ;;
        --truth)
            TRUTH="$2"; shift 2 ;;
        --out-prefix)
            OUT_PREFIX="$2"; shift 2 ;;
        --chunks-dir)
            CHUNKS_DIR="$2"; shift 2 ;;
        --chr)
            CHR_ARG="$2"; shift 2 ;;
        --concat-force)
            CONCAT_FORCE=true; shift ;;
        --)
            shift
            EXTRA_ARGS+=("$@")
            break ;;
        *)
            EXTRA_ARGS+=("$1"); shift ;;
    esac
done

if [[ -n "${IMPUTED}" && -n "${CHUNKS_DIR}" ]]; then
    echo "[ERROR] --imputed and --chunks-dir are mutually exclusive." >&2
    exit 1
fi
if [[ -z "${IMPUTED}" && -z "${CHUNKS_DIR}" ]]; then
    usage
    exit 1
fi
if [[ -z "${TRUTH}" ]]; then
    usage
    exit 1
fi

if [[ -n "${IMPUTED}" && ! -f "${IMPUTED}" ]]; then
    echo "[ERROR] Imputed VCF not found: ${IMPUTED}" >&2
    exit 1
fi
if [[ -n "${CHUNKS_DIR}" && ! -d "${CHUNKS_DIR}" ]]; then
    echo "[ERROR] Chunks directory not found: ${CHUNKS_DIR}" >&2
    exit 1
fi
if [[ ! -f "${TRUTH}" ]]; then
    echo "[ERROR] Truth VCF not found: ${TRUTH}" >&2
    exit 1
fi

if [[ -n "${CHUNKS_DIR}" ]]; then
    CHUNKS_DIR="$(cd "${CHUNKS_DIR}" && pwd)"
fi

if [[ -z "${OUT_PREFIX}" ]]; then
    if [[ -n "${CHUNKS_DIR}" ]]; then
        OUTPUT_DIR_GUESS=""
        if candidate_dir="$(cd "${CHUNKS_DIR}/../.." 2>/dev/null && pwd)"; then
            [[ -f "${candidate_dir}/run_manifest.tsv" ]] && OUTPUT_DIR_GUESS="${candidate_dir}"
        fi
        if [[ -n "${OUTPUT_DIR_GUESS}" ]]; then
            OUT_PREFIX="${OUTPUT_DIR_GUESS}/eval/dosage_eval"
        else
            OUT_PREFIX="${CHUNKS_DIR}/dosage_eval"
        fi
    else
        base="$(basename "${IMPUTED}")"
        dir="$(cd "$(dirname "${IMPUTED}")" && pwd)"
        case "${base}" in
            *.vcf.gz) OUT_PREFIX="${dir}/${base%.vcf.gz}" ;;
            *.vcf) OUT_PREFIX="${dir}/${base%.vcf}" ;;
            *.bcf) OUT_PREFIX="${dir}/${base%.bcf}" ;;
            *) OUT_PREFIX="${dir}/${base}" ;;
        esac
    fi
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
if [[ -n "${CHUNKS_DIR}" && ! -f "${CONCAT_SCRIPT}" ]]; then
    echo "[ERROR] Missing helper: ${CONCAT_SCRIPT}" >&2
    exit 1
fi

OUT_DIR="$(cd "$(dirname "${OUT_PREFIX}")" && pwd)"
LOG_DIR="${OUT_DIR}/slurm"
mkdir -p "${LOG_DIR}"

# If already running under SLURM, just execute directly (no re-submit).
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    if [[ -n "${CHUNKS_DIR}" ]]; then
        concat_cmd=(bash "${CONCAT_SCRIPT}" --chunks-dir "${CHUNKS_DIR}")
        [[ -n "${CHR_ARG}" ]] && concat_cmd+=(--chr "${CHR_ARG}")
        [[ "${CONCAT_FORCE}" == "true" ]] && concat_cmd+=(--force)
        echo "[INFO] Running inside SLURM job ${SLURM_JOB_ID}: ${concat_cmd[*]}"
        IMPUTED="$("${concat_cmd[@]}")"
        echo "[INFO] Concatenated imputed VCF: ${IMPUTED}"
    fi
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

if [[ -n "${CHUNKS_DIR}" ]]; then
    # Two chained steps (concat, then evaluate) don't fit cleanly through --wrap's
    # single command string, so generate a small script and submit that instead.
    concat_args=(--chunks-dir "${CHUNKS_DIR}")
    [[ -n "${CHR_ARG}" ]] && concat_args+=(--chr "${CHR_ARG}")
    [[ "${CONCAT_FORCE}" == "true" ]] && concat_args+=(--force)
    concat_args_quoted="$(printf " %q" "${concat_args[@]}")"

    dosage_args=(--truth "${TRUTH}" --out-prefix "${OUT_PREFIX}")
    if [[ "${#EXTRA_ARGS[@]}" -gt 0 ]]; then
        dosage_args+=("${EXTRA_ARGS[@]}")
    fi
    dosage_args_quoted="$(printf " %q" "${dosage_args[@]}")"

    CHAIN_SCRIPT="${LOG_DIR}/dosage_r2_concat_chain_$(date +%Y%m%d_%H%M%S).sh"
    cat <<EOF > "${CHAIN_SCRIPT}"
#!/bin/bash
set -euo pipefail
IMPUTED_PATH="\$(bash "${CONCAT_SCRIPT}"${concat_args_quoted})"
echo "[INFO] Concatenated imputed VCF: \${IMPUTED_PATH}" >&2
exec bash "${DOSAGE_SCRIPT}" --imputed "\${IMPUTED_PATH}"${dosage_args_quoted}
EOF
    chmod +x "${CHAIN_SCRIPT}"

    echo "[INFO] Submitting: sbatch ${sbatch_args[*]} ${CHAIN_SCRIPT} (concat + dosage_r2 chain)"
    sbatch "${sbatch_args[@]}" "${CHAIN_SCRIPT}"
else
    cmd=(bash "${DOSAGE_SCRIPT}" --imputed "${IMPUTED}" --truth "${TRUTH}" --out-prefix "${OUT_PREFIX}")
    if [[ "${#EXTRA_ARGS[@]}" -gt 0 ]]; then
        cmd+=("${EXTRA_ARGS[@]}")
    fi
    cmd_str="$(printf " %q" "${cmd[@]}")"

    echo "[INFO] Submitting: sbatch ${sbatch_args[*]} --wrap=\"${cmd_str}\""
    sbatch "${sbatch_args[@]}" --wrap="${cmd_str}"
fi
