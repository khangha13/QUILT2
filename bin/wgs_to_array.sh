#!/bin/bash
#SBATCH --job-name=wgs2array
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# Submit with: sbatch bin/wgs_to_array.sh /path/input.vcf.gz [output.vcf.gz]

# Lightweight SLURM wrapper around utils/wgs_to_array_vcf.py.
# Assumes input VCFs are biallelic SNPs. Only GT is recoded; other fields unchanged.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
PY_SCRIPT="${ROOT_DIR}/utils/wgs_to_array_vcf.py"
ENV_FILE="${ROOT_DIR}/config/environment.sh"
ENV_TEMPLATE="${ROOT_DIR}/config/environment.template.sh"
LIB="${ROOT_DIR}/lib/functions.sh"

usage() {
    cat <<'EOF'
Usage:
  sbatch bin/wgs_to_array.sh /path/input.vcf.gz [output.vcf.gz]

Notes:
  - Activates the QUILT2 conda env if configured (QUILT2_CONDA_ENV).
  - Output defaults to <input_basename>_array.vcf.gz in the same directory if not provided.
  - Adds a tabix index if bcftools is available.
EOF
}

INPUT="${1:-}"
OUTPUT="${2:-}"

if [[ -z "${INPUT}" ]]; then
    usage
    exit 1
fi

if [[ ! -f "${INPUT}" ]]; then
    echo "[ERROR] Input VCF not found: ${INPUT}" >&2
    exit 1
fi

# Derive default output path if not supplied.
if [[ -z "${OUTPUT}" ]]; then
    base="$(basename "${INPUT}")"
    dir="$(cd "$(dirname "${INPUT}")" && pwd)"
    case "${base}" in
        *.vcf.gz) OUTPUT="${dir}/${base%.vcf.gz}_array.vcf.gz" ;;
        *.vcf) OUTPUT="${dir}/${base%.vcf}_array.vcf.gz" ;;
        *.bcf) OUTPUT="${dir}/${base%.bcf}_array.vcf.gz" ;;
        *) OUTPUT="${dir}/${base}.array.vcf.gz" ;;
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

# Source shared helpers when available.
if [[ -f "${LIB}" ]]; then
    # shellcheck source=/dev/null
    source "${LIB}"
else
    # Minimal fallbacks
    log_info() { echo "[INFO] $*" >&2; }
    log_warn() { echo "[WARN] $*" >&2; }
    log_error() { echo "[ERROR] $*" >&2; }
    run_cmd() { "$@"; }
    load_quilt_env() { return 0; }
    ensure_bcftools() { command -v bcftools >/dev/null 2>&1; }
fi

log_info "Input VCF : ${INPUT}"
log_info "Output VCF: ${OUTPUT}"

# Try to activate conda env (from lib/functions.sh).
if ! load_quilt_env; then
    log_warn "Proceeding without activating conda env; ensure pysam is available."
fi

# Verify python and helper script.
PYTHON_BIN="${PYTHON_BIN:-python}"
if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    log_error "Python not found (tried ${PYTHON_BIN}); set PYTHON_BIN if needed."
    exit 1
fi
if [[ ! -f "${PY_SCRIPT}" ]]; then
    log_error "Helper script missing: ${PY_SCRIPT}"
    exit 1
fi

run_cmd "${PYTHON_BIN}" "${PY_SCRIPT}" -i "${INPUT}" -o "${OUTPUT}"

# Optional indexing if bcftools exists.
if command -v bcftools >/dev/null 2>&1; then
    run_cmd bcftools index -f -c "${OUTPUT}"
else
    log_warn "bcftools not found; output not indexed."
fi

log_info "Done. Output: ${OUTPUT}"
