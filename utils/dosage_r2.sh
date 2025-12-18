#!/bin/bash
# Compute concordance and dosage r2 between an imputed VCF (recoded to array)
# and a truth VCF. The imputed VCF is forced through utils/wgs_to_array_vcf.py
# to align AT/CG coding with the array format; truth is assumed already on the
# desired scale. Dosages are derived from GT; DS/GP tags are no longer used.
# Outputs include per-variant metrics plus a per-site/per-sample concordance
# matrix in Parquet (unless disabled).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
LIB_PATH="${ROOT_DIR}/lib/functions.sh"

if [[ ! -f "${LIB_PATH}" ]]; then
    echo "[ERROR] Missing helper library: ${LIB_PATH}" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "${LIB_PATH}"

usage() {
    cat <<'EOF'
Usage: dosage_r2.sh --imputed <imputed.vcf.gz> --truth <truth.vcf.gz> --out-prefix <path> [options]

Required:
  --imputed PATH        Imputed VCF/BCF (array-coded or will be recoded)
  --truth PATH          Truth VCF/BCF (GT used as reference)
  --out-prefix PREFIX   Output prefix for metrics/plots (e.g., results/dosage_eval)

Options:
  --samples FILE        File with sample IDs to evaluate (one per line). Defaults to intersection of VCF samples.
  --region STR          Region string (e.g., chr1:1-1e6) to limit evaluation.
  --use-vcfpp           Additionally run vcfppR::vcfcomp/vcfplot on harmonized VCFs (requires vcfppR).
  --no-parquet          Skip writing per-site/per-sample concordance Parquet (default: write it).
  --no-biallelic-only   Do not restrict to biallelic SNPs (default is to restrict).
  --no-plots            Skip plotting; only metrics TSVs are produced.
  --keep-temp           Do not delete the temporary working directory.
  --help                Show this message and exit.

Outputs (PREFIX.*):
  metrics.tsv           Per-variant metrics (r2, concordance, maf, counts)
  summary.tsv           Overall summary metrics
  maf_bins.tsv          Metrics summarized by MAF bins
  concordance.parquet   Wide per-site/per-sample concordance (0/1/NA)
  r2_hist.png           Distribution of r2 (if plots enabled)
  concordance_hist.png  Distribution of concordance (if plots enabled)
  r2_vs_maf.png         Scatter/hex plot of r2 vs MAF (if plots enabled)
EOF
}

IMPUTED_VCF=""
TRUTH_VCF=""
OUT_PREFIX=""
SAMPLE_FILE=""
REGION=""
BIALLELIC_ONLY=true
RUN_PLOTS=true
KEEP_TEMP=false
USE_VCFPP=false
WRITE_PARQUET=true

while [[ $# -gt 0 ]]; do
    case "$1" in
        --imputed)
            IMPUTED_VCF="$2"; shift 2 ;;
        --truth)
            TRUTH_VCF="$2"; shift 2 ;;
        --out-prefix)
            OUT_PREFIX="$2"; shift 2 ;;
        --samples)
            SAMPLE_FILE="$2"; shift 2 ;;
        --region)
            REGION="$2"; shift 2 ;;
        --use-vcfpp)
            USE_VCFPP=true; shift ;;
        --no-parquet)
            WRITE_PARQUET=false; shift ;;
        --no-biallelic-only)
            BIALLELIC_ONLY=false; shift ;;
        --no-plots)
            RUN_PLOTS=false; shift ;;
        --keep-temp)
            KEEP_TEMP=true; shift ;;
        --help|-h)
            usage; exit 0 ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ -z "${IMPUTED_VCF}" || -z "${TRUTH_VCF}" ]]; then
    usage
    exit 1
fi

if [[ -z "${OUT_PREFIX}" ]]; then
    base="$(basename "${IMPUTED_VCF}")"
    dir="$(cd "$(dirname "${IMPUTED_VCF}")" && pwd)"
    case "${base}" in
        *.vcf.gz) OUT_PREFIX="${dir}/${base%.vcf.gz}" ;;
        *.vcf) OUT_PREFIX="${dir}/${base%.vcf}" ;;
        *.bcf) OUT_PREFIX="${dir}/${base%.bcf}" ;;
        *) OUT_PREFIX="${dir}/${base}" ;;
    esac
    log_info "Defaulting --out-prefix to ${OUT_PREFIX}"
fi

if [[ ! -f "${IMPUTED_VCF}" ]]; then
    log_error "Imputed VCF not found: ${IMPUTED_VCF}"
    exit 1
fi
if [[ ! -f "${TRUTH_VCF}" ]]; then
    log_error "Truth VCF not found: ${TRUTH_VCF}"
    exit 1
fi

OUT_DIR="$(cd "$(dirname "${OUT_PREFIX}")" && pwd)"
OUT_BASENAME="$(basename "${OUT_PREFIX}")"
OUT_PREFIX="${OUT_DIR}/${OUT_BASENAME}"
mkdir -p "${OUT_DIR}"

TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/dosage_r2.XXXXXX")"
cleanup() {
    if [[ "${KEEP_TEMP}" == "true" ]]; then
        log_info "Temporary files kept at: ${TMP_DIR}"
    else
        rm -rf "${TMP_DIR}"
    fi
}
trap cleanup EXIT

log_info "Temporary working directory: ${TMP_DIR}"

# Activate 'rplot' conda environment (similar to step1d)
CONDA_ENV="${CONDA_ENV:-rplot}"

if command -v module >/dev/null 2>&1; then
    if module load miniforge/25.3.0-3 >/dev/null 2>&1; then
        log_info "Loaded miniforge/25.3.0-3 module"
    else
        log_warn "Failed to load miniforge/25.3.0-3 module"
    fi
    if module load bcftools >/dev/null 2>&1; then
        log_info "Loaded bcftools module"
    else
        log_warn "Failed to load bcftools module"
    fi
fi

conda_activated=false
if [[ -n "${ROOTMINIFORGE:-}" && -f "${ROOTMINIFORGE}/etc/profile.d/conda.sh" ]]; then
    # shellcheck source=/dev/null
    source "${ROOTMINIFORGE}/etc/profile.d/conda.sh"
    if conda activate "${CONDA_ENV}" >/dev/null 2>&1; then
        log_info "Conda environment '${CONDA_ENV}' activated"
        conda_activated=true
    else
        log_warn "Failed to activate conda environment: ${CONDA_ENV}"
    fi
elif command -v conda >/dev/null 2>&1; then
    conda_base="$(conda info --base 2>/dev/null || true)"
    if [[ -n "${conda_base}" && -f "${conda_base}/etc/profile.d/conda.sh" ]]; then
        # shellcheck source=/dev/null
        source "${conda_base}/etc/profile.d/conda.sh"
        if conda activate "${CONDA_ENV}" >/dev/null 2>&1; then
            log_info "Conda environment '${CONDA_ENV}' activated"
            conda_activated=true
        else
            log_warn "Failed to activate conda environment: ${CONDA_ENV}"
        fi
    fi
fi

if [[ "${conda_activated}" != "true" ]]; then
    log_warn "Continuing without activating conda env; ensure bcftools, python, and R are available."
fi

ensure_bcftools || exit 1
require_cmd Rscript || exit 1

maybe_index() {
    local vcf="$1"
    if [[ "${vcf}" =~ \.vcf\.gz$ || "${vcf}" =~ \.bcf$ ]]; then
        run_cmd bcftools index -f -c "${vcf}"
    fi
}

maybe_index "${IMPUTED_VCF}"
maybe_index "${TRUTH_VCF}"

# Recode imputed VCF to array-style genotypes to align with array output.
PYTHON_BIN="${PYTHON_BIN:-python}"
ARRAY_SCRIPT="${ROOT_DIR}/utils/wgs_to_array_vcf.py"
if [[ ! -f "${ARRAY_SCRIPT}" ]]; then
    log_error "Recoder not found: ${ARRAY_SCRIPT}"
    exit 1
fi
if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    log_error "Python not found (set PYTHON_BIN if needed)"
    exit 1
fi

IMPUTED_ARRAY="${TMP_DIR}/imputed.array.vcf.gz"

log_info "Recoding imputed VCF to array scale"
run_cmd "${PYTHON_BIN}" "${ARRAY_SCRIPT}" -i "${IMPUTED_VCF}" -o "${IMPUTED_ARRAY}"

run_cmd bcftools index -f -c "${IMPUTED_ARRAY}"
maybe_index "${TRUTH_VCF}"

# Derive sample set
SAMPLE_SET="${SAMPLE_FILE}"
if [[ -n "${SAMPLE_FILE}" ]]; then
    if [[ ! -f "${SAMPLE_FILE}" ]]; then
        log_error "Sample file not found: ${SAMPLE_FILE}"
        exit 1
    fi
else
    log_info "Deriving common sample set from both VCFs"
    bcftools query -l "${IMPUTED_VCF}" | sort > "${TMP_DIR}/imputed.samples"
    bcftools query -l "${TRUTH_VCF}" | sort > "${TMP_DIR}/truth.samples"
    comm -12 "${TMP_DIR}/imputed.samples" "${TMP_DIR}/truth.samples" > "${TMP_DIR}/common.samples"
    if [[ ! -s "${TMP_DIR}/common.samples" ]]; then
        log_error "No overlapping samples found between VCFs"
        exit 1
    fi
    SAMPLE_SET="${TMP_DIR}/common.samples"
fi

# Intersect sites (allele-aware) and restrict to region if provided
log_info "Intersecting sites between imputed and truth VCFs"
ISEC_DIR="${TMP_DIR}/isec"
mkdir -p "${ISEC_DIR}"
run_cmd bcftools isec -c all -p "${ISEC_DIR}" -Oz "${IMPUTED_ARRAY}" "${TRUTH_VCF}"
SITE_VCF="${ISEC_DIR}/0002.vcf.gz"
if [[ ! -s "${SITE_VCF}" ]]; then
    log_error "Intersection produced no common sites"
    exit 1
fi
SITE_LIST="${TMP_DIR}/sites.txt"
bcftools query -f '%CHROM\t%POS\t%POS\n' "${SITE_VCF}" > "${SITE_LIST}"

REGION_ARGS=()
if [[ -n "${REGION}" ]]; then
    REGION_ARGS=( -r "${REGION}" )
fi

BIALLELIC_ARGS=()
if [[ "${BIALLELIC_ONLY}" == "true" ]]; then
    BIALLELIC_ARGS=( -m2 -M2 -v snps )
fi

log_info "Extracting harmonized VCFs"
IMPUTED_COMMON="${TMP_DIR}/imputed.common.vcf.gz"
TRUTH_COMMON="${TMP_DIR}/truth.common.vcf.gz"

run_cmd bcftools view -R "${SITE_LIST}" -S "${SAMPLE_SET}" "${REGION_ARGS[@]}" "${BIALLELIC_ARGS[@]}" -Oz -o "${IMPUTED_COMMON}" "${IMPUTED_ARRAY}"
run_cmd bcftools view -R "${SITE_LIST}" -S "${SAMPLE_SET}" "${REGION_ARGS[@]}" "${BIALLELIC_ARGS[@]}" -Oz -o "${TRUTH_COMMON}" "${TRUTH_VCF}"
run_cmd bcftools index -f -c "${IMPUTED_COMMON}"
run_cmd bcftools index -f -c "${TRUTH_COMMON}"

# Build header
readarray -t SAMPLE_ORDER < "${SAMPLE_SET}"
HEADER="CHROM\tPOS\tREF\tALT\tID"
for s in "${SAMPLE_ORDER[@]}"; do
    HEADER="${HEADER}\t${s}"
done

IMPUTED_GT_TSV="${TMP_DIR}/imputed_gt.tsv"
TRUTH_GT_TSV="${TMP_DIR}/truth_gt.tsv"

log_info "Extracting imputed genotypes (array-coded GT)"
echo -e "${HEADER}" > "${IMPUTED_GT_TSV}"
if ! bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID[\t%GT]\n" -S "${SAMPLE_SET}" "${IMPUTED_COMMON}" >> "${IMPUTED_GT_TSV}"; then
    log_error "Failed to extract GT from imputed VCF"
    exit 1
fi

log_info "Extracting truth genotypes (GT)"
echo -e "${HEADER}" > "${TRUTH_GT_TSV}"
if ! bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID[\t%GT]\n" -S "${SAMPLE_SET}" "${TRUTH_COMMON}" >> "${TRUTH_GT_TSV}"; then
    log_error "Failed to extract GT from truth VCF"
    exit 1
fi

R_HELPER="${SCRIPT_DIR}/dosage_r2.R"
if [[ ! -f "${R_HELPER}" ]]; then
    log_error "Missing R helper script: ${R_HELPER}"
    exit 1
fi

R_ARGS=(
    "${R_HELPER}"
    "--imputed-gt" "${IMPUTED_GT_TSV}"
    "--truth-gt" "${TRUTH_GT_TSV}"
    "--samples" "${SAMPLE_SET}"
    "--out-prefix" "${OUT_PREFIX}"
    "--write-parquet" "${WRITE_PARQUET}"
)

if [[ "${USE_VCFPP}" == "true" ]]; then
    R_ARGS+=( "--use-vcfpp" "--vcfpp-imputed" "${IMPUTED_COMMON}" "--vcfpp-truth" "${TRUTH_COMMON}" )
fi

if [[ "${RUN_PLOTS}" == "true" ]]; then
    R_ARGS+=( "--plots" )
fi

log_info "Running R metrics helper"
run_cmd Rscript "${R_ARGS[@]}"

log_info "Dosage r2 and concordance metrics written to ${OUT_PREFIX}.*"
