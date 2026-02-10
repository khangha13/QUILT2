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

Environment:
  - Loads module "${MINIFORGE_MODULE:-miniforge}" and activates conda env "${CONDA_ENV:-myenv_py310}".
  - Override with MINIFORGE_MODULE and CONDA_ENV environment variables.

Contig Name Handling:
  The canonical contig style is "ChrNN" (e.g., Chr01, Chr02, Chr17).
  The script automatically detects and normalizes contig naming conventions.
  Supported styles: "Chr01" (ChrNN), "chr1" (chrN), "1" (N).
  Any VCF using a non-canonical style is automatically renamed to ChrNN format.

Position-Only Intersection:
  The script uses position-only matching (CHROM + POS) to find common variants.
  REF/ALT columns are allowed to differ between imputed and truth VCFs because the
  AT/CG recoding (wgs_to_array_vcf.py) makes dosages comparable regardless of which
  physical alleles appear as REF/ALT (array TOP-strand vs reference-genome alleles).

Outputs (PREFIX.*):
  metrics.tsv                          Per-variant metrics (r2, concordance, maf, counts)
  summary.tsv                          Overall summary metrics
  maf_bins.tsv                         Metrics summarized by MAF bins
  concordance.parquet                  Wide per-site/per-sample concordance (0/1/NA)
  imputed.array_coded.common.vcf.gz   Translated imputed VCF at common positions
  truth.common.vcf.gz                 Truth VCF at common positions
  imputed.AB_format.tsv               Imputed genotypes in A/B format (A/A, A/B, B/B)
  truth.AB_format.tsv                 Truth genotypes in A/B format (A/A, A/B, B/B)
  r2_hist.png                          Distribution of r2 (if plots enabled)
  concordance_hist.png                 Distribution of concordance (if plots enabled)
  r2_vs_maf.png                        Scatter/hex plot of r2 vs MAF (if plots enabled)
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

# Activate conda environment (explicit, no fallback).
CONDA_ENV="${CONDA_ENV:-myenv_py310}"
MINIFORGE_MODULE="${MINIFORGE_MODULE:-miniforge}"
BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"

if command -v module >/dev/null 2>&1; then
    if module load "${MINIFORGE_MODULE}" >/dev/null 2>&1; then
        log_info "Loaded ${MINIFORGE_MODULE} module"
    else
        log_error "Failed to load module: ${MINIFORGE_MODULE}"
        exit 1
    fi
else
    log_error "module command not found; cannot load ${MINIFORGE_MODULE}"
    exit 1
fi

if [[ -z "${ROOTMINIFORGE:-}" ]]; then
    log_error "ROOTMINIFORGE is not set; cannot source conda.sh"
    exit 1
fi
if [[ ! -f "${ROOTMINIFORGE}/etc/profile.d/conda.sh" ]]; then
    log_error "conda.sh not found at ${ROOTMINIFORGE}/etc/profile.d/conda.sh"
    exit 1
fi

# shellcheck source=/dev/null
source "${ROOTMINIFORGE}/etc/profile.d/conda.sh"
if conda activate "${CONDA_ENV}" >/dev/null 2>&1; then
    log_info "Conda environment '${CONDA_ENV}' activated"
else
    log_error "Failed to activate conda environment: ${CONDA_ENV}"
    exit 1
fi

ensure_rscript() {
    if command -v Rscript >/dev/null 2>&1; then
        return 0
    fi

    # Auto-fix for missing Rscript: install R + required packages into the active conda env.
    # This matches the recommended manual fix:
    #   conda activate myenv_py310
    #   conda install -y r-base r-data.table r-ggplot2 r-arrow
    if ! command -v conda >/dev/null 2>&1; then
        log_error "Rscript not found, and conda is not available to install it."
        return 1
    fi

    log_warn "Rscript not found; installing R into conda env '${CONDA_ENV}'"
    if ! conda install -y r-base r-data.table r-ggplot2 r-arrow; then
        log_error "Failed to install R packages into conda env '${CONDA_ENV}'."
        return 1
    fi

    if ! command -v Rscript >/dev/null 2>&1; then
        log_error "Rscript still not found after conda install in env '${CONDA_ENV}'."
        return 1
    fi
}

ensure_bcftools || exit 1
ensure_rscript || exit 1
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
PYTHON_BIN="${PYTHON_BIN:-python3}"
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

# =============================================================================
# CONTIG NAME DETECTION AND NORMALIZATION
# =============================================================================
# Detect contig naming conventions in each VCF and create rename maps if needed.
# Common patterns: "Chr01", "chr1", "1", "NC_..." etc.

detect_contig_style() {
    local vcf="$1"
    local first_contig
    first_contig="$(bcftools view -H "${vcf}" 2>/dev/null | head -1 | cut -f1)"
    if [[ -z "${first_contig}" ]]; then
        # Fallback to header contigs
        first_contig="$(bcftools view -h "${vcf}" 2>/dev/null | grep '^##contig=<ID=' | head -1 | sed 's/.*ID=\([^,>]*\).*/\1/')"
    fi
    
    if [[ "${first_contig}" =~ ^Chr[0-9]+ ]]; then
        echo "ChrNN"  # e.g., Chr01, Chr02
    elif [[ "${first_contig}" =~ ^chr[0-9]+ ]]; then
        echo "chrN"   # e.g., chr1, chr2
    elif [[ "${first_contig}" =~ ^[0-9]+$ ]]; then
        echo "N"      # e.g., 1, 2
    else
        echo "other"  # Unknown pattern
    fi
}

# Build a contig rename map file: OLD_NAME<TAB>NEW_NAME
# Converts source style to target style for chromosomes 1-17 (apple) plus common others
build_contig_rename_map() {
    local source_style="$1"
    local target_style="$2"
    local map_file="$3"
    
    : > "${map_file}"
    
    for i in $(seq 1 22); do  # Cover up to 22 for generality (human) + apple 1-17
        local src_name tgt_name
        case "${source_style}" in
            ChrNN) src_name="$(printf 'Chr%02d' "${i}")" ;;
            chrN)  src_name="chr${i}" ;;
            N)     src_name="${i}" ;;
            *)     src_name="${i}" ;;
        esac
        case "${target_style}" in
            ChrNN) tgt_name="$(printf 'Chr%02d' "${i}")" ;;
            chrN)  tgt_name="chr${i}" ;;
            N)     tgt_name="${i}" ;;
            *)     tgt_name="${i}" ;;
        esac
        if [[ "${src_name}" != "${tgt_name}" ]]; then
            printf '%s\t%s\n' "${src_name}" "${tgt_name}" >> "${map_file}"
        fi
    done
    
    # Also handle X, Y, M/MT if present
    for special in X Y M MT; do
        local src_name tgt_name
        case "${source_style}" in
            ChrNN) src_name="Chr${special}" ;;
            chrN)  src_name="chr${special}" ;;
            N)     src_name="${special}" ;;
            *)     src_name="${special}" ;;
        esac
        case "${target_style}" in
            ChrNN) tgt_name="Chr${special}" ;;
            chrN)  tgt_name="chr${special}" ;;
            N)     tgt_name="${special}" ;;
            *)     tgt_name="${special}" ;;
        esac
        if [[ "${src_name}" != "${tgt_name}" ]]; then
            printf '%s\t%s\n' "${src_name}" "${tgt_name}" >> "${map_file}"
        fi
    done
}

# Normalize contig names in a VCF to match target style
normalize_contigs() {
    local input_vcf="$1"
    local output_vcf="$2"
    local rename_map="$3"
    
    if [[ ! -s "${rename_map}" ]]; then
        # No renaming needed; just copy/link
        cp "${input_vcf}" "${output_vcf}"
    else
        log_info "Renaming contigs using map: ${rename_map}"
        run_cmd bcftools annotate --rename-chrs "${rename_map}" "${input_vcf}" -Oz -o "${output_vcf}"
    fi
    run_cmd bcftools index -f -c "${output_vcf}"
}

log_info "Detecting contig naming conventions"
IMPUTED_CONTIG_STYLE="$(detect_contig_style "${IMPUTED_ARRAY}")"
TRUTH_CONTIG_STYLE="$(detect_contig_style "${TRUTH_VCF}")"
log_info "Imputed VCF contig style: ${IMPUTED_CONTIG_STYLE}"
log_info "Truth VCF contig style: ${TRUTH_CONTIG_STYLE}"

# Canonical contig style is ChrNN (e.g., Chr01, Chr02). Normalize any deviant VCFs.
CANONICAL_STYLE="ChrNN"

IMPUTED_NORMALIZED="${IMPUTED_ARRAY}"
if [[ "${IMPUTED_CONTIG_STYLE}" != "${CANONICAL_STYLE}" ]]; then
    log_warn "Imputed VCF uses '${IMPUTED_CONTIG_STYLE}' contigs. Normalizing to canonical '${CANONICAL_STYLE}'."
    CONTIG_RENAME_MAP="${TMP_DIR}/imputed_contig_rename.map"
    build_contig_rename_map "${IMPUTED_CONTIG_STYLE}" "${CANONICAL_STYLE}" "${CONTIG_RENAME_MAP}"
    if [[ -s "${CONTIG_RENAME_MAP}" ]]; then
        log_info "Imputed contig rename map (first 10 entries):"
        head -10 "${CONTIG_RENAME_MAP}" >&2 || true
        IMPUTED_NORMALIZED="${TMP_DIR}/imputed.array.renamed.vcf.gz"
        normalize_contigs "${IMPUTED_ARRAY}" "${IMPUTED_NORMALIZED}" "${CONTIG_RENAME_MAP}"
    else
        log_warn "No rename mappings generated for imputed; proceeding with original contig names."
    fi
fi

TRUTH_NORMALIZED="${TRUTH_VCF}"
if [[ "${TRUTH_CONTIG_STYLE}" != "${CANONICAL_STYLE}" ]]; then
    log_warn "Truth VCF uses '${TRUTH_CONTIG_STYLE}' contigs. Normalizing to canonical '${CANONICAL_STYLE}'."
    CONTIG_RENAME_MAP="${TMP_DIR}/truth_contig_rename.map"
    build_contig_rename_map "${TRUTH_CONTIG_STYLE}" "${CANONICAL_STYLE}" "${CONTIG_RENAME_MAP}"
    if [[ -s "${CONTIG_RENAME_MAP}" ]]; then
        log_info "Truth contig rename map (first 10 entries):"
        head -10 "${CONTIG_RENAME_MAP}" >&2 || true
        TRUTH_NORMALIZED="${TMP_DIR}/truth.renamed.vcf.gz"
        normalize_contigs "${TRUTH_VCF}" "${TRUTH_NORMALIZED}" "${CONTIG_RENAME_MAP}"
    else
        log_warn "No rename mappings generated for truth; proceeding with original contig names."
    fi
fi

# =============================================================================
# POSITION-ONLY INTERSECTION (ignoring REF/ALT)
# =============================================================================
# This approach finds overlapping positions first, then filters each VCF to those
# positions. REF/ALT may still differ at shared positions.
#
# NOTE: Potential future optimization - `bcftools isec -n=2 -c all` can perform
# position-only matching directly and may be faster for large files. Example:
#   bcftools isec -n=2 -c all imputed.vcf.gz truth.vcf.gz -p isec_output/
# The current manual approach (bcftools query + comm) is kept for explicit control
# and detailed logging of the intersection process.

log_info "Extracting positions from each VCF (position-only intersection)"
IMPUTED_POS="${TMP_DIR}/imputed.pos.txt"
TRUTH_POS="${TMP_DIR}/truth.pos.txt"

bcftools query -f '%CHROM\t%POS\n' "${IMPUTED_NORMALIZED}" | sort -u > "${IMPUTED_POS}"
bcftools query -f '%CHROM\t%POS\n' "${TRUTH_NORMALIZED}" | sort -u > "${TRUTH_POS}"

log_info "Finding common positions"
COMMON_POS="${TMP_DIR}/common.pos.txt"
comm -12 "${IMPUTED_POS}" "${TRUTH_POS}" > "${COMMON_POS}"

COMMON_POS_COUNT="$(wc -l < "${COMMON_POS}" | tr -d ' ')"
IMPUTED_POS_COUNT="$(wc -l < "${IMPUTED_POS}" | tr -d ' ')"
TRUTH_POS_COUNT="$(wc -l < "${TRUTH_POS}" | tr -d ' ')"

log_info "Position counts: imputed=${IMPUTED_POS_COUNT}, truth=${TRUTH_POS_COUNT}, common=${COMMON_POS_COUNT}"

if [[ "${COMMON_POS_COUNT}" -eq 0 ]]; then
    log_error "No overlapping positions found between VCFs."
    log_error "This may indicate completely different contig names or genomic regions."
    log_error "Sample imputed positions:"
    head -5 "${IMPUTED_POS}" >&2 || true
    log_error "Sample truth positions:"
    head -5 "${TRUTH_POS}" >&2 || true
    exit 1
fi

# Convert to targets format (CHROM<TAB>POS) for bcftools -T
# Using -T instead of -R avoids pulling extra records at multi-allelic positions
SITE_LIST="${TMP_DIR}/sites.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2}' "${COMMON_POS}" > "${SITE_LIST}"

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

REGION_ARGS=()
if [[ -n "${REGION}" ]]; then
    REGION_ARGS=( -r "${REGION}" )
fi

BIALLELIC_ARGS=()
if [[ "${BIALLELIC_ONLY}" == "true" ]]; then
    BIALLELIC_ARGS=( -m2 -M2 -v snps )
fi

# =============================================================================
# FILTER TO COMMON POSITIONS (position-only matching)
# =============================================================================
# The AT/CG recoding (wgs_to_array_vcf.py) makes dosages comparable even when
# REF/ALT columns differ between imputed and truth VCFs. The array truth VCF
# uses TOP-strand allele coding (e.g., REF=A ALT=G) while the imputed VCF
# preserves reference-genome alleles (e.g., REF=T ALT=C). Since the AT/CG
# grouping is strand-invariant (A<->T both map to group 0, C<->G both map to
# group 1), the GT dosage (0/0->0, 0/1->1, 1/1->2) is identical regardless
# of which physical alleles are listed as REF/ALT.
#
# Therefore we match on (CHROM, POS) only and skip REF/ALT filtering.

log_info "Filtering VCFs to ${COMMON_POS_COUNT} common positions"
IMPUTED_COMMON="${TMP_DIR}/imputed.common.vcf.gz"
TRUTH_COMMON="${TMP_DIR}/truth.common.vcf.gz"

# Use -T (targets) for exact position matching; avoids extra records from -R region overlaps
run_cmd bcftools view -T "${SITE_LIST}" -S "${SAMPLE_SET}" "${REGION_ARGS[@]}" "${BIALLELIC_ARGS[@]}" -Oz -o "${IMPUTED_COMMON}" "${IMPUTED_NORMALIZED}"
run_cmd bcftools view -T "${SITE_LIST}" -S "${SAMPLE_SET}" "${REGION_ARGS[@]}" "${BIALLELIC_ARGS[@]}" -Oz -o "${TRUTH_COMMON}" "${TRUTH_NORMALIZED}"
run_cmd bcftools index -f -c "${IMPUTED_COMMON}"
run_cmd bcftools index -f -c "${TRUTH_COMMON}"

# Sanity check: both VCFs should have variants and positions should overlap
IMPUTED_COMMON_COUNT="$(bcftools view -H "${IMPUTED_COMMON}" | wc -l | tr -d ' ')"
TRUTH_COMMON_COUNT="$(bcftools view -H "${TRUTH_COMMON}" | wc -l | tr -d ' ')"
log_info "Variants after filtering: imputed=${IMPUTED_COMMON_COUNT}, truth=${TRUTH_COMMON_COUNT}"

if [[ "${IMPUTED_COMMON_COUNT}" -eq 0 || "${TRUTH_COMMON_COUNT}" -eq 0 ]]; then
    log_error "No variants remaining after filtering. Check VCF content and contig naming."
    exit 1
fi

# Report position-level overlap (informational)
bcftools query -f '%CHROM\t%POS\n' "${IMPUTED_COMMON}" | sort -u > "${TMP_DIR}/imputed.common.pos"
bcftools query -f '%CHROM\t%POS\n' "${TRUTH_COMMON}" | sort -u > "${TMP_DIR}/truth.common.pos"
SHARED_POS="$(comm -12 "${TMP_DIR}/imputed.common.pos" "${TMP_DIR}/truth.common.pos" | wc -l | tr -d ' ')"
log_info "Shared positions in filtered VCFs: ${SHARED_POS}"

# =============================================================================
# SAVE TRANSLATED VCFs FOR FUTURE CONVENIENCE
# =============================================================================
# Copy the position-matched VCFs to the output prefix for reuse / inspection.

IMPUTED_OUT_VCF="${OUT_PREFIX}.imputed.array_coded.common.vcf.gz"
TRUTH_OUT_VCF="${OUT_PREFIX}.truth.common.vcf.gz"

log_info "Saving translated VCFs for future convenience"
cp "${IMPUTED_COMMON}" "${IMPUTED_OUT_VCF}"
cp "${IMPUTED_COMMON}.csi" "${IMPUTED_OUT_VCF}.csi" 2>/dev/null || bcftools index -f -c "${IMPUTED_OUT_VCF}"
cp "${TRUTH_COMMON}" "${TRUTH_OUT_VCF}"
cp "${TRUTH_COMMON}.csi" "${TRUTH_OUT_VCF}.csi" 2>/dev/null || bcftools index -f -c "${TRUTH_OUT_VCF}"
log_info "Saved: ${IMPUTED_OUT_VCF}"
log_info "Saved: ${TRUTH_OUT_VCF}"

# Build header (shared by A/B conversion and GT extraction below)
readarray -t SAMPLE_ORDER < "${SAMPLE_SET}"
HEADER="CHROM\tPOS\tREF\tALT\tID"
for s in "${SAMPLE_ORDER[@]}"; do
    HEADER="${HEADER}\t${s}"
done

# =============================================================================
# A/B FORMAT CONVERSION
# =============================================================================
# Produce TSV files with A/B genotypes for both imputed and truth VCFs.
# In the AT/CG coding scheme:
#   Allele index 0  =  AT group  →  A
#   Allele index 1  =  CG group  →  B
# So: 0/0 → A/A,  0/1 or 1/0 → A/B,  1/1 → B/B,  ./. → ./.
# Phased genotypes (|) are also handled: 0|0 → A|A, 0|1 → A|B, etc.

gt_to_ab() {
    # Convert GT indices to A/B notation in a TSV stream (columns 6+ are GT fields).
    # First line (header) is passed through unchanged.
    awk -F'\t' 'BEGIN{OFS="\t"}
    NR==1 {print; next}
    {
        for (i=6; i<=NF; i++) {
            gt = $i
            # Handle phased (|) and unphased (/)
            if      (gt == "0/0" || gt == "0|0") $i = "A/A"
            else if (gt == "0/1" || gt == "0|1") $i = "A/B"
            else if (gt == "1/0" || gt == "1|0") $i = "B/A"
            else if (gt == "1/1" || gt == "1|1") $i = "B/B"
            else $i = "./."
        }
        print
    }'
}

log_info "Generating A/B format genotype files"
IMPUTED_AB_TSV="${OUT_PREFIX}.imputed.AB_format.tsv"
TRUTH_AB_TSV="${OUT_PREFIX}.truth.AB_format.tsv"

# Extract GT as indices, then pipe through gt_to_ab to convert
{
    echo -e "${HEADER}"
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID[\t%GT]\n" -S "${SAMPLE_SET}" "${IMPUTED_COMMON}"
} | gt_to_ab > "${IMPUTED_AB_TSV}"

{
    echo -e "${HEADER}"
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID[\t%GT]\n" -S "${SAMPLE_SET}" "${TRUTH_COMMON}"
} | gt_to_ab > "${TRUTH_AB_TSV}"

IMPUTED_AB_VARIANTS="$(( $(wc -l < "${IMPUTED_AB_TSV}" | tr -d ' ') - 1 ))"
TRUTH_AB_VARIANTS="$(( $(wc -l < "${TRUTH_AB_TSV}" | tr -d ' ') - 1 ))"
log_info "A/B format files: imputed=${IMPUTED_AB_VARIANTS} variants, truth=${TRUTH_AB_VARIANTS} variants"
log_info "Saved: ${IMPUTED_AB_TSV}"
log_info "Saved: ${TRUTH_AB_TSV}"

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
