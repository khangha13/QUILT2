#!/bin/bash
# =============================================================================
# dosage_r2.sh  –  Imputation quality evaluation (concordance & dosage r²)
# =============================================================================
#
# Compares an imputed VCF against a truth (array) VCF to compute per-variant
# dosage r², genotype concordance, and MAF-binned summaries.
#
# Pipeline overview:
#   1. Normalize contig names to canonical ChrNN format
#   2. Find overlapping positions (CHROM + POS) between imputed and truth VCFs
#   3. Deduplicate multi-allelic positions (bcftools norm -d snps)
#   4. Remove strand-ambiguous (A/T, T/A, C/G, G/C) loci
#   5. Translate both VCFs to A/B format:
#        - Imputed: decode GT to nucleotides via REF/ALT, group by AT/CG
#        - Truth:   simple index mapping (0/0→A/A, 0/1→A/B, 1/1→B/B)
#   6. Feed A/B genotype TSVs to R for r², concordance, and plots
#
# Dosages are derived from GT fields only; DS/GP tags are not used.
# Outputs include per-variant metrics, MAF-bin summaries, a per-site/per-sample
# concordance matrix in Parquet (unless disabled), and diagnostic plots.

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

# =============================================================================
# USAGE
# =============================================================================

usage() {
    cat <<'EOF'
Usage: dosage_r2.sh --imputed <imputed.vcf.gz> --truth <truth.vcf.gz> --out-prefix <path> [options]

Required:
  --imputed PATH        Imputed VCF/BCF (raw WGS-style genotypes)
  --truth PATH          Truth VCF/BCF (array genotypes; GT used as reference)
  --out-prefix PREFIX   Output prefix for metrics/plots (e.g., results/dosage_eval)

Options:
  --samples FILE        File with sample IDs to evaluate (one per line).
                        Defaults to the intersection of samples in both VCFs.
  --region STR          Region string (e.g., Chr01:1-1e6) to limit evaluation.
  --use-vcfpp           Additionally run vcfppR comparison on unambiguous VCFs.
  --no-parquet          Skip writing the concordance Parquet file.
  --no-biallelic-only   Do not restrict to biallelic SNPs (default: restrict).
  --no-plots            Skip plotting; only produce metric TSVs.
  --keep-temp           Do not delete the temporary working directory.
  --help                Show this message and exit.

Environment:
  MINIFORGE_MODULE      Module name for miniforge (default: miniforge)
  CONDA_ENV             Conda environment to activate (default: myenv_py310)
  BCFTOOLS_MODULE       Module name for bcftools (default: bcftools/1.18-gcc-12.3.0)

Contig Name Handling:
  The canonical contig style is "ChrNN" (e.g., Chr01, Chr02, Chr17).
  Supported input styles: "Chr01" (ChrNN), "chr1" (chrN), "1" (N).
  Any VCF using a non-canonical style is automatically renamed to ChrNN.

Processing Steps:
  1. Position-only intersection finds common (CHROM, POS) between VCFs.
  2. Multi-allelic positions (decomposed into biallelic rows) are deduplicated.
  3. Strand-ambiguous loci (REF/ALT = A/T, T/A, C/G, G/C) are removed because
     they cannot be reliably assigned to the A or B allele class.
  4. Both VCFs are translated to A/B format:
       - Imputed: nucleotide pairs are grouped by AT/CG:
           Both in {A,T} → A/A | Mixed {A,T}×{C,G} → A/B | Both in {C,G} → B/B
       - Truth: allele indices are mapped directly:
           0/0 → A/A | 0/1 → A/B | 1/1 → B/B
  5. A/B dosages (A/A=0, A/B=1, B/B=2) are compared to compute r² and concordance.

Outputs (PREFIX.*):
  IMPUTED_overlapped_only.vcf.gz              Imputed VCF at common positions (deduped)
  TRUTH_overlapped_only.vcf.gz                Truth VCF at common positions (deduped)
  IMPUTED_overlapped_unambiguous_only.vcf.gz  Imputed VCF after removing ambiguous loci
  TRUTH_overlapped_unambiguous_only.vcf.gz    Truth VCF after removing ambiguous loci
  ambiguous_loci_removed.tsv                  Strand-ambiguous positions that were removed
  imputed.AB_format.tsv                       Imputed genotypes in A/B format
  truth.AB_format.tsv                         Truth genotypes in A/B format
  translation_exceptions.tsv                  Unexpected GTs found during translation (if any)
  metrics.tsv                                 Per-variant r², concordance, MAF
  summary.tsv                                 Overall summary statistics
  maf_bins.tsv                                Metrics aggregated by MAF bins
  concordance.parquet                         Per-site per-sample concordance (0/1/NA)
  r2_hist.png                                 Distribution of r² (if plots enabled)
  concordance_hist.png                        Distribution of concordance (if plots enabled)
  r2_vs_maf.png                               r² vs MAF heatmap (if plots enabled)
EOF
}

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

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
        --imputed)        IMPUTED_VCF="$2"; shift 2 ;;
        --truth)          TRUTH_VCF="$2"; shift 2 ;;
        --out-prefix)     OUT_PREFIX="$2"; shift 2 ;;
        --samples)        SAMPLE_FILE="$2"; shift 2 ;;
        --region)         REGION="$2"; shift 2 ;;
        --use-vcfpp)      USE_VCFPP=true; shift ;;
        --no-parquet)     WRITE_PARQUET=false; shift ;;
        --no-biallelic-only) BIALLELIC_ONLY=false; shift ;;
        --no-plots)       RUN_PLOTS=false; shift ;;
        --keep-temp)      KEEP_TEMP=true; shift ;;
        --help|-h)        usage; exit 0 ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# =============================================================================
# INPUT VALIDATION
# =============================================================================

if [[ -z "${IMPUTED_VCF}" || -z "${TRUTH_VCF}" ]]; then
    usage
    exit 1
fi

# Derive default output prefix from imputed VCF name if not supplied.
if [[ -z "${OUT_PREFIX}" ]]; then
    base="$(basename "${IMPUTED_VCF}")"
    dir="$(cd "$(dirname "${IMPUTED_VCF}")" && pwd)"
    case "${base}" in
        *.vcf.gz) OUT_PREFIX="${dir}/${base%.vcf.gz}" ;;
        *.vcf)    OUT_PREFIX="${dir}/${base%.vcf}" ;;
        *.bcf)    OUT_PREFIX="${dir}/${base%.bcf}" ;;
        *)        OUT_PREFIX="${dir}/${base}" ;;
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

# =============================================================================
# TEMPORARY DIRECTORY + CLEANUP
# =============================================================================

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

# =============================================================================
# ENVIRONMENT SETUP (conda, bcftools, R)
# =============================================================================

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

# Auto-install R if missing from the conda environment.
ensure_rscript() {
    if command -v Rscript >/dev/null 2>&1; then
        return 0
    fi
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
ensure_rscript  || exit 1
require_cmd Rscript || exit 1

# Utility: index a VCF/BCF if it looks like a compressed file.
maybe_index() {
    local vcf="$1"
    if [[ "${vcf}" =~ \.vcf\.gz$ || "${vcf}" =~ \.bcf$ ]]; then
        run_cmd bcftools index -f -c "${vcf}"
    fi
}

maybe_index "${IMPUTED_VCF}"
maybe_index "${TRUTH_VCF}"

# =============================================================================
# STEP 1: CONTIG NAME DETECTION AND NORMALIZATION
# =============================================================================
# The canonical contig style for this pipeline is ChrNN (e.g., Chr01, Chr02).
# Both VCFs are inspected: if a VCF uses a different convention (e.g., "chr1",
# "1"), its contigs are renamed to ChrNN using bcftools annotate --rename-chrs.
# This ensures position-level matching works even if the two VCFs use different
# naming schemes.

detect_contig_style() {
    # Determine the contig naming convention used in a VCF.
    # Returns one of: "ChrNN", "chrN", "N", "other".
    local vcf="$1"
    local first_contig
    first_contig="$(bcftools view -H "${vcf}" 2>/dev/null | head -1 | cut -f1)"
    if [[ -z "${first_contig}" ]]; then
        first_contig="$(bcftools view -h "${vcf}" 2>/dev/null \
            | grep '^##contig=<ID=' | head -1 \
            | sed 's/.*ID=\([^,>]*\).*/\1/')"
    fi

    if   [[ "${first_contig}" =~ ^Chr[0-9]+ ]]; then echo "ChrNN"
    elif [[ "${first_contig}" =~ ^chr[0-9]+ ]]; then echo "chrN"
    elif [[ "${first_contig}" =~ ^[0-9]+$ ]];   then echo "N"
    else echo "other"
    fi
}

build_contig_rename_map() {
    # Build a TSV rename map (OLD<TAB>NEW) from source_style to target_style.
    # Covers chromosomes 1-22 for generality plus X, Y, M, MT.
    local source_style="$1"
    local target_style="$2"
    local map_file="$3"

    : > "${map_file}"

    for i in $(seq 1 22); do
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

    for special in X Y M MT; do
        local src_name tgt_name
        case "${source_style}" in
            ChrNN) src_name="Chr${special}" ;; chrN) src_name="chr${special}" ;;
            N)     src_name="${special}" ;;     *)    src_name="${special}" ;;
        esac
        case "${target_style}" in
            ChrNN) tgt_name="Chr${special}" ;; chrN) tgt_name="chr${special}" ;;
            N)     tgt_name="${special}" ;;     *)    tgt_name="${special}" ;;
        esac
        if [[ "${src_name}" != "${tgt_name}" ]]; then
            printf '%s\t%s\n' "${src_name}" "${tgt_name}" >> "${map_file}"
        fi
    done
}

normalize_contigs() {
    # Rename contigs in a VCF using a rename map produced by build_contig_rename_map.
    local input_vcf="$1"
    local output_vcf="$2"
    local rename_map="$3"

    if [[ ! -s "${rename_map}" ]]; then
        cp "${input_vcf}" "${output_vcf}"
    else
        log_info "Renaming contigs using map: ${rename_map}"
        run_cmd bcftools annotate --rename-chrs "${rename_map}" "${input_vcf}" -Oz -o "${output_vcf}"
    fi
    run_cmd bcftools index -f -c "${output_vcf}"
}

log_info "Detecting contig naming conventions"
IMPUTED_CONTIG_STYLE="$(detect_contig_style "${IMPUTED_VCF}")"
TRUTH_CONTIG_STYLE="$(detect_contig_style "${TRUTH_VCF}")"
log_info "Imputed VCF contig style: ${IMPUTED_CONTIG_STYLE}"
log_info "Truth VCF contig style: ${TRUTH_CONTIG_STYLE}"

CANONICAL_STYLE="ChrNN"

# Normalize imputed VCF contigs if needed.
IMPUTED_NORMALIZED="${IMPUTED_VCF}"
if [[ "${IMPUTED_CONTIG_STYLE}" != "${CANONICAL_STYLE}" ]]; then
    log_warn "Imputed VCF uses '${IMPUTED_CONTIG_STYLE}' contigs. Normalizing to '${CANONICAL_STYLE}'."
    CONTIG_MAP="${TMP_DIR}/imputed_contig_rename.map"
    build_contig_rename_map "${IMPUTED_CONTIG_STYLE}" "${CANONICAL_STYLE}" "${CONTIG_MAP}"
    if [[ -s "${CONTIG_MAP}" ]]; then
        log_info "Imputed contig rename map (first 10 entries):"
        head -10 "${CONTIG_MAP}" >&2 || true
        IMPUTED_NORMALIZED="${TMP_DIR}/imputed.renamed.vcf.gz"
        normalize_contigs "${IMPUTED_VCF}" "${IMPUTED_NORMALIZED}" "${CONTIG_MAP}"
    else
        log_warn "No rename mappings generated for imputed; proceeding with original names."
    fi
fi

# Normalize truth VCF contigs if needed.
TRUTH_NORMALIZED="${TRUTH_VCF}"
if [[ "${TRUTH_CONTIG_STYLE}" != "${CANONICAL_STYLE}" ]]; then
    log_warn "Truth VCF uses '${TRUTH_CONTIG_STYLE}' contigs. Normalizing to '${CANONICAL_STYLE}'."
    CONTIG_MAP="${TMP_DIR}/truth_contig_rename.map"
    build_contig_rename_map "${TRUTH_CONTIG_STYLE}" "${CANONICAL_STYLE}" "${CONTIG_MAP}"
    if [[ -s "${CONTIG_MAP}" ]]; then
        log_info "Truth contig rename map (first 10 entries):"
        head -10 "${CONTIG_MAP}" >&2 || true
        TRUTH_NORMALIZED="${TMP_DIR}/truth.renamed.vcf.gz"
        normalize_contigs "${TRUTH_VCF}" "${TRUTH_NORMALIZED}" "${CONTIG_MAP}"
    else
        log_warn "No rename mappings generated for truth; proceeding with original names."
    fi
fi

# =============================================================================
# STEP 2: POSITION-ONLY INTERSECTION
# =============================================================================
# Extract unique (CHROM, POS) from each VCF and find the overlap using comm.
# NOTE: A future optimization could use `bcftools isec -n=2 -c all` for
# position-only matching, which may be faster on very large files. The current
# manual approach is retained for explicit control and detailed logging.

log_info "Extracting positions from each VCF (position-only intersection)"
IMPUTED_POS="${TMP_DIR}/imputed.pos.txt"
TRUTH_POS="${TMP_DIR}/truth.pos.txt"

bcftools query -f '%CHROM\t%POS\n' "${IMPUTED_NORMALIZED}" | sort -u > "${IMPUTED_POS}"
bcftools query -f '%CHROM\t%POS\n' "${TRUTH_NORMALIZED}"   | sort -u > "${TRUTH_POS}"

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
    log_error "Sample imputed positions:" ; head -5 "${IMPUTED_POS}" >&2 || true
    log_error "Sample truth positions:"   ; head -5 "${TRUTH_POS}" >&2 || true
    exit 1
fi

# Convert to targets format (CHROM<TAB>POS) for bcftools -T.
SITE_LIST="${TMP_DIR}/sites.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2}' "${COMMON_POS}" > "${SITE_LIST}"

# =============================================================================
# SAMPLE SET DERIVATION
# =============================================================================
# Use the user-supplied sample file, or derive the intersection of samples
# present in both VCFs.

SAMPLE_SET="${SAMPLE_FILE}"
if [[ -n "${SAMPLE_FILE}" ]]; then
    if [[ ! -f "${SAMPLE_FILE}" ]]; then
        log_error "Sample file not found: ${SAMPLE_FILE}"
        exit 1
    fi
else
    log_info "Deriving common sample set from both VCFs"
    bcftools query -l "${IMPUTED_VCF}" | sort > "${TMP_DIR}/imputed.samples"
    bcftools query -l "${TRUTH_VCF}"   | sort > "${TMP_DIR}/truth.samples"
    comm -12 "${TMP_DIR}/imputed.samples" "${TMP_DIR}/truth.samples" > "${TMP_DIR}/common.samples"
    if [[ ! -s "${TMP_DIR}/common.samples" ]]; then
        log_error "No overlapping samples found between VCFs"
        exit 1
    fi
    SAMPLE_SET="${TMP_DIR}/common.samples"
fi

SAMPLE_COUNT="$(wc -l < "${SAMPLE_SET}" | tr -d ' ')"
log_info "Evaluating ${SAMPLE_COUNT} samples"

REGION_ARGS=()
if [[ -n "${REGION}" ]]; then
    REGION_ARGS=( -r "${REGION}" )
fi

BIALLELIC_ARGS=()
if [[ "${BIALLELIC_ONLY}" == "true" ]]; then
    BIALLELIC_ARGS=( -m2 -M2 -v snps )
fi

# =============================================================================
# STEP 3: FILTER TO COMMON POSITIONS + DEDUPLICATE
# =============================================================================
# Filter both VCFs to overlapping positions using -T (targets) for exact
# position matching, then deduplicate with `bcftools norm -d snps` to remove
# multi-allelic sites that were decomposed into multiple biallelic rows at the
# same position (e.g., pos 1000: A→G and A→T as separate records).

log_info "Filtering VCFs to ${COMMON_POS_COUNT} common positions"

IMPUTED_PREDEDUP="${TMP_DIR}/imputed.prededup.vcf.gz"
TRUTH_PREDEDUP="${TMP_DIR}/truth.prededup.vcf.gz"

run_cmd bcftools view -T "${SITE_LIST}" -S "${SAMPLE_SET}" "${REGION_ARGS[@]}" "${BIALLELIC_ARGS[@]}" \
    -Oz -o "${IMPUTED_PREDEDUP}" "${IMPUTED_NORMALIZED}"
run_cmd bcftools view -T "${SITE_LIST}" -S "${SAMPLE_SET}" "${REGION_ARGS[@]}" "${BIALLELIC_ARGS[@]}" \
    -Oz -o "${TRUTH_PREDEDUP}" "${TRUTH_NORMALIZED}"
run_cmd bcftools index -f -c "${IMPUTED_PREDEDUP}"
run_cmd bcftools index -f -c "${TRUTH_PREDEDUP}"

IMPUTED_PREDEDUP_N="$(bcftools view -H "${IMPUTED_PREDEDUP}" | wc -l | tr -d ' ')"
TRUTH_PREDEDUP_N="$(bcftools view -H "${TRUTH_PREDEDUP}" | wc -l | tr -d ' ')"

# Deduplicate: keep one SNP per position.
IMPUTED_OVERLAPPED="${TMP_DIR}/imputed.overlapped.vcf.gz"
TRUTH_OVERLAPPED="${TMP_DIR}/truth.overlapped.vcf.gz"

log_info "Deduplicating multi-allelic positions (imputed=${IMPUTED_PREDEDUP_N}, truth=${TRUTH_PREDEDUP_N})"
run_cmd bcftools norm -d snps "${IMPUTED_PREDEDUP}" -Oz -o "${IMPUTED_OVERLAPPED}"
run_cmd bcftools norm -d snps "${TRUTH_PREDEDUP}"   -Oz -o "${TRUTH_OVERLAPPED}"
run_cmd bcftools index -f -c "${IMPUTED_OVERLAPPED}"
run_cmd bcftools index -f -c "${TRUTH_OVERLAPPED}"

IMPUTED_OVERLAPPED_N="$(bcftools view -H "${IMPUTED_OVERLAPPED}" | wc -l | tr -d ' ')"
TRUTH_OVERLAPPED_N="$(bcftools view -H "${TRUTH_OVERLAPPED}" | wc -l | tr -d ' ')"
log_info "Variants after dedup: imputed=${IMPUTED_OVERLAPPED_N}, truth=${TRUTH_OVERLAPPED_N}"

if [[ "${IMPUTED_PREDEDUP_N}" -ne "${IMPUTED_OVERLAPPED_N}" ]]; then
    log_warn "Removed $((IMPUTED_PREDEDUP_N - IMPUTED_OVERLAPPED_N)) duplicate position(s) from imputed VCF"
fi
if [[ "${TRUTH_PREDEDUP_N}" -ne "${TRUTH_OVERLAPPED_N}" ]]; then
    log_warn "Removed $((TRUTH_PREDEDUP_N - TRUTH_OVERLAPPED_N)) duplicate position(s) from truth VCF"
fi

if [[ "${IMPUTED_OVERLAPPED_N}" -eq 0 || "${TRUTH_OVERLAPPED_N}" -eq 0 ]]; then
    log_error "No variants remaining after filtering. Check VCF content and contig naming."
    exit 1
fi

# Save overlapped-only VCFs to output prefix for inspection / reuse.
IMPUTED_OVERLAP_OUT="${OUT_PREFIX}.IMPUTED_overlapped_only.vcf.gz"
TRUTH_OVERLAP_OUT="${OUT_PREFIX}.TRUTH_overlapped_only.vcf.gz"

log_info "Saving overlapped VCFs"
cp "${IMPUTED_OVERLAPPED}" "${IMPUTED_OVERLAP_OUT}"
cp "${IMPUTED_OVERLAPPED}.csi" "${IMPUTED_OVERLAP_OUT}.csi" 2>/dev/null || bcftools index -f -c "${IMPUTED_OVERLAP_OUT}"
cp "${TRUTH_OVERLAPPED}" "${TRUTH_OVERLAP_OUT}"
cp "${TRUTH_OVERLAPPED}.csi" "${TRUTH_OVERLAP_OUT}.csi" 2>/dev/null || bcftools index -f -c "${TRUTH_OVERLAP_OUT}"
log_info "Saved: ${IMPUTED_OVERLAP_OUT}"
log_info "Saved: ${TRUTH_OVERLAP_OUT}"

# =============================================================================
# STEP 4: REMOVE STRAND-AMBIGUOUS LOCI
# =============================================================================
# Strand-ambiguous SNPs have REF/ALT that are complementary pairs:
#   A/T, T/A  (A and T are complements)
#   C/G, G/C  (C and G are complements)
# For these loci it is impossible to determine which strand the array is
# reporting on, so the A/B assignment is unreliable. We remove them from BOTH
# VCFs and record which positions were dropped.

log_info "Identifying strand-ambiguous loci in imputed VCF"

# Extract all variants from the imputed overlapped VCF and flag ambiguous ones.
# Ambiguous = REF/ALT is one of: A/T, T/A, C/G, G/C.
AMBIGUOUS_LOCI="${TMP_DIR}/ambiguous_loci.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${IMPUTED_OVERLAPPED}" \
    | awk -F'\t' 'BEGIN{OFS="\t"}
        ($3=="A" && $4=="T") || ($3=="T" && $4=="A") ||
        ($3=="C" && $4=="G") || ($3=="G" && $4=="C") {print}
      ' > "${AMBIGUOUS_LOCI}"

AMBIG_COUNT="$(wc -l < "${AMBIGUOUS_LOCI}" | tr -d ' ')"
log_info "Found ${AMBIG_COUNT} strand-ambiguous position(s) to remove"

# Write the report to the output prefix (always, even if 0 rows).
AMBIG_REPORT="${OUT_PREFIX}.ambiguous_loci_removed.tsv"
{
    echo -e "CHROM\tPOS\tREF\tALT"
    cat "${AMBIGUOUS_LOCI}"
} > "${AMBIG_REPORT}"
log_info "Ambiguous loci report: ${AMBIG_REPORT}"

# Build a targets exclusion file (CHROM<TAB>POS) for bcftools -T ^file.
AMBIG_POSITIONS="${TMP_DIR}/ambiguous_positions.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2}' "${AMBIGUOUS_LOCI}" > "${AMBIG_POSITIONS}"

# Filter both VCFs: exclude ambiguous positions.
IMPUTED_UNAMBIG="${TMP_DIR}/imputed.unambiguous.vcf.gz"
TRUTH_UNAMBIG="${TMP_DIR}/truth.unambiguous.vcf.gz"

if [[ "${AMBIG_COUNT}" -gt 0 ]]; then
    log_info "Removing ${AMBIG_COUNT} ambiguous loci from both VCFs"
    run_cmd bcftools view -T ^"${AMBIG_POSITIONS}" "${IMPUTED_OVERLAPPED}" -Oz -o "${IMPUTED_UNAMBIG}"
    run_cmd bcftools view -T ^"${AMBIG_POSITIONS}" "${TRUTH_OVERLAPPED}"   -Oz -o "${TRUTH_UNAMBIG}"
else
    log_info "No ambiguous loci to remove; copying overlapped VCFs as-is"
    cp "${IMPUTED_OVERLAPPED}" "${IMPUTED_UNAMBIG}"
    cp "${TRUTH_OVERLAPPED}"   "${TRUTH_UNAMBIG}"
fi
run_cmd bcftools index -f -c "${IMPUTED_UNAMBIG}"
run_cmd bcftools index -f -c "${TRUTH_UNAMBIG}"

IMPUTED_UNAMBIG_N="$(bcftools view -H "${IMPUTED_UNAMBIG}" | wc -l | tr -d ' ')"
TRUTH_UNAMBIG_N="$(bcftools view -H "${TRUTH_UNAMBIG}" | wc -l | tr -d ' ')"
log_info "Variants after ambiguous removal: imputed=${IMPUTED_UNAMBIG_N}, truth=${TRUTH_UNAMBIG_N}"

if [[ "${IMPUTED_UNAMBIG_N}" -eq 0 || "${TRUTH_UNAMBIG_N}" -eq 0 ]]; then
    log_error "No variants remaining after ambiguous loci removal."
    exit 1
fi

# Save unambiguous VCFs to output prefix.
IMPUTED_UNAMBIG_OUT="${OUT_PREFIX}.IMPUTED_overlapped_unambiguous_only.vcf.gz"
TRUTH_UNAMBIG_OUT="${OUT_PREFIX}.TRUTH_overlapped_unambiguous_only.vcf.gz"

log_info "Saving unambiguous VCFs"
cp "${IMPUTED_UNAMBIG}" "${IMPUTED_UNAMBIG_OUT}"
cp "${IMPUTED_UNAMBIG}.csi" "${IMPUTED_UNAMBIG_OUT}.csi" 2>/dev/null || bcftools index -f -c "${IMPUTED_UNAMBIG_OUT}"
cp "${TRUTH_UNAMBIG}" "${TRUTH_UNAMBIG_OUT}"
cp "${TRUTH_UNAMBIG}.csi" "${TRUTH_UNAMBIG_OUT}.csi" 2>/dev/null || bcftools index -f -c "${TRUTH_UNAMBIG_OUT}"
log_info "Saved: ${IMPUTED_UNAMBIG_OUT}"
log_info "Saved: ${TRUTH_UNAMBIG_OUT}"

# =============================================================================
# STEP 5: TRANSLATE BOTH VCFs TO A/B FORMAT
# =============================================================================
# Both VCFs are converted to a common A/B genotype representation:
#
# --- IMPUTED VCF (nucleotide-based translation) ---
# For each variant, the GT allele indices (0 = REF, 1 = ALT) are decoded back
# to actual nucleotides using the REF and ALT columns. Each nucleotide is then
# classified into one of two groups:
#   Group A (AT):  A or T
#   Group B (CG):  C or G
#
# The diploid genotype is then assigned:
#   Both alleles in {A, T}            → A/A  (e.g., AA, AT, TA, TT)
#   One allele {A,T} + one {C,G}     → A/B  (e.g., AC, AG, TC, TG, CA, GA, CT, GT)
#   Both alleles in {C, G}           → B/B  (e.g., CC, CG, GC, GG)
#   Missing (./.)                    → ./.
#   Anything else                    → recorded as exception
#
# Since strand-ambiguous loci (A/T, T/A, C/G, G/C) were already removed in
# Step 4, every remaining variant has one allele from {A,T} and one from {C,G}.
# This means:
#   - Homozygous REF: if REF∈{A,T} → A/A; if REF∈{C,G} → B/B
#   - Heterozygous:   always A/B (one from each group)
#   - Homozygous ALT: if ALT∈{A,T} → A/A; if ALT∈{C,G} → B/B
#
# --- TRUTH VCF (index-based translation) ---
# Array genotypes use TOP-strand allele convention where:
#   Allele index 0 (REF) = A allele
#   Allele index 1 (ALT) = B allele
# Therefore:
#   0/0 → A/A   |   0/1 or 1/0 → A/B   |   1/1 → B/B
#   ./.  → ./.
#   Anything else → recorded as exception

# Build the TSV header (shared by both A/B files).
readarray -t SAMPLE_ORDER < "${SAMPLE_SET}"
HEADER="CHROM\tPOS\tREF\tALT\tID"
for s in "${SAMPLE_ORDER[@]}"; do
    HEADER="${HEADER}\t${s}"
done

# --- Imputed A/B translation ---
# The awk script receives `bcftools query` output with columns:
#   $1=CHROM  $2=POS  $3=REF  $4=ALT  $5=ID  $6..=GT per sample
# For each GT field, it decodes allele indices to nucleotides, classifies them,
# and outputs the A/B genotype. Exceptions (unexpected genotypes) are written
# to a separate file via a redirect.

log_info "Translating imputed VCF to A/B format"
IMPUTED_AB_TSV="${OUT_PREFIX}.imputed.AB_format.tsv"
IMPUTED_EXCEPTIONS="${TMP_DIR}/imputed_exceptions.tsv"

{
    echo -e "${HEADER}"
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID[\t%GT]\n" \
        -S "${SAMPLE_SET}" "${IMPUTED_UNAMBIG}"
} | awk -F'\t' -v exc_file="${IMPUTED_EXCEPTIONS}" '
BEGIN { OFS = "\t" }

# Helper: classify a nucleotide into AT (="A") or CG (="B"), or "?" if unknown.
function nuc_group(base) {
    if (base == "A" || base == "T") return "A"
    if (base == "C" || base == "G") return "B"
    return "?"
}

# Pass the header line through, replacing REF/ALT with A/B.
NR == 1 {
    $3 = "REF"; $4 = "ALT"   # keep column names as-is for clarity
    print
    next
}

# Data lines.
{
    ref = $3   # REF nucleotide
    alt = $4   # ALT nucleotide

    # Replace REF/ALT columns with A and B to reflect the new coding.
    $3 = "A"
    $4 = "B"

    # Translate each sample GT (columns 6 onward).
    for (i = 6; i <= NF; i++) {
        gt = $i

        # Handle missing genotypes.
        if (gt == "./." || gt == ".|.") {
            $i = "./."
            continue
        }

        # Split GT on / or | to get the two allele indices.
        sep = "/"
        n = split(gt, idx, "/")
        if (n != 2) {
            n = split(gt, idx, "|")
            sep = "|"
        }
        if (n != 2) {
            # Unexpected GT format — record as exception.
            print $1, $2, ref, alt, "sample_col=" i, gt >> exc_file
            $i = "./."
            continue
        }

        # Decode allele indices to actual nucleotides.
        # Index 0 = REF, index 1 = ALT. Anything else is unexpected.
        nuc1 = ""
        nuc2 = ""
        if      (idx[1] == "0") nuc1 = ref
        else if (idx[1] == "1") nuc1 = alt
        else if (idx[1] == ".") { $i = "./."; continue }
        else { print $1, $2, ref, alt, "sample_col=" i, gt >> exc_file; $i = "./."; continue }

        if      (idx[2] == "0") nuc2 = ref
        else if (idx[2] == "1") nuc2 = alt
        else if (idx[2] == ".") { $i = "./."; continue }
        else { print $1, $2, ref, alt, "sample_col=" i, gt >> exc_file; $i = "./."; continue }

        # Classify each nucleotide into group A (AT) or B (CG).
        g1 = nuc_group(nuc1)
        g2 = nuc_group(nuc2)

        if (g1 == "?" || g2 == "?") {
            # Non-ACGT nucleotide — exception.
            print $1, $2, ref, alt, "sample_col=" i, gt, nuc1, nuc2 >> exc_file
            $i = "./."
            continue
        }

        $i = g1 sep g2   # e.g., "A/B", "B/A", "A/A", "B/B"
    }

    print
}
' > "${IMPUTED_AB_TSV}"

IMPUTED_AB_N="$(( $(wc -l < "${IMPUTED_AB_TSV}" | tr -d ' ') - 1 ))"
log_info "Imputed A/B format: ${IMPUTED_AB_N} variants written to ${IMPUTED_AB_TSV}"

# --- Truth A/B translation ---
# Simple index-based mapping: 0/0→A/A, 0/1→A/B, 1/0→A/B, 1/1→B/B, ./.→./.
# Any other GT string is recorded as an exception.

log_info "Translating truth VCF to A/B format"
TRUTH_AB_TSV="${OUT_PREFIX}.truth.AB_format.tsv"
TRUTH_EXCEPTIONS="${TMP_DIR}/truth_exceptions.tsv"

{
    echo -e "${HEADER}"
    bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID[\t%GT]\n" \
        -S "${SAMPLE_SET}" "${TRUTH_UNAMBIG}"
} | awk -F'\t' -v exc_file="${TRUTH_EXCEPTIONS}" '
BEGIN { OFS = "\t" }

# Pass the header line through, replacing REF/ALT with A/B.
NR == 1 {
    $3 = "REF"; $4 = "ALT"
    print
    next
}

# Data lines.
{
    ref_orig = $3
    alt_orig = $4

    # Replace REF/ALT columns with A and B.
    $3 = "A"
    $4 = "B"

    # Translate each sample GT (columns 6 onward).
    for (i = 6; i <= NF; i++) {
        gt = $i

        # Standard mappings.
        if      (gt == "0/0" || gt == "0|0") $i = "A/A"
        else if (gt == "0/1" || gt == "0|1") $i = "A/B"
        else if (gt == "1/0" || gt == "1|0") $i = "A/B"
        else if (gt == "1/1" || gt == "1|1") $i = "B/B"
        else if (gt == "./." || gt == ".|.") $i = "./."
        else {
            # Unexpected GT — record as exception.
            print $1, $2, ref_orig, alt_orig, "sample_col=" i, gt >> exc_file
            $i = "./."
        }
    }

    print
}
' > "${TRUTH_AB_TSV}"

TRUTH_AB_N="$(( $(wc -l < "${TRUTH_AB_TSV}" | tr -d ' ') - 1 ))"
log_info "Truth A/B format: ${TRUTH_AB_N} variants written to ${TRUTH_AB_TSV}"

# --- Exception reporting ---
# Merge any exceptions from imputed and truth translations into a single report.
# Only created if exceptions exist.
EXCEPTION_REPORT="${OUT_PREFIX}.translation_exceptions.tsv"
EXCEPTION_COUNT=0

if [[ -s "${IMPUTED_EXCEPTIONS}" || -s "${TRUTH_EXCEPTIONS}" ]]; then
    {
        echo -e "SOURCE\tCHROM\tPOS\tREF\tALT\tCONTEXT\tGT\tEXTRA"
        if [[ -s "${IMPUTED_EXCEPTIONS}" ]]; then
            awk -F'\t' 'BEGIN{OFS="\t"} {print "imputed", $0}' "${IMPUTED_EXCEPTIONS}"
        fi
        if [[ -s "${TRUTH_EXCEPTIONS}" ]]; then
            awk -F'\t' 'BEGIN{OFS="\t"} {print "truth", $0}' "${TRUTH_EXCEPTIONS}"
        fi
    } > "${EXCEPTION_REPORT}"

    EXCEPTION_COUNT="$(( $(wc -l < "${EXCEPTION_REPORT}" | tr -d ' ') - 1 ))"
    log_warn "${EXCEPTION_COUNT} translation exception(s) recorded in: ${EXCEPTION_REPORT}"
else
    log_info "No translation exceptions found"
fi

# =============================================================================
# STEP 6: COMPUTE METRICS VIA R
# =============================================================================
# Feed the A/B format TSVs to the R helper script which:
#   - Converts A/B genotypes to dosages (A/A=0, A/B=1, B/B=2)
#   - Computes per-variant r², concordance, and MAF
#   - Produces MAF-bin summaries
#   - Writes a concordance parquet and diagnostic plots

R_HELPER="${SCRIPT_DIR}/dosage_r2.R"
if [[ ! -f "${R_HELPER}" ]]; then
    log_error "Missing R helper script: ${R_HELPER}"
    exit 1
fi

R_ARGS=(
    "${R_HELPER}"
    "--imputed-gt" "${IMPUTED_AB_TSV}"
    "--truth-gt"   "${TRUTH_AB_TSV}"
    "--samples"    "${SAMPLE_SET}"
    "--out-prefix" "${OUT_PREFIX}"
    "--write-parquet" "${WRITE_PARQUET}"
)

if [[ "${USE_VCFPP}" == "true" ]]; then
    R_ARGS+=( "--use-vcfpp"
              "--vcfpp-imputed" "${IMPUTED_UNAMBIG}"
              "--vcfpp-truth"   "${TRUTH_UNAMBIG}" )
fi

if [[ "${RUN_PLOTS}" == "true" ]]; then
    R_ARGS+=( "--plots" )
fi

log_info "Running R metrics helper"
run_cmd Rscript "${R_ARGS[@]}"

log_info "Dosage r² and concordance metrics written to ${OUT_PREFIX}.*"
