#!/bin/bash
# =============================================================================
# test_concordance_check_with_array_validation.sh  –  Quarto-ready sample concordance matcher
# =============================================================================
#
# Matches samples from a nucleotide VCF/BCF (--vcf1) against samples from a
# truth/array VCF/BCF (--truth/--vcf) using the same overlap and imputed-side
# A/B harmonisation strategy as modules/evaluate/dosage_r2.sh. The truth VCF
# is expected to be array-coded such that GT index 0 means A and GT index 1
# means B; REF/ALT are retained for QC/filtering, not for truth GT decoding.
#
# Outputs are Quarto-friendly TSVs plus intermediate VCF/TSV artifacts for
# download and audit.

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
Usage:
  test_concordance_check_with_array_validation.sh --vcf1 <query.vcf.gz|query.bcf> --truth <truth.vcf.gz|truth.bcf> [options]
  test_concordance_check_with_array_validation.sh --vcf1 <query.vcf.gz|query.bcf> --vcf   <truth.vcf.gz|truth.bcf> [options]

Required:
  --vcf1 PATH             Query VCF/BCF with standard nucleotide alleles.
  --truth PATH            Truth/array VCF/BCF using array A/B genotype coding
                          (GT index 0 -> A, 1 -> B; REF/ALT kept for QC).
  --vcf PATH              Alias for --truth.

Optional:
  --out-prefix PREFIX     Output prefix. Defaults to the basename of --vcf1 in
                          its source directory.
  --samples-vcf1 FILE     Sample IDs to retain from --vcf1, one per line.
  --samples-truth FILE    Sample IDs to retain from --truth, one per line.
  --region STR            Restrict processing to a region, e.g. Chr01:1-5000000.
  --keep-temp             Keep the temporary working directory.
  --force                 Remove existing outputs and recompute everything.
  --help                  Show this message and exit.

Primary outputs:
  <prefix>.best_matches.tsv
  <prefix>.pairwise_concordance.tsv
  <prefix>.pipeline_audit.tsv
  <prefix>.output_manifest.tsv

Intermediate outputs retained by default:
  <prefix>.VCF1_overlapped_only.vcf.gz
  <prefix>.TRUTH_overlapped_only.vcf.gz
  <prefix>.VCF1_overlapped_unambiguous_only.vcf.gz
  <prefix>.TRUTH_overlapped_unambiguous_only.vcf.gz
  <prefix>.vcf1.AB_format.tsv
  <prefix>.truth.AB_format.tsv
  <prefix>.vcf1.translation_exceptions.tsv
  <prefix>.truth.translation_exceptions.tsv
  <prefix>.duplicates_removed.vcf1.tsv
  <prefix>.duplicates_removed.truth.tsv
  <prefix>.ambiguous_loci_removed.tsv
EOF
}

VCF1=""
TRUTH_VCF=""
TRUTH_ALIAS=""
OUT_PREFIX=""
VCF1_SAMPLE_FILE=""
TRUTH_SAMPLE_FILE=""
REGION=""
KEEP_TEMP=false
FORCE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --vcf1)           VCF1="$2"; shift 2 ;;
        --truth)          TRUTH_VCF="$2"; shift 2 ;;
        --vcf)            TRUTH_ALIAS="$2"; shift 2 ;;
        --out-prefix)     OUT_PREFIX="$2"; shift 2 ;;
        --samples-vcf1)   VCF1_SAMPLE_FILE="$2"; shift 2 ;;
        --samples-truth)  TRUTH_SAMPLE_FILE="$2"; shift 2 ;;
        --region)         REGION="$2"; shift 2 ;;
        --keep-temp)      KEEP_TEMP=true; shift ;;
        --force)          FORCE=true; shift ;;
        --help|-h)        usage; exit 0 ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ -n "${TRUTH_VCF}" && -n "${TRUTH_ALIAS}" && "${TRUTH_VCF}" != "${TRUTH_ALIAS}" ]]; then
    log_error "--truth and --vcf were both supplied with different paths"
    exit 1
fi
if [[ -z "${TRUTH_VCF}" ]]; then
    TRUTH_VCF="${TRUTH_ALIAS}"
fi

if [[ -z "${VCF1}" || -z "${TRUTH_VCF}" ]]; then
    usage
    exit 1
fi

if [[ ! -f "${VCF1}" ]]; then
    log_error "VCF1 input not found: ${VCF1}"
    exit 1
fi
if [[ ! -f "${TRUTH_VCF}" ]]; then
    log_error "Truth input not found: ${TRUTH_VCF}"
    exit 1
fi

if [[ -z "${OUT_PREFIX}" ]]; then
    base="$(basename "${VCF1}")"
    dir="$(cd "$(dirname "${VCF1}")" && pwd)"
    case "${base}" in
        *.vcf.gz) OUT_PREFIX="${dir}/${base%.vcf.gz}" ;;
        *.vcf)    OUT_PREFIX="${dir}/${base%.vcf}" ;;
        *.bcf)    OUT_PREFIX="${dir}/${base%.bcf}" ;;
        *)        OUT_PREFIX="${dir}/${base}" ;;
    esac
    log_info "Defaulting --out-prefix to ${OUT_PREFIX}"
fi

OUT_DIR="$(cd "$(dirname "${OUT_PREFIX}")" && pwd)"
OUT_BASENAME="$(basename "${OUT_PREFIX}")"
OUT_PREFIX="${OUT_DIR}/${OUT_BASENAME}"
mkdir -p "${OUT_DIR}"

TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/concordance_check.XXXXXX")"
cleanup() {
    if [[ "${KEEP_TEMP}" == "true" ]]; then
        log_info "Temporary files kept at: ${TMP_DIR}"
    else
        rm -rf "${TMP_DIR}"
    fi
}
trap cleanup EXIT

log_info "Temporary working directory: ${TMP_DIR}"

if ! load_quilt_env; then
    log_warn "Environment bootstrap failed; continuing with the current PATH"
fi

ensure_rscript() {
    if command -v Rscript >/dev/null 2>&1; then
        return 0
    fi
    if ! command -v conda >/dev/null 2>&1; then
        log_error "Rscript not found, and conda is not available to install it."
        return 1
    fi
    log_warn "Rscript not found; installing R packages into the active conda environment"
    if ! conda install -y r-base r-data.table; then
        log_error "Failed to install Rscript/data.table with conda."
        return 1
    fi
    if ! command -v Rscript >/dev/null 2>&1; then
        log_error "Rscript still not found after conda installation."
        return 1
    fi
}

ensure_r_package() {
    local pkg="$1"
    shift
    if Rscript -e "quit(status = !requireNamespace('${pkg}', quietly = TRUE))" >/dev/null 2>&1; then
        return 0
    fi
    if ! command -v conda >/dev/null 2>&1; then
        log_error "R package '${pkg}' is missing, and conda is unavailable."
        return 1
    fi
    log_warn "R package '${pkg}' not found; installing ${*}"
    if ! conda install -y "$@"; then
        log_error "Failed to install required R package(s): $*"
        return 1
    fi
    if ! Rscript -e "quit(status = !requireNamespace('${pkg}', quietly = TRUE))" >/dev/null 2>&1; then
        log_error "R package '${pkg}' still unavailable after installation."
        return 1
    fi
}

require_cmd awk || exit 1
require_cmd comm || exit 1
require_cmd sort || exit 1
ensure_bcftools || exit 1
ensure_rscript || exit 1
ensure_r_package "data.table" r-data.table || exit 1

maybe_index() {
    local vcf="$1"
    if [[ "${vcf}" =~ \.vcf\.gz$ || "${vcf}" =~ \.bcf$ ]]; then
        run_cmd bcftools index -f -c "${vcf}"
    fi
}

count_variants() {
    local vcf="$1"
    bcftools view -H "${vcf}" | wc -l | tr -d ' '
}

count_table_rows() {
    local path="$1"
    if [[ ! -s "${path}" ]]; then
        echo "0"
        return 0
    fi
    local rows
    rows="$(wc -l < "${path}" | tr -d ' ')"
    if [[ "${rows}" -le 1 ]]; then
        echo "0"
    else
        echo "$((rows - 1))"
    fi
}

step_done() {
    [[ "${FORCE}" == "true" ]] && return 1
    for f in "$@"; do
        [[ -s "${f}" ]] || return 1
    done
    return 0
}

PIPELINE_AUDIT_LOG="${OUT_PREFIX}.pipeline_audit.tsv"
BEST_MATCHES_TSV="${OUT_PREFIX}.best_matches.tsv"
PAIRWISE_TSV="${OUT_PREFIX}.pairwise_concordance.tsv"
OUTPUT_MANIFEST="${OUT_PREFIX}.output_manifest.tsv"

VCF1_OVERLAP_OUT="${OUT_PREFIX}.VCF1_overlapped_only.vcf.gz"
TRUTH_OVERLAP_OUT="${OUT_PREFIX}.TRUTH_overlapped_only.vcf.gz"
VCF1_DUP_REPORT="${OUT_PREFIX}.duplicates_removed.vcf1.tsv"
TRUTH_DUP_REPORT="${OUT_PREFIX}.duplicates_removed.truth.tsv"

VCF1_UNAMBIG_OUT="${OUT_PREFIX}.VCF1_overlapped_unambiguous_only.vcf.gz"
TRUTH_UNAMBIG_OUT="${OUT_PREFIX}.TRUTH_overlapped_unambiguous_only.vcf.gz"
AMBIG_REPORT="${OUT_PREFIX}.ambiguous_loci_removed.tsv"

VCF1_AB_TSV="${OUT_PREFIX}.vcf1.AB_format.tsv"
TRUTH_AB_TSV="${OUT_PREFIX}.truth.AB_format.tsv"
VCF1_EXCEPTION_REPORT="${OUT_PREFIX}.vcf1.translation_exceptions.tsv"
TRUTH_EXCEPTION_REPORT="${OUT_PREFIX}.truth.translation_exceptions.tsv"

if [[ "${FORCE}" == "true" ]]; then
    log_warn "--force: deleting existing concordance utility outputs"
    rm -f "${PIPELINE_AUDIT_LOG}" "${BEST_MATCHES_TSV}" "${PAIRWISE_TSV}" "${OUTPUT_MANIFEST}"
    rm -f "${VCF1_OVERLAP_OUT}" "${VCF1_OVERLAP_OUT}.csi" "${TRUTH_OVERLAP_OUT}" "${TRUTH_OVERLAP_OUT}.csi"
    rm -f "${VCF1_UNAMBIG_OUT}" "${VCF1_UNAMBIG_OUT}.csi" "${TRUTH_UNAMBIG_OUT}" "${TRUTH_UNAMBIG_OUT}.csi"
    rm -f "${VCF1_DUP_REPORT}" "${TRUTH_DUP_REPORT}" "${AMBIG_REPORT}"
    rm -f "${VCF1_AB_TSV}" "${TRUTH_AB_TSV}"
    rm -f "${VCF1_EXCEPTION_REPORT}" "${TRUTH_EXCEPTION_REPORT}"
fi

printf "step\ttimestamp\tcount_before\tcount_after\tcount_removed\tnote\n" > "${PIPELINE_AUDIT_LOG}"

log_audit_step() {
    local step="$1"
    local before="$2"
    local after="$3"
    local note="$4"
    local removed="NA"
    if [[ "${before}" =~ ^[0-9]+$ && "${after}" =~ ^[0-9]+$ ]]; then
        removed="$((before - after))"
    fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "${step}" "$(date -Iseconds)" "${before}" "${after}" "${removed}" "${note}" >> "${PIPELINE_AUDIT_LOG}"
}

prepare_sample_set() {
    local label="$1"
    local input_vcf="$2"
    local requested_file="$3"
    local resolved_file="$4"

    local available_file="${TMP_DIR}/${label}.available.samples"
    local available_sorted="${TMP_DIR}/${label}.available.sorted"
    local cleaned_request="${TMP_DIR}/${label}.requested.cleaned"
    local cleaned_sorted="${TMP_DIR}/${label}.requested.sorted"
    local missing_file="${TMP_DIR}/${label}.requested.missing"

    bcftools query -l "${input_vcf}" > "${available_file}"
    if [[ ! -s "${available_file}" ]]; then
        log_error "No samples found in ${label} VCF: ${input_vcf}"
        exit 1
    fi

    sort "${available_file}" > "${available_sorted}"
    local available_count
    available_count="$(wc -l < "${available_file}" | tr -d ' ')"

    if [[ -n "${requested_file}" ]]; then
        if [[ ! -f "${requested_file}" ]]; then
            log_error "Requested sample file for ${label} not found: ${requested_file}"
            exit 1
        fi
        awk 'NF > 0 && !seen[$0]++ { print $0 }' "${requested_file}" > "${cleaned_request}"
        if [[ ! -s "${cleaned_request}" ]]; then
            log_error "Requested sample file for ${label} is empty after removing blank/duplicate entries: ${requested_file}"
            exit 1
        fi
        sort "${cleaned_request}" > "${cleaned_sorted}"
        comm -23 "${cleaned_sorted}" "${available_sorted}" > "${missing_file}" || true
        if [[ -s "${missing_file}" ]]; then
            log_error "Requested ${label} samples were not found in ${input_vcf}. Example:"
            head -5 "${missing_file}" >&2 || true
            exit 1
        fi
        cp "${cleaned_request}" "${resolved_file}"
    else
        cp "${available_file}" "${resolved_file}"
    fi

    if [[ ! -s "${resolved_file}" ]]; then
        log_error "Resolved ${label} sample set is empty"
        exit 1
    fi

    local kept_count
    kept_count="$(wc -l < "${resolved_file}" | tr -d ' ')"
    log_info "${label} samples to evaluate: ${kept_count}"
    log_audit_step "sample_subset_${label}" "${available_count}" "${kept_count}" "sample subset prepared for ${label}"
}

detect_contig_style() {
    local vcf="$1"
    local first_contig
    first_contig="$(bcftools view -H "${vcf}" 2>/dev/null | head -1 | cut -f1)"
    if [[ -z "${first_contig}" ]]; then
        first_contig="$(bcftools view -h "${vcf}" 2>/dev/null \
            | grep '^##contig=<ID=' | head -1 \
            | sed 's/.*ID=\([^,>]*\).*/\1/')"
    fi

    if [[ "${first_contig}" =~ ^Chr[0-9]+ ]]; then
        echo "ChrNN"
    elif [[ "${first_contig}" =~ ^chr[0-9]+ ]]; then
        echo "chrN"
    elif [[ "${first_contig}" =~ ^[0-9]+$ ]]; then
        echo "N"
    else
        echo "other"
    fi
}

build_contig_rename_map() {
    local source_style="$1"
    local target_style="$2"
    local map_file="$3"

    : > "${map_file}"

    for i in $(seq 1 17); do
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
            printf "%s\t%s\n" "${src_name}" "${tgt_name}" >> "${map_file}"
        fi
    done
}

normalize_contigs() {
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

restrict_to_apple_contigs() {
    local input_vcf="$1"
    local output_vcf="$2"

    run_cmd bcftools view -t "${APPLE_CONTIG_LIST}" "${input_vcf}" -Oz -o "${output_vcf}"
    run_cmd bcftools index -f -c "${output_vcf}"
}

report_removed_duplicates() {
    local prededup_vcf="$1"
    local deduped_vcf="$2"
    local report="$3"
    local label="$4"
    local pre_list="${TMP_DIR}/${label}.prededup.variants.txt"
    local post_list="${TMP_DIR}/${label}.deduped.variants.txt"

    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${prededup_vcf}" | sort > "${pre_list}"
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${deduped_vcf}" | sort > "${post_list}"

    {
        echo -e "CHROM\tPOS\tREF\tALT"
        comm -23 "${pre_list}" "${post_list}"
    } > "${report}"
}

assert_matching_positions() {
    local left_vcf="$1"
    local right_vcf="$2"
    local left_label="$3"
    local right_label="$4"
    local left_pos="${TMP_DIR}/${left_label}.positions.tsv"
    local right_pos="${TMP_DIR}/${right_label}.positions.tsv"
    local left_only="${TMP_DIR}/${left_label}.only.positions.tsv"
    local right_only="${TMP_DIR}/${right_label}.only.positions.tsv"

    bcftools query -f '%CHROM\t%POS\n' "${left_vcf}" | sort -u > "${left_pos}"
    bcftools query -f '%CHROM\t%POS\n' "${right_vcf}" | sort -u > "${right_pos}"
    comm -23 "${left_pos}" "${right_pos}" > "${left_only}" || true
    comm -13 "${left_pos}" "${right_pos}" > "${right_only}" || true

    if [[ -s "${left_only}" || -s "${right_only}" ]]; then
        log_error "Post-dedup position sets diverged between ${left_label} and ${right_label}"
        if [[ -s "${left_only}" ]]; then
            log_error "Example ${left_label}-only position: $(head -1 "${left_only}")"
        fi
        if [[ -s "${right_only}" ]]; then
            log_error "Example ${right_label}-only position: $(head -1 "${right_only}")"
        fi
        exit 1
    fi
}

build_ab_header() {
    local sample_file="$1"
    local header="CHROM\tPOS\tREF\tALT\tID"
    local sample
    while IFS= read -r sample || [[ -n "${sample}" ]]; do
        header="${header}\t${sample}"
    done < "${sample_file}"
    printf "%b" "${header}"
}

finalize_exception_report() {
    local tmp_file="$1"
    local out_file="$2"
    {
        echo -e "CHROM\tPOS\tREF\tALT\tCONTEXT\tGT\tEXTRA"
        if [[ -s "${tmp_file}" ]]; then
            awk -F'\t' 'BEGIN { OFS = "\t" }
                {
                    extra = ""
                    if (NF >= 7) {
                        extra = $7
                        for (i = 8; i <= NF; i++) extra = extra ";" $i
                    }
                    print $1, $2, $3, $4, $5, $6, extra
                }' "${tmp_file}"
        fi
    } > "${out_file}"
}

translate_query_vcf_to_ab() {
    local input_vcf="$1"
    local sample_file="$2"
    local output_tsv="$3"
    local exception_tmp="$4"
    local header
    header="$(build_ab_header "${sample_file}")"

    {
        echo -e "${header}"
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID[\t%GT]\n" -S "${sample_file}" "${input_vcf}"
    } | awk -F'\t' -v exc_file="${exception_tmp}" '
BEGIN { OFS = "\t" }

function nuc_group(base) {
    if (base == "A" || base == "T") return "A"
    if (base == "C" || base == "G") return "B"
    return "?"
}

NR == 1 { print; next }

{
    ref = $3
    alt = $4

    for (i = 6; i <= NF; i++) {
        gt = $i

        if (gt == "./." || gt == ".|.") {
            $i = "./."
            continue
        }

        sep = "/"
        n = split(gt, idx, "/")
        if (n != 2) {
            n = split(gt, idx, "|")
            sep = "|"
        }
        if (n != 2) {
            print $1, $2, ref, alt, "sample_col=" i, gt >> exc_file
            $i = "./."
            continue
        }

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

        g1 = nuc_group(nuc1)
        g2 = nuc_group(nuc2)

        if (g1 == "?" || g2 == "?") {
            print $1, $2, ref, alt, "sample_col=" i, gt, nuc1, nuc2 >> exc_file
            $i = "./."
            continue
        }

        $i = g1 sep g2
    }

    print
}' > "${output_tsv}"
}

decode_truth_vcf_ab() {
    local input_vcf="$1"
    local sample_file="$2"
    local output_tsv="$3"
    local exception_tmp="$4"
    local header
    header="$(build_ab_header "${sample_file}")"

    {
        echo -e "${header}"
        bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%ID[\t%GT]\n" -S "${sample_file}" "${input_vcf}"
    } | awk -F'\t' -v exc_file="${exception_tmp}" '
BEGIN { OFS = "\t" }

NR == 1 { print; next }

{
    ref = $3
    alt = $4

    for (i = 6; i <= NF; i++) {
        gt = $i

        if (gt == "./." || gt == ".|.") {
            $i = "./."
            continue
        }

        sep = "/"
        n = split(gt, idx, "/")
        if (n != 2) {
            n = split(gt, idx, "|")
            sep = "|"
        }
        if (n != 2) {
            print $1, $2, ref, alt, "sample_col=" i, gt >> exc_file
            $i = "./."
            continue
        }

        a1 = ""
        a2 = ""
        if      (idx[1] == "0") a1 = "A"
        else if (idx[1] == "1") a1 = "B"
        else if (idx[1] == ".") { $i = "./."; continue }
        else { print $1, $2, ref, alt, "sample_col=" i, gt >> exc_file; $i = "./."; continue }

        if      (idx[2] == "0") a2 = "A"
        else if (idx[2] == "1") a2 = "B"
        else if (idx[2] == ".") { $i = "./."; continue }
        else { print $1, $2, ref, alt, "sample_col=" i, gt >> exc_file; $i = "./."; continue }

        $i = a1 sep a2
    }

    print
}' > "${output_tsv}"
}

write_output_manifest() {
    printf "file_type\tpath\tdescription\n" > "${OUTPUT_MANIFEST}"

    add_manifest_row() {
        local file_type="$1"
        local path="$2"
        local description="$3"
        printf "%s\t%s\t%s\n" "${file_type}" "${path}" "${description}" >> "${OUTPUT_MANIFEST}"
    }

    add_manifest_row "best_matches" "${BEST_MATCHES_TSV}" "Best and runner-up truth sample for each VCF1 sample"
    add_manifest_row "pairwise_concordance" "${PAIRWISE_TSV}" "Tidy all-vs-all sample concordance table"
    add_manifest_row "pipeline_audit" "${PIPELINE_AUDIT_LOG}" "Per-step counts for QC and Quarto reporting"
    add_manifest_row "vcf1_ab_tsv" "${VCF1_AB_TSV}" "VCF1 genotypes translated into A/B format"
    add_manifest_row "truth_ab_tsv" "${TRUTH_AB_TSV}" "Truth genotypes decoded from array A/B genotype indices"
    add_manifest_row "vcf1_translation_exceptions" "${VCF1_EXCEPTION_REPORT}" "VCF1 translation exceptions encountered during A/B conversion"
    add_manifest_row "truth_translation_exceptions" "${TRUTH_EXCEPTION_REPORT}" "Truth genotype decode exceptions encountered while converting array A/B indices"
    add_manifest_row "duplicates_removed_vcf1" "${VCF1_DUP_REPORT}" "Duplicate loci removed from the VCF1 overlap set"
    add_manifest_row "duplicates_removed_truth" "${TRUTH_DUP_REPORT}" "Duplicate loci removed from the truth overlap set"
    add_manifest_row "ambiguous_loci_removed" "${AMBIG_REPORT}" "Strand-ambiguous loci removed from both inputs"
    add_manifest_row "vcf1_overlap_vcf" "${VCF1_OVERLAP_OUT}" "VCF1 filtered to overlapping positions and deduplicated loci"
    add_manifest_row "truth_overlap_vcf" "${TRUTH_OVERLAP_OUT}" "Truth filtered to overlapping positions and deduplicated loci"
    add_manifest_row "vcf1_unambiguous_vcf" "${VCF1_UNAMBIG_OUT}" "VCF1 overlap set after removing strand-ambiguous loci"
    add_manifest_row "truth_unambiguous_vcf" "${TRUTH_UNAMBIG_OUT}" "Truth overlap set after removing strand-ambiguous loci"
    add_manifest_row "output_manifest" "${OUTPUT_MANIFEST}" "Manifest describing all concordance utility outputs"
}

maybe_index "${VCF1}"
maybe_index "${TRUTH_VCF}"

VCF1_SAMPLE_SET="${TMP_DIR}/vcf1.samples.txt"
TRUTH_SAMPLE_SET="${TMP_DIR}/truth.samples.txt"
prepare_sample_set "vcf1" "${VCF1}" "${VCF1_SAMPLE_FILE}" "${VCF1_SAMPLE_SET}"
prepare_sample_set "truth" "${TRUTH_VCF}" "${TRUTH_SAMPLE_FILE}" "${TRUTH_SAMPLE_SET}"

REGION_ARGS=()
if [[ -n "${REGION}" ]]; then
    REGION_ARGS=( -r "${REGION}" )
fi

BIALLELIC_ARGS=( -m2 -M2 -v snps )
APPLE_CONTIG_LIST="Chr01,Chr02,Chr03,Chr04,Chr05,Chr06,Chr07,Chr08,Chr09,Chr10,Chr11,Chr12,Chr13,Chr14,Chr15,Chr16,Chr17"

VCF1_CONTIG_STYLE="$(detect_contig_style "${VCF1}")"
TRUTH_CONTIG_STYLE="$(detect_contig_style "${TRUTH_VCF}")"
CANONICAL_STYLE="ChrNN"

VCF1_NORMALIZED="${VCF1}"
if [[ "${VCF1_CONTIG_STYLE}" != "${CANONICAL_STYLE}" ]]; then
    VCF1_CONTIG_MAP="${TMP_DIR}/vcf1_contig_rename.map"
    build_contig_rename_map "${VCF1_CONTIG_STYLE}" "${CANONICAL_STYLE}" "${VCF1_CONTIG_MAP}"
    if [[ -s "${VCF1_CONTIG_MAP}" ]]; then
        VCF1_NORMALIZED="${TMP_DIR}/vcf1.normalized.vcf.gz"
        normalize_contigs "${VCF1}" "${VCF1_NORMALIZED}" "${VCF1_CONTIG_MAP}"
    else
        log_warn "No VCF1 contig rename mappings were generated; using the original file"
    fi
fi

TRUTH_NORMALIZED="${TRUTH_VCF}"
if [[ "${TRUTH_CONTIG_STYLE}" != "${CANONICAL_STYLE}" ]]; then
    TRUTH_CONTIG_MAP="${TMP_DIR}/truth_contig_rename.map"
    build_contig_rename_map "${TRUTH_CONTIG_STYLE}" "${CANONICAL_STYLE}" "${TRUTH_CONTIG_MAP}"
    if [[ -s "${TRUTH_CONTIG_MAP}" ]]; then
        TRUTH_NORMALIZED="${TMP_DIR}/truth.normalized.vcf.gz"
        normalize_contigs "${TRUTH_VCF}" "${TRUTH_NORMALIZED}" "${TRUTH_CONTIG_MAP}"
    else
        log_warn "No truth contig rename mappings were generated; using the original file"
    fi
fi

VCF1_APPLE_ONLY="${TMP_DIR}/vcf1.apple_only.vcf.gz"
TRUTH_APPLE_ONLY="${TMP_DIR}/truth.apple_only.vcf.gz"
log_info "Restricting concordance inputs to canonical apple chromosomes (${APPLE_CONTIG_LIST})"
restrict_to_apple_contigs "${VCF1_NORMALIZED}" "${VCF1_APPLE_ONLY}"
restrict_to_apple_contigs "${TRUTH_NORMALIZED}" "${TRUTH_APPLE_ONLY}"

if step_done "${VCF1_OVERLAP_OUT}" "${TRUTH_OVERLAP_OUT}" "${VCF1_DUP_REPORT}" "${TRUTH_DUP_REPORT}"; then
    log_info "[SKIP] Overlap/dedup step outputs already exist"
    VCF1_OVERLAPPED_N="$(count_variants "${VCF1_OVERLAP_OUT}")"
    TRUTH_OVERLAPPED_N="$(count_variants "${TRUTH_OVERLAP_OUT}")"
    log_audit_step "overlap_dedup_vcf1" "${VCF1_OVERLAPPED_N}" "${VCF1_OVERLAPPED_N}" "reused existing overlapping VCF1 output"
    log_audit_step "overlap_dedup_truth" "${TRUTH_OVERLAPPED_N}" "${TRUTH_OVERLAPPED_N}" "reused existing overlapping truth output"
else
    log_info "Filtering both inputs before position intersection"
    VCF1_FILTERED="${TMP_DIR}/vcf1.filtered.vcf.gz"
    TRUTH_FILTERED="${TMP_DIR}/truth.filtered.vcf.gz"
    run_cmd bcftools view -S "${VCF1_SAMPLE_SET}" "${REGION_ARGS[@]}" "${BIALLELIC_ARGS[@]}" \
        -Oz -o "${VCF1_FILTERED}" "${VCF1_APPLE_ONLY}"
    run_cmd bcftools view -S "${TRUTH_SAMPLE_SET}" "${REGION_ARGS[@]}" "${BIALLELIC_ARGS[@]}" \
        -Oz -o "${TRUTH_FILTERED}" "${TRUTH_APPLE_ONLY}"
    run_cmd bcftools index -f -c "${VCF1_FILTERED}"
    run_cmd bcftools index -f -c "${TRUTH_FILTERED}"

    log_info "Extracting common positions between VCF1 and truth"
    VCF1_POS="${TMP_DIR}/vcf1.positions.tsv"
    TRUTH_POS="${TMP_DIR}/truth.positions.tsv"
    COMMON_POS="${TMP_DIR}/common.positions.tsv"
    SITE_LIST="${TMP_DIR}/common.sites.tsv"

    bcftools query -f '%CHROM\t%POS\n' "${VCF1_FILTERED}" | sort -u > "${VCF1_POS}"
    bcftools query -f '%CHROM\t%POS\n' "${TRUTH_FILTERED}" | sort -u > "${TRUTH_POS}"
    comm -12 "${VCF1_POS}" "${TRUTH_POS}" > "${COMMON_POS}"

    COMMON_POS_COUNT="$(wc -l < "${COMMON_POS}" | tr -d ' ')"
    if [[ "${COMMON_POS_COUNT}" -eq 0 ]]; then
        log_error "No overlapping positions found between VCF1 and truth"
        exit 1
    fi
    awk -F'\t' 'BEGIN { OFS = "\t" } { print $1, $2 }' "${COMMON_POS}" > "${SITE_LIST}"

    VCF1_PREDEDUP="${TMP_DIR}/vcf1.prededup.vcf.gz"
    TRUTH_PREDEDUP="${TMP_DIR}/truth.prededup.vcf.gz"
    VCF1_OVERLAPPED_TMP="${TMP_DIR}/vcf1.overlapped.vcf.gz"
    TRUTH_OVERLAPPED_TMP="${TMP_DIR}/truth.overlapped.vcf.gz"

    run_cmd bcftools view -T "${SITE_LIST}" -Oz -o "${VCF1_PREDEDUP}" "${VCF1_FILTERED}"
    run_cmd bcftools view -T "${SITE_LIST}" -Oz -o "${TRUTH_PREDEDUP}" "${TRUTH_FILTERED}"
    run_cmd bcftools index -f -c "${VCF1_PREDEDUP}"
    run_cmd bcftools index -f -c "${TRUTH_PREDEDUP}"

    VCF1_PREDEDUP_N="$(count_variants "${VCF1_PREDEDUP}")"
    TRUTH_PREDEDUP_N="$(count_variants "${TRUTH_PREDEDUP}")"

    run_cmd bcftools norm -d snps "${VCF1_PREDEDUP}" -Oz -o "${VCF1_OVERLAPPED_TMP}"
    run_cmd bcftools norm -d snps "${TRUTH_PREDEDUP}" -Oz -o "${TRUTH_OVERLAPPED_TMP}"
    run_cmd bcftools index -f -c "${VCF1_OVERLAPPED_TMP}"
    run_cmd bcftools index -f -c "${TRUTH_OVERLAPPED_TMP}"

    VCF1_OVERLAPPED_N="$(count_variants "${VCF1_OVERLAPPED_TMP}")"
    TRUTH_OVERLAPPED_N="$(count_variants "${TRUTH_OVERLAPPED_TMP}")"
    assert_matching_positions "${VCF1_OVERLAPPED_TMP}" "${TRUTH_OVERLAPPED_TMP}" "vcf1_dedup" "truth_dedup"

    if [[ "${VCF1_PREDEDUP_N}" -ne "${VCF1_OVERLAPPED_N}" ]]; then
        report_removed_duplicates "${VCF1_PREDEDUP}" "${VCF1_OVERLAPPED_TMP}" "${VCF1_DUP_REPORT}" "vcf1"
    else
        echo -e "CHROM\tPOS\tREF\tALT" > "${VCF1_DUP_REPORT}"
    fi
    if [[ "${TRUTH_PREDEDUP_N}" -ne "${TRUTH_OVERLAPPED_N}" ]]; then
        report_removed_duplicates "${TRUTH_PREDEDUP}" "${TRUTH_OVERLAPPED_TMP}" "${TRUTH_DUP_REPORT}" "truth"
    else
        echo -e "CHROM\tPOS\tREF\tALT" > "${TRUTH_DUP_REPORT}"
    fi

    if [[ "${VCF1_OVERLAPPED_N}" -eq 0 || "${TRUTH_OVERLAPPED_N}" -eq 0 ]]; then
        log_error "No variants remain after overlap filtering and deduplication"
        exit 1
    fi

    cp "${VCF1_OVERLAPPED_TMP}" "${VCF1_OVERLAP_OUT}"
    cp "${TRUTH_OVERLAPPED_TMP}" "${TRUTH_OVERLAP_OUT}"
    run_cmd bcftools index -f -c "${VCF1_OVERLAP_OUT}"
    run_cmd bcftools index -f -c "${TRUTH_OVERLAP_OUT}"

    log_audit_step "overlap_dedup_vcf1" "${VCF1_PREDEDUP_N}" "${VCF1_OVERLAPPED_N}" "positions intersected and duplicate loci removed for VCF1"
    log_audit_step "overlap_dedup_truth" "${TRUTH_PREDEDUP_N}" "${TRUTH_OVERLAPPED_N}" "positions intersected and duplicate loci removed for truth"
fi

if step_done "${VCF1_UNAMBIG_OUT}" "${TRUTH_UNAMBIG_OUT}" "${AMBIG_REPORT}"; then
    log_info "[SKIP] Strand ambiguity filtering outputs already exist"
    VCF1_UNAMBIG_N="$(count_variants "${VCF1_UNAMBIG_OUT}")"
    TRUTH_UNAMBIG_N="$(count_variants "${TRUTH_UNAMBIG_OUT}")"
    log_audit_step "remove_ambiguous_vcf1" "${VCF1_UNAMBIG_N}" "${VCF1_UNAMBIG_N}" "reused existing unambiguous VCF1 output"
    log_audit_step "remove_ambiguous_truth" "${TRUTH_UNAMBIG_N}" "${TRUTH_UNAMBIG_N}" "reused existing unambiguous truth output"
else
    log_info "Identifying strand-ambiguous loci"
    AMBIGUOUS_LOCI="${TMP_DIR}/ambiguous_loci.tsv"
    {
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${VCF1_OVERLAP_OUT}" \
            | awk -F'\t' 'BEGIN { OFS = "\t" }
                ($3=="A" && $4=="T") || ($3=="T" && $4=="A") ||
                ($3=="C" && $4=="G") || ($3=="G" && $4=="C") { print $1, $2, $3, $4, "vcf1" }'
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${TRUTH_OVERLAP_OUT}" \
            | awk -F'\t' 'BEGIN { OFS = "\t" }
                ($3=="A" && $4=="T") || ($3=="T" && $4=="A") ||
                ($3=="C" && $4=="G") || ($3=="G" && $4=="C") { print $1, $2, $3, $4, "truth" }'
    } | sort -u > "${AMBIGUOUS_LOCI}"

    {
        echo -e "CHROM\tPOS\tREF\tALT\tSOURCE"
        cat "${AMBIGUOUS_LOCI}"
    } > "${AMBIG_REPORT}"

    AMBIG_POSITIONS="${TMP_DIR}/ambiguous_positions.tsv"
    awk -F'\t' 'BEGIN { OFS = "\t" } { print $1, $2 }' "${AMBIGUOUS_LOCI}" | sort -u > "${AMBIG_POSITIONS}"

    VCF1_OVERLAPPED_N="$(count_variants "${VCF1_OVERLAP_OUT}")"
    TRUTH_OVERLAPPED_N="$(count_variants "${TRUTH_OVERLAP_OUT}")"

    if [[ -s "${AMBIG_POSITIONS}" ]]; then
        run_cmd bcftools view -T ^"${AMBIG_POSITIONS}" "${VCF1_OVERLAP_OUT}" -Oz -o "${VCF1_UNAMBIG_OUT}"
        run_cmd bcftools view -T ^"${AMBIG_POSITIONS}" "${TRUTH_OVERLAP_OUT}" -Oz -o "${TRUTH_UNAMBIG_OUT}"
    else
        cp "${VCF1_OVERLAP_OUT}" "${VCF1_UNAMBIG_OUT}"
        cp "${TRUTH_OVERLAP_OUT}" "${TRUTH_UNAMBIG_OUT}"
    fi
    run_cmd bcftools index -f -c "${VCF1_UNAMBIG_OUT}"
    run_cmd bcftools index -f -c "${TRUTH_UNAMBIG_OUT}"

    VCF1_UNAMBIG_N="$(count_variants "${VCF1_UNAMBIG_OUT}")"
    TRUTH_UNAMBIG_N="$(count_variants "${TRUTH_UNAMBIG_OUT}")"

    if [[ "${VCF1_UNAMBIG_N}" -eq 0 || "${TRUTH_UNAMBIG_N}" -eq 0 ]]; then
        log_error "No variants remain after removing strand-ambiguous loci"
        exit 1
    fi

    log_audit_step "remove_ambiguous_vcf1" "${VCF1_OVERLAPPED_N}" "${VCF1_UNAMBIG_N}" "strand-ambiguous loci removed from VCF1 overlap set"
    log_audit_step "remove_ambiguous_truth" "${TRUTH_OVERLAPPED_N}" "${TRUTH_UNAMBIG_N}" "strand-ambiguous loci removed from truth overlap set"
fi

if step_done "${VCF1_AB_TSV}" "${TRUTH_AB_TSV}" "${VCF1_EXCEPTION_REPORT}" "${TRUTH_EXCEPTION_REPORT}"; then
    log_info "[SKIP] A/B preparation outputs already exist"
    VCF1_AB_N="$(count_table_rows "${VCF1_AB_TSV}")"
    TRUTH_AB_N="$(count_table_rows "${TRUTH_AB_TSV}")"
    log_audit_step "ab_translate_vcf1" "${VCF1_AB_N}" "${VCF1_AB_N}" "reused existing VCF1 A/B translation output"
    log_audit_step "ab_decode_truth" "${TRUTH_AB_N}" "${TRUTH_AB_N}" "reused existing truth A/B decode output"
else
    log_info "Preparing A/B genotype tables"
    VCF1_EXCEPTION_TMP="${TMP_DIR}/vcf1.translation_exceptions.tmp.tsv"
    TRUTH_EXCEPTION_TMP="${TMP_DIR}/truth.translation_exceptions.tmp.tsv"
    : > "${VCF1_EXCEPTION_TMP}"
    : > "${TRUTH_EXCEPTION_TMP}"

    translate_query_vcf_to_ab "${VCF1_UNAMBIG_OUT}" "${VCF1_SAMPLE_SET}" "${VCF1_AB_TSV}" "${VCF1_EXCEPTION_TMP}"
    decode_truth_vcf_ab "${TRUTH_UNAMBIG_OUT}" "${TRUTH_SAMPLE_SET}" "${TRUTH_AB_TSV}" "${TRUTH_EXCEPTION_TMP}"

    finalize_exception_report "${VCF1_EXCEPTION_TMP}" "${VCF1_EXCEPTION_REPORT}"
    finalize_exception_report "${TRUTH_EXCEPTION_TMP}" "${TRUTH_EXCEPTION_REPORT}"

    VCF1_UNAMBIG_N="$(count_variants "${VCF1_UNAMBIG_OUT}")"
    TRUTH_UNAMBIG_N="$(count_variants "${TRUTH_UNAMBIG_OUT}")"
    VCF1_AB_N="$(count_table_rows "${VCF1_AB_TSV}")"
    TRUTH_AB_N="$(count_table_rows "${TRUTH_AB_TSV}")"

    log_audit_step "ab_translate_vcf1" "${VCF1_UNAMBIG_N}" "${VCF1_AB_N}" "VCF1 translated to A/B genotype space"
    log_audit_step "ab_decode_truth" "${TRUTH_UNAMBIG_N}" "${TRUTH_AB_N}" "truth decoded from array A/B genotype indices"
fi

PAIRWISE_R_SCRIPT="${TMP_DIR}/pairwise_concordance.R"
cat > "${PAIRWISE_R_SCRIPT}" <<'RSCRIPT'
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  opts <- list(vcf1_gt = NULL, truth_gt = NULL, pairwise_out = NULL, best_out = NULL)
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!key %in% c("--vcf1-gt", "--truth-gt", "--pairwise-out", "--best-out")) {
      stop("Unknown argument: ", key)
    }
    if (i == length(args)) stop(key, " requires a value")
    val <- args[[i + 1]]
    switch(key,
      "--vcf1-gt" = opts$vcf1_gt <- val,
      "--truth-gt" = opts$truth_gt <- val,
      "--pairwise-out" = opts$pairwise_out <- val,
      "--best-out" = opts$best_out <- val
    )
    i <- i + 2
  }
  missing <- names(opts)[vapply(opts, is.null, logical(1))]
  if (length(missing) > 0) stop("Missing required arguments: ", paste(missing, collapse = ", "))
  opts
}

opts <- parse_args(args)

load_gt_table <- function(path) {
  dt <- fread(path, na.strings = c(".", "NA"), showProgress = FALSE)
  if (ncol(dt) < 6) stop("Expected metadata + at least one sample column in ", path)
  meta <- dt[, c("CHROM", "POS", "REF", "ALT", "ID"), with = FALSE]
  mat <- as.matrix(dt[, -(1:5)])
  list(meta = meta, mat = mat, samples = colnames(dt)[-(1:5)])
}

variant_key <- function(meta_dt) {
  paste(meta_dt$CHROM, meta_dt$POS, sep = ":")
}

align_to_truth <- function(truth_meta, vcf1_meta, vcf1_mat) {
  tk <- variant_key(truth_meta)
  vk <- variant_key(vcf1_meta)

  if (anyDuplicated(tk)) stop("Truth table has duplicated variants by (CHROM,POS)")
  if (anyDuplicated(vk)) stop("VCF1 table has duplicated variants by (CHROM,POS)")

  idx <- match(tk, vk)
  if (any(is.na(idx))) {
    missing_keys <- tk[is.na(idx)]
    stop("VCF1 table is missing ", length(missing_keys), " truth variants by (CHROM,POS). Example missing key: ", missing_keys[[1]])
  }

  list(meta = truth_meta, mat = vcf1_mat[idx, , drop = FALSE])
}

gt_to_dosage <- function(gt_mat) {
  flat <- gsub("\\|", "/", gt_mat)
  dosage <- rep(NA_real_, length(flat))
  dosage[flat == "A/A"] <- 0
  dosage[flat == "A/B" | flat == "B/A"] <- 1
  dosage[flat == "B/B"] <- 2
  dim(dosage) <- dim(gt_mat)
  dosage
}

same_or_both_na <- function(x, y) {
  (is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)
}

truth_gt <- load_gt_table(opts$truth_gt)
vcf1_gt <- load_gt_table(opts$vcf1_gt)
aligned <- align_to_truth(truth_gt$meta, vcf1_gt$meta, vcf1_gt$mat)
vcf1_gt$meta <- aligned$meta
vcf1_gt$mat <- aligned$mat

vcf1_dosage <- gt_to_dosage(vcf1_gt$mat)
truth_dosage <- gt_to_dosage(truth_gt$mat)

vcf1_samples <- vcf1_gt$samples
truth_samples <- truth_gt$samples

if (length(vcf1_samples) == 0 || length(truth_samples) == 0) {
  stop("Both VCF1 and truth must contain at least one sample")
}

pairwise_list <- vector("list", length(vcf1_samples) * length(truth_samples))
k <- 1L
for (i in seq_along(vcf1_samples)) {
  q <- vcf1_dosage[, i]
  for (j in seq_along(truth_samples)) {
    t <- truth_dosage[, j]
    keep <- !(is.na(q) | is.na(t))
    n_compared <- sum(keep)
    mismatch_count <- if (n_compared > 0) sum(q[keep] != t[keep]) else NA_integer_
    concordance <- if (n_compared > 0) 1 - (mismatch_count / n_compared) else NA_real_
    pairwise_list[[k]] <- data.table(
      vcf1_sample = vcf1_samples[[i]],
      truth_sample = truth_samples[[j]],
      concordance = concordance,
      n_compared = n_compared,
      mismatch_count = mismatch_count
    )
    k <- k + 1L
  }
}

pairwise_dt <- rbindlist(pairwise_list)
if (!any(pairwise_dt$n_compared > 0, na.rm = TRUE)) {
  stop("Zero comparable genotypes across all VCF1/truth sample pairs")
}

pairwise_dt[, sort_concordance := fifelse(is.na(concordance), -Inf, concordance)]
pairwise_dt[, sort_n_compared := fifelse(is.na(n_compared), -1L, n_compared)]
pairwise_dt[, sort_mismatch := fifelse(is.na(mismatch_count), .Machine$integer.max, mismatch_count)]

setorder(pairwise_dt, vcf1_sample, -sort_concordance, -sort_n_compared, sort_mismatch, truth_sample)
pairwise_dt[, match_rank_within_vcf1 := seq_len(.N), by = vcf1_sample]

pairwise_out <- pairwise_dt[, .(
  vcf1_sample,
  truth_sample,
  concordance,
  n_compared,
  mismatch_count,
  match_rank_within_vcf1
)]
fwrite(pairwise_out, opts$pairwise_out, sep = "\t")

best_dt <- pairwise_dt[match_rank_within_vcf1 == 1L, .(
  vcf1_sample,
  best_truth_sample = truth_sample,
  best_concordance = concordance,
  best_n_compared = n_compared,
  best_mismatch_count = mismatch_count
)]

second_dt <- pairwise_dt[match_rank_within_vcf1 == 2L, .(
  vcf1_sample,
  second_truth_sample = truth_sample,
  second_concordance = concordance,
  second_n_compared = n_compared
)]

tie_dt <- pairwise_dt[, {
  top <- .SD[1L]
  .(
    tie_flag = .SD[
      same_or_both_na(concordance, top$concordance[[1]]) &
      same_or_both_na(n_compared, top$n_compared[[1]]) &
      same_or_both_na(mismatch_count, top$mismatch_count[[1]])
    , .N] > 1L
  )
}, by = vcf1_sample]

best_dt <- merge(best_dt, second_dt, by = "vcf1_sample", all.x = TRUE, sort = FALSE)
best_dt <- merge(best_dt, tie_dt, by = "vcf1_sample", all.x = TRUE, sort = FALSE)
best_dt[, concordance_delta := fifelse(
  is.na(best_concordance) | is.na(second_concordance),
  NA_real_,
  best_concordance - second_concordance
)]

setcolorder(best_dt, c(
  "vcf1_sample",
  "best_truth_sample",
  "best_concordance",
  "best_n_compared",
  "best_mismatch_count",
  "second_truth_sample",
  "second_concordance",
  "second_n_compared",
  "concordance_delta",
  "tie_flag"
))
fwrite(best_dt, opts$best_out, sep = "\t")
RSCRIPT

if step_done "${PAIRWISE_TSV}" "${BEST_MATCHES_TSV}"; then
    log_info "[SKIP] Pairwise concordance outputs already exist"
    PAIRWISE_ROWS="$(count_table_rows "${PAIRWISE_TSV}")"
    BEST_ROWS="$(count_table_rows "${BEST_MATCHES_TSV}")"
    log_audit_step "pairwise_concordance" "${PAIRWISE_ROWS}" "${PAIRWISE_ROWS}" "reused existing pairwise concordance output"
    log_audit_step "best_match_summary" "${BEST_ROWS}" "${BEST_ROWS}" "reused existing best-match summary output"
else
    run_cmd Rscript "${PAIRWISE_R_SCRIPT}" \
        --vcf1-gt "${VCF1_AB_TSV}" \
        --truth-gt "${TRUTH_AB_TSV}" \
        --pairwise-out "${PAIRWISE_TSV}" \
        --best-out "${BEST_MATCHES_TSV}"

    PAIRWISE_ROWS="$(count_table_rows "${PAIRWISE_TSV}")"
    BEST_ROWS="$(count_table_rows "${BEST_MATCHES_TSV}")"
    VCF1_SAMPLE_COUNT="$(wc -l < "${VCF1_SAMPLE_SET}" | tr -d ' ')"
    TRUTH_SAMPLE_COUNT="$(wc -l < "${TRUTH_SAMPLE_SET}" | tr -d ' ')"
    SAMPLE_PAIR_COUNT="$((VCF1_SAMPLE_COUNT * TRUTH_SAMPLE_COUNT))"

    log_audit_step "pairwise_concordance" "${SAMPLE_PAIR_COUNT}" "${PAIRWISE_ROWS}" "all VCF1/truth sample pairs scored for genotype concordance"
    log_audit_step "best_match_summary" "${VCF1_SAMPLE_COUNT}" "${BEST_ROWS}" "best and runner-up truth matches selected for each VCF1 sample"
fi

write_output_manifest
log_info "Best matches written to ${BEST_MATCHES_TSV}"
log_info "Pairwise concordance written to ${PAIRWISE_TSV}"
log_info "Output manifest written to ${OUTPUT_MANIFEST}"
