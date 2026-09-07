#!/usr/bin/env bash
# Prepare fixed HortQG array-position masks and run the unchanged WGS evaluator.
#
# This is a Bunya-only, manuscript-specific utility. It deliberately keeps all
# paths and run selections in the configuration block below instead of adding
# another interface to the QUILT2 pipeline.
# Despite its historical name, the unchanged WGS evaluator compares QUILT2
# FORMAT/GT hard-call ALT counts with truth GT; it does not read FORMAT/DS.
#
# Submit every invocation as a modest Slurm coordinator job; preprocessing is
# intentionally refused on a login node. Suggested sequence:
#   1. DRY_RUN=true,  EXECUTION_SCOPE="full"    # validate all 60 analyses
#   2. DRY_RUN=false, EXECUTION_SCOPE="pilot"   # Chr01 RosBREED pilot
#   3. DRY_RUN=false, EXECUTION_SCOPE="one_run" # five masks for one run
#   4. DRY_RUN=false, EXECUTION_SCOPE="full"    # all remaining analyses
#   5. Rerun step 4 after Slurm completion to validate and stage TSV outputs.

set -euo pipefail
export LC_ALL=C

# =============================================================================
# Hardcoded Bunya configuration
# =============================================================================

DRY_RUN=true
FORCE=false
RESUBMIT_INCOMPLETE=false
EXECUTION_SCOPE="full" # pilot | one_run | full

QUILT2_ROOT="/scratch/project/bigdata_apple/QUILT2"
IMPUTED_BASE="/scratch/project/bigdata_apple/downsampling/NCBI_downsampling_imputed"
WGS_TRUTH_DIR="/QRISdata/Q8367/WGS_Reference_Panel/NCBI_truth_set/7.Consolidated_VCF"
REFERENCE_FASTA="/QRISdata/Q8367/Reference_Genome/GDDH13_1-1_formatted.fasta"

# Copy the two HortQG source files into this directory before running on Bunya.
OUTPUT_ROOT="/scratch/project/bigdata_apple/downsampling/HortQG_WGS_mask_sensitivity"
HORTQG_SOURCE_DIR="${OUTPUT_ROOT}/source"
HORTQG_VCF="${HORTQG_SOURCE_DIR}/HortQG_database__483724variants__48individuals.vcf"
HORTQG_METADATA="${HORTQG_SOURCE_DIR}/HortQG_database__48individuals_metadata.tsv"

DOSAGE_R2="${QUILT2_ROOT}/bin/dosage_r2_sbatch.sh"
CONCAT_IMPUTED="${QUILT2_ROOT}/modules/evaluate/concat_imputed.sh"
WGS_EVALUATOR="${QUILT2_ROOT}/modules/evaluate/dosage_r2_wgs.sh"
WGS_EVALUATOR_R="${QUILT2_ROOT}/modules/evaluate/dosage_r2_wgs.R"

# These hashes pin the local HortQG snapshot used to define the five masks.
EXPECTED_HORTQG_VCF_SHA256="eb16e8a4eabbeb94be2bed7cd13532101a7025d6ee1fec443f6d6029e65e0107"
EXPECTED_HORTQG_METADATA_SHA256="c91d8e70a997010a2743084b2e9909cf9fd6c394af6171cb3bad23c31278a8eb"

PILOT_RUN_ID="unbiased_filtered_liao"
PILOT_MASK_ID="rosbreed"
PILOT_CHROMOSOME="Chr01"
ONE_RUN_ID="unbiased_filtered_liao"

# Mask columns: short ID, report label, metadata dataset, expected observed SNPs,
# expected canonical Chr01-Chr17 positions, expected canonical-mask SHA-256.
mask_rows() {
    cat <<'EOF'
fruitbreedomics	FruitBreedomics	FruitBreedomics_001	322330	268652	ba1c939c2010c56875e987f15d1d44692f99d4a396a5afbb346e7ff99490a282
grove	Grove	Grove_Dataset_001	319203	268871	d46204ad5119f6f5ebe2bcba7ea7ac7329ad272e0a7350d3b19ae648839af22d
howard	Howard	Howard_Dataset_001	323926	268926	4a04bb42ddfd7960a119446819b21631f3b3fec7285638ca6f04356afb5ef6d0
migicovsky	Migicovsky	Migicovsky_2021_Dataset_001	438622	383523	89801993ef418f5f376bfb42792c369406cd6208da1b5be013fb258bc876ab52
rosbreed	RosBREED	RosBREED_Dataset_005	12136	9606	b3d67a67c110bce0ca7a9e4eb6d4c1b7c0f85c2d809b7a39075f048cf50a2142
EOF
}

# Run columns: ID, run type, treatment, panel, chunks directory.
run_rows() {
    cat <<EOF
biased_filtered_liao	Biased	Filtered	Liao	${IMPUTED_BASE}/filter/NCBI_WGS_Liao/chunks/imputed
biased_filtered_ncbi	Biased	Filtered	NCBI	${IMPUTED_BASE}/filter/NCBI_WGS_NCBI/chunks/imputed
biased_filtered_combined	Biased	Filtered	Combined	${IMPUTED_BASE}/filter/NCBI_WGS_Combined/chunks/imputed
biased_no_filter_liao	Biased	No_filter	Liao	${IMPUTED_BASE}/no_filter/NCBI_WGS_no_filter_Liao/chunks/imputed
biased_no_filter_ncbi	Biased	No_filter	NCBI	${IMPUTED_BASE}/no_filter/NCBI_WGS_no_filter_NCBI/chunks/imputed
biased_no_filter_combined	Biased	No_filter	Combined	${IMPUTED_BASE}/no_filter/NCBI_WGS_no_filter_Combined/chunks/imputed
unbiased_filtered_liao	Unbiased	Filtered	Liao	${IMPUTED_BASE}/filter/NCBI_WGS_Liao_holdout/chunks/imputed
unbiased_filtered_ncbi	Unbiased	Filtered	NCBI	${IMPUTED_BASE}/filter/NCBI_WGS_NCBI_holdout/chunks/imputed
unbiased_filtered_combined	Unbiased	Filtered	Combined	${IMPUTED_BASE}/filter/NCBI_WGS_Combined_holdout/chunks/imputed
unbiased_no_filter_liao	Unbiased	No_filter	Liao	${IMPUTED_BASE}/no_filter/NCBI_WGS_no_filter_Liao_holdout/chunks/imputed
unbiased_no_filter_ncbi	Unbiased	No_filter	NCBI	${IMPUTED_BASE}/no_filter/NCBI_WGS_no_filter_NCBI_holdout/chunks/imputed
unbiased_no_filter_combined	Unbiased	No_filter	Combined	${IMPUTED_BASE}/no_filter/NCBI_WGS_no_filter_Combined_holdout/chunks/imputed
EOF
}

# =============================================================================
# Generic helpers
# =============================================================================

log() { printf '[INFO] %s\n' "$*" >&2; }
warn() { printf '[WARN] %s\n' "$*" >&2; }
die() { printf '[ERROR] %s\n' "$*" >&2; exit 1; }

require_file() { [[ -f "$1" ]] || die "File not found: $1"; }
require_dir() { [[ -d "$1" ]] || die "Directory not found: $1"; }
require_command() { command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"; }

validate_boolean() {
    [[ "$2" == "true" || "$2" == "false" ]] || die "$1 must be true or false."
}

sha256_file() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | awk '{print $1}'
    else
        shasum -a 256 "$1" | awk '{print $1}'
    fi
}

sha256_text_file() {
    sha256_file "$1"
}

evaluator_bundle_sha256() {
    local payload="${TMP_WORK}/evaluator_bundle_sha256.tsv" path
    : > "${payload}"
    for path in "${DOSAGE_R2}" "${WGS_EVALUATOR}" "${WGS_EVALUATOR_R}" "${CONCAT_IMPUTED}"; do
        printf '%s\t%s\n' "$(basename "${path}")" "$(sha256_file "${path}")" >> "${payload}"
    done
    sha256_text_file "${payload}"
}

md5_file() {
    if command -v md5sum >/dev/null 2>&1; then
        md5sum "$1" | awk '{print $1}'
    else
        md5 -q "$1"
    fi
}

file_stamp() {
    local path="$1" metadata
    if metadata="$(stat -c '%s:%Y' "${path}" 2>/dev/null)"; then
        printf '%s:%s\n' "${path}" "${metadata}"
    else
        metadata="$(stat -f '%z:%m' "${path}")"
        printf '%s:%s\n' "${path}" "${metadata}"
    fi
}

canonical_chr() {
    local value="$1" number
    value="${value#chr}"
    value="${value#Chr}"
    [[ "${value}" =~ ^[0-9]+$ ]] || return 1
    number=$((10#${value}))
    (( number >= 0 && number <= 99 )) || return 1
    printf 'Chr%02d\n' "${number}"
}

make_rename_map() {
    local input="$1" expected="$2" output="$3" source canonical matches=0
    : > "${output}"
    while IFS= read -r source; do
        canonical="$(canonical_chr "${source}" 2>/dev/null || true)"
        [[ "${canonical}" == "${expected}" ]] || continue
        printf '%s\t%s\n' "${source}" "${expected}" >> "${output}"
        matches=$((matches + 1))
    done < <(bcftools view -h "${input}" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p')
    (( matches == 1 )) || die "Expected exactly one contig mapping to ${expected} in ${input}; found ${matches}."
}

install_generated_file() {
    local source="$1" destination="$2"
    mkdir -p "$(dirname "${destination}")"
    if [[ -f "${destination}" ]]; then
        if cmp -s "${source}" "${destination}"; then
            rm -f "${source}"
            return
        fi
        [[ "${FORCE}" == "true" ]] || die "Generated file differs from existing file; set FORCE=true to replace: ${destination}"
    fi
    mv -f "${source}" "${destination}"
}

print_command() {
    printf '[DRY-RUN]'
    printf ' %q' "$@"
    printf '\n'
}

is_selected_run() {
    local run_id="$1"
    case "${EXECUTION_SCOPE}" in
        pilot) [[ "${run_id}" == "${PILOT_RUN_ID}" ]] ;;
        one_run) [[ "${run_id}" == "${ONE_RUN_ID}" ]] ;;
        full) return 0 ;;
    esac
}

is_selected_mask() {
    local mask_id="$1"
    [[ "${EXECUTION_SCOPE}" != "pilot" || "${mask_id}" == "${PILOT_MASK_ID}" ]]
}

scope_suffix() {
    if [[ "${EXECUTION_SCOPE}" == "pilot" ]]; then
        printf 'pilot_%s\n' "${PILOT_CHROMOSOME}"
    else
        printf 'full\n'
    fi
}

TMP_WORK="$(mktemp -d "${TMPDIR:-/tmp}/hortqg_wgs_mask_sensitivity.XXXXXX")"
COORDINATOR_LOCK=""
REPAIRED_VCF_RESULT=""
CONCATENATED_RESULT=""
MASKED_IMPUTED_RESULT=""
STAGED_METRICS_RESULT=""
cleanup() {
    [[ -z "${COORDINATOR_LOCK}" ]] || rmdir "${COORDINATOR_LOCK}" 2>/dev/null || true
    rm -rf "${TMP_WORK}"
}
trap cleanup EXIT

acquire_coordinator_lock() {
    local lock_path="${OUTPUT_ROOT}/.wgs_mask_sensitivity_coordinator.lock"
    if ! mkdir "${lock_path}" 2>/dev/null; then
        die "Another live utility coordinator may be using ${OUTPUT_ROOT}. After confirming no coordinator is running, remove the stale empty lock directory: ${lock_path}"
    fi
    COORDINATOR_LOCK="${lock_path}"
}

# =============================================================================
# Preflight and source repair
# =============================================================================

preflight() {
    (( $# == 0 )) || die "This fixed utility accepts no command-line arguments; edit its configuration block instead."
    [[ "$(uname -s)" == "Linux" ]] || die "This utility is Bunya/Linux-only and must not be executed on local macOS."
    validate_boolean DRY_RUN "${DRY_RUN}"
    validate_boolean FORCE "${FORCE}"
    validate_boolean RESUBMIT_INCOMPLETE "${RESUBMIT_INCOMPLETE}"
    case "${EXECUTION_SCOPE}" in
        pilot|one_run|full) ;;
        *) die "EXECUTION_SCOPE must be pilot, one_run, or full." ;;
    esac
    if [[ "${FORCE}" == "true" && "${RESUBMIT_INCOMPLETE}" == "true" ]]; then
        warn "FORCE and RESUBMIT_INCOMPLETE are both enabled; confirm prior Slurm jobs were cancelled before running."
    fi
    [[ -n "${SLURM_JOB_ID:-}" ]] || {
        die "Submit this utility with sbatch; mask construction and VCF preprocessing must not run on a Bunya login node."
    }

    if ! command -v bcftools >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
        module load "${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}" || die "Could not load the Bunya bcftools module."
    fi
    for command in awk bcftools bgzip cmp env grep sed sort; do
        require_command "${command}"
    done
    if ! command -v sha256sum >/dev/null 2>&1; then
        require_command shasum
    fi
    if ! command -v md5sum >/dev/null 2>&1; then
        require_command md5
    fi
    [[ "${DRY_RUN}" == "true" ]] || require_command sbatch

    require_file "${HORTQG_VCF}"
    require_file "${HORTQG_METADATA}"
    require_file "${REFERENCE_FASTA}"
    require_file "${REFERENCE_FASTA}.fai"
    require_file "${DOSAGE_R2}"
    require_file "${CONCAT_IMPUTED}"
    require_file "${WGS_EVALUATOR}"
    require_file "${WGS_EVALUATOR_R}"
    require_dir "${WGS_TRUTH_DIR}"

    while IFS=$'\t' read -r _ _ _ _ chunks_dir; do
        require_dir "${chunks_dir}"
    done < <(run_rows)

    local chr truth_vcf
    for number in $(seq 1 17); do
        printf -v chr 'Chr%02d' "${number}"
        truth_vcf="${WGS_TRUTH_DIR}/${chr}_consolidated.vcf.gz"
        require_file "${truth_vcf}"
        [[ -f "${truth_vcf}.csi" || -f "${truth_vcf}.tbi" ]] || die "Truth VCF is not indexed: ${truth_vcf}"
    done

    local observed_vcf_hash observed_metadata_hash
    observed_vcf_hash="$(sha256_file "${HORTQG_VCF}")"
    observed_metadata_hash="$(sha256_file "${HORTQG_METADATA}")"
    [[ "${observed_vcf_hash}" == "${EXPECTED_HORTQG_VCF_SHA256}" ]] || {
        die "HortQG VCF SHA-256 differs from the fixed source snapshot: ${observed_vcf_hash}"
    }
    [[ "${observed_metadata_hash}" == "${EXPECTED_HORTQG_METADATA_SHA256}" ]] || {
        die "HortQG metadata SHA-256 differs from the fixed source snapshot: ${observed_metadata_hash}"
    }
}

repair_hortqg_vcf() {
    local cache_root="$1" repaired_dir="${cache_root}/source" repaired_vcf
    local temporary_vcf="${TMP_WORK}/hortqg.repaired.vcf.gz"
    repaired_vcf="${repaired_dir}/HortQG_database__483724variants__48individuals.repaired.vcf.gz"
    mkdir -p "${repaired_dir}"

    if [[ -f "${repaired_vcf}" && ( -f "${repaired_vcf}.csi" || -f "${repaired_vcf}.tbi" ) ]]; then
        REPAIRED_VCF_RESULT="${repaired_vcf}"
        return
    fi

    awk '
      /^##FORMAT=<ID=GT,/ { gt++; if (gt > 1) next }
      /^##INFO=<ID=PR,/ { pr++; if (pr > 1) next }
      { print }
    ' "${HORTQG_VCF}" | bgzip -c > "${temporary_vcf}"
    bcftools index -f -c "${temporary_vcf}"

    [[ "$(bcftools view -h "${temporary_vcf}" | grep -c '^##FORMAT=<ID=GT,')" == "1" ]] || die "Repaired HortQG VCF does not contain exactly one FORMAT/GT header."
    [[ "$(bcftools view -h "${temporary_vcf}" | grep -c '^##INFO=<ID=PR,')" == "1" ]] || die "Repaired HortQG VCF does not contain exactly one INFO/PR header."

    mv -f "${temporary_vcf}" "${repaired_vcf}"
    mv -f "${temporary_vcf}.csi" "${repaired_vcf}.csi"
    REPAIRED_VCF_RESULT="${repaired_vcf}"
}

# =============================================================================
# Fixed mask construction
# =============================================================================

build_canonical_masks() {
    local cache_root="$1" repaired_vcf="$2" masks_dir="${cache_root}/masks"
    local base_manifest="${TMP_WORK}/mask_manifest.base.tsv"
    local all_samples="${TMP_WORK}/hortqg.samples.txt"
    mkdir -p "${masks_dir}"
    bcftools query -l "${repaired_vcf}" > "${all_samples}"

    printf 'scope_id\tscope_label\tdataset_name\tsource_snp_count\tcanonical_position_count\tcanonical_mask_sha256\tsample_count\tsamples_sha256\tsamples_file\tcanonical_mask\tsource_vcf_sha256\tmetadata_sha256\tsource_exclusion_count\tsource_exclusions_file\n' > "${base_manifest}"

    local mask_id scope_label dataset_name expected_source expected_canonical expected_hash
    while IFS=$'\t' read -r mask_id scope_label dataset_name expected_source expected_canonical expected_hash; do
        local samples_tmp="${TMP_WORK}/${mask_id}.samples.txt"
        local observed_records="${TMP_WORK}/${mask_id}.observed.tsv"
        local canonical_raw="${TMP_WORK}/${mask_id}.canonical.raw.tsv"
        local canonical_sorted="${TMP_WORK}/${mask_id}.canonical.positions.tsv"
        local source_exclusions_tmp="${TMP_WORK}/${mask_id}.source_exclusions.tmp.tsv"
        local source_exclusions_sorted="${TMP_WORK}/${mask_id}.source_exclusions.tsv"
        local sample_csv sample_count samples_hash source_count canonical_raw_count canonical_count canonical_hash
        local observed_count source_exclusion_count
        local samples_path="${masks_dir}/${mask_id}.samples.txt"
        local canonical_path="${masks_dir}/${mask_id}.canonical.positions.tsv"
        local source_exclusions_path="${masks_dir}/${mask_id}.source_exclusions.tsv"

        awk -F '\t' -v wanted="${dataset_name}" '
          NR == 1 {
            for (i = 1; i <= NF; i++) {
              if ($i == "individual") sample_col = i
              if ($i == "dataset_name") dataset_col = i
            }
            if (!sample_col || !dataset_col) exit 2
            next
          }
          {
            n = split($dataset_col, labels, /[[:space:]]*;[[:space:]]*/)
            for (i = 1; i <= n; i++) {
              if (labels[i] == wanted) {
                print $sample_col
                break
              }
            }
          }
        ' "${HORTQG_METADATA}" | sort -u > "${samples_tmp}"
        [[ -s "${samples_tmp}" ]] || die "No metadata samples found for ${dataset_name}."
        while IFS= read -r sample; do
            grep -Fqx "${sample}" "${all_samples}" || die "Metadata sample is absent from HortQG VCF: ${sample}"
        done < "${samples_tmp}"
        sample_csv="$(paste -sd, "${samples_tmp}")"
        sample_count="$(wc -l < "${samples_tmp}" | tr -d ' ')"
        samples_hash="$(sha256_file "${samples_tmp}")"

        bcftools view -s "${sample_csv}" -Ou "${repaired_vcf}" \
            | bcftools view -i 'N_PASS(GT!="mis")>0' -Ou \
            | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' > "${observed_records}"

        source_count="$(awk -F '\t' '
          $3 ~ /^[ACGTN]$/ && $4 ~ /^[ACGTN]$/ && $4 != "." { n++ }
          END { print n + 0 }
        ' "${observed_records}")"
        [[ "${source_count}" == "${expected_source}" ]] || {
            die "${mask_id} source count drift: expected ${expected_source}, found ${source_count}."
        }

        awk -F '\t' '
          $1 ~ /^[0-9]+$/ && $1 + 0 >= 1 && $1 + 0 <= 17 &&
          $3 ~ /^[ACGTN]$/ && $4 ~ /^[ACGTN]$/ && $4 != "." {
            printf "Chr%02d\t%s\n", $1 + 0, $2
          }
        ' "${observed_records}" > "${canonical_raw}"
        awk -F '\t' 'BEGIN { OFS = FS }
          $4 == "." {
            print $1, $2, $3, $4, "missing_alt"
            next
          }
          $3 !~ /^[ACGTN]$/ || $4 !~ /^[ACGTN]$/ {
            print $1, $2, $3, $4, "non_snp_or_non_biallelic_source_record"
            next
          }
          $1 !~ /^[0-9]+$/ || $1 + 0 < 1 || $1 + 0 > 17 {
            print $1, $2, $3, $4, "noncanonical_source_chromosome"
          }
        ' "${observed_records}" > "${source_exclusions_tmp}"
        {
            printf 'CHROM\tPOS\tREF\tALT\treason\n'
            sort -k1,1V -k2,2n -k5,5 "${source_exclusions_tmp}"
        } > "${source_exclusions_sorted}"
        canonical_raw_count="$(wc -l < "${canonical_raw}" | tr -d ' ')"
        observed_count="$(wc -l < "${observed_records}" | tr -d ' ')"
        source_exclusion_count="$(wc -l < "${source_exclusions_tmp}" | tr -d ' ')"
        sort -k1,1 -k2,2n -u "${canonical_raw}" > "${canonical_sorted}"
        canonical_count="$(wc -l < "${canonical_sorted}" | tr -d ' ')"
        (( canonical_raw_count + source_exclusion_count == observed_count )) || {
            die "Source exclusion accounting does not reconcile for ${mask_id}."
        }
        [[ "${canonical_raw_count}" == "${canonical_count}" ]] || die "${mask_id} contains duplicate canonical CHROM:POS rows."
        [[ "${canonical_count}" == "${expected_canonical}" ]] || {
            die "${mask_id} canonical count drift: expected ${expected_canonical}, found ${canonical_count}."
        }
        canonical_hash="$(sha256_text_file "${canonical_sorted}")"
        [[ "${canonical_hash}" == "${expected_hash}" ]] || {
            die "${mask_id} canonical mask SHA-256 drift: expected ${expected_hash}, found ${canonical_hash}."
        }

        install_generated_file "${samples_tmp}" "${samples_path}"
        install_generated_file "${canonical_sorted}" "${canonical_path}"
        install_generated_file "${source_exclusions_sorted}" "${source_exclusions_path}"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${mask_id}" "${scope_label}" "${dataset_name}" "${source_count}" \
            "${canonical_count}" "${canonical_hash}" "${sample_count}" "${samples_hash}" \
            "${samples_path}" "${canonical_path}" "${EXPECTED_HORTQG_VCF_SHA256}" \
            "${EXPECTED_HORTQG_METADATA_SHA256}" "${source_exclusion_count}" \
            "${source_exclusions_path}" >> "${base_manifest}"
        log "Validated ${scope_label}: ${source_count} observed SNPs; ${canonical_count} canonical Chr01-Chr17 positions; ${source_exclusion_count} source records excluded."
    done < <(mask_rows)

    install_generated_file "${base_manifest}" "${masks_dir}/mask_manifest.base.tsv"
}

# =============================================================================
# WGS-eligible masks and masked truth copies
# =============================================================================

truth_source_signature() {
    local output="$1" chr truth_vcf
    : > "${output}"
    printf 'truth-normalization-v2-post-normalization-duplicates\n' >> "${output}"
    bcftools --version | sed -n '1p' >> "${output}"
    file_stamp "${REFERENCE_FASTA}" >> "${output}"
    file_stamp "${REFERENCE_FASTA}.fai" >> "${output}"
    for number in $(seq 1 17); do
        printf -v chr 'Chr%02d' "${number}"
        truth_vcf="${WGS_TRUTH_DIR}/${chr}_consolidated.vcf.gz"
        file_stamp "${truth_vcf}" >> "${output}"
        [[ ! -f "${truth_vcf}.csi" ]] || file_stamp "${truth_vcf}.csi" >> "${output}"
        [[ ! -f "${truth_vcf}.tbi" ]] || file_stamp "${truth_vcf}.tbi" >> "${output}"
    done
}

build_masked_truth() {
    local cache_root="$1" masks_dir="${cache_root}/masks" truth_root="${cache_root}/truth"
    local base_manifest="${masks_dir}/mask_manifest.base.tsv"
    local final_manifest="${TMP_WORK}/mask_manifest.tsv"
    local truth_signature_file="${TMP_WORK}/truth_source_signature.txt"
    local truth_signature
    truth_source_signature "${truth_signature_file}"
    truth_signature="$(sha256_text_file "${truth_signature_file}")"

    printf 'scope_id\tscope_label\tdataset_name\tsource_snp_count\tcanonical_position_count\twgs_eligible_position_count\tmissing_from_wgs_truth_count\tmissing_alt_wgs_truth_count\tnon_biallelic_wgs_truth_count\tnormalization_excluded_wgs_truth_count\tduplicate_wgs_truth_position_count\tcanonical_mask_sha256\tmask_sha256\tsample_count\tsamples_sha256\tsamples_file\tcanonical_mask\tevaluation_mask\texclusions_file\twgs_truth_record_exclusion_count\twgs_truth_record_exclusions_file\twgs_truth_record_exclusions_sha256\ttruth_dir\tsource_vcf_sha256\tmetadata_sha256\ttruth_source_signature\tcoordinate_mapping_assumption\tsource_exclusion_count\tsource_exclusions_file\tmasked_truth_sha256_manifest\tmasked_truth_content_sha256\n' > "${final_manifest}"

    local mask_id scope_label dataset_name source_count canonical_count canonical_hash
    local sample_count samples_hash samples_path canonical_path _source_vcf_hash _metadata_hash
    local source_exclusion_count source_exclusions_file
    while IFS=$'\t' read -r mask_id scope_label dataset_name source_count canonical_count canonical_hash sample_count samples_hash samples_path canonical_path _source_vcf_hash _metadata_hash source_exclusion_count source_exclusions_file; do
        [[ "${mask_id}" != "scope_id" ]] || continue
        local raw_truth="${TMP_WORK}/${mask_id}.truth.raw.tsv"
        local normalized_truth="${TMP_WORK}/${mask_id}.truth.normalized.tsv"
        local truth_events="${TMP_WORK}/${mask_id}.truth.events.tsv"
        local eligible_tmp="${TMP_WORK}/${mask_id}.eligible.tmp.tsv"
        local eligible_sorted="${TMP_WORK}/${mask_id}.positions.tsv"
        local exclusions_tmp="${TMP_WORK}/${mask_id}.exclusions.tmp.tsv"
        local exclusions_sorted="${TMP_WORK}/${mask_id}.exclusions.tsv"
        local record_exclusion_events="${TMP_WORK}/${mask_id}.truth_record_exclusion.events.tsv"
        local record_exclusions_tmp="${TMP_WORK}/${mask_id}.truth_record_exclusions.tmp.tsv"
        local record_exclusions_sorted="${TMP_WORK}/${mask_id}.truth_record_exclusions.tsv"
        local evaluation_mask="${masks_dir}/${mask_id}.positions.tsv"
        local exclusions_file="${masks_dir}/${mask_id}.exclusions.tsv"
        local record_exclusions_file="${masks_dir}/${mask_id}.truth_record_exclusions.tsv"
        local truth_dir="${truth_root}/${mask_id}"
        local truth_checksum_manifest="${truth_dir}/masked_truth_sha256.tsv"
        local chr truth_vcf rename_map
        : > "${raw_truth}"
        : > "${normalized_truth}"

        for number in $(seq 1 17); do
            printf -v chr 'Chr%02d' "${number}"
            truth_vcf="${WGS_TRUTH_DIR}/${chr}_consolidated.vcf.gz"
            rename_map="${TMP_WORK}/${mask_id}.${chr}.truth.rename.tsv"
            local normalized_candidate="${TMP_WORK}/${mask_id}.${chr}.truth.normalized.vcf.gz"
            make_rename_map "${truth_vcf}" "${chr}" "${rename_map}"
            bcftools annotate --rename-chrs "${rename_map}" -Ou "${truth_vcf}" \
                | bcftools view -T "${canonical_path}" -Ou \
                | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' >> "${raw_truth}"
            bcftools annotate --rename-chrs "${rename_map}" -Ou "${truth_vcf}" \
                | bcftools view -T "${canonical_path}" -m2 -M2 -v snps -Ou \
                | bcftools norm -f "${REFERENCE_FASTA}" -Ou \
                | bcftools norm -d exact -Oz -o "${normalized_candidate}"
            bcftools index -f -c "${normalized_candidate}"
            bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "${normalized_candidate}" >> "${normalized_truth}"
        done

        awk 'BEGIN { OFS = "\t" } { print "raw", $0 }' "${raw_truth}" > "${truth_events}"
        awk 'BEGIN { OFS = "\t" } { print "normalized", $0 }' "${normalized_truth}" >> "${truth_events}"
        awk -F '\t' -v eligible="${eligible_tmp}" -v excluded="${exclusions_tmp}" '
          NR == FNR { mask[$1 SUBSEP $2] = $1 FS $2; next }
          {
            key = $2 SUBSEP $3
            if ($1 == "raw") {
              raw_total[key]++
              if ($5 == ".") raw_missing_alt[key]++
              if ($4 ~ /^[ACGTN]$/ && $5 ~ /^[ACGTN]$/ && $5 != ".") raw_valid[key]++
            } else if ($1 == "normalized") {
              normalized_count[key]++
            }
          }
          END {
            for (key in mask) {
              if (normalized_count[key] == 1) print mask[key] > eligible
              else if (normalized_count[key] > 1) print mask[key], "duplicate_wgs_truth_position" > excluded
              else if (!raw_total[key]) print mask[key], "absent_from_wgs_truth" > excluded
              else if (raw_missing_alt[key] == raw_total[key]) print mask[key], "missing_alt_wgs_truth" > excluded
              else if (!raw_valid[key]) print mask[key], "non_biallelic_wgs_truth" > excluded
              else print mask[key], "excluded_during_wgs_truth_normalization" > excluded
            }
          }
        ' "${canonical_path}" "${truth_events}"
        sort -k1,1 -k2,2n -u "${eligible_tmp}" > "${eligible_sorted}"
        {
            printf 'CHROM\tPOS\treason\n'
            sort -k1,1 -k2,2n -k3,3 "${exclusions_tmp}"
        } > "${exclusions_sorted}"

        awk 'BEGIN { OFS = "\t" } { print "normalized", $0 }' "${normalized_truth}" > "${record_exclusion_events}"
        awk 'BEGIN { OFS = "\t" } NR > 1 { print "position", $1, $2, $3 }' "${exclusions_sorted}" >> "${record_exclusion_events}"
        awk 'BEGIN { OFS = "\t" } { print "raw", $0 }' "${raw_truth}" >> "${record_exclusion_events}"
        awk -F '\t' 'BEGIN { OFS = FS }
          $1 == "normalized" {
            normalized[$2 SUBSEP $3 SUBSEP $4 SUBSEP $5]++
            next
          }
          $1 == "position" {
            position_reason[$2 SUBSEP $3] = $4
            next
          }
          $1 == "raw" {
            position_key = $2 SUBSEP $3
            exact_key = position_key SUBSEP $4 SUBSEP $5
            reason = ""
            if ($5 == ".") reason = "missing_alt_wgs_truth_record"
            else if ($4 !~ /^[ACGTN]$/ || $5 !~ /^[ACGTN]$/) reason = "non_snp_or_non_biallelic_wgs_truth_record"
            else if (position_reason[position_key] == "duplicate_wgs_truth_position") reason = "duplicate_wgs_truth_position"
            else if (position_reason[position_key] == "excluded_during_wgs_truth_normalization") reason = "excluded_during_wgs_truth_normalization"
            else {
              retained_seen[exact_key]++
              if (retained_seen[exact_key] > normalized[exact_key]) reason = "collapsed_exact_duplicate_wgs_truth_record"
            }
            if (reason != "") print $2, $3, $4, $5, reason
          }
        ' "${record_exclusion_events}" > "${record_exclusions_tmp}"
        {
            printf 'CHROM\tPOS\tREF\tALT\treason\n'
            sort -k1,1 -k2,2n -k3,3 -k4,4 -k5,5 "${record_exclusions_tmp}"
        } > "${record_exclusions_sorted}"

        local eligible_count missing_count missing_alt_count non_biallelic_count normalization_excluded_count
        local duplicate_count evaluation_hash
        eligible_count="$(wc -l < "${eligible_sorted}" | tr -d ' ')"
        missing_count="$(awk -F '\t' '$3=="absent_from_wgs_truth"{n++} END{print n+0}' "${exclusions_sorted}")"
        missing_alt_count="$(awk -F '\t' '$3=="missing_alt_wgs_truth"{n++} END{print n+0}' "${exclusions_sorted}")"
        non_biallelic_count="$(awk -F '\t' '$3=="non_biallelic_wgs_truth"{n++} END{print n+0}' "${exclusions_sorted}")"
        normalization_excluded_count="$(awk -F '\t' '$3=="excluded_during_wgs_truth_normalization"{n++} END{print n+0}' "${exclusions_sorted}")"
        duplicate_count="$(awk -F '\t' '$3=="duplicate_wgs_truth_position"{n++} END{print n+0}' "${exclusions_sorted}")"
        (( eligible_count + missing_count + missing_alt_count + non_biallelic_count + normalization_excluded_count + duplicate_count == canonical_count )) || {
            die "WGS eligibility accounting does not reconcile for ${mask_id}."
        }
        evaluation_hash="$(sha256_text_file "${eligible_sorted}")"
        install_generated_file "${eligible_sorted}" "${evaluation_mask}"
        install_generated_file "${exclusions_sorted}" "${exclusions_file}"
        install_generated_file "${record_exclusions_sorted}" "${record_exclusions_file}"
        local record_exclusion_count record_exclusions_hash
        record_exclusion_count="$(awk 'NR>1{n++} END{print n+0}' "${record_exclusions_file}")"
        record_exclusions_hash="$(sha256_file "${record_exclusions_file}")"

        local signature_payload="${TMP_WORK}/${mask_id}.truth.signature.payload"
        local truth_cache_signature marker observed_marker=""
        printf '%s\n%s\n%s\n' "truth-cache-v2-post-normalization-duplicates" \
            "${evaluation_hash}" "${truth_signature}" > "${signature_payload}"
        truth_cache_signature="$(sha256_text_file "${signature_payload}")"
        marker="${truth_dir}/.complete_signature"
        [[ ! -f "${marker}" ]] || observed_marker="$(<"${marker}")"
        if [[ -n "${observed_marker}" && "${observed_marker}" != "${truth_cache_signature}" && "${FORCE}" != "true" ]]; then
            die "Masked-truth cache is stale for ${mask_id}; set FORCE=true to rebuild: ${truth_dir}"
        fi

        local cache_complete=true chr_count output_vcf output_count expected_truth_files=0
        [[ -s "${truth_checksum_manifest}" ]] || cache_complete=false
        if [[ -s "${truth_checksum_manifest}" && "$(head -n 1 "${truth_checksum_manifest}")" != $'file\tsha256\trecord_count' ]]; then
            cache_complete=false
        fi
        for number in $(seq 1 17); do
            printf -v chr 'Chr%02d' "${number}"
            chr_count="$(awk -F '\t' -v wanted="${chr}" '$1==wanted{n++} END{print n+0}' "${evaluation_mask}")"
            output_vcf="${truth_dir}/${chr}_consolidated.vcf.gz"
            if (( chr_count == 0 )); then
                [[ ! -e "${output_vcf}" && ! -e "${output_vcf}.csi" && ! -e "${output_vcf}.tbi" ]] || cache_complete=false
            elif [[ ! -s "${output_vcf}" || !( -f "${output_vcf}.csi" || -f "${output_vcf}.tbi" ) ]]; then
                cache_complete=false
            else
                expected_truth_files=$((expected_truth_files + 1))
                output_count="$(bcftools view -H "${output_vcf}" | wc -l | tr -d ' ')"
                [[ "${output_count}" == "${chr_count}" ]] || cache_complete=false
                if [[ -s "${truth_checksum_manifest}" ]]; then
                    local checksum_name checksum_hash checksum_count
                    checksum_name="$(basename "${output_vcf}")"
                    checksum_hash="$(awk -F '\t' -v wanted="${checksum_name}" '$1==wanted{print $2}' "${truth_checksum_manifest}")"
                    checksum_count="$(awk -F '\t' -v wanted="${checksum_name}" '$1==wanted{print $3}' "${truth_checksum_manifest}")"
                    [[ "${checksum_hash}" == "$(sha256_file "${output_vcf}")" && "${checksum_count}" == "${output_count}" ]] || cache_complete=false
                fi
            fi
        done
        if [[ -s "${truth_checksum_manifest}" ]]; then
            [[ "$(awk 'NR>1{n++} END{print n+0}' "${truth_checksum_manifest}")" == "${expected_truth_files}" ]] || cache_complete=false
        fi
        if [[ "${observed_marker}" == "${truth_cache_signature}" && "${cache_complete}" != "true" && "${FORCE}" != "true" ]]; then
            die "Masked-truth cache is incomplete for ${mask_id}; set FORCE=true to rebuild: ${truth_dir}"
        fi

        if [[ "${observed_marker}" != "${truth_cache_signature}" || "${cache_complete}" != "true" || "${FORCE}" == "true" ]]; then
            mkdir -p "${truth_dir}"
            for number in $(seq 1 17); do
                printf -v chr 'Chr%02d' "${number}"
                local output_tmp normalized_candidate
                chr_count="$(awk -F '\t' -v wanted="${chr}" '$1==wanted{n++} END{print n+0}' "${evaluation_mask}")"
                output_vcf="${truth_dir}/${chr}_consolidated.vcf.gz"
                rm -f "${output_vcf}" "${output_vcf}.csi" "${output_vcf}.tbi"
                if (( chr_count == 0 )); then
                    continue
                fi
                normalized_candidate="${TMP_WORK}/${mask_id}.${chr}.truth.normalized.vcf.gz"
                output_tmp="${TMP_WORK}/${mask_id}.${chr}.truth.vcf.gz"
                bcftools view -T "${evaluation_mask}" -Oz -o "${output_tmp}" "${normalized_candidate}"
                bcftools index -f -c "${output_tmp}"
                output_count="$(bcftools view -H "${output_tmp}" | wc -l | tr -d ' ')"
                [[ "${output_count}" == "${chr_count}" ]] || {
                    die "Masked truth count changed during normalization for ${mask_id}/${chr}: expected ${chr_count}, found ${output_count}."
                }
                mv -f "${output_tmp}" "${output_vcf}"
                mv -f "${output_tmp}.csi" "${output_vcf}.csi"
            done
            local checksum_tmp="${TMP_WORK}/${mask_id}.masked_truth_sha256.tsv"
            printf 'file\tsha256\trecord_count\n' > "${checksum_tmp}"
            for number in $(seq 1 17); do
                printf -v chr 'Chr%02d' "${number}"
                output_vcf="${truth_dir}/${chr}_consolidated.vcf.gz"
                [[ -f "${output_vcf}" ]] || continue
                output_count="$(bcftools view -H "${output_vcf}" | wc -l | tr -d ' ')"
                printf '%s\t%s\t%s\n' "$(basename "${output_vcf}")" \
                    "$(sha256_file "${output_vcf}")" "${output_count}" >> "${checksum_tmp}"
            done
            mv -f "${checksum_tmp}" "${truth_checksum_manifest}"
            printf '%s\n' "${truth_cache_signature}" > "${marker}.tmp.$$"
            mv -f "${marker}.tmp.$$" "${marker}"
        fi

        require_file "${truth_checksum_manifest}"
        local truth_content_hash
        truth_content_hash="$(sha256_file "${truth_checksum_manifest}")"

        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "${mask_id}" "${scope_label}" "${dataset_name}" "${source_count}" \
            "${canonical_count}" "${eligible_count}" "${missing_count}" "${missing_alt_count}" "${non_biallelic_count}" \
            "${normalization_excluded_count}" "${duplicate_count}" "${canonical_hash}" "${evaluation_hash}" "${sample_count}" "${samples_hash}" \
            "${samples_path}" "${canonical_path}" "${evaluation_mask}" "${exclusions_file}" \
            "${record_exclusion_count}" "${record_exclusions_file}" "${record_exclusions_hash}" \
            "${truth_dir}" "${EXPECTED_HORTQG_VCF_SHA256}" "${EXPECTED_HORTQG_METADATA_SHA256}" \
            "${truth_signature}" "numeric_1-17_to_Chr01-Chr17_no_declared_source_assembly" \
            "${source_exclusion_count}" "${source_exclusions_file}" \
            "${truth_checksum_manifest}" "${truth_content_hash}" >> "${final_manifest}"
        log "${scope_label}: ${eligible_count} WGS-eligible positions; excluded ${missing_count} absent, ${missing_alt_count} ALT-missing, ${non_biallelic_count} non-biallelic, ${normalization_excluded_count} during normalization, and ${duplicate_count} post-normalization duplicate-position sites."
    done < "${base_manifest}"

    install_generated_file "${final_manifest}" "${masks_dir}/mask_manifest.tsv"
}

# =============================================================================
# Masked imputed inputs and exact-match preflight
# =============================================================================

prepare_concatenated_run() {
    local cache_root="$1" run_id="$2" chunks_dir="$3" baseline_signature="$4"
    local concat_dir="${cache_root}/concat/${run_id}" concat_args
    local output_vcf="${concat_dir}/imputed.all_chroms.vcf.gz"
    local marker="${concat_dir}/.source_run_signature" observed_marker=""
    [[ ! -f "${marker}" ]] || observed_marker="$(<"${marker}")"
    if [[ "${observed_marker}" == "${baseline_signature}" && -s "${output_vcf}" && \
          ( -f "${output_vcf}.csi" || -f "${output_vcf}.tbi" ) && "${FORCE}" != "true" ]]; then
        bcftools view -h "${output_vcf}" >/dev/null || die "Concatenated VCF cannot be opened; set FORCE=true to rebuild: ${output_vcf}"
        bcftools index -n "${output_vcf}" >/dev/null || die "Concatenated VCF index is invalid; set FORCE=true to rebuild: ${output_vcf}"
        CONCATENATED_RESULT="${output_vcf}"
        return
    fi
    if [[ "${FORCE}" != "true" && ( -e "${output_vcf}" || -e "${output_vcf}.csi" || \
          -e "${output_vcf}.tbi" || -e "${marker}" ) ]]; then
        die "Concatenated utility cache is stale, incomplete, or unversioned for ${run_id}; set FORCE=true to rebuild: ${concat_dir}"
    fi
    rm -f "${marker}"
    concat_args=(--chunks-dir "${chunks_dir}" --out-dir "${concat_dir}" --force)
    local produced produced_log="${TMP_WORK}/${run_id}.concat.stdout"
    bash "${CONCAT_IMPUTED}" "${concat_args[@]}" > "${produced_log}"
    produced="$(tail -n 1 "${produced_log}")"
    [[ "${produced}" == "${output_vcf}" ]] || die "Unexpected concatenated VCF path for ${run_id}: ${produced}"
    [[ -s "${produced}" ]] || die "Concatenation did not create a nonempty VCF for ${run_id}: ${produced}"
    [[ -f "${produced}.csi" || -f "${produced}.tbi" ]] || die "Concatenation did not create an index for ${run_id}: ${produced}"
    bcftools view -h "${produced}" >/dev/null || die "Concatenated VCF cannot be opened for ${run_id}: ${produced}"
    bcftools index -n "${produced}" >/dev/null || die "Concatenated VCF index is invalid for ${run_id}: ${produced}"
    mkdir -p "${concat_dir}"
    printf '%s\n' "${baseline_signature}" > "${marker}.tmp.$$"
    mv -f "${marker}.tmp.$$" "${marker}"
    CONCATENATED_RESULT="${produced}"
}

prepare_masked_imputed() {
    local cache_root="$1" run_id="$2" concatenated="$3" mask_id="$4" evaluation_mask="$5"
    local mask_hash="$6" baseline_signature="$7" output_dir="${cache_root}/imputed/${run_id}/${mask_id}"
    local output_vcf="${output_dir}/imputed.masked.vcf.gz" marker="${output_dir}/.complete_signature"
    local payload="${TMP_WORK}/${run_id}.${mask_id}.imputed.signature.payload" signature existing=""
    printf '%s\n' "masked-imputed-cache-v2" > "${payload}"
    file_stamp "${concatenated}" >> "${payload}"
    file_stamp "${REFERENCE_FASTA}" >> "${payload}"
    file_stamp "${REFERENCE_FASTA}.fai" >> "${payload}"
    printf '%s\n%s\n' "${mask_hash}" "${baseline_signature}" >> "${payload}"
    bcftools --version | sed -n '1p' >> "${payload}"
    signature="$(sha256_text_file "${payload}")"
    [[ ! -f "${marker}" ]] || existing="$(<"${marker}")"

    if [[ "${existing}" == "${signature}" && -s "${output_vcf}" && ( -f "${output_vcf}.csi" || -f "${output_vcf}.tbi" ) && "${FORCE}" != "true" ]]; then
        bcftools view -h "${output_vcf}" >/dev/null || die "Masked-imputed VCF cannot be opened; set FORCE=true to rebuild: ${output_vcf}"
        bcftools index -n "${output_vcf}" >/dev/null || die "Masked-imputed VCF index is invalid; set FORCE=true to rebuild: ${output_vcf}"
        MASKED_IMPUTED_RESULT="${output_vcf}"
        return
    fi
    if [[ "${FORCE}" != "true" && ( -e "${output_vcf}" || -e "${output_vcf}.csi" || \
          -e "${output_vcf}.tbi" || -e "${marker}" ) ]]; then
        die "Masked-imputed cache is stale, incomplete, or unversioned; set FORCE=true to rebuild: ${output_dir}"
    fi

    mkdir -p "${output_dir}"
    rm -f "${marker}"
    local rename_map="${TMP_WORK}/${run_id}.imputed.rename.tsv" source canonical
    : > "${rename_map}"
    while IFS= read -r source; do
        canonical="$(canonical_chr "${source}" 2>/dev/null || true)"
        [[ "${canonical}" =~ ^Chr(0[1-9]|1[0-7])$ ]] || continue
        printf '%s\t%s\n' "${source}" "${canonical}" >> "${rename_map}"
    done < <(bcftools view -h "${concatenated}" | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p')
    [[ "$(cut -f2 "${rename_map}" | sort -u | wc -l | tr -d ' ')" == "$(wc -l < "${rename_map}" | tr -d ' ')" ]] || {
        die "Multiple imputed contigs map to one canonical chromosome: ${concatenated}"
    }

    local output_tmp="${TMP_WORK}/${run_id}.${mask_id}.imputed.vcf.gz"
    bcftools annotate --rename-chrs "${rename_map}" -Ou "${concatenated}" \
        | bcftools view -T "${evaluation_mask}" -m2 -M2 -v snps -Ou \
        | bcftools norm -f "${REFERENCE_FASTA}" -Ou \
        | bcftools norm -d exact -Oz -o "${output_tmp}"
    bcftools index -f -c "${output_tmp}"
    rm -f "${output_vcf}.csi" "${output_vcf}.tbi"
    mv -f "${output_tmp}" "${output_vcf}"
    mv -f "${output_tmp}.csi" "${output_vcf}.csi"
    bcftools view -h "${output_vcf}" >/dev/null || die "New masked-imputed VCF cannot be opened: ${output_vcf}"
    bcftools index -n "${output_vcf}" >/dev/null || die "New masked-imputed VCF index is invalid: ${output_vcf}"
    printf '%s\n' "${signature}" > "${marker}.tmp.$$"
    mv -f "${marker}.tmp.$$" "${marker}"
    MASKED_IMPUTED_RESULT="${output_vcf}"
}

preflight_exact_matches() {
    local masked_imputed="$1" truth_dir="$2" output="$3" chromosome_limit="${4:-}"
    local evaluation_mask="$5" chr truth_vcf count mask_count total=0
    printf 'chromosome\texact_match_variant_count\n' > "${output}"
    for number in $(seq 1 17); do
        printf -v chr 'Chr%02d' "${number}"
        [[ -z "${chromosome_limit}" || "${chr}" == "${chromosome_limit}" ]] || continue
        truth_vcf="${truth_dir}/${chr}_consolidated.vcf.gz"
        mask_count="$(awk -F '\t' -v wanted="${chr}" '$1==wanted{n++} END{print n+0}' "${evaluation_mask}")"
        if (( mask_count == 0 )); then
            [[ ! -e "${truth_vcf}" && ! -e "${truth_vcf}.csi" && ! -e "${truth_vcf}.tbi" ]] || {
                die "Unexpected masked-truth file for an empty mask chromosome: ${truth_vcf}"
            }
            continue
        fi
        require_file "${truth_vcf}"
        [[ -f "${truth_vcf}.csi" || -f "${truth_vcf}.tbi" ]] || die "Masked truth VCF is not indexed: ${truth_vcf}"
        count="$(
            bcftools isec -c none -n=2 -w1 -r "${chr}" -Ou "${masked_imputed}" "${truth_vcf}" \
                | bcftools view -H \
                | wc -l \
                | tr -d ' '
        )"
        printf '%s\t%s\n' "${chr}" "${count}" >> "${output}"
        total=$((total + count))
    done
    (( total > 0 )) || die "No exact WGS/imputed matches remain for ${masked_imputed}."
}

# =============================================================================
# Submission, collection, and manifests
# =============================================================================

validate_and_summarise_metrics() {
    local metrics="$1" summary_out="$2"
    local counts="${TMP_WORK}/metrics.counts.$$"
    local expected_header=$'sample\tr_overall\tr2_overall\tn_variants\tr_maf_[0.0,0.1)\tr2_maf_[0.0,0.1)\tn_maf_[0.0,0.1)\tr_maf_[0.1,0.2)\tr2_maf_[0.1,0.2)\tn_maf_[0.1,0.2)\tr_maf_[0.2,0.3)\tr2_maf_[0.2,0.3)\tn_maf_[0.2,0.3)\tr_maf_[0.3,0.4)\tr2_maf_[0.3,0.4)\tn_maf_[0.3,0.4)\tr_maf_[0.4,0.5]\tr2_maf_[0.4,0.5]\tn_maf_[0.4,0.5]'
    [[ "$(head -n 1 "${metrics}")" == "${expected_header}" ]] || {
        die "Per-sample metrics header does not match the unchanged WGS schema: ${metrics}"
    }
    awk -F '\t' '
      NR == 1 {
        for (i = 1; i <= NF; i++) {
          if ($i == "sample") sample_col = i
          if ($i == "r2_overall") r2_col = i
          if ($i == "n_variants") n_col = i
          if ($i ~ /^r2_maf_\[/) {
            maf_r2[++n_maf_r2] = i
            maf_r2_name[n_maf_r2] = $i
          }
          if ($i ~ /^n_maf_\[/) {
            maf_n[++n_maf_n] = i
            maf_n_name[n_maf_n] = $i
          }
        }
        if (!sample_col || !r2_col || !n_col || n_maf_r2 != 5 || n_maf_n != 5) exit 2
        for (j = 1; j <= n_maf_r2; j++) {
          expected_count = maf_r2_name[j]
          sub(/^r2_/, "n_", expected_count)
          found = 0
          for (k = 1; k <= n_maf_n; k++) {
            if (maf_n_name[k] == expected_count) found = 1
          }
          if (!found) exit 2
        }
        next
      }
      {
        if ($sample_col == "" || seen[$sample_col]++) exit 3
        if ($r2_col !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ || $r2_col < 0 || $r2_col > 1) exit 4
        if ($n_col !~ /^[0-9]+$/) exit 5
        for (j = 1; j <= n_maf_r2; j++) {
          value = $(maf_r2[j])
          if (value !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/ || value < 0 || value > 1) exit 4
        }
        for (j = 1; j <= n_maf_n; j++) {
          if ($(maf_n[j]) !~ /^[0-9]+$/) exit 5
        }
        print $n_col
      }
    ' "${metrics}" > "${counts}" || die "Invalid per-sample metrics TSV: ${metrics}"
    [[ -s "${counts}" ]] || die "Metrics TSV contains no samples: ${metrics}"
    sort -n "${counts}" -o "${counts}"
    awk '
      { values[NR] = $1 }
      END {
        if (NR % 2) median = values[(NR + 1) / 2]
        else median = (values[NR / 2] + values[NR / 2 + 1]) / 2
        printf "%d\t%d\t%.10g\t%d\n", NR, values[1], median, values[NR]
      }
    ' "${counts}" > "${summary_out}"
}

assert_metrics_samples() {
    local metrics="$1" expected_samples="$2" label="$3"
    local observed="${TMP_WORK}/$(basename "${metrics}").samples.$$"
    local expected_sorted="${TMP_WORK}/$(basename "${expected_samples}").expected_samples.$$"
    awk -F '\t' 'NR > 1 { print $1 }' "${metrics}" | sort -u > "${observed}"
    awk 'NF { print $1 }' "${expected_samples}" | sort -u > "${expected_sorted}"
    cmp -s "${observed}" "${expected_sorted}" || die "Sample IDs differ from the genome-wide baseline for ${label}."
}

manifest_value() {
    local manifest="$1" key="$2"
    awk -F '\t' -v wanted="${key}" '$1 == wanted { print $2; exit }' "${manifest}"
}

genomewide_result_dir() {
    local chunks_dir="$1"
    [[ "${chunks_dir}" == */chunks/imputed ]] || die "Cannot derive genome-wide result path from chunks directory: ${chunks_dir}"
    printf '%s/eval/dosage_eval_wgs\n' "${chunks_dir%/chunks/imputed}"
}

validate_genomewide_baseline() {
    local run_id="$1" chunks_dir="$2" output="$3"
    local baseline_dir baseline_manifest baseline_metrics baseline_samples baseline_signature
    baseline_dir="$(genomewide_result_dir "${chunks_dir}")"
    baseline_manifest="${baseline_dir}/run_manifest.tsv"
    baseline_metrics="${baseline_dir}/per_sample_metrics.tsv"
    baseline_samples="${baseline_dir}/intermediate/common_samples.txt"
    require_file "${baseline_dir}/.complete"
    require_file "${baseline_dir}/.evaluation_mode"
    require_file "${baseline_dir}/.run_signature"
    require_file "${baseline_manifest}"
    require_file "${baseline_metrics}"
    require_file "${baseline_samples}"
    [[ "$(<"${baseline_dir}/.evaluation_mode")" == "wgs" ]] || die "Genome-wide baseline is not WGS mode: ${baseline_dir}"
    [[ "$(manifest_value "${baseline_manifest}" truth_mode)" == "wgs" ]] || die "Genome-wide truth mode mismatch: ${baseline_manifest}"
    [[ "$(manifest_value "${baseline_manifest}" output_schema)" == "wgs-gt-isec-v5" ]] || die "Genome-wide evaluator schema mismatch: ${baseline_manifest}"
    [[ "$(manifest_value "${baseline_manifest}" comparison_field)" == "GT" ]] || die "Genome-wide comparison field is not GT: ${baseline_manifest}"
    [[ "$(manifest_value "${baseline_manifest}" genotype_encoding)" == "ALT_COUNT_0_1_2" ]] || die "Genome-wide genotype encoding mismatch: ${baseline_manifest}"
    local baseline_chromosomes expected_chromosomes observed_chromosomes
    baseline_chromosomes="$(manifest_value "${baseline_manifest}" chromosomes)"
    expected_chromosomes="${TMP_WORK}/${run_id}.expected_baseline_chromosomes.txt"
    observed_chromosomes="${TMP_WORK}/${run_id}.observed_baseline_chromosomes.txt"
    : > "${expected_chromosomes}"
    for number in $(seq 1 17); do
        local required_chr
        printf -v required_chr 'Chr%02d' "${number}"
        printf '%s\n' "${required_chr}" >> "${expected_chromosomes}"
    done
    printf '%s\n' "${baseline_chromosomes}" | tr ',' '\n' | awk 'NF' | sort > "${observed_chromosomes}"
    cmp -s "${expected_chromosomes}" "${observed_chromosomes}" || {
        die "Genome-wide baseline must cover exactly Chr01-Chr17 once each: ${baseline_manifest}"
    }

    local probe_parent="${TMP_WORK}/genomewide_signature" probe_dir="${TMP_WORK}/genomewide_signature/${run_id}"
    mkdir -p "${probe_parent}"
    bash "${WGS_EVALUATOR}" --chunks-dir "${chunks_dir}" \
        --truth-dataset-dir "${WGS_TRUTH_DIR}" --reference-fasta "${REFERENCE_FASTA}" \
        --out-prefix "${probe_dir}" --prepare-only
    require_file "${probe_dir}/.run_signature"
    baseline_signature="$(<"${baseline_dir}/.run_signature")"
    [[ "$(<"${probe_dir}/.run_signature")" == "${baseline_signature}" ]] || {
        die "Genome-wide baseline is stale or differs from current chunks/truth/reference/filter settings: ${baseline_dir}"
    }
    [[ "$(manifest_value "${baseline_manifest}" run_signature)" == "${baseline_signature}" ]] || {
        die "Genome-wide run manifest signature mismatch: ${baseline_manifest}"
    }

    local baseline_summary="${TMP_WORK}/${run_id}.genomewide.metrics_summary.tsv"
    validate_and_summarise_metrics "${baseline_metrics}" "${baseline_summary}"
    assert_metrics_samples "${baseline_metrics}" "${baseline_samples}" "${run_id} genome-wide result"
    printf '%s\t%s\t%s\t%s\n' "${baseline_dir}" "${baseline_signature}" \
        "$(md5_file "${baseline_metrics}")" "${baseline_samples}" > "${output}"
}

validate_completed_result() {
    local run_id="$1" mask_id="$2" masked_imputed="$3" truth_dir="$4"
    local result_dir="$5" chr_csv="$6" baseline_samples="$7"
    bash "${WGS_EVALUATOR}" --imputed "${masked_imputed}" \
        --truth-dataset-dir "${truth_dir}" --reference-fasta "${REFERENCE_FASTA}" \
        --chr "${chr_csv}" --out-prefix "${result_dir}" --prepare-only
    require_file "${result_dir}/.complete"
    require_file "${result_dir}/.run_signature"
    require_file "${result_dir}/run_manifest.tsv"
    require_file "${result_dir}/intermediate/common_samples.txt"
    [[ "$(manifest_value "${result_dir}/run_manifest.tsv" run_signature)" == "$(<"${result_dir}/.run_signature")" ]] || {
        die "Completed evaluator run-manifest signature mismatch for ${run_id}/${mask_id}."
    }
    cmp -s <(sort -u "${result_dir}/intermediate/common_samples.txt") <(sort -u "${baseline_samples}") || {
        die "Evaluator sample IDs differ from the genome-wide baseline for ${run_id}/${mask_id}."
    }
    assert_metrics_samples "${result_dir}/per_sample_metrics.tsv" "${baseline_samples}" "${run_id}/${mask_id}"
}

stage_metrics() {
    local run_type="$1" treatment="$2" panel="$3" mask_id="$4" metrics="$5"
    local run_dir treatment_dir filename destination temporary
    run_dir="$(printf '%s' "${run_type}" | tr '[:upper:]' '[:lower:]')_run"
    if [[ "${treatment}" == "Filtered" ]]; then
        treatment_dir="filter_treatment"
    else
        treatment_dir="no_filter_treatment"
    fi
    filename="${run_type}_${treatment}_WGS_2x_${panel}__scope-${mask_id}.per_sample_metrics.tsv"
    destination="${OUTPUT_ROOT}/staging/${run_dir}/${treatment_dir}/${filename}"
    mkdir -p "$(dirname "${destination}")"
    if [[ -f "${destination}" ]]; then
        if cmp -s "${metrics}" "${destination}"; then
            STAGED_METRICS_RESULT="${destination}"
            return
        fi
        [[ "${FORCE}" == "true" ]] || die "Staged metrics differ from the validated evaluator output; set FORCE=true to replace: ${destination}"
    fi
    temporary="${destination}.tmp.$$"
    cp "${metrics}" "${temporary}"
    mv -f "${temporary}" "${destination}"
    cmp -s "${metrics}" "${destination}" || die "Staged metrics differ from evaluator output: ${destination}"
    STAGED_METRICS_RESULT="${destination}"
}

extract_receipt_ids() {
    local receipt="$1" array_id finalizer_id
    if [[ "$(head -n 1 "${receipt}")" == $'run_id\tmask_id\tarray_job_id\tfinalizer_job_id\t'* ]]; then
        array_id="$(awk -F '\t' 'NR==2{print $3}' "${receipt}")"
        finalizer_id="$(awk -F '\t' 'NR==2{print $4}' "${receipt}")"
    else
        array_id="$(sed -n 's/.*Submitted WGS chromosome array job \([0-9][0-9]*\).*/\1/p' "${receipt}" | tail -n 1)"
        finalizer_id="$(sed -n 's/.*Submitted WGS finalizer job \([0-9][0-9]*\).*/\1/p' "${receipt}" | tail -n 1)"
    fi
    printf '%s\t%s\n' "${array_id:-NA}" "${finalizer_id:-NA}"
}

validate_submission_receipt() {
    local receipt="$1" expected_evaluator_hash="$2" expected_mask_hash="$3"
    local expected_baseline_signature="$4" result_dir="$5"
    local receipt_run receipt_mask _array _finalizer _submitted receipt_evaluator
    local receipt_mask_hash receipt_baseline_signature receipt_run_signature
    require_file "${receipt}"
    [[ "$(head -n 1 "${receipt}")" == $'run_id\tmask_id\tarray_job_id\tfinalizer_job_id\tsubmitted_at\tevaluator_sha256\tmask_sha256\tgenomewide_run_signature\tmasked_run_signature' ]] || {
        die "Submission receipt lacks the required provenance fields: ${receipt}"
    }
    IFS=$'\t' read -r receipt_run receipt_mask _array _finalizer _submitted receipt_evaluator \
        receipt_mask_hash receipt_baseline_signature receipt_run_signature < <(sed -n '2p' "${receipt}")
    [[ -n "${receipt_run}" && -n "${receipt_mask}" ]] || die "Submission receipt has no data row: ${receipt}"
    [[ "${receipt_evaluator}" == "${expected_evaluator_hash}" ]] || die "Evaluator hash differs from the submission receipt: ${receipt}"
    [[ "${receipt_mask_hash}" == "${expected_mask_hash}" ]] || die "Mask hash differs from the submission receipt: ${receipt}"
    [[ "${receipt_baseline_signature}" == "${expected_baseline_signature}" ]] || die "Genome-wide baseline signature differs from the submission receipt: ${receipt}"
    [[ "${receipt_run_signature}" == "$(<"${result_dir}/.run_signature")" ]] || die "Masked evaluator signature differs from the submission receipt: ${receipt}"
}

write_treatment_manifests() {
    local global_manifest="$1" run_type run_dir treatment treatment_dir subset destination
    for run_type in Biased Unbiased; do
        run_dir="$(printf '%s' "${run_type}" | tr '[:upper:]' '[:lower:]')_run"
        for treatment in Filtered No_filter; do
            [[ "${treatment}" == "Filtered" ]] && treatment_dir="filter_treatment" || treatment_dir="no_filter_treatment"
            destination="${OUTPUT_ROOT}/staging/${run_dir}/${treatment_dir}/wgs_array_mask_sensitivity_manifest.tsv"
            subset="${TMP_WORK}/${run_type}.${treatment}.manifest.tsv"
            awk -F '\t' -v OFS='\t' -v wanted_run="${run_type}" -v wanted_treatment="${treatment}" '
              NR == 1 || ($1 == wanted_run && $2 == wanted_treatment)
            ' "${global_manifest}" > "${subset}"
            mkdir -p "$(dirname "${destination}")"
            cp "${subset}" "${destination}.tmp.$$"
            mv -f "${destination}.tmp.$$" "${destination}"
        done
    done
}

dry_run_commands() {
    local masks_root="$1" mask_manifest="${masks_root}/masks/mask_manifest.tsv"
    local run_id run_type treatment panel chunks_dir
    local mask_id scope_label _dataset source_count canonical_count eligible_count
    local _missing _missing_alt _non_biallelic _normalization_excluded _duplicates _canonical_hash mask_hash _rest
    local planned=0
    while IFS=$'\t' read -r run_id run_type treatment panel chunks_dir; do
        local baseline_dir
        baseline_dir="$(genomewide_result_dir "${chunks_dir}")"
        require_file "${baseline_dir}/.complete"
        require_file "${baseline_dir}/.run_signature"
        require_file "${baseline_dir}/run_manifest.tsv"
        require_file "${baseline_dir}/per_sample_metrics.tsv"
        while IFS=$'\t' read -r mask_id scope_label _dataset source_count canonical_count eligible_count _missing _missing_alt _non_biallelic _normalization_excluded _duplicates _canonical_hash mask_hash _rest; do
            [[ "${mask_id}" != "scope_id" ]] || continue
            local masked_imputed="${OUTPUT_ROOT}/cache/imputed/${run_id}/${mask_id}/imputed.masked.vcf.gz"
            local truth_dir="${OUTPUT_ROOT}/cache/truth/${mask_id}"
            local result_dir="${OUTPUT_ROOT}/results/${run_type}/${treatment}/${panel}/${mask_id}"
            local command=(env -u SLURM_JOB_ID -u SLURM_ARRAY_JOB_ID -u SLURM_ARRAY_TASK_ID
                bash "${DOSAGE_R2}" --truth-mode wgs --imputed "${masked_imputed}"
                --truth-dataset-dir "${truth_dir}" --reference-fasta "${REFERENCE_FASTA}"
                --out-prefix "${result_dir}")
            [[ "${FORCE}" != "true" ]] || command+=(-- --force)
            print_command "${command[@]}"
            planned=$((planned + 1))
        done < "${mask_manifest}"
    done < <(run_rows)
    [[ "${planned}" == "60" ]] || die "Dry-run expected 60 analyses; found ${planned}."
    log "Dry-run validated five masks, 12 run paths, and ${planned} planned evaluations. No persistent files or jobs were created."
}

run_analyses() {
    local cache_root="$1" mask_manifest="${cache_root}/masks/mask_manifest.tsv"
    local manifest_tmp="${TMP_WORK}/analysis_manifest.tsv"
    local receipts_root="${OUTPUT_ROOT}/receipts" logs_root="${OUTPUT_ROOT}/submission_logs"
    mkdir -p "${receipts_root}" "${logs_root}" "${OUTPUT_ROOT}/results" "${OUTPUT_ROOT}/staging"

    local environment_file="${QUILT2_ROOT}/config/environment.sh"
    local environment_template="${QUILT2_ROOT}/config/environment.template.sh"
    if [[ -f "${environment_file}" ]]; then
        # shellcheck source=/dev/null
        source "${environment_file}"
    elif [[ -f "${environment_template}" ]]; then
        # shellcheck source=/dev/null
        source "${environment_template}"
    fi
    local truth_filter_enabled truth_min_gq truth_min_dp evaluator_sha256
    truth_filter_enabled="${QUILT2_WGS_TRUTH_FILTER_ENABLED:-true}"
    truth_min_gq="${QUILT2_WGS_TRUTH_MIN_GQ:-60}"
    truth_min_dp="${QUILT2_WGS_TRUTH_MIN_DP:-10}"
    # Composite hash of the unchanged wrapper, WGS shell/R evaluator, and concat helper.
    evaluator_sha256="$(evaluator_bundle_sha256)"

    printf 'run_type\ttreatment\treference_panel\tscope_id\tscope_label\tsource_snp_count\tcanonical_position_count\twgs_eligible_position_count\texact_match_variant_count\tmask_sha256\tmask_membership_sha256\tmasked_truth_content_sha256\tcoordinate_mapping_assumption\ttruth_filter_enabled\ttruth_min_gq\ttruth_min_dp\tevaluator_sha256\tgenomewide_result_dir\tgenomewide_run_signature\tgenomewide_metrics_md5\timputed_chunks_dir\tmasked_imputed_vcf\tmasked_truth_dir\tresult_dir\tstatus\tarray_job_id\tfinalizer_job_id\tmetrics_tsv\tmetrics_md5\tn_samples\tmin_usable_variants\tmedian_usable_variants\tmax_usable_variants\n' > "${manifest_tmp}"

    local run_id run_type treatment panel chunks_dir concatenated=""
    while IFS=$'\t' read -r run_id run_type treatment panel chunks_dir; do
        local run_selected=false
        local baseline_dir baseline_signature="NA" baseline_metrics_md5="NA" baseline_samples="NA"
        is_selected_run "${run_id}" && run_selected=true
        baseline_dir="$(genomewide_result_dir "${chunks_dir}")"
        if [[ "${run_selected}" == "true" ]]; then
            local baseline_info="${TMP_WORK}/${run_id}.genomewide_baseline.tsv"
            validate_genomewide_baseline "${run_id}" "${chunks_dir}" "${baseline_info}"
            IFS=$'\t' read -r baseline_dir baseline_signature baseline_metrics_md5 baseline_samples < "${baseline_info}"
            log "Preparing concatenated imputed VCF for ${run_id}."
            prepare_concatenated_run "${cache_root}" "${run_id}" "${chunks_dir}" "${baseline_signature}"
            concatenated="${CONCATENATED_RESULT}"
            require_file "${concatenated}"
            [[ -f "${concatenated}.csi" || -f "${concatenated}.tbi" ]] || die "Concatenated VCF is not indexed: ${concatenated}"
        fi

        local mask_id scope_label dataset_name source_count canonical_count eligible_count
        local missing_count missing_alt_count non_biallelic_count normalization_excluded_count duplicate_count
        local canonical_hash mask_hash
        local sample_count samples_hash samples_path canonical_path evaluation_mask exclusions_file
        local record_exclusion_count record_exclusions_file record_exclusions_hash truth_dir
        local source_vcf_hash metadata_hash truth_signature coordinate_assumption
        local source_exclusion_count source_exclusions_file
        local truth_checksum_manifest truth_content_hash
        while IFS=$'\t' read -r mask_id scope_label dataset_name source_count canonical_count eligible_count missing_count missing_alt_count non_biallelic_count normalization_excluded_count duplicate_count canonical_hash mask_hash sample_count samples_hash samples_path canonical_path evaluation_mask exclusions_file record_exclusion_count record_exclusions_file record_exclusions_hash truth_dir source_vcf_hash metadata_hash truth_signature coordinate_assumption source_exclusion_count source_exclusions_file truth_checksum_manifest truth_content_hash; do
            [[ "${mask_id}" != "scope_id" ]] || continue
            local selected=false scope_tag="full" imputed_cache_root="${cache_root}"
            if [[ "${run_selected}" == "true" ]] && is_selected_mask "${mask_id}"; then selected=true; fi
            local masked_imputed="${cache_root}/imputed/${run_id}/${mask_id}/imputed.masked.vcf.gz"
            local result_dir="${OUTPUT_ROOT}/results/${run_type}/${treatment}/${panel}/${mask_id}"
            local receipt="${receipts_root}/${run_id}__${mask_id}.tsv"
            local attempt="${receipts_root}/${run_id}__${mask_id}.attempt.tsv"
            local submission_log="${logs_root}/${run_id}__${mask_id}.log"
            local exact_file="${cache_root}/imputed/${run_id}/${mask_id}/exact_matches_by_chromosome.tsv"
            if [[ "${selected}" == "true" && "${EXECUTION_SCOPE}" == "pilot" ]]; then
                scope_tag="$(scope_suffix)"
                imputed_cache_root="${cache_root}/${scope_tag}"
                masked_imputed="${imputed_cache_root}/imputed/${run_id}/${mask_id}/imputed.masked.vcf.gz"
                result_dir="${OUTPUT_ROOT}/results/${scope_tag}/${run_type}/${treatment}/${panel}/${mask_id}"
                receipt="${receipts_root}/${scope_tag}__${run_id}__${mask_id}.tsv"
                attempt="${receipts_root}/${scope_tag}__${run_id}__${mask_id}.attempt.tsv"
                submission_log="${logs_root}/${scope_tag}__${run_id}__${mask_id}.log"
                exact_file="${imputed_cache_root}/imputed/${run_id}/${mask_id}/exact_matches_by_chromosome.tsv"
            fi

            local status="planned" array_job_id="NA" finalizer_job_id="NA" staged_metrics="NA" metrics_md5="NA"
            local exact_count="NA" n_samples="NA" min_n="NA" median_n="NA" max_n="NA"
            local complete_file="${result_dir}/.complete" metrics="${result_dir}/per_sample_metrics.tsv"

            if [[ "${selected}" == "true" ]]; then
                local prepared
                prepare_masked_imputed "${imputed_cache_root}" "${run_id}" "${concatenated}" "${mask_id}" "${evaluation_mask}" "${mask_hash}" "${baseline_signature}"
                prepared="${MASKED_IMPUTED_RESULT}"
                [[ "${prepared}" == "${masked_imputed}" ]] || die "Internal masked-imputed path mismatch: ${prepared} != ${masked_imputed}"
                mkdir -p "$(dirname "${exact_file}")"
                local chromosome_limit=""
                [[ "${scope_tag}" != "full" ]] && chromosome_limit="${PILOT_CHROMOSOME}"
                preflight_exact_matches "${masked_imputed}" "${truth_dir}" "${exact_file}.tmp.$$" "${chromosome_limit}" "${evaluation_mask}"
                mv -f "${exact_file}.tmp.$$" "${exact_file}"
                exact_count="$(awk -F '\t' 'NR>1{s+=$2} END{print s+0}' "${exact_file}")"
                local chr_csv
                local zero_match_chromosomes
                zero_match_chromosomes="$(awk -F '\t' 'NR>1 && $2==0{print $1}' "${exact_file}" | paste -sd,)"
                [[ -z "${zero_match_chromosomes}" ]] || {
                    die "Mask-bearing chromosome(s) have no exact imputed/truth allele matches for ${run_id}/${mask_id}: ${zero_match_chromosomes}"
                }
                chr_csv="$(awk -F '\t' 'NR>1{print $1}' "${exact_file}" | paste -sd,)"
                [[ -n "${chr_csv}" ]] || die "No exact-match chromosomes found for ${run_id}/${mask_id}."

                if [[ -f "${complete_file}" && ! -s "${metrics}" && "${FORCE}" != "true" ]]; then
                    die "Completed result is missing per-sample metrics; inspect or set FORCE=true: ${result_dir}"
                elif [[ -f "${complete_file}" && -s "${metrics}" && "${FORCE}" != "true" ]]; then
                    [[ -f "${receipt}" ]] || {
                        die "Completed result has no utility submission receipt and cannot be staged automatically: ${result_dir}"
                    }
                    validate_completed_result "${run_id}" "${mask_id}" "${masked_imputed}" "${truth_dir}" \
                        "${result_dir}" "${chr_csv}" "${baseline_samples}"
                    validate_submission_receipt "${receipt}" "${evaluator_sha256}" "${mask_hash}" \
                        "${baseline_signature}" "${result_dir}"
                    local metrics_summary="${TMP_WORK}/${run_id}.${mask_id}.metrics_summary.tsv"
                    validate_and_summarise_metrics "${metrics}" "${metrics_summary}"
                    IFS=$'\t' read -r n_samples min_n median_n max_n < "${metrics_summary}"
                    require_file "${result_dir}/qc/filter_summary.tsv"
                    local evaluator_exact_count
                    evaluator_exact_count="$(awk -F '\t' '$1=="variants" && $2=="ALL" && $3=="exact_matches"{print $4}' "${result_dir}/qc/filter_summary.tsv")"
                    [[ "${evaluator_exact_count}" =~ ^[0-9]+$ ]] || die "Could not read exact-match count from ${result_dir}/qc/filter_summary.tsv"
                    [[ "${exact_count}" == "${evaluator_exact_count}" ]] || {
                        die "Exact-match audit differs for ${run_id}/${mask_id}: preflight ${exact_count}, evaluator ${evaluator_exact_count}."
                    }
                    if [[ "${scope_tag}" == "full" ]]; then
                        stage_metrics "${run_type}" "${treatment}" "${panel}" "${mask_id}" "${metrics}"
                        staged_metrics="${STAGED_METRICS_RESULT}"
                        metrics_md5="$(md5_file "${staged_metrics}")"
                        status="staged"
                    else
                        status="pilot_complete"
                    fi
                    if [[ -f "${receipt}" ]]; then
                        IFS=$'\t' read -r array_job_id finalizer_job_id < <(extract_receipt_ids "${receipt}")
                    fi
                elif [[ -f "${receipt}" && ! -f "${complete_file}" && "${RESUBMIT_INCOMPLETE}" != "true" ]]; then
                    IFS=$'\t' read -r array_job_id finalizer_job_id < <(extract_receipt_ids "${receipt}")
                    status="submitted_or_incomplete"
                    warn "Skipping receipt-backed incomplete analysis: ${run_id}/${mask_id}. Check Slurm and ${submission_log}."
                elif [[ -f "${attempt}" && ! -f "${receipt}" && ! -f "${complete_file}" && "${RESUBMIT_INCOMPLETE}" != "true" ]]; then
                    status="ambiguous_submission_attempt"
                    warn "Skipping an analysis with an ambiguous prior submission attempt: ${run_id}/${mask_id}. Confirm/cancel any jobs, then set RESUBMIT_INCOMPLETE=true if needed."
                else
                    [[ ! -f "${receipt}" && ! -f "${attempt}" ]] || {
                        warn "Explicitly resubmitting ${run_id}/${mask_id}; confirm prior Slurm jobs were cancelled."
                    }
                    local command=(env -u SLURM_JOB_ID -u SLURM_ARRAY_JOB_ID -u SLURM_ARRAY_TASK_ID
                        bash "${DOSAGE_R2}" --truth-mode wgs --imputed "${masked_imputed}"
                        --truth-dataset-dir "${truth_dir}" --reference-fasta "${REFERENCE_FASTA}"
                        --chr "${chr_csv}" --out-prefix "${result_dir}")
                    [[ "${FORCE}" != "true" ]] || command+=(-- --force)
                    log "Submitting ${run_id}/${mask_id} (${exact_count} structural exact matches across ${chr_csv})."
                    mkdir -p "$(dirname "${submission_log}")"
                    local attempted_at="$(date -u +%Y-%m-%dT%H:%M:%SZ)"
                    if [[ -f "${attempt}" ]]; then
                        cp "${attempt}" "${attempt}.tmp.$$"
                    else
                        printf 'run_id\tmask_id\tattempted_at\tcoordinator_job_id\tresult_dir\tevaluator_sha256\tmask_sha256\tgenomewide_run_signature\n' > "${attempt}.tmp.$$"
                    fi
                    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                        "${run_id}" "${mask_id}" "${attempted_at}" "${SLURM_JOB_ID}" \
                        "${result_dir}" "${evaluator_sha256}" "${mask_hash}" "${baseline_signature}" >> "${attempt}.tmp.$$"
                    mv -f "${attempt}.tmp.$$" "${attempt}"
                    if "${command[@]}" 2>&1 | tee "${submission_log}.tmp.$$"; then
                        mv -f "${submission_log}.tmp.$$" "${submission_log}"
                    else
                        local rc=${PIPESTATUS[0]}
                        mv -f "${submission_log}.tmp.$$" "${submission_log}"
                        die "Submission failed or became ambiguous for ${run_id}/${mask_id} with exit status ${rc}; the attempt record prevents automatic resubmission."
                    fi
                    IFS=$'\t' read -r array_job_id finalizer_job_id < <(extract_receipt_ids "${submission_log}")
                    [[ "${array_job_id}" != "NA" && "${finalizer_job_id}" != "NA" ]] || {
                        die "Could not parse submitted job IDs from ${submission_log}; the attempt record prevents automatic resubmission."
                    }
                    require_file "${result_dir}/.run_signature"
                    {
                        printf 'run_id\tmask_id\tarray_job_id\tfinalizer_job_id\tsubmitted_at\tevaluator_sha256\tmask_sha256\tgenomewide_run_signature\tmasked_run_signature\n'
                        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                            "${run_id}" "${mask_id}" "${array_job_id}" "${finalizer_job_id}" \
                            "${attempted_at}" "${evaluator_sha256}" "${mask_hash}" \
                            "${baseline_signature}" "$(<"${result_dir}/.run_signature")"
                    } > "${receipt}.tmp.$$"
                    mv -f "${receipt}.tmp.$$" "${receipt}"
                    status="submitted"
                fi
            elif [[ -f "${complete_file}" ]]; then
                status="complete_not_selected"
            elif [[ -f "${receipt}" ]]; then
                IFS=$'\t' read -r array_job_id finalizer_job_id < <(extract_receipt_ids "${receipt}")
                status="submitted_or_incomplete"
            elif [[ -f "${attempt}" ]]; then
                status="ambiguous_submission_attempt"
            fi

            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${run_type}" "${treatment}" "${panel}" "${mask_id}" "${scope_label}" \
                "${source_count}" "${canonical_count}" "${eligible_count}" "${exact_count}" "${mask_hash}" \
                "${samples_hash}" "${truth_content_hash}" "${coordinate_assumption}" "${truth_filter_enabled}" "${truth_min_gq}" "${truth_min_dp}" \
                "${evaluator_sha256}" "${baseline_dir}" "${baseline_signature}" "${baseline_metrics_md5}" \
                "${chunks_dir}" "${masked_imputed}" "${truth_dir}" "${result_dir}" "${status}" \
                "${array_job_id}" "${finalizer_job_id}" "${staged_metrics}" "${metrics_md5}" "${n_samples}" "${min_n}" \
                "${median_n}" "${max_n}" >> "${manifest_tmp}"
        done < "${mask_manifest}"
    done < <(run_rows)

    local manifest="${OUTPUT_ROOT}/analysis_manifest.tsv"
    cp "${manifest_tmp}" "${manifest}.tmp.$$"
    mv -f "${manifest}.tmp.$$" "${manifest}"
    cp "${manifest}" "${OUTPUT_ROOT}/staging/analysis_manifest.tsv.tmp.$$"
    mv -f "${OUTPUT_ROOT}/staging/analysis_manifest.tsv.tmp.$$" "${OUTPUT_ROOT}/staging/analysis_manifest.tsv"
    write_treatment_manifests "${manifest}"
    log "Updated analysis manifest: ${manifest}"
}

main() {
    preflight "$@"

    local cache_root
    if [[ "${DRY_RUN}" == "true" ]]; then
        cache_root="${TMP_WORK}/dry_run_cache"
    else
        cache_root="${OUTPUT_ROOT}/cache"
        mkdir -p "${cache_root}"
        acquire_coordinator_lock
    fi

    local repaired_vcf
    repair_hortqg_vcf "${cache_root}"
    repaired_vcf="${REPAIRED_VCF_RESULT}"
    build_canonical_masks "${cache_root}" "${repaired_vcf}"
    build_masked_truth "${cache_root}"

    if [[ "${DRY_RUN}" == "true" ]]; then
        dry_run_commands "${cache_root}"
        return
    fi

    run_analyses "${cache_root}"
}

main "$@"
