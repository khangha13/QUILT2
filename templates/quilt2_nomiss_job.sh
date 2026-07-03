#!/bin/bash
# QUILT2 Phase 1 array worker: filter panel VCF by phased rate (per chromosome)
set -euo pipefail

# The orchestrator substitutes this placeholder with an absolute path.
QUILT2_ROOT="${QUILT2_ROOT:-__QUILT2_ROOT__}"

if [ ! -d "${QUILT2_ROOT}" ]; then
    echo "[quilt2_nomiss_job] ERROR: QUILT2_ROOT does not exist: ${QUILT2_ROOT}" >&2
    exit 1
fi

source "${QUILT2_ROOT}/lib/functions.sh"

if [ "$#" -lt 6 ]; then
    log_error "Usage: quilt2_nomiss_job.sh <WORK_DIR> <REFERENCE_PANEL_DIR> <PANEL_STANDARDISED_DIR> <PANEL_NOMISS_DIR> <MIN_VALID_GT_RATE> <CHR_MANIFEST> [BCFTOOLS_MODULE] [QUILT2_CONDA_ENV] [FAIL_FLAG]"
    exit 1
fi

WORK_DIR="$1"
REFERENCE_PANEL_DIR="$2"
PANEL_STANDARDISED_DIR="$3"
PANEL_NOMISS_DIR="$4"
MIN_VALID_GT_RATE="$5"
CHR_MANIFEST="$6"
BCFTOOLS_MODULE="${7:-${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}}"
QUILT2_CONDA_ENV="${8:-${QUILT2_CONDA_ENV:-quilt2}}"
FAIL_FLAG="${9:-${NOMISS_FAIL_FLAG:-}}"
STANDARDISE_NAME="${STANDARDISE_NAME:-false}"
STANDARDISE_NAME_FORCE="${STANDARDISE_NAME_FORCE:-false}"
STANDARDISE_SUFFIX="${STANDARDISE_SUFFIX:-_chr}"
REFERENCE_FASTA="${REFERENCE_FASTA:-}"
REFERENCE_FASTA_INDEX="${REFERENCE_FASTA_INDEX:-${REFERENCE_FASTA}.fai}"

# Propagate DRY_RUN if set by orchestrator
DRY_RUN="${DRY_RUN:-false}"

# Export flags consumed by helpers
export REMOVE_MISSING="${REMOVE_MISSING:-false}"
export MIN_VALID_GT_RATE
export PANEL_OUT_DIR="${PANEL_NOMISS_DIR}"
export MISSING_REPORT="${MISSING_REPORT:-${PANEL_NOMISS_DIR%/}/missing_sites_removed.tsv}"
export CHUNK_FILE="" # ensure normalize_panel_vcf contig warning is meaningful

# Guardrails
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    log_error "SLURM_ARRAY_TASK_ID is not set; this script must run as a SLURM array job."
    exit 1
fi

if [ ! -f "${CHR_MANIFEST}" ]; then
    log_error "Chromosome manifest not found: ${CHR_MANIFEST}"
    exit 1
fi

readarray -t CHR_LIST < "${CHR_MANIFEST}"
chr_count="${#CHR_LIST[@]}"
if ! [[ "${SLURM_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [ "${SLURM_ARRAY_TASK_ID}" -lt 0 ] || [ "${SLURM_ARRAY_TASK_ID}" -ge "${chr_count}" ]; then
    log_error "Array index ${SLURM_ARRAY_TASK_ID} out of range 0..$((chr_count-1))"
    exit 1
fi

CHR="${CHR_LIST[$SLURM_ARRAY_TASK_ID]}"

if [[ -n "${FAIL_FLAG}" ]]; then
    on_error() {
        log_error "Phase 1 filtering failed for ${CHR}; setting failure flag ${FAIL_FLAG}"
        touch "${FAIL_FLAG}" || true
    }
    trap on_error ERR
fi

# Ensure output directories exist
mkdir -p "${PANEL_STANDARDISED_DIR}" "${PANEL_NOMISS_DIR}"

# Tooling
export BCFTOOLS_MODULE QUILT2_CONDA_ENV
load_quilt_env || exit 1
ensure_bcftools || exit 1

# Check if VCF contigs are already in ChrNN format (no rename needed)
vcf_has_chr_contigs() {
    local src="$1"
    # Get first data contig from VCF
    local first_contig
    first_contig="$(bcftools view -H "${src}" 2>/dev/null | head -1 | cut -f1)"
    if [[ -z "${first_contig}" ]]; then
        # Empty VCF or can't read; check header contigs instead
        first_contig="$(bcftools view -h "${src}" 2>/dev/null | grep '^##contig=<ID=' | head -1 | sed 's/.*ID=\([^,>]*\).*/\1/')"
    fi
    # Return true if contig starts with "Chr" (case-insensitive)
    [[ "${first_contig,,}" == chr* ]]
}

# Ensure VCF has a tabix index (required for bcftools to handle undefined contigs)
ensure_vcf_indexed() {
    local src="$1"
    
    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "+ tabix -f -p vcf \"${src}\""
        return 0
    fi
    
    # Check if .tbi index exists and is newer than the VCF
    if [[ -f "${src}.tbi" && "${src}.tbi" -nt "${src}" ]]; then
        log_info "Tabix index already exists for ${src}"
        return 0
    fi
    
    # Create tabix index (handles undefined contigs in header)
    log_info "Creating tabix index for ${src}"
    if tabix -f -p vcf "${src}" 2>/dev/null; then
        return 0
    else
        log_warn "tabix indexing failed for ${src}; bcftools may fail on undefined contigs"
        return 1
    fi
}

standardize_panel_vcf() {
    local chr="$1"
    local src_dir="$2"
    local out_dir="$3"
    local suffix="$4"
    local force="$5"

    local src
    src="$(pick_panel_vcf "${src_dir}" "${chr}")"
    if [[ -z "${src}" ]]; then
        log_error "No source VCF found to standardise for ${chr} in ${src_dir}"
        return 1
    fi

    mkdir -p "${out_dir}"
    local dest="${out_dir%/}/${chr}${suffix}.vcf.gz"

    if [[ -f "${dest}" && "${force}" != "true" ]]; then
        if bcftools index -f -c "${dest}" >/dev/null 2>&1; then
            log_info "Standardised VCF already exists for ${chr}: ${dest}"
            echo "${dest}"
            return 0
        else
            log_warn "Existing standardised VCF for ${chr} is not sortable/indexable; rebuilding with --standardise-name-force semantics."
            force="true"
        fi
    fi

    # Check if contigs already have Chr prefix — skip rename if so
    if vcf_has_chr_contigs "${src}"; then
        log_info "VCF already has Chr-prefixed contigs; skipping rename for ${chr}"
        # Just copy/link and index
        if [[ "${DRY_RUN}" == "true" ]]; then
            echo "+ cp \"${src}\" \"${dest}\""
            echo "+ bcftools index -f -c \"${dest}\""
            echo "${dest}"
            return 0
        fi
        cp "${src}" "${dest}"
        bcftools index -f -c "${dest}"
        echo "${dest}"
        return 0
    fi

    # Ensure source VCF has tabix index (allows bcftools to handle undefined contigs)
    ensure_vcf_indexed "${src}" || true

    # Build a simple rename map 1..17 -> Chr01..Chr17 (matches apple panel); already-Chr contigs remain unchanged.
    local rename_map
    rename_map="$(mktemp)"
    trap 'rm -f "${rename_map}"' RETURN
    : > "${rename_map}"
    for i in $(seq 1 17); do
        printf "%d\tChr%02d\n" "${i}" "${i}" >> "${rename_map}"
    done

    if [[ "${DRY_RUN}" == "true" ]]; then
        cat <<EOF
+ bcftools annotate --rename-chrs "${rename_map}" "${src}" | bcftools sort -Oz -o "${dest}"
+ bcftools index -f -c "${dest}"
EOF
        echo "${dest}"
        return 0
    fi

    log_info "Standardising ${chr}: ${src} -> ${dest}"
    tmp_sorted="${dest}.sorted.tmp.vcf.gz"
    rm -f "${tmp_sorted}" "${dest}"
    if ! bcftools annotate --rename-chrs "${rename_map}" "${src}" | bcftools sort -Oz -o "${tmp_sorted}"; then
        log_error "Failed to standardise+sort ${chr}"
        return 1
    fi
    mv "${tmp_sorted}" "${dest}"
    if ! bcftools index -f -c "${dest}"; then
        log_error "Indexing failed for standardised VCF ${dest}"
        return 1
    fi
    echo "${dest}"
}

log_info "Phase 1: panel prep for ${CHR} (standardise=${STANDARDISE_NAME}, remove_missing=${REMOVE_MISSING}, min valid GT rate ${MIN_VALID_GT_RATE})"
panel_source_dir="${REFERENCE_PANEL_DIR}"
if [[ "${STANDARDISE_NAME}" == "true" ]]; then
    std_vcf="$(standardize_panel_vcf "${CHR}" "${REFERENCE_PANEL_DIR}" "${PANEL_STANDARDISED_DIR}" "${STANDARDISE_SUFFIX}" "${STANDARDISE_NAME_FORCE}")" || exit 1
    panel_source_dir="$(cd "$(dirname "${std_vcf}")" && pwd)"
fi
cleaned_vcf=""
if [[ "${REMOVE_MISSING}" == "true" ]]; then
    cleaned_vcf="$(normalize_panel_vcf "${CHR}" "${panel_source_dir}")"
else
    # No filtering requested; use the standardised (or original) panel as the output for downstream steps.
    cleaned_vcf="${panel_source_dir%/}/${CHR}${STANDARDISE_SUFFIX}.vcf.gz"
    if [[ ! -f "${cleaned_vcf}" ]]; then
        cleaned_vcf="$(pick_panel_vcf "${panel_source_dir}" "${CHR}")"
    fi
fi

# If remove-missing not requested, still ensure the source panel is indexed.
if [[ "${REMOVE_MISSING}" != "true" ]]; then
    if [[ -n "${cleaned_vcf}" && -f "${cleaned_vcf}" && "${DRY_RUN}" != "true" ]]; then
        if [[ "${cleaned_vcf}" =~ \.vcf\.gz$ ]]; then
            bcftools index -f -c "${cleaned_vcf}"
        fi
    fi
fi

if [[ "${DRY_RUN}" != "true" ]]; then
    if [[ -z "${cleaned_vcf}" || ! -s "${cleaned_vcf}" ]]; then
        log_error "Panel VCF missing or empty for ${CHR}: ${cleaned_vcf:-<empty>}"
        exit 1
    fi
    if [[ ! -f "${cleaned_vcf}.csi" && ! -f "${cleaned_vcf}.tbi" ]]; then
        log_error "Index for panel VCF missing for ${CHR}: ${cleaned_vcf}(.csi|.tbi)"
        exit 1
    fi
fi

log_info "Phase 1: completed ${CHR} → ${cleaned_vcf}"
