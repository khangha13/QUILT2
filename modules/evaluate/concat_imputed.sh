#!/bin/bash
# =============================================================================
# concat_imputed.sh – Concatenate per-chunk QUILT2 imputed VCFs into
#                      per-chromosome and genome-wide VCFs.
# =============================================================================
#
# run_quilt2.sh (Phase 2) writes one VCF per chunk under:
#   <chunks-dir>/<chr>/quilt2.diploid.<chr>.<start>-<end>.vcf.gz
#
# This script concatenates those into:
#   <out-dir>/<chr>/imputed.<chr>.vcf.gz     (per chromosome)
#   <out-dir>/imputed.all_chroms.vcf.gz      (genome-wide, only if 2+ chromosomes)
#
# Chunk order/completeness is resolved from the run_quilt2.sh chunk manifest
# when available (preferred: exact chunk boundaries, and missing chunks are
# detected before concatenating). Falls back to numeric-sorting chunk
# filenames by their start coordinate when no manifest is found — chunk start
# coordinates are not zero-padded, so a plain lexical/glob sort would not
# match genomic order.
#
# Adjacent chunks commonly overlap (e.g. QUILT::quilt_chunk_map()-derived
# chunks share a small boundary buffer), and each chunk's VCF contains every
# variant across its full regionStart-regionEnd, including that overlap. This
# script trims each chunk to end just before the next chunk's start (based on
# manifest/filename boundaries, not a guessed overlap size) before
# concatenating, so positions stay non-decreasing and `bcftools index` on the
# result succeeds.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Script lives at <repo>/modules/evaluate/concat_imputed.sh, so climb two levels.
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
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
Usage: concat_imputed.sh --chunks-dir <dir> [options]

Required:
  --chunks-dir PATH   Directory of per-chunk imputed VCFs, i.e. the
                       OUTPUT_DIR/chunks/imputed layout produced by run_quilt2.sh:
                       <chunks-dir>/<chr>/quilt2.diploid.<chr>.<start>-<end>.vcf.gz

Options:
  --manifest PATH     run_quilt2.sh chunk manifest (chunk_id|chr|start|end|buffer
                       lines, one per chunk). If omitted, auto-detected by reading
                       <chunks-dir>/../../run_manifest.tsv's "chunk_manifest" field
                       (the standard OUTPUT_DIR/chunks/imputed layout). Falls back
                       to filename-based ordering if no manifest is found or given.
  --chr LIST          Comma/space-separated chromosomes to include. Default: every
                       chromosome present in the manifest, or (in the no-manifest
                       fallback) every subdirectory of --chunks-dir.
  --out-dir PATH      Where to write concatenated VCFs. Default: same as --chunks-dir.
  --force             Re-concatenate even if outputs already exist.
  --help              Show this message and exit.

Outputs:
  <out-dir>/<chr>/imputed.<chr>.vcf.gz     One per chromosome (+ .csi index)
  <out-dir>/imputed.all_chroms.vcf.gz      Genome-wide, only if 2+ chromosomes
                                            were concatenated (+ .csi index)

Prints exactly one path to stdout on success: the genome-wide VCF if produced,
otherwise the single per-chromosome VCF. All logging goes to stderr, so this
script is safe to capture with $(...).
EOF
}

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

CHUNKS_DIR=""
MANIFEST_PATH=""
CHR_ARG=""
OUT_DIR=""
FORCE=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --chunks-dir) CHUNKS_DIR="$2"; shift 2 ;;
        --manifest)   MANIFEST_PATH="$2"; shift 2 ;;
        --chr)        CHR_ARG="$2"; shift 2 ;;
        --out-dir)    OUT_DIR="$2"; shift 2 ;;
        --force)      FORCE=true; shift ;;
        --help|-h)    usage; exit 0 ;;
        *)
            log_error "Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

if [[ -z "${CHUNKS_DIR}" ]]; then
    usage
    exit 1
fi
if [[ ! -d "${CHUNKS_DIR}" ]]; then
    log_error "Chunks directory not found: ${CHUNKS_DIR}"
    exit 1
fi
CHUNKS_DIR="$(cd "${CHUNKS_DIR}" && pwd)"

if [[ -z "${OUT_DIR}" ]]; then
    OUT_DIR="${CHUNKS_DIR}"
fi
mkdir -p "${OUT_DIR}"
OUT_DIR="$(cd "${OUT_DIR}" && pwd)"

ensure_bcftools || exit 1

TMP_DIR="$(mktemp -d "${TMPDIR:-/tmp}/concat_imputed.XXXXXX")"
cleanup() { rm -rf "${TMP_DIR}"; }
trap cleanup EXIT

# =============================================================================
# MANIFEST RESOLUTION
# =============================================================================

if [[ -z "${MANIFEST_PATH}" ]]; then
    OUTPUT_DIR_GUESS=""
    if candidate_dir="$(cd "${CHUNKS_DIR}/../.." 2>/dev/null && pwd)"; then
        OUTPUT_DIR_GUESS="${candidate_dir}"
    fi
    if [[ -n "${OUTPUT_DIR_GUESS}" && -f "${OUTPUT_DIR_GUESS}/run_manifest.tsv" ]]; then
        AUTO_RUN_MANIFEST="${OUTPUT_DIR_GUESS}/run_manifest.tsv"
        detected="$(awk -F'\t' '$1=="chunk_manifest"{print $2}' "${AUTO_RUN_MANIFEST}")"
        if [[ -n "${detected}" && -f "${detected}" ]]; then
            MANIFEST_PATH="${detected}"
            log_info "Auto-detected chunk manifest via ${AUTO_RUN_MANIFEST}: ${MANIFEST_PATH}"
        fi
    fi
fi

if [[ -n "${MANIFEST_PATH}" && ! -f "${MANIFEST_PATH}" ]]; then
    log_error "Chunk manifest not found: ${MANIFEST_PATH}"
    exit 1
fi

USE_MANIFEST=false
if [[ -n "${MANIFEST_PATH}" ]]; then
    USE_MANIFEST=true
    log_info "Using chunk manifest: ${MANIFEST_PATH}"
else
    log_warn "No chunk manifest found or given; falling back to filename-based chunk ordering."
    log_warn "This cannot detect missing chunks; correctness depends on run_quilt2.sh's naming convention."
fi

# =============================================================================
# CHROMOSOME LIST RESOLUTION
# =============================================================================

declare -a CHR_LIST=()
if [[ -n "${CHR_ARG}" ]]; then
    IFS=', ' read -r -a CHR_LIST <<< "${CHR_ARG}"
elif [[ "${USE_MANIFEST}" == "true" ]]; then
    while IFS= read -r chr; do
        CHR_LIST+=( "${chr}" )
    done < <(awk -F'|' '{print $2}' "${MANIFEST_PATH}" | awk '!seen[$0]++')
else
    while IFS= read -r dir; do
        CHR_LIST+=( "$(basename "${dir}")" )
    done < <(find "${CHUNKS_DIR}" -maxdepth 1 -mindepth 1 -type d | sort)
fi

if [[ "${#CHR_LIST[@]}" -eq 0 ]]; then
    log_error "No chromosomes resolved (check --chunks-dir / --chr / --manifest)."
    exit 1
fi

log_info "Chromosomes in scope: ${CHR_LIST[*]}"

# =============================================================================
# CHUNK ORDERING
# =============================================================================

# Writes the ordered list of chunks for one chromosome to stdout as
# "start\tend\tpath" lines (start/end kept so the caller can trim overlaps).
# Exit status: 0 = ok (rows printed), 1 = no chunks found for this chromosome,
# 2 = manifest lists chunks that are missing on disk (hard error, already logged).
build_chr_filelist() {
    local chr="$1"
    local chr_dir="${CHUNKS_DIR%/}/${chr}"

    if [[ "${USE_MANIFEST}" == "true" ]]; then
        local -a lines=()
        local -a missing=()
        local chunk_id m_chr start end buffer f
        while IFS='|' read -r chunk_id m_chr start end buffer; do
            [[ "${m_chr}" == "${chr}" ]] || continue
            f="${chr_dir}/quilt2.diploid.${chr}.${start}-${end}.vcf.gz"
            if [[ -s "${f}" ]]; then
                lines+=( "${start}"$'\t'"${end}"$'\t'"${f}" )
            else
                missing+=( "${f}" )
            fi
        done < "${MANIFEST_PATH}"

        if [[ "${#lines[@]}" -eq 0 && "${#missing[@]}" -eq 0 ]]; then
            return 1
        fi
        if [[ "${#missing[@]}" -gt 0 ]]; then
            log_error "Missing ${#missing[@]} chunk VCF(s) for ${chr} (expected per manifest):"
            printf '  %s\n' "${missing[@]}" >&2
            return 2
        fi
        printf '%s\n' "${lines[@]}"
        return 0
    fi

    if [[ ! -d "${chr_dir}" ]]; then
        return 1
    fi
    local -a files=()
    while IFS= read -r f; do
        files+=( "${f}" )
    done < <(find "${chr_dir}" -maxdepth 1 -type f -name "quilt2.diploid.${chr}.*.vcf.gz")
    if [[ "${#files[@]}" -eq 0 ]]; then
        return 1
    fi
    for f in "${files[@]}"; do
        local start end
        read -r start end < <(basename "${f}" | awk -F'[.-]' '{print $4, $5}')
        printf '%s\t%s\t%s\n' "${start}" "${end}" "${f}"
    done | sort -k1,1n
    return 0
}

# Given ordered "start\tend\tpath" chunk rows for one chromosome, trims each
# chunk's end to stop just before the next chunk's start (no-op if chunks
# don't actually overlap), writing the resulting (possibly-trimmed) VCF paths
# to stdout in order, ready for `bcftools concat --naive -f`.
trim_overlaps() {
    local chr="$1"
    local rows_file="$2"
    local -a starts=() ends=() paths=()
    local start end path

    while IFS=$'\t' read -r start end path; do
        starts+=( "${start}" )
        ends+=( "${end}" )
        paths+=( "${path}" )
    done < "${rows_file}"

    local n="${#paths[@]}"
    local i trimmed_end
    for (( i = 0; i < n; i++ )); do
        trimmed_end="${ends[i]}"
        if (( i + 1 < n )) && (( starts[i+1] - 1 < trimmed_end )); then
            trimmed_end=$(( starts[i+1] - 1 ))
        fi

        if (( trimmed_end < starts[i] )); then
            log_error "Overlap trimming left an empty region for ${chr} chunk ${starts[i]}-${ends[i]}; chunk manifest may be malformed."
            return 1
        fi

        if (( trimmed_end < ends[i] )); then
            local trimmed_path="${TMP_DIR}/${chr}.chunk_$(printf '%03d' "${i}").trimmed.vcf.gz"
            log_info "Trimming overlap for ${chr} chunk ${starts[i]}-${ends[i]} -> ${starts[i]}-${trimmed_end} (next chunk starts at $(( trimmed_end + 1 )))"
            bcftools view -t "${chr}:${starts[i]}-${trimmed_end}" -Oz -o "${trimmed_path}" "${paths[i]}"
            echo "${trimmed_path}"
        else
            echo "${paths[i]}"
        fi
    done
    return 0
}

# =============================================================================
# PER-CHROMOSOME CONCAT
# =============================================================================

declare -a CHR_OUTPUTS=()

for chr in "${CHR_LIST[@]}"; do
    chr_out_dir="${OUT_DIR%/}/${chr}"
    chr_output="${chr_out_dir}/imputed.${chr}.vcf.gz"

    if [[ "${FORCE}" != "true" && -s "${chr_output}" && ( -f "${chr_output}.csi" || -f "${chr_output}.tbi" ) ]]; then
        log_info "[SKIP] Per-chromosome VCF already exists for ${chr}: ${chr_output}"
        CHR_OUTPUTS+=( "${chr_output}" )
        continue
    fi

    rows_file="${TMP_DIR}/${chr}.rows.tsv"
    status=0
    build_chr_filelist "${chr}" > "${rows_file}" || status=$?
    if [[ "${status}" -eq 1 ]]; then
        log_warn "No chunk VCFs found for ${chr}; skipping."
        continue
    elif [[ "${status}" -ne 0 ]]; then
        exit 1
    fi

    n_chunks="$(wc -l < "${rows_file}" | tr -d ' ')"
    filelist_file="${TMP_DIR}/${chr}.filelist"
    trim_overlaps "${chr}" "${rows_file}" > "${filelist_file}" || exit 1

    mkdir -p "${chr_out_dir}"
    log_info "Concatenating ${n_chunks} chunk(s) for ${chr} -> ${chr_output}"
    bcftools concat --naive -Oz -f "${filelist_file}" -o "${chr_output}"
    bcftools index -f -c "${chr_output}"

    CHR_OUTPUTS+=( "${chr_output}" )
done

if [[ "${#CHR_OUTPUTS[@]}" -eq 0 ]]; then
    log_error "No per-chromosome VCFs were produced; nothing to concatenate."
    exit 1
fi

# =============================================================================
# GENOME-WIDE CONCAT
# =============================================================================

if [[ "${#CHR_OUTPUTS[@]}" -eq 1 ]]; then
    FINAL_VCF="${CHR_OUTPUTS[0]}"
else
    FINAL_VCF="${OUT_DIR%/}/imputed.all_chroms.vcf.gz"
    if [[ "${FORCE}" != "true" && -s "${FINAL_VCF}" && ( -f "${FINAL_VCF}.csi" || -f "${FINAL_VCF}.tbi" ) ]]; then
        log_info "[SKIP] Genome-wide VCF already exists: ${FINAL_VCF}"
    else
        genome_filelist="${TMP_DIR}/genome.filelist"
        printf '%s\n' "${CHR_OUTPUTS[@]}" > "${genome_filelist}"
        log_info "Concatenating ${#CHR_OUTPUTS[@]} chromosome VCFs -> ${FINAL_VCF}"
        # Not --naive here: each per-chromosome VCF's header only declares its
        # own contig (inherited from QUILT2.R's single-chromosome chunk output),
        # so headers differ across files. Regular concat merges/unions headers;
        # --naive requires byte-identical headers and rejects this.
        bcftools concat -Oz -f "${genome_filelist}" -o "${FINAL_VCF}"
        bcftools index -f -c "${FINAL_VCF}"
    fi
fi

log_info "Done. Final VCF: ${FINAL_VCF}"
echo "${FINAL_VCF}"
