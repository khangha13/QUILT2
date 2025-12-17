#!/bin/bash
# QUILT2 array worker: processes one chunk per SLURM array task
set -euo pipefail

# The orchestrator substitutes this placeholder with an absolute path.
QUILT2_ROOT="${QUILT2_ROOT:-__QUILT2_ROOT__}"

if [ ! -d "${QUILT2_ROOT}" ]; then
    echo "[quilt2_job] ERROR: QUILT2_ROOT does not exist: ${QUILT2_ROOT}" >&2
    exit 1
fi

source "${QUILT2_ROOT}/lib/functions.sh"

if [ "$#" -lt 17 ]; then
    log_error "Usage: quilt2_job.sh <WORK_DIR> <CHUNK_MANIFEST> <REFERENCE_PANEL_DIR> <GENETIC_MAP> <GENETIC_MAP_IS_DIR> <QUILT2_PREP_SCRIPT> <QUILT2_RUN_SCRIPT> <BAMLIST> <NGEN> <REMOVE_MISSING> <MIN_VALID_GT_RATE> <PREP_ONLY> <IMPUTE_ONLY> <OUTPUT_DIR> <PANEL_OUT_DIR> <RDATA_DIR> <TMP_DIR> [BCFTOOLS_MODULE] [QUILT2_CONDA_ENV]"
    exit 1
fi

WORK_DIR="$1"
CHUNK_MANIFEST="$2"
REFERENCE_PANEL_DIR="$3"
GENETIC_MAP_INPUT="$4"
GENETIC_MAP_IS_DIR="$5"
QUILT2_PREP_SCRIPT="$6"
QUILT2_RUN_SCRIPT="$7"
BAMLIST="$8"
NGEN="$9"
REMOVE_MISSING="${10}"
MIN_VALID_GT_RATE="${11}"
PREP_ONLY="${12}"
IMPUTE_ONLY="${13}"
OUTPUT_DIR="${14}"
PANEL_OUT_DIR="${15}"
RDATA_DIR="${16}"
TMP_DIR="${17}"
BCFTOOLS_MODULE="${18:-${BCFTOOLS_MODULE:-bcftools/1.18-GCC-12.3.0}}"
QUILT2_CONDA_ENV="${19:-${QUILT2_CONDA_ENV:-quilt2}}"

# Optional evaluation inputs (may be blank)
TRUTH_VCF="${20:-}"
EVAL_OUTPUT_DIR="${21:-}"

# Export flags for helper functions
export REMOVE_MISSING MIN_VALID_GT_RATE PANEL_OUT_DIR
export MISSING_REPORT="${PANEL_OUT_DIR%/}/missing_sites_removed.tsv"
export CHUNK_FILE="" # manifest-driven; no single chunk file

# Validate array context
if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    log_error "SLURM_ARRAY_TASK_ID is not set; this script must run as a SLURM array job."
    exit 1
fi

if [ ! -f "${CHUNK_MANIFEST}" ]; then
    log_error "Chunk manifest not found: ${CHUNK_MANIFEST}"
    exit 1
fi

CHUNKS=()
while IFS= read -r line; do
    CHUNKS+=( "${line}" )
done < "${CHUNK_MANIFEST}"
chunk_count="${#CHUNKS[@]}"
if ! [[ "${SLURM_ARRAY_TASK_ID}" =~ ^[0-9]+$ ]] || [ "${SLURM_ARRAY_TASK_ID}" -lt 0 ] || [ "${SLURM_ARRAY_TASK_ID}" -ge "${chunk_count}" ]; then
    log_error "Array index ${SLURM_ARRAY_TASK_ID} out of range 0..$((chunk_count-1))"
    exit 1
fi

CHUNK="${CHUNKS[$SLURM_ARRAY_TASK_ID]}"
IFS='|' read -r CHUNK_ID CHR START END BUFFER <<< "${CHUNK}"

log_info "Array task ${SLURM_ARRAY_TASK_ID}/${chunk_count} processing chunk ${CHUNK_ID} (${CHR}:${START}-${END}, buffer=${BUFFER})"

# Prepare environment and toolchain
export BCFTOOLS_MODULE QUILT2_CONDA_ENV
load_quilt_env
ensure_bcftools || exit 1
require_cmd Rscript || exit 1

# Resolve genetic map for this chromosome
if [[ "${GENETIC_MAP_IS_DIR}" == "true" ]]; then
    GENETIC_MAP_DIR="${GENETIC_MAP_INPUT}"
else
    GENETIC_MAP_DIR="$(dirname "${GENETIC_MAP_INPUT}")"
fi
GENETIC_MAP_FILE="$(resolve_genetic_map "${CHR}" "${GENETIC_MAP_IS_DIR}" "${GENETIC_MAP_INPUT}" "${GENETIC_MAP_DIR}")" || exit 1

# Ensure output dirs exist
mkdir -p "${OUTPUT_DIR}" "${PANEL_OUT_DIR}" "${RDATA_DIR}" "${TMP_DIR}"

# Panel VCF selection / optional missingness filtering
panel_vcf="$(normalize_panel_vcf "${CHR}" "${REFERENCE_PANEL_DIR}")" || exit 1

# Validate BAM list (always required for impute; not needed only when PREP_ONLY)
if [[ "${PREP_ONLY}" != "true" ]]; then
    if [[ -z "${BAMLIST}" || ! -f "${BAMLIST}" ]]; then
        log_error "BAM list is required for imputation; not found: ${BAMLIST:-<empty>}"
        exit 1
    fi
fi

# Prepare reference for the chunk unless in impute-only mode
prepare_reference_chunk() {
    local chr="$1"
    local start="$2"
    local end="$3"
    local buffer="$4"
    local panel="$5"

    local prepared_file="${RDATA_DIR%/}/QUILT_prepared_reference.${chr}.${start}.${end}.RData"
    if [[ "${IMPUTE_ONLY}" == "true" ]]; then
        echo "${prepared_file}"
        return 0
    fi

    if [[ -f "${prepared_file}" ]]; then
        log_info "Prepared reference already exists for ${chr}:${start}-${end}; skipping."
        echo "${prepared_file}"
        return 0
    fi

    log_info "Preparing reference for ${chr}:${start}-${end}"
    local cmd=(Rscript "${QUILT2_PREP_SCRIPT}"
        "--genetic_map_file=${GENETIC_MAP_FILE}"
        "--reference_vcf_file=${panel}"
        "--chr=${chr}"
        "--regionStart=${start}"
        "--regionEnd=${end}"
        "--nGen=${NGEN}"
        "--buffer=${buffer}"
        "--outputdir=${OUTPUT_DIR}"
    )

    run_cmd "${cmd[@]}"
    echo "${prepared_file}"
}

impute_chunk() {
    local chr="$1"
    local start="$2"
    local end="$3"
    local buffer="$4"
    local prepared_file="$5"

    if [[ "${DRY_RUN:-false}" != "true" && ! -f "${prepared_file}" ]]; then
        log_error "Prepared reference missing for ${chr}:${start}-${end}: ${prepared_file}"
        return 1
    fi

    local output_vcf="${OUTPUT_DIR%/}/quilt2.diploid.${chr}.${start}-${end}.vcf.gz"
    if [[ -f "${output_vcf}" && "${PREP_ONLY}" == "false" ]]; then
        log_info "Imputation output exists for ${chr}:${start}-${end}; skipping."
        return 0
    fi

    log_info "Running QUILT2 (diploid SNP) for ${chr}:${start}-${end}"
    local cmd=(Rscript "${QUILT2_RUN_SCRIPT}"
        "--prepared_reference_filename=${prepared_file}"
        "--bamlist=${BAMLIST}"
        "--chr=${chr}"
        "--regionStart=${start}"
        "--regionEnd=${end}"
        "--nGen=${NGEN}"
        "--buffer=${buffer}"
        "--output_filename=${output_vcf}"
    )

    run_cmd "${cmd[@]}"
    echo "${output_vcf}"
}

prepared_file="$(prepare_reference_chunk "${CHR}" "${START}" "${END}" "${BUFFER}" "${panel_vcf}")" || exit 1

if [[ "${PREP_ONLY}" != "true" ]]; then
    imputed_vcf="$(impute_chunk "${CHR}" "${START}" "${END}" "${BUFFER}" "${prepared_file}")" || exit 1

    # Optional per-chunk evaluation (only if truth provided)
    if [[ -n "${TRUTH_VCF}" && -f "${TRUTH_VCF}" ]]; then
        eval_out_dir="${EVAL_OUTPUT_DIR:-${OUTPUT_DIR%/}/eval}"
        mkdir -p "${eval_out_dir}"
        base="$(basename "${imputed_vcf%.vcf.gz}")"
        out_prefix="${eval_out_dir%/}/${base}"
        log_info "Evaluating ${imputed_vcf} against truth ${TRUTH_VCF}"
        if [[ "${DRY_RUN:-false}" == "true" ]]; then
            cat <<EOF
+ Rscript -e "library(vcfppR); res <- vcfcomp(test='${imputed_vcf}', truth='${TRUTH_VCF}', stats='r2', formats=c('DS','GT')); write.table(res, file='${out_prefix}.r2.tsv', sep='\t', quote=FALSE, row.names=FALSE); pdf('${out_prefix}.r2.pdf'); vcfplot(res, col=2, cex=1.2, lwd=2, type='b'); dev.off()"
EOF
        else
            Rscript - "${imputed_vcf}" "${TRUTH_VCF}" "${out_prefix}" <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
imputed <- args[[1]]
truth <- args[[2]]
out_prefix <- args[[3]]
if (!requireNamespace("vcfppR", quietly = TRUE)) {
  stop("vcfppR is required for evaluation. Install it before running evaluation.")
}
library(vcfppR)
res <- vcfcomp(test = imputed, truth = truth,
               stats = "r2",
               region = NULL,
               formats = c("DS","GT"))
utils::write.table(res, file = paste0(out_prefix, ".r2.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
grDevices::pdf(paste0(out_prefix, ".r2.pdf"))
vcfplot(res, col = 2, cex = 1.2, lwd = 2, type = "b")
grDevices::dev.off()
RSCRIPT
        fi
    fi
fi

log_info "Completed chunk ${CHUNK_ID} (${CHR}:${START}-${END})"
