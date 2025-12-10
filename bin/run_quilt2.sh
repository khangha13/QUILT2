#!/bin/bash
# QUILT2 orchestrator (SLURM array, mirrors Step1C pattern)
set -euo pipefail

ORCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QUILT2_ROOT="$(cd "${ORCH_DIR}/.." && pwd)"
export QUILT2_ROOT

source "${QUILT2_ROOT}/lib/functions.sh"
source "${QUILT2_ROOT}/config/quilt2_config.sh"
load_quilt_env || true

usage() {
    cat <<'EOF'
Usage:
  bash bin/run_quilt2.sh -i <work_dir> --genetic-map <map|dir> [options]

Required:
  -i, --input-dir PATH         Working folder (contains 7.Consolidated_VCF or 8.Imputated_VCF_BEAGLE)
  --genetic-map PATH           Genetic map file or directory with per-chromosome maps

Chunk specification (one of):
  --auto-chunk-map             Use QUILT::quilt_chunk_map (requires R pkg QUILT)
  --chunk-file PATH            TSV with chunk definitions (chr start end [buffer] OR chunk chr start end [buffer])
  --region-start N --region-end N   Single region applied to all chromosomes

Core options:
  --chr LIST                   Comma/space list of chromosomes (default apple Chr01-17)
  --buffer N                   Buffer bp (default 500000)
  --n-gen N                    nGen passed to QUILT2 (default 100)
  --bamlist PATH               BAM list (defaults to <work_dir>/bamlist.txt or bamlist.1.0.txt)
  --reference-panel-dir PATH   Panel VCF dir (defaults to 8.Imputated_VCF_BEAGLE then 7.Consolidated_VCF then work_dir)
  --quilt2-home PATH           Directory containing QUILT2.R and QUILT2_prepare_reference.R
  --quilt2-prepare-script PATH Override QUILT2_prepare_reference.R
  --quilt2-run-script PATH     Override QUILT2.R
  --remove-missing             Filter panel variants by phased rate
  --min-phased-rate FLOAT      Fraction with phased genotypes to keep (default 0.95)
  --prepare-only               Run prepare phase only
  --impute-only                Skip prepare; assumes prepared reference exists
  --dry-run                    Print commands without executing
  --standardise-name           Create/use Chr01-style renamed VCFs (numeric -> ChrNN)
  --standardise-name-force     Force re-run of renaming even if outputs exist
  --truth-vcf PATH             Optional truth VCF for evaluation
  --eval-output PATH           Optional evaluation output directory

SLURM config (env-driven; see config/quilt2_config.sh):
  QUILT2_ACCOUNT, QUILT2_PARTITION, QUILT2_QOS, QUILT2_NODES, QUILT2_NTASKS,
  QUILT2_CPUS_PER_TASK, QUILT2_MEMORY, QUILT2_TIME_LIMIT, QUILT2_ARRAY_MAX
EOF
}

DEFAULT_CHROMS=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17)
DEFAULT_BUFFER=500000
DEFAULT_NGEN=100

ORIG_ARGS=("$@")
INPUT_DIR=""
REFERENCE_PANEL_DIR=""
GENETIC_MAP_FILE="${QUILT2_GENETIC_MAP:-}"
REFERENCE_FASTA="${PIPELINE_REFERENCE_FASTA:-}"
BAMLIST=""
CHROM_ARG=""
REGION_START=1
REGION_END=""
BUFFER="${DEFAULT_BUFFER}"
NGEN="${DEFAULT_NGEN}"
CHUNK_FILE=""
AUTO_CHUNK_MAP="false"
QUILT2_HOME=""
QUILT2_PREP_SCRIPT="${QUILT2_PREP_SCRIPT:-}"
QUILT2_RUN_SCRIPT="${QUILT2_RUN_SCRIPT:-}"
TRUTH_VCF=""
EVAL_OUTPUT_DIR=""
REMOVE_MISSING="false"
MIN_PHASED_RATE="0.95"
PREP_ONLY="false"
IMPUTE_ONLY="false"
DRY_RUN="false"
STANDARDISE_NAME="false"
STANDARDISE_NAME_FORCE="false"
BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"
QUILT2_CONDA_ENV="${QUILT2_CONDA_ENV:-quilt2}"

SUBMIT_SELF="true"
while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input-dir) INPUT_DIR="$2"; shift 2 ;;
        --reference-panel-dir) REFERENCE_PANEL_DIR="$2"; shift 2 ;;
        --genetic-map) GENETIC_MAP_FILE="$2"; shift 2 ;;
        --reference-fasta) REFERENCE_FASTA="$2"; shift 2 ;;
        --bamlist) BAMLIST="$2"; shift 2 ;;
        --chr) CHROM_ARG="$2"; shift 2 ;;
        --region-start) REGION_START="$2"; shift 2 ;;
        --region-end) REGION_END="$2"; shift 2 ;;
        --buffer) BUFFER="$2"; shift 2 ;;
        --n-gen) NGEN="$2"; shift 2 ;;
        --chunk-file) CHUNK_FILE="$2"; shift 2 ;;
        --auto-chunk-map) AUTO_CHUNK_MAP="true"; shift ;;
        --quilt2-home) QUILT2_HOME="$2"; shift 2 ;;
        --quilt2-prepare-script) QUILT2_PREP_SCRIPT="$2"; shift 2 ;;
        --quilt2-run-script) QUILT2_RUN_SCRIPT="$2"; shift 2 ;;
        --truth-vcf) TRUTH_VCF="$2"; shift 2 ;;
        --eval-output) EVAL_OUTPUT_DIR="$2"; shift 2 ;;
        --remove-missing) REMOVE_MISSING="true"; shift ;;
        --min-phased-rate) MIN_PHASED_RATE="$2"; shift 2 ;;
        --prepare-only) PREP_ONLY="true"; shift ;;
        --impute-only) IMPUTE_ONLY="true"; shift ;;
        --dry-run) DRY_RUN="true"; shift ;;
        --standardise-name) STANDARDISE_NAME="true"; shift ;;
        --standardise-name-force) STANDARDISE_NAME_FORCE="true"; shift ;;
        --submit-self) SUBMIT_SELF="$2"; shift 2 ;;
        --no-submit|--submit-self=false) SUBMIT_SELF="false"; shift ;;
        -h|--help) usage; exit 0 ;;
        --) shift; break ;;
        -*)
            log_error "Unknown option: $1"
            usage; exit 1 ;;
        *) break ;;
    esac
done

if [[ -z "${INPUT_DIR}" ]]; then
    log_error "Missing --input-dir"
    usage; exit 1
fi

if [[ "${PREP_ONLY}" == "true" && "${IMPUTE_ONLY}" == "true" ]]; then
    log_error "--prepare-only and --impute-only cannot be used together."
    exit 1
fi

WORK_DIR="$(cd "${INPUT_DIR}" && pwd)"
OUTPUT_DIR="${WORK_DIR%/}/quilt2_output"
PANEL_OUT_DIR="${OUTPUT_DIR}/panel"
RDATA_DIR="${OUTPUT_DIR}/RData"
TMP_DIR="${OUTPUT_DIR}/tmp"
SLURM_DIR="${WORK_DIR%/}/quilt2_slurm"
mkdir -p "${OUTPUT_DIR}" "${PANEL_OUT_DIR}" "${RDATA_DIR}" "${TMP_DIR}" "${SLURM_DIR}"
MISSING_REPORT="${PANEL_OUT_DIR}/missing_sites_removed.tsv"
NOMISS_FAIL_FLAG="${SLURM_DIR}/quilt2_nomiss_failed.flag"

if [[ "${SUBMIT_SELF}" == "true" && -z "${SLURM_JOB_ID:-}" ]]; then
    MASTER_SCRIPT="${SLURM_DIR}/quilt2_master_$(date +%Y%m%d_%H%M%S).sh"
    args_quoted="$(printf " %q" "${ORIG_ARGS[@]}")"
    {
    cat <<EOF
#!/bin/bash
#SBATCH --job-name=Q2_MASTER
#SBATCH --output=${SLURM_DIR}/quilt2_master_%j.output
#SBATCH --error=${SLURM_DIR}/quilt2_master_%j.error

export QUILT2_ROOT="${QUILT2_ROOT}"

bash "${QUILT2_ROOT}/bin/run_quilt2.sh"${args_quoted}
EOF
    } > "${MASTER_SCRIPT}"
    chmod +x "${MASTER_SCRIPT}"

    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "[DRY-RUN] Would submit master: sbatch ${MASTER_SCRIPT}"
        exit 0
    fi

    master_job_id="$(sbatch "${MASTER_SCRIPT}" | awk '{print $4}')"
    if [[ -z "${master_job_id}" ]]; then
        log_error "Failed to submit master job."
        exit 1
    fi
    echo "${master_job_id}" > "${SLURM_DIR}/quilt2_master_job_id.txt"
    log_info "Submitted master job ${master_job_id}"
    log_info "SLURM script: ${MASTER_SCRIPT}"
    log_info "SLURM logs:   ${SLURM_DIR}/quilt2_master_%j.(output|error)"
    exit 0
fi

if [[ -z "${REFERENCE_PANEL_DIR}" ]]; then
    for candidate in "${WORK_DIR}/8.Imputated_VCF_BEAGLE" "${WORK_DIR}/7.Consolidated_VCF" "${WORK_DIR}"; do
        if [[ -d "${candidate}" ]]; then
            REFERENCE_PANEL_DIR="${candidate}"
            break
        fi
    done
fi
if [[ -z "${REFERENCE_PANEL_DIR}" ]]; then
    log_error "Unable to resolve reference panel directory. Set --reference-panel-dir."
    exit 1
fi
REFERENCE_PANEL_DIR="$(cd "${REFERENCE_PANEL_DIR}" && pwd)"

if [[ -z "${GENETIC_MAP_FILE}" ]]; then
    log_error "Genetic map is required (--genetic-map or QUILT2_GENETIC_MAP)."
    exit 1
fi

# Determine if genetic map is a directory
GENETIC_MAP_IS_DIR="false"
if [[ -d "${GENETIC_MAP_FILE}" ]]; then
    GENETIC_MAP_IS_DIR="true"
    GENETIC_MAP_FILE="$(cd "${GENETIC_MAP_FILE}" && pwd)"
elif [[ -f "${GENETIC_MAP_FILE}" ]]; then
    GENETIC_MAP_FILE="$(cd "$(dirname "${GENETIC_MAP_FILE}")" && pwd)/$(basename "${GENETIC_MAP_FILE}")"
else
    log_error "Genetic map not found: ${GENETIC_MAP_FILE}"
    exit 1
fi

if [[ -n "${TRUTH_VCF}" && ! -f "${TRUTH_VCF}" ]]; then
    log_error "Truth VCF not found: ${TRUTH_VCF}"
    exit 1
fi

if [[ -z "${BAMLIST}" && "${PREP_ONLY}" != "true" ]]; then
    for candidate in "${WORK_DIR}/bamlist.txt" "${WORK_DIR}/bamlist.1.0.txt" "${WORK_DIR}/bamlist.tsv"; do
        if [[ -f "${candidate}" ]]; then
            BAMLIST="${candidate}"
            break
        fi
    done
fi
if [[ -z "${BAMLIST}" && "${PREP_ONLY}" != "true" ]]; then
    log_error "BAM list not found. Provide --bamlist or place bamlist.txt in the working folder (required for impute-only and full runs)."
    exit 1
fi

wait_for_slurm_job() {
    local job_id="$1" desc="$2" fail_flag="$3"
    log_info "Waiting for ${desc} job ${job_id} to finish..."
    while true; do
        if command -v squeue >/dev/null 2>&1; then
            queue_out="$(squeue -h -j "${job_id}" 2>/dev/null || true)"
            if [[ -n "${queue_out}" ]]; then
                sleep 20
                continue
            fi
        else
            # No scheduler query available; exit loop and rely on output/flag checks.
            sleep 20
        fi
        break
    done
    if [[ -n "${fail_flag}" && -f "${fail_flag}" ]]; then
        log_error "${desc} failed; see ${fail_flag}"
        exit 1
    fi
}

# SLURM config (shared by both phases)
config="$(get_quilt2_config)"
CFG_ACCOUNT="" CFG_PARTITION="" CFG_QOS="" CFG_NODES="" CFG_NTASKS="" CFG_CPUS="" CFG_MEMORY="" CFG_TIME="" CFG_ARRAY_MAX=""
while IFS='=' read -r k v; do
    case "${k}" in
        ACCOUNT) CFG_ACCOUNT="${v}" ;;
        PARTITION) CFG_PARTITION="${v}" ;;
        QOS) CFG_QOS="${v}" ;;
        NODES) CFG_NODES="${v}" ;;
        NTASKS) CFG_NTASKS="${v}" ;;
        CPUS) CFG_CPUS="${v}" ;;
        MEMORY) CFG_MEMORY="${v}" ;;
        TIME) CFG_TIME="${v}" ;;
        ARRAY_MAX) CFG_ARRAY_MAX="${v}" ;;
    esac
done <<< "${config}"

CHR_LIST=("${DEFAULT_CHROMS[@]}")
if [[ -n "${CHROM_ARG}" ]]; then
    IFS=', ' read -r -a CHR_LIST <<< "${CHROM_ARG}"
fi

# Phase 1: panel prep (standardise and/or remove-missing) as a SLURM array (per chromosome)
RUN_PHASE1="false"
if [[ "${REMOVE_MISSING}" == "true" || "${STANDARDISE_NAME}" == "true" ]]; then
    RUN_PHASE1="true"
fi

if [[ "${RUN_PHASE1}" == "true" ]]; then
    log_info "Phase 1: submitting panel prep array over ${#CHR_LIST[@]} chromosomes (standardise=${STANDARDISE_NAME}, remove_missing=${REMOVE_MISSING})"
    rm -f "${NOMISS_FAIL_FLAG}"

    configured_array_limit="${CFG_ARRAY_MAX:-0}"
    PHASE1_CHR_LIST=( "${CHR_LIST[@]}" )
    if [[ "${configured_array_limit}" -gt 0 && "${#PHASE1_CHR_LIST[@]}" -gt "${configured_array_limit}" ]]; then
        log_warn "Phase 1 chromosome count (${#PHASE1_CHR_LIST[@]}) exceeds array cap (${configured_array_limit}); truncating manifest."
        PHASE1_CHR_LIST=( "${PHASE1_CHR_LIST[@]:0:${configured_array_limit}}" )
    fi

    nomiss_manifest="${SLURM_DIR}/quilt2_nomiss_chr_$(date +%Y%m%d_%H%M%S).txt"
    : > "${nomiss_manifest}"
    for chr in "${PHASE1_CHR_LIST[@]}"; do
        echo "${chr}" >> "${nomiss_manifest}"
    done

    nomiss_array_max=$(( ${#PHASE1_CHR_LIST[@]} - 1 ))

    NOMISS_TEMPLATE="${QUILT2_ROOT}/templates/quilt2_nomiss_job.sh"
    if [[ ! -f "${NOMISS_TEMPLATE}" ]]; then
        log_error "Phase 1 template not found: ${NOMISS_TEMPLATE}"
        exit 1
    fi

    NOMISS_SCRIPT="${SLURM_DIR}/quilt2_nomiss_$(date +%Y%m%d_%H%M%S).sh"
    {
    cat <<EOF
#!/bin/bash
#SBATCH --job-name=Q2NOMISS
EOF
    [[ -n "${CFG_ACCOUNT}" ]]   && echo "#SBATCH --account=${CFG_ACCOUNT}"
    [[ -n "${CFG_PARTITION}" ]] && echo "#SBATCH --partition=${CFG_PARTITION}"
    [[ -n "${CFG_QOS}" ]]       && echo "#SBATCH --qos=${CFG_QOS}"
    cat <<EOF
#SBATCH --nodes=${CFG_NODES}
#SBATCH --ntasks=${CFG_NTASKS}
#SBATCH --cpus-per-task=${CFG_CPUS}
#SBATCH --mem=${CFG_MEMORY}
#SBATCH --time=${CFG_TIME}
#SBATCH --array=0-${nomiss_array_max}
#SBATCH --output=${SLURM_DIR}/quilt2_nomiss_%A_%a.output
#SBATCH --error=${SLURM_DIR}/quilt2_nomiss_%A_%a.error

export QUILT2_ROOT="${QUILT2_ROOT}"
export DRY_RUN="${DRY_RUN}"
export NOMISS_FAIL_FLAG="${NOMISS_FAIL_FLAG}"
export MISSING_REPORT="${MISSING_REPORT}"
export STANDARDISE_NAME="${STANDARDISE_NAME}"
export STANDARDISE_NAME_FORCE="${STANDARDISE_NAME_FORCE}"
export STANDARDISE_SUFFIX="_chr"

bash "${NOMISS_TEMPLATE}" \
  "${WORK_DIR}" \
  "${REFERENCE_PANEL_DIR}" \
  "${PANEL_OUT_DIR}" \
  "${MIN_PHASED_RATE}" \
  "${nomiss_manifest}" \
  "${BCFTOOLS_MODULE}" \
  "${QUILT2_CONDA_ENV}" \
  "${NOMISS_FAIL_FLAG}"
EOF
    } > "${NOMISS_SCRIPT}"
    chmod +x "${NOMISS_SCRIPT}"

    log_info "Generated Phase 1 SLURM script: ${NOMISS_SCRIPT}"
    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "[DRY-RUN] Would submit Phase 1: sbatch ${NOMISS_SCRIPT}"
    else
        nomiss_job_id="$(sbatch "${NOMISS_SCRIPT}" | awk '{print $4}')"
        if [[ -z "${nomiss_job_id}" ]]; then
            log_error "Failed to submit Phase 1 SLURM job."
            exit 1
        fi
        echo "${nomiss_job_id}" > "${SLURM_DIR}/quilt2_nomiss_job_id.txt"
        log_info "Submitted Phase 1 array job ${nomiss_job_id} (script: ${NOMISS_SCRIPT})"
        wait_for_slurm_job "${nomiss_job_id}" "Phase 1 remove-missing" "${NOMISS_FAIL_FLAG}"
    fi

    # Validate Phase 1 outputs before proceeding
    if [[ "${DRY_RUN}" != "true" ]]; then
        for chr in "${PHASE1_CHR_LIST[@]}"; do
            if [[ "${STANDARDISE_NAME}" == "true" ]]; then
                std="${PANEL_OUT_DIR%/}/${chr}_chr.vcf.gz"
                if [[ ! -s "${std}" ]]; then
                    log_error "Phase 1 standardised VCF missing for ${chr}: ${std}"
                    exit 1
                fi
                if [[ ! -f "${std}.csi" && ! -f "${std}.tbi" ]]; then
                    log_error "Index missing for ${chr}: ${std}(.csi|.tbi)"
                    exit 1
                fi
            fi
            if [[ "${REMOVE_MISSING}" == "true" ]]; then
                cleaned="${PANEL_OUT_DIR%/}/quilt.nomiss.${chr}.vcf.gz"
                if [[ ! -s "${cleaned}" ]]; then
                    log_error "Phase 1 filtered VCF missing for ${chr}: ${cleaned}"
                    exit 1
                fi
                if [[ ! -f "${cleaned}.csi" && ! -f "${cleaned}.tbi" ]]; then
                    log_error "Index missing for ${chr}: ${cleaned}(.csi|.tbi)"
                    exit 1
                fi
            fi
        done
        log_info "Phase 1 outputs verified for all chromosomes."
    fi
fi

# Resolve QUILT2 scripts
resolve_quilt2_script() {
    local target="$1"
    local explicit="$2"
    local filename="$3"
    if [[ -n "${explicit}" ]]; then
        echo "${explicit}"
        return
    fi
    if [[ -n "${QUILT2_HOME:-}" && -f "${QUILT2_HOME%/}/${filename}" ]]; then
        echo "${QUILT2_HOME%/}/${filename}"
        return
    fi
    if command -v "${filename}" >/dev/null 2>&1; then
        command -v "${filename}"
        return
    fi
    log_error "Cannot locate ${target} (${filename}). Set --${target//_/-}-script or --quilt2-home."
    exit 1
}

QUILT2_PREP_SCRIPT="$(resolve_quilt2_script quilt2_prepare "${QUILT2_PREP_SCRIPT}" "QUILT2_prepare_reference.R")"
QUILT2_RUN_SCRIPT="$(resolve_quilt2_script quilt2_run "${QUILT2_RUN_SCRIPT}" "QUILT2.R")"

# Chunk assembly (auto, file, or region)
declare -a CHUNKS=()

add_chunk() {
    local chunk_id="$1" chr="$2" start="$3" end="$4" buffer="$5"
    if [[ -z "${chr}" || -z "${start}" || -z "${end}" ]]; then
        log_error "Invalid chunk definition: chr=${chr} start=${start} end=${end}"
        exit 1
    fi
    if ! [[ "${start}" =~ ^[0-9]+$ && "${end}" =~ ^[0-9]+$ ]]; then
        log_error "Chunk ${chunk_id} has non-numeric coordinates: start=${start} end=${end}"
        exit 1
    fi
    CHUNKS+=("${chunk_id}|${chr}|${start}|${end}|${buffer}")
}

if [[ -z "${CHUNK_FILE}" && "${AUTO_CHUNK_MAP}" == "true" ]]; then
    require_cmd Rscript || exit 1
    if ! Rscript -e "quit(status = !requireNamespace('QUILT', quietly = TRUE))"; then
        log_error "--auto-chunk-map requested but R package 'QUILT' is not installed."
        exit 1
    fi
    CHUNK_FILE="${TMP_DIR%/}/quilt_auto_chunks.tsv"
    log_info "Auto-deriving chunks with QUILT::quilt_chunk_map into ${CHUNK_FILE}"

    chr_map_pairs=""
    for chr in "${CHR_LIST[@]}"; do
        map_file="$(resolve_genetic_map "${chr}" "${GENETIC_MAP_IS_DIR}" "${GENETIC_MAP_FILE}" "${GENETIC_MAP_FILE}")" || exit 1
        if [[ -n "${chr_map_pairs}" ]]; then
            chr_map_pairs="${chr_map_pairs};"
        fi
        chr_map_pairs="${chr_map_pairs}${chr}:${map_file}"
    done

    r_cmd=$(cat <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
out_file <- args[[1]]
chr_map_pairs <- args[[2]]
suppressMessages(library(QUILT))
lines <- c()
pairs <- strsplit(chr_map_pairs, ";")[[1]]
for (pair in pairs) {
  parts <- strsplit(pair, ":")[[1]]
  chr <- parts[1]
  map_file <- paste(parts[-1], collapse = ":")
  dat <- QUILT::quilt_chunk_map(chr, map_file)
  if (!all(c("chunk","chr","region") %in% names(dat))) {
    stop("quilt_chunk_map output missing required columns for chr ", chr)
  }
  for (i in seq_len(nrow(dat))) {
    reg <- dat$region[i]
    reg_parts <- unlist(strsplit(reg, "[:-]"))
    if (length(reg_parts) != 3) stop("Unexpected region format: ", reg)
    start_int <- as.integer(as.numeric(reg_parts[2]))
    end_int <- as.integer(as.numeric(reg_parts[3]))
    lines <- c(lines, sprintf("%s\t%s\t%d\t%d", dat$chunk[i], reg_parts[1], start_int, end_int))
  }
}
writeLines(lines, out_file)
RSCRIPT
)
    run_cmd Rscript - "${CHUNK_FILE}" "${chr_map_pairs}" <<< "${r_cmd}"
fi

if [[ -n "${CHUNK_FILE}" ]]; then
    if [[ ! -f "${CHUNK_FILE}" ]]; then
        log_error "Chunk file not found: ${CHUNK_FILE}"
        exit 1
    fi
    while read -r c1 c2 c3 c4 c5; do
        [[ -z "${c1}" ]] && continue
        [[ "${c1}" =~ ^# ]] && continue
        lc1="$(echo "${c1}" | tr '[:upper:]' '[:lower:]')"
        lc2="$(echo "${c2}" | tr '[:upper:]' '[:lower:]')"
        if [[ "${lc1}" == "chr" || "${lc1}" == "chrom" || "${lc1}" == "chunk" ]]; then
            continue
        fi
        if [[ "${lc1}" == chr* ]]; then
            add_chunk "${c1}_${c2}_${c3}" "${c1}" "${c2}" "${c3}" "${c4:-${BUFFER}}"
        elif [[ "${lc2}" == chr* ]]; then
            add_chunk "${c1}" "${c2}" "${c3}" "${c4}" "${c5:-${BUFFER}}"
        elif [[ "${c4}" =~ ^[0-9]+$ && "${c4}" -gt 1000 ]]; then
            add_chunk "${c1}" "${c2}" "${c3}" "${c4}" "${c5:-${BUFFER}}"
        else
            add_chunk "${c1}_${c2}_${c3}" "${c1}" "${c2}" "${c3}" "${c4:-${BUFFER}}"
        fi
    done < "${CHUNK_FILE}"
else
    if [[ -z "${REGION_END}" ]]; then
        log_error "No region specification provided. Use --auto-chunk-map or --chunk-file or --region-end."
        exit 1
    fi
    for chr in "${CHR_LIST[@]}"; do
        add_chunk "${chr}_${REGION_START}_${REGION_END}" "${chr}" "${REGION_START}" "${REGION_END}" "${BUFFER}"
    done
fi

if [[ "${#CHUNKS[@]}" -eq 0 ]]; then
    log_error "No chunks assembled; aborting."
    exit 1
fi

# Persist manifest
MANIFEST_FILE="${TMP_DIR%/}/quilt2_chunks_$(date +%Y%m%d_%H%M%S).txt"
create_chunk_manifest "${MANIFEST_FILE}" "${CHUNKS[@]}" || exit 1

array_max=$(( ${#CHUNKS[@]} - 1 ))
configured_array_limit="${CFG_ARRAY_MAX:-0}"
if [[ "${configured_array_limit}" -gt 0 && "${array_max}" -ge "${configured_array_limit}" ]]; then
    log_warn "Chunk count (${#CHUNKS[@]}) exceeds configured array limit (${configured_array_limit}); truncating."
    array_max=$((configured_array_limit - 1))
fi

SLURM_SCRIPT="${SLURM_DIR}/quilt2_array_$(date +%Y%m%d_%H%M%S).sh"
TEMPLATE="${QUILT2_ROOT}/templates/quilt2_job.sh"
PHASE2_PANEL_DIR="${REFERENCE_PANEL_DIR}"
if [[ "${RUN_PHASE1}" == "true" ]]; then
    PHASE2_PANEL_DIR="${PANEL_OUT_DIR}"
fi

if [[ ! -f "${TEMPLATE}" ]]; then
    log_error "Template not found: ${TEMPLATE}"
    exit 1
fi

{
cat <<EOF
#!/bin/bash
#SBATCH --job-name=QUILT2
EOF
[[ -n "${CFG_ACCOUNT}" ]]   && echo "#SBATCH --account=${CFG_ACCOUNT}"
[[ -n "${CFG_PARTITION}" ]] && echo "#SBATCH --partition=${CFG_PARTITION}"
[[ -n "${CFG_QOS}" ]]       && echo "#SBATCH --qos=${CFG_QOS}"
cat <<EOF
#SBATCH --nodes=${CFG_NODES}
#SBATCH --ntasks=${CFG_NTASKS}
#SBATCH --cpus-per-task=${CFG_CPUS}
#SBATCH --mem=${CFG_MEMORY}
#SBATCH --time=${CFG_TIME}
#SBATCH --array=0-${array_max}
#SBATCH --output=${SLURM_DIR}/quilt2_%A_%a.output
#SBATCH --error=${SLURM_DIR}/quilt2_%A_%a.error

export QUILT2_ROOT="${QUILT2_ROOT}"
export DRY_RUN="${DRY_RUN}"

bash "${TEMPLATE}" \
  "${WORK_DIR}" \
  "${MANIFEST_FILE}" \
  "${PHASE2_PANEL_DIR}" \
  "${GENETIC_MAP_FILE}" \
  "${GENETIC_MAP_IS_DIR}" \
  "${QUILT2_PREP_SCRIPT}" \
  "${QUILT2_RUN_SCRIPT}" \
  "${BAMLIST}" \
  "${NGEN}" \
  "${REMOVE_MISSING}" \
  "${MIN_PHASED_RATE}" \
  "${PREP_ONLY}" \
  "${IMPUTE_ONLY}" \
  "${OUTPUT_DIR}" \
  "${PANEL_OUT_DIR}" \
  "${RDATA_DIR}" \
  "${TMP_DIR}" \
  "${BCFTOOLS_MODULE}" \
  "${QUILT2_CONDA_ENV}" \
  "${TRUTH_VCF}" \
  "${EVAL_OUTPUT_DIR}"
EOF
} > "${SLURM_SCRIPT}"
chmod +x "${SLURM_SCRIPT}"

log_info "Generated SLURM script: ${SLURM_SCRIPT}"
log_info "Chunks: ${#CHUNKS[@]} (array 0-${array_max})"

if [[ "${DRY_RUN}" == "true" ]]; then
    echo "[DRY-RUN] Would submit: sbatch ${SLURM_SCRIPT}"
    exit 0
fi

job_id="$(sbatch "${SLURM_SCRIPT}" | awk '{print $4}')"
if [[ -z "${job_id}" ]]; then
    log_error "Failed to submit SLURM job."
    exit 1
fi

echo "${job_id}" > "${SLURM_DIR}/quilt2_job_id.txt"
log_info "Submitted QUILT2 array job ${job_id} (script: ${SLURM_SCRIPT})"
