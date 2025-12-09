#!/bin/bash
# QUILT2 orchestrator (SLURM array, mirrors Step1C pattern)
set -euo pipefail

ORCH_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
QUILT2_ROOT="$(cd "${ORCH_DIR}/.." && pwd)"
export QUILT2_ROOT

source "${QUILT2_ROOT}/lib/functions.sh"
source "${QUILT2_ROOT}/config/quilt2_config.sh"

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
BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"
QUILT2_CONDA_ENV="${QUILT2_CONDA_ENV:-quilt2}"

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
mkdir -p "${OUTPUT_DIR}" "${PANEL_OUT_DIR}" "${RDATA_DIR}" "${TMP_DIR}"
MISSING_REPORT="${PANEL_OUT_DIR}/missing_sites_removed.tsv"

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

CHR_LIST=("${DEFAULT_CHROMS[@]}")
if [[ -n "${CHROM_ARG}" ]]; then
    IFS=', ' read -r -a CHR_LIST <<< "${CHROM_ARG}"
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

# Build SLURM script from template
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

array_max=$(( ${#CHUNKS[@]} - 1 ))
configured_array_limit="${CFG_ARRAY_MAX:-0}"
if [[ "${configured_array_limit}" -gt 0 && "${array_max}" -ge "${configured_array_limit}" ]]; then
    log_warn "Chunk count (${#CHUNKS[@]}) exceeds configured array limit (${configured_array_limit}); truncating."
    array_max=$((configured_array_limit - 1))
fi

SLURM_DIR="${WORK_DIR%/}/quilt2_slurm"
mkdir -p "${SLURM_DIR}"
SLURM_SCRIPT="${SLURM_DIR}/quilt2_array_$(date +%Y%m%d_%H%M%S).sh"
TEMPLATE="${QUILT2_ROOT}/templates/quilt2_job.sh"

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
  "${REFERENCE_PANEL_DIR}" \
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
