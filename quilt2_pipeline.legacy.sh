#!/bin/bash
set -euo pipefail

# QUILT2 pipeline wrapper for apple data. Uses the provided working folder as
# the hub for inputs/outputs, reusing GATK pipeline defaults where possible.
#
# NOTE: This pipeline currently runs interactively (no SLURM orchestration).
#       SLURM array parallelization will be added once the pipeline enters
#       production. All outputs are written to WORK_DIR/quilt2_output/.
#
# QUILT2 does NOT require a reference FASTA - it uses the reference panel VCF
# directly for imputation. See: https://github.com/rwdavies/QUILT/blob/master/README_QUILT2.org

# Reference FASTA used only for chromosome validation via .fai (not passed to QUILT2)
DEFAULT_APPLE_REF="/QRISdata/Q8367/Reference_Genome/GDDH13_1-1_formatted.fasta"
DEFAULT_CHROMS=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17)
DEFAULT_BUFFER=500000
DEFAULT_NGEN=100

log_info() { echo "[INFO] $*" >&2; }
log_warn() { echo "[WARN] $*" >&2; }
log_error() { echo "[ERROR] $*" >&2; }

usage() {
    cat <<'EOF'
Usage:
  bash quilt2_pipeline.sh -i <work_dir> --genetic-map <map.txt> --auto-chunk-map [options]
  bash quilt2_pipeline.sh -i <work_dir> --genetic-map <map.txt> --chunk-file <chunks.tsv> [options]
  bash quilt2_pipeline.sh -i <work_dir> --genetic-map <map.txt> --region-start <start> --region-end <end> [options]

Required:
  -i, --input-dir PATH        Working folder (all outputs stay here).
  --genetic-map PATH          Genetic map file OR directory containing per-chromosome maps.
                              If a directory, expects files like: Chr01.map, genetic_map.Chr01.txt, etc.
                              Can also set via QUILT2_GENETIC_MAP environment variable.

Region specification (one of the following):
  --auto-chunk-map            Derive chunks automatically using QUILT::quilt_chunk_map (requires QUILT R pkg).
  --chunk-file PATH           TSV with either: chr start end [buffer] OR chunk chr start end [buffer].
  --region-start N            Region start (bp); defaults to 1. Requires --region-end if no chunk-file/auto-chunk.
  --region-end N              Region end (bp). Required if neither --chunk-file nor --auto-chunk-map is provided.

Core options:
  --chr LIST                  Comma/space list of chromosomes (default apple Chr01-17).
  --buffer N                  Buffer in bp passed to QUILT2 (default 500000).
  --n-gen N                   nGen passed to QUILT2 (default 100).
  --bamlist PATH              BAM list (defaults to <work_dir>/bamlist.txt or bamlist.1.0.txt).

Reference/panel options:
  --reference-fasta PATH      Reference FASTA for chromosome validation only (not used by QUILT2).
                              Defaults to PIPELINE_REFERENCE_FASTA or apple GDDH13 path.
  --reference-panel-dir PATH  Directory with phased panel VCFs; defaults to
                              <work_dir>/8.Imputated_VCF_BEAGLE, then 7.Consolidated_VCF, then <work_dir>.
  --quilt2-home PATH          Directory containing QUILT2.R and QUILT2_prepare_reference.R.
  --quilt2-prepare-script PATH  Override path to QUILT2_prepare_reference.R.
  --quilt2-run-script PATH      Override path to QUILT2.R.
  --truth-vcf PATH            Optional truth VCF for evaluation (vcfppR); contigs must match outputs.
  --eval-output PATH          Optional evaluation output directory (default: <work_dir>/quilt2_output/eval).
  --remove-missing            Filter variants: keep only those with ≥95% phased genotypes (default: false).
                              Use --min-phased-rate to adjust the threshold.
  --min-phased-rate FLOAT     Minimum fraction of samples with phased genotypes (0.0-1.0, default: 0.95).
                              Variants below this threshold are removed when --remove-missing is set.
  --prepare-only              Stop after QUILT2_prepare_reference.R.
  --impute-only               Assume prepared reference exists; skip prepare step.

Other:
  --dry-run                   Print commands without executing.
  Evaluation uses Rscript + vcfppR; ensure those are available in your environment.
  -h, --help                  Show this message.

Chunk TSV hint: a header is ignored. Columns may be:
  chr start end [buffer]           (buffer falls back to --buffer)
  chunk chr start end [buffer]     (chunk name is only used for logging)
EOF
}

INPUT_DIR=""
REFERENCE_PANEL_DIR=""
GENETIC_MAP_FILE="${QUILT2_GENETIC_MAP:-}"
REFERENCE_FASTA="${PIPELINE_REFERENCE_FASTA:-${DEFAULT_APPLE_REF}}"
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
QUILT2_CONDA_ENV="${QUILT2_CONDA_ENV:-quilt2}"
MISSING_REPORT=""
REMOVE_MISSING="false"
MIN_PHASED_RATE="0.95"
PREP_ONLY="false"
IMPUTE_ONLY="false"
DRY_RUN="false"
BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--input-dir)
            INPUT_DIR="$2"; shift 2 ;;
        --reference-panel-dir)
            REFERENCE_PANEL_DIR="$2"; shift 2 ;;
        --genetic-map)
            GENETIC_MAP_FILE="$2"; shift 2 ;;
        --reference-fasta)
            REFERENCE_FASTA="$2"; shift 2 ;;
        --bamlist)
            BAMLIST="$2"; shift 2 ;;
        --chr)
            CHROM_ARG="$2"; shift 2 ;;
        --region-start)
            REGION_START="$2"; shift 2 ;;
        --region-end)
            REGION_END="$2"; shift 2 ;;
        --buffer)
            BUFFER="$2"; shift 2 ;;
        --n-gen)
            NGEN="$2"; shift 2 ;;
        --chunk-file)
            CHUNK_FILE="$2"; shift 2 ;;
        --auto-chunk-map)
            AUTO_CHUNK_MAP="true"; shift ;;
        --quilt2-home)
            QUILT2_HOME="$2"; shift 2 ;;
        --quilt2-prepare-script)
            QUILT2_PREP_SCRIPT="$2"; shift 2 ;;
        --quilt2-run-script)
            QUILT2_RUN_SCRIPT="$2"; shift 2 ;;
        --truth-vcf)
            TRUTH_VCF="$2"; shift 2 ;;
        --eval-output)
            EVAL_OUTPUT_DIR="$2"; shift 2 ;;
        --remove-missing)
            REMOVE_MISSING="true"; shift ;;
        --min-phased-rate)
            MIN_PHASED_RATE="$2"; shift 2 ;;
        --prepare-only)
            PREP_ONLY="true"; shift ;;
        --impute-only)
            IMPUTE_ONLY="true"; shift ;;
        --dry-run)
            DRY_RUN="true"; shift ;;
        -h|--help)
            usage; exit 0 ;;
        --)
            shift; break ;;
        -*)
            log_error "Unknown option: $1"
            usage
            exit 1 ;;
        *)
            break ;;
    esac
done

if [[ -z "${INPUT_DIR}" ]]; then
    log_error "Missing --input-dir"
    usage
    exit 1
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
if [[ -z "${EVAL_OUTPUT_DIR}" ]]; then
    EVAL_OUTPUT_DIR="${OUTPUT_DIR}/eval"
fi
MISSING_REPORT="${PANEL_OUT_DIR}/missing_sites_removed.tsv"

mkdir -p "${OUTPUT_DIR}" "${PANEL_OUT_DIR}" "${RDATA_DIR}" "${TMP_DIR}"

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
    log_error "Genetic map is required (--genetic-map or set QUILT2_GENETIC_MAP)."
    exit 1
fi

# Determine if genetic map is a file or directory
GENETIC_MAP_IS_DIR="false"
if [[ -d "${GENETIC_MAP_FILE}" ]]; then
    GENETIC_MAP_IS_DIR="true"
    GENETIC_MAP_DIR="$(cd "${GENETIC_MAP_FILE}" && pwd)"
    log_info "Genetic map directory: ${GENETIC_MAP_DIR}"
elif [[ -f "${GENETIC_MAP_FILE}" ]]; then
    GENETIC_MAP_DIR=""
    log_info "Genetic map file: ${GENETIC_MAP_FILE}"
else
    log_error "Genetic map not found: ${GENETIC_MAP_FILE} (expected file or directory)"
    exit 1
fi

# REFERENCE_FASTA is optional - used only for chromosome validation via .fai
# QUILT2 does not require it (uses reference panel VCF directly)

if [[ -n "${TRUTH_VCF}" && ! -f "${TRUTH_VCF}" ]]; then
    log_error "Truth VCF not found: ${TRUTH_VCF}"
    exit 1
fi

if [[ -z "${BAMLIST}" ]]; then
    for candidate in "${WORK_DIR}/bamlist.txt" "${WORK_DIR}/bamlist.1.0.txt" "${WORK_DIR}/bamlist.tsv"; do
        if [[ -f "${candidate}" ]]; then
            BAMLIST="${candidate}"
            break
        fi
    done
fi

if [[ "${IMPUTE_ONLY}" == "false" && "${PREP_ONLY}" == "false" ]] || [[ "${IMPUTE_ONLY}" == "true" ]]; then
    if [[ -z "${BAMLIST}" ]]; then
        log_error "BAM list not found. Provide --bamlist or place bamlist.txt in the working folder."
        exit 1
    fi
    if [[ ! -f "${BAMLIST}" ]]; then
        log_error "BAM list file does not exist: ${BAMLIST}"
        exit 1
    fi
fi

require_cmd() {
    local cmd="$1"
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        log_error "Command not found: ${cmd}"
        exit 1
    fi
}

ensure_bcftools() {
    if command -v bcftools >/dev/null 2>&1; then
        return
    fi

    if command -v module >/dev/null 2>&1; then
        if [[ -n "${BCFTOOLS_MODULE}" ]]; then
            if module load "${BCFTOOLS_MODULE}" >/dev/null 2>&1; then
                log_info "Loaded ${BCFTOOLS_MODULE} module for bcftools"
            else
                log_warn "Failed to load ${BCFTOOLS_MODULE}; bcftools still unavailable"
            fi
        fi
    else
        log_warn "module command not found; cannot auto-load bcftools."
    fi

    if ! command -v bcftools >/dev/null 2>&1; then
        log_error "bcftools not found in PATH. Install it or load module ${BCFTOOLS_MODULE}."
        exit 1
    fi
}

load_quilt_env() {
    local conda_env="${QUILT2_CONDA_ENV}"
    local conda_sourced="false"

    if command -v module >/dev/null 2>&1; then
        if module load miniforge/25.3.0-3 >/dev/null 2>&1; then
            log_info "Loaded miniforge/25.3.0-3 module"
        else
            log_error "Failed to load miniforge/25.3.0-3 module"
            exit 1
        fi
    else
        log_warn "module command not found; skipping module load."
    fi

    if [ -n "${ROOTMINIFORGE:-}" ] && [ -f "${ROOTMINIFORGE}/etc/profile.d/conda.sh" ]; then
        # shellcheck source=/dev/null
        source "${ROOTMINIFORGE}/etc/profile.d/conda.sh"
        conda_sourced="true"
    elif command -v conda >/dev/null 2>&1; then
        local conda_base
        conda_base="$(conda info --base 2>/dev/null || true)"
        if [ -n "${conda_base}" ] && [ -f "${conda_base}/etc/profile.d/conda.sh" ]; then
            # shellcheck source=/dev/null
            source "${conda_base}/etc/profile.d/conda.sh"
            conda_sourced="true"
        fi
    fi

    if [ "${conda_sourced}" != "true" ]; then
        log_warn "Conda init script not found; assuming required tools are already on PATH"
        return
    fi

    if conda activate "${conda_env}" >/dev/null 2>&1; then
        log_info "Activated conda environment: ${conda_env}"
    else
        log_error "Failed to activate conda environment: ${conda_env}"
        exit 1
    fi
}

if command -v module >/dev/null 2>&1; then
    module purge
fi

load_quilt_env

ensure_bcftools

for tool in Rscript bcftools; do
    require_cmd "${tool}"
done

resolve_quilt2_script() {
    local target="$1"
    local explicit_path="$2"
    local filename="$3"

    if [[ -n "${explicit_path}" ]]; then
        echo "${explicit_path}"
        return
    fi

    if [[ -n "${QUILT2_HOME:-}" ]] && [[ -f "${QUILT2_HOME%/}/${filename}" ]]; then
        echo "${QUILT2_HOME%/}/${filename}"
        return
    fi

    if command -v "${filename}" >/dev/null 2>&1; then
        command -v "${filename}"
        return
    fi

    log_error "Cannot locate ${target}. Set --${target//_/-}-script or --quilt2-home."
    exit 1
}

QUILT2_PREP_SCRIPT="$(resolve_quilt2_script quilt2_prepare "${QUILT2_PREP_SCRIPT}" "QUILT2_prepare_reference.R")"
QUILT2_RUN_SCRIPT="$(resolve_quilt2_script quilt2_run "${QUILT2_RUN_SCRIPT}" "QUILT2.R")"

if [[ ! -f "${QUILT2_PREP_SCRIPT}" ]]; then
    log_error "QUILT2_prepare_reference.R not found at ${QUILT2_PREP_SCRIPT}"
    exit 1
fi
if [[ ! -f "${QUILT2_RUN_SCRIPT}" ]]; then
    log_error "QUILT2.R not found at ${QUILT2_RUN_SCRIPT}"
    exit 1
fi

# Resolve genetic map file for a specific chromosome.
# If GENETIC_MAP_IS_DIR=true, searches for chr-specific file in directory.
# Otherwise returns the single GENETIC_MAP_FILE.
resolve_genetic_map() {
    local chr="$1"
    if [[ "${GENETIC_MAP_IS_DIR}" == "true" ]]; then
        # Search for per-chromosome map file with common naming patterns
        local -a patterns=(
            "${GENETIC_MAP_DIR}/${chr}.map"
            "${GENETIC_MAP_DIR}/${chr}.txt"
            "${GENETIC_MAP_DIR}/${chr}.txt.gz"
            "${GENETIC_MAP_DIR}/genetic_map.${chr}.txt"
            "${GENETIC_MAP_DIR}/genetic_map.${chr}.txt.gz"
            "${GENETIC_MAP_DIR}/genetic_map_${chr}.txt"
            "${GENETIC_MAP_DIR}/genetic_map_${chr}.txt.gz"
            "${GENETIC_MAP_DIR}/${chr}_genetic_map.txt"
            "${GENETIC_MAP_DIR}/${chr}_genetic_map.txt.gz"
        )
        # Also try glob for any file containing the chromosome name
        local -a globbed=()
        while IFS= read -r path; do
            globbed+=( "${path}" )
        done < <(find "${GENETIC_MAP_DIR}" -maxdepth 1 -type f \( -name "*${chr}*.txt" -o -name "*${chr}*.txt.gz" -o -name "*${chr}*.map" \) 2>/dev/null | head -1)
        patterns+=( "${globbed[@]}" )

        for candidate in "${patterns[@]}"; do
            if [[ -f "${candidate}" ]]; then
                echo "${candidate}"
                return 0
            fi
        done
        log_error "Genetic map for ${chr} not found in ${GENETIC_MAP_DIR}"
        exit 1
    else
        echo "${GENETIC_MAP_FILE}"
    fi
}

readarray -t CHR_LIST <<< "$(printf "%s\n" "${DEFAULT_CHROMS[@]}")"
if [[ -n "${CHROM_ARG}" ]]; then
    # Accept comma or space separated values
    IFS=', ' read -r -a CHR_LIST <<< "${CHROM_ARG}"
fi

# Validate chromosomes against reference .fai if available (per GATK §5.4 lessons)
# REFERENCE_FASTA is optional; only used for this validation, not passed to QUILT2
if [[ -n "${REFERENCE_FASTA}" && -f "${REFERENCE_FASTA}.fai" ]]; then
    for chr in "${CHR_LIST[@]}"; do
        if ! awk -v c="${chr}" '$1==c{found=1; exit} END{exit !found}' "${REFERENCE_FASTA}.fai"; then
            log_warn "Chromosome '${chr}' not found in ${REFERENCE_FASTA}.fai; ensure naming matches reference."
        fi
    done
elif [[ -n "${REFERENCE_FASTA}" ]]; then
    log_warn "Reference FASTA index (.fai) not found at ${REFERENCE_FASTA}.fai; skipping chromosome validation."
fi

# Define run_cmd early so it's available for auto-chunk-map and later stages
run_cmd() {
    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "+ $*"
    else
        "$@"
    fi
}

declare -a CHUNKS=()

add_chunk() {
    local chunk_id="$1"
    local chr="$2"
    local start="$3"
    local end="$4"
    local buffer="$5"

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
    # Derive chunks using QUILT::quilt_chunk_map and the provided genetic map.
    if ! Rscript -e "quit(status = !requireNamespace('QUILT', quietly = TRUE))"; then
        log_error "--auto-chunk-map requested but the R package 'QUILT' is not installed."
        exit 1
    fi
    CHUNK_FILE="${TMP_DIR%/}/quilt_auto_chunks.tsv"
    log_info "Auto-deriving chunks with QUILT::quilt_chunk_map into ${CHUNK_FILE}"

    # Build a comma-separated list of chr:mapfile pairs for the R script
    chr_map_pairs=""
    for chr in "${CHR_LIST[@]}"; do
        chr_map="$(resolve_genetic_map "${chr}")"
        if [[ -n "${chr_map_pairs}" ]]; then
            chr_map_pairs="${chr_map_pairs};"
        fi
        chr_map_pairs="${chr_map_pairs}${chr}:${chr_map}"
    done

    r_cmd=$(cat <<'RSCRIPT'
args <- commandArgs(trailingOnly = TRUE)
out_file <- args[[1]]
chr_map_pairs <- args[[2]]
suppressMessages(library(QUILT))
lines <- c()
# Parse chr:mapfile pairs separated by ;
pairs <- strsplit(chr_map_pairs, ";")[[1]]
for (pair in pairs) {
  parts <- strsplit(pair, ":")[[1]]
  chr <- parts[1]
  map_file <- paste(parts[-1], collapse = ":")  # Handle paths with colons
  dat <- QUILT::quilt_chunk_map(chr, map_file)
  if (!all(c("chunk","chr","region") %in% names(dat))) {
    stop("quilt_chunk_map output missing required columns for chr ", chr)
  }
  for (i in seq_len(nrow(dat))) {
    reg <- dat$region[i]
    # region format like "chr:start-end"
    reg_parts <- unlist(strsplit(reg, "[:-]"))
    if (length(reg_parts) != 3) stop("Unexpected region format: ", reg)
    # Convert coordinates to integers and format explicitly to avoid scientific notation
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
        # Detect format based on column patterns:
        # Format A: "chr start end [buffer]" - c1 is chromosome, c2/c3 are coordinates
        # Format B: "chunk chr start end [buffer]" - c1 is chunk ID, c2 is chromosome, c3/c4 are coordinates
        #
        # Heuristic: if c1 starts with 'chr' or c2 is not a valid coordinate (not purely numeric or too small),
        # treat as Format A. Otherwise, if c4 looks like a coordinate (large number), treat as Format B.
        if [[ "${lc1}" == chr* ]]; then
            # Format A: chr start end [buffer] (chromosome starts with 'chr')
            add_chunk "${c1}_${c2}_${c3}" "${c1}" "${c2}" "${c3}" "${c4:-${BUFFER}}"
        elif [[ "${lc2}" == chr* ]]; then
            # Format B: chunk chr start end [buffer] (chromosome starts with 'chr')
            add_chunk "${c1}" "${c2}" "${c3}" "${c4}" "${c5:-${BUFFER}}"
        elif [[ "${c4}" =~ ^[0-9]+$ && "${c4}" -gt 1000 ]]; then
            # Format B: chunk chr start end [buffer] (c4 looks like a coordinate - numeric and large)
            # This handles numeric chromosome names like "1", "2", etc.
            add_chunk "${c1}" "${c2}" "${c3}" "${c4}" "${c5:-${BUFFER}}"
        else
            # Fallback: assume chr start end [buffer] with non-standard chr naming
            add_chunk "${c1}_${c2}_${c3}" "${c1}" "${c2}" "${c3}" "${c4:-${BUFFER}}"
        fi
    done < "${CHUNK_FILE}"
else
    if [[ -z "${REGION_END}" ]]; then
        log_error "No region specification provided. Use one of:"
        log_error "  --auto-chunk-map           (derive chunks from genetic map)"
        log_error "  --chunk-file <file.tsv>    (provide explicit chunks)"
        log_error "  --region-end <bp>          (use same region for all chromosomes)"
        exit 1
    fi
    for chr in "${CHR_LIST[@]}"; do
        add_chunk "${chr}_${REGION_START}_${REGION_END}" "${chr}" "${REGION_START}" "${REGION_END}" "${BUFFER}"
    done
fi

run_evaluation() {
    local truth_vcf="$1"; shift
    local -a imputed_vcfs=("$@")

    if [[ -z "${truth_vcf}" ]]; then
        return
    fi
    if [[ "${#imputed_vcfs[@]}" -eq 0 ]]; then
        log_warn "No imputed VCFs to evaluate."
        return
    fi

    mkdir -p "${EVAL_OUTPUT_DIR}"

    for vcf in "${imputed_vcfs[@]}"; do
        [[ -f "${vcf}" ]] || { log_warn "Skipping evaluation, missing VCF: ${vcf}"; continue; }
        local base
        base="$(basename "${vcf%.vcf.gz}")"
        local out_prefix="${EVAL_OUTPUT_DIR%/}/${base}"
        log_info "Evaluating ${vcf} against truth ${truth_vcf}"
        if [[ "${DRY_RUN}" == "true" ]]; then
            cat <<EOF
+ Rscript -e "library(vcfppR); res <- vcfcomp(test='${vcf}', truth='${truth_vcf}', stats='r2', formats=c('DS','GT')); write.table(res, file='${out_prefix}.r2.tsv', sep='\t', quote=FALSE, row.names=FALSE); pdf('${out_prefix}.r2.pdf'); vcfplot(res, col=2, cex=1.2, lwd=2, type='b'); dev.off()"
EOF
            continue
        fi

        Rscript - "${vcf}" "${truth_vcf}" "${out_prefix}" <<'RSCRIPT'
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
    done
}

# Pick the best panel VCF for a chromosome from REFERENCE_PANEL_DIR.
# Priority: *phased.vcf.gz > *beagle.vcf.gz > first match.
# Returns empty string if no suitable VCF found.
pick_panel_vcf() {
    local chr="$1"
    local best_phased=""
    local best_beagle=""
    local best_any=""

    local -a candidates=(
        "${REFERENCE_PANEL_DIR}/apple_panel.refpol.${chr}.vcf.gz"
        "${REFERENCE_PANEL_DIR}/panel.snps.clean__${chr}.vcf.gz"
        "${REFERENCE_PANEL_DIR}/${chr}.vcf.gz"
    )
    local -a globbed=()
    # Use patterns with separators to avoid matching chr1 with chr10, chr11, etc.
    # Match: chr_*, chr.*, chr-* patterns
    while IFS= read -r path; do
        globbed+=( "${path}" )
    done < <(find "${REFERENCE_PANEL_DIR}" -maxdepth 1 \( -type f -o -type l \) \( -name "${chr}_*.vcf.gz" -o -name "${chr}.*.vcf.gz" -o -name "${chr}-*.vcf.gz" \) 2>/dev/null | sort)
    candidates+=( "${globbed[@]}" )

    for candidate in "${candidates[@]}"; do
        [[ -s "${candidate}" ]] || continue
        local base
        base="$(basename "${candidate}")"
        if [[ "${base}" =~ phased\.vcf\.gz$ ]]; then
            best_phased="${candidate}"
            break
        elif [[ -z "${best_beagle}" && "${base}" =~ beagle\.vcf\.gz$ ]]; then
            best_beagle="${candidate}"
        elif [[ -z "${best_any}" ]]; then
            best_any="${candidate}"
        fi
    done

    if [[ -n "${best_phased}" ]]; then
        echo "${best_phased}"
    elif [[ -n "${best_beagle}" ]]; then
        echo "${best_beagle}"
    elif [[ -n "${best_any}" ]]; then
        echo "${best_any}"
    else
        echo ""
    fi
}

normalize_panel_vcf() {
    local chr="$1"
    local invcf=""
    invcf="$(pick_panel_vcf "${chr}")"

    if [[ -z "${invcf}" ]]; then
        log_error "Reference panel VCF for ${chr} not found in ${REFERENCE_PANEL_DIR}."
        exit 1
    fi

    log_info "Using Step1C panel VCF for ${chr}: ${invcf}"
    if [[ "${DRY_RUN}" != "true" ]]; then
        if [[ "${invcf}" =~ \.vcf\.gz$ ]]; then
            if ! bcftools index -f -c "${invcf}"; then
                log_error "Failed to index panel VCF: ${invcf}"
                exit 1
            fi
        fi
    else
        echo "+ bcftools index -f -c \"${invcf}\""
    fi

    if [[ "${REMOVE_MISSING}" == "false" ]]; then
        local first_contig=""
        first_contig=$(gzip -dc "${invcf}" 2>/dev/null | awk '!/^#/ {print $1; exit}')
        if [[ -n "${first_contig}" && -z "${CHUNK_FILE}" ]]; then
            log_warn "Detected contig name '${first_contig}' in ${invcf}. Ensure --chr/--region values match this naming."
        fi
        echo "${invcf}"
        return 0
    fi

    local cleaned_vcf="${PANEL_OUT_DIR}/quilt.nomiss.${chr}.vcf.gz"

    # Filter expression: keep variants where ≥MIN_PHASED_RATE of samples have phased genotypes (contain "|")
    # This removes variants with too many missing or unphased genotypes
    local filter_expr='COUNT(GT~"[|]") >= '"${MIN_PHASED_RATE}"' * N_SAMPLES'

    if [[ "${DRY_RUN}" == "true" ]]; then
        cat <<EOF
+ bcftools view -H "${invcf}" | wc -l    # total variants
+ bcftools view -e '${filter_expr}' -H "${invcf}" | wc -l    # removed (below ${MIN_PHASED_RATE} phased rate)
+ bcftools view -i '${filter_expr}' "${invcf}" -Oz -o "${cleaned_vcf}"
+ bcftools index -f -c "${cleaned_vcf}"
EOF
        echo "${cleaned_vcf}"
        return 0
    fi

    if [[ ! -f "${cleaned_vcf}" ]]; then
        local total removed kept
        total=$(bcftools view -H "${invcf}" | wc -l | awk '{print $1}')
        removed=$(bcftools view -e "${filter_expr}" -H "${invcf}" | wc -l | awk '{print $1}')
        if ! bcftools view -i "${filter_expr}" "${invcf}" -Oz -o "${cleaned_vcf}"; then
            log_error "bcftools filter (min phased rate ${MIN_PHASED_RATE}) failed for ${invcf}"
            exit 1
        fi
        run_cmd bcftools index -f -c "${cleaned_vcf}"
        kept=$((total - removed))
        if [[ ! -f "${MISSING_REPORT}" ]]; then
            echo -e "chromosome\tinput_vcf\ttotal\tremoved\tkept\tmin_phased_rate\toutput_vcf" > "${MISSING_REPORT}"
        fi
        echo -e "${chr}\t${invcf}\t${total}\t${removed}\t${kept}\t${MIN_PHASED_RATE}\t${cleaned_vcf}" >> "${MISSING_REPORT}"
        log_info "Filtered ${chr}: ${removed}/${total} variants removed (below ${MIN_PHASED_RATE} phased rate); ${kept} kept. Output: ${cleaned_vcf}"
    else
        log_info "Cleaned panel already exists for ${chr}: ${cleaned_vcf}"
    fi

    local first_contig=""
    first_contig=$(gzip -dc "${cleaned_vcf}" 2>/dev/null | awk '!/^#/ {print $1; exit}')
    if [[ -n "${first_contig}" && -z "${CHUNK_FILE}" ]]; then
        log_warn "Detected contig name '${first_contig}' in ${cleaned_vcf}. Ensure --chr/--region values match this naming."
    fi

    echo "${cleaned_vcf}"
}

prepare_reference_chunk() {
    local chr="$1"
    local start="$2"
    local end="$3"
    local buffer="$4"
    local panel_vcf="$5"

    local prepared_file="${RDATA_DIR}/QUILT_prepared_reference.${chr}.${start}.${end}.RData"
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
    local chr_map_file
    chr_map_file="$(resolve_genetic_map "${chr}")"
    local cmd=(Rscript "${QUILT2_PREP_SCRIPT}"
        "--genetic_map_file=${chr_map_file}"
        "--reference_vcf_file=${panel_vcf}"
        "--chr=${chr}"
        "--regionStart=${start}"
        "--regionEnd=${end}"
        "--nGen=${NGEN}"
        "--buffer=${buffer}"
        "--outputdir=${OUTPUT_DIR}"
    )

    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "+ ${cmd[*]}"
    else
        "${cmd[@]}"
    fi

    echo "${prepared_file}"
}

impute_chunk() {
    local chr="$1"
    local start="$2"
    local end="$3"
    local buffer="$4"
    local prepared_file="$5"

    if [[ "${DRY_RUN}" != "true" && ! -f "${prepared_file}" ]]; then
        log_error "Prepared reference missing for ${chr}:${start}-${end}: ${prepared_file}"
        exit 1
    fi

    local output_vcf="${OUTPUT_DIR}/quilt2.diploid.${chr}.${start}-${end}.vcf.gz"
    if [[ -f "${output_vcf}" && "${PREP_ONLY}" == "false" ]]; then
        log_info "Imputation output exists for ${chr}:${start}-${end}; skipping."
        return
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

    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "+ ${cmd[*]}"
    else
        "${cmd[@]}"
    fi

    OUTPUT_VCFS+=("${output_vcf}")
}

log_info "Working directory: ${WORK_DIR}"
log_info "Panel dir: ${REFERENCE_PANEL_DIR}"
log_info "Output dir: ${OUTPUT_DIR}"
if [[ "${GENETIC_MAP_IS_DIR}" == "true" ]]; then
    log_info "Genetic map: ${GENETIC_MAP_DIR} (per-chromosome)"
else
    log_info "Genetic map: ${GENETIC_MAP_FILE}"
fi
if [[ -n "${REFERENCE_FASTA}" ]]; then
    log_info "Reference FASTA (validation only): ${REFERENCE_FASTA}"
fi
if [[ -n "${BAMLIST}" ]]; then
    log_info "BAM list: ${BAMLIST}"
fi
if [[ -n "${TRUTH_VCF}" ]]; then
    log_info "Truth VCF for evaluation: ${TRUTH_VCF}"
fi

declare -A PANEL_CACHE=()
declare -a OUTPUT_VCFS=()

# -----------------------------------------------------------------------------
# Main processing loop
# NOTE: Currently runs interactively (sequential). SLURM array parallelization
#       will be added once the pipeline enters production. All intermediate and
#       output files are written to ${OUTPUT_DIR} (under WORK_DIR/quilt2_output/).
# -----------------------------------------------------------------------------
for chunk in "${CHUNKS[@]}"; do
    IFS='|' read -r chunk_id chr start end buffer <<< "${chunk}"
    log_info "---- Chunk ${chunk_id} (${chr}:${start}-${end}, buffer=${buffer}) ----"

    prepared_file="${RDATA_DIR}/QUILT_prepared_reference.${chr}.${start}.${end}.RData"

    if [[ "${IMPUTE_ONLY}" != "true" ]]; then
        if [[ -z "${PANEL_CACHE[${chr}]+x}" ]]; then
            PANEL_CACHE["${chr}"]="$(normalize_panel_vcf "${chr}")"
        fi
        panel_vcf="${PANEL_CACHE[${chr}]}"
        prepared_file="$(prepare_reference_chunk "${chr}" "${start}" "${end}" "${buffer}" "${panel_vcf}")"
    fi

    if [[ "${PREP_ONLY}" == "true" ]]; then
        continue
    fi

    impute_chunk "${chr}" "${start}" "${end}" "${buffer}" "${prepared_file}"
done

run_evaluation "${TRUTH_VCF}" "${OUTPUT_VCFS[@]}"

log_info "QUILT2 pipeline finished. Outputs in ${OUTPUT_DIR}"
