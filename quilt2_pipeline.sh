#!/bin/bash
set -euo pipefail

# QUILT2 pipeline wrapper for apple data. Uses the provided working folder as
# the hub for inputs/outputs, reusing GATK pipeline defaults where possible.

DEFAULT_APPLE_REF="/QRISdata/Q8367/Reference_Genome/GDDH13_1-1_formatted.fasta"
DEFAULT_CHROMS=(Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17)
DEFAULT_BUFFER=500000
DEFAULT_NGEN=100

log_info() { echo "[INFO] $*"; }
log_warn() { echo "[WARN] $*"; }
log_error() { echo "[ERROR] $*" >&2; }

usage() {
    cat <<'EOF'
Usage:
  bash quilt2_pipeline.sh -i <work_dir> --genetic-map <map.txt> --region-start <start> --region-end <end> [options]
  bash quilt2_pipeline.sh -i <work_dir> --chunk-file <chunks.tsv> [options]

Required:
  -i, --input-dir PATH        Working folder (all outputs stay here).
  --genetic-map PATH          Genetic map file (or set QUILT2_GENETIC_MAP).
  --region-start N            Region start (bp) unless --chunk-file is provided.
  --region-end N              Region end (bp) unless --chunk-file is provided.

Core options:
  --chunk-file PATH           TSV with either: chr start end [buffer] OR chunk chr start end [buffer].
  --chr LIST                  Comma/space list of chromosomes (default apple Chr01-17).
  --buffer N                  Buffer in bp passed to QUILT2 (default 500000).
  --n-gen N                   nGen passed to QUILT2 (default 100).
  --bamlist PATH              BAM list (defaults to <work_dir>/bamlist.txt or bamlist.1.0.txt).

Reference/panel options:
  --reference-fasta PATH      Reference FASTA (default PIPELINE_REFERENCE_FASTA or apple GDDH13 path).
  --reference-panel-dir PATH  Directory with phased panel VCFs; defaults to
                              <work_dir>/8.Imputated_VCF_BEAGLE, then 7.Consolidated_VCF, then <work_dir>.
  --quilt2-home PATH          Directory containing QUILT2.R and QUILT2_prepare_reference.R.
  --quilt2-prepare-script PATH  Override path to QUILT2_prepare_reference.R.
  --quilt2-run-script PATH      Override path to QUILT2.R.
  --truth-vcf PATH            Optional truth VCF for evaluation (vcfppR); contigs must match outputs.
  --eval-output PATH          Optional evaluation output directory (default: <work_dir>/quilt2_output/eval).
  --remove-missing            Drop variants where all genotypes are missing (default: false).
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
QUILT2_HOME=""
QUILT2_PREP_SCRIPT="${QUILT2_PREP_SCRIPT:-}"
QUILT2_RUN_SCRIPT="${QUILT2_RUN_SCRIPT:-}"
TRUTH_VCF=""
EVAL_OUTPUT_DIR=""
QUILT2_CONDA_ENV="${QUILT2_CONDA_ENV:-quilt2}"
MISSING_REPORT=""
REMOVE_MISSING="false"
PREP_ONLY="false"
IMPUTE_ONLY="false"
DRY_RUN="false"

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

if [[ -z "${REFERENCE_FASTA}" ]]; then
    log_error "Reference FASTA is required (--reference-fasta or PIPELINE_REFERENCE_FASTA)."
    exit 1
fi

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

load_quilt_env

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

if [[ ! -f "${GENETIC_MAP_FILE}" ]]; then
    log_error "Genetic map not found: ${GENETIC_MAP_FILE}"
    exit 1
fi
if [[ ! -f "${REFERENCE_FASTA}" ]]; then
    log_error "Reference FASTA not found: ${REFERENCE_FASTA}"
    exit 1
fi

readarray -t CHR_LIST <<< "$(printf "%s\n" "${DEFAULT_CHROMS[@]}")"
if [[ -n "${CHROM_ARG}" ]]; then
    # Accept comma or space separated values
    IFS=', ' read -r -a CHR_LIST <<< "${CHROM_ARG}"
fi

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
        else
            add_chunk "${c1:-chunk}" "${c1}" "${c2}" "${c3}" "${c4:-${BUFFER}}"
        fi
    done < "${CHUNK_FILE}"
else
    if [[ -z "${REGION_END}" ]]; then
        log_error "Provide --region-end or a --chunk-file."
        exit 1
    fi
    for chr in "${CHR_LIST[@]}"; do
        add_chunk "${chr}_${REGION_START}_${REGION_END}" "${chr}" "${REGION_START}" "${REGION_END}" "${BUFFER}"
    done
fi

run_cmd() {
    if [[ "${DRY_RUN}" == "true" ]]; then
        echo "+ $*"
    else
        "$@"
    fi
}

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

normalize_panel_vcf() {
    local chr="$1"
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
        while IFS= read -r path; do
            globbed+=( "${path}" )
        done < <(find "${REFERENCE_PANEL_DIR}" -maxdepth 1 -type f -name "${chr}*.vcf.gz" | sort)
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

    local invcf=""
    invcf="$(pick_panel_vcf "${chr}")"

    if [[ -z "${invcf}" ]]; then
        log_error "Reference panel VCF for ${chr} not found in ${REFERENCE_PANEL_DIR}."
        exit 1
    fi

    log_info "Using Step1C panel VCF for ${chr}: ${invcf}"
    if [[ "${DRY_RUN}" != "true" ]]; then
        if [[ "${invcf}" =~ \.vcf\.gz$ ]]; then
            run_cmd bcftools index -f -c "${invcf}"
        fi
    else
        echo "+ bcftools index -f -c \"${invcf}\""
    fi

    if [[ "${REMOVE_MISSING}" == "false" ]]; then
        local first_contig=""
        first_contig=$(zcat "${invcf}" 2>/dev/null | awk '!/^#/ {print $1; exit}')
        if [[ -n "${first_contig}" && -z "${CHUNK_FILE}" ]]; then
            log_warn "Detected contig name '${first_contig}' in ${invcf}. Ensure --chr/--region values match this naming."
        fi
        echo "${invcf}"
        return 0
    fi

    local cleaned_vcf="${PANEL_OUT_DIR}/quilt.nomiss.${chr}.vcf.gz"

    if [[ "${DRY_RUN}" == "true" ]]; then
        cat <<EOF
+ bcftools view -H "${invcf}" | wc -l    # total variants
+ bcftools view -e 'COUNT(GT!=".|.")>0' -H "${invcf}" | wc -l    # removed (all missing)
+ bcftools view -i 'COUNT(GT!=".|.")>0' "${invcf}" -Oz -o "${cleaned_vcf}"
+ bcftools index -f -c "${cleaned_vcf}"
EOF
        echo "${cleaned_vcf}"
        return 0
    fi

    if [[ ! -f "${cleaned_vcf}" ]]; then
        local total missing kept
        total=$(bcftools view -H "${invcf}" | wc -l | awk '{print $1}')
        missing=$(bcftools view -e 'COUNT(GT!=".|.")>0' -H "${invcf}" | wc -l | awk '{print $1}')
        if ! bcftools view -i 'COUNT(GT!=".|.")>0' "${invcf}" -Oz -o "${cleaned_vcf}"; then
            log_error "bcftools filter (remove all-missing genotypes) failed for ${invcf}"
            exit 1
        fi
        run_cmd bcftools index -f -c "${cleaned_vcf}"
        kept=$((total - missing))
        if [[ ! -f "${MISSING_REPORT}" ]]; then
            echo -e "chromosome\tinput_vcf\tremoved_all_missing\tkept_variants\toutput_vcf" > "${MISSING_REPORT}"
        fi
        echo -e "${chr}\t${invcf}\t${missing}\t${kept}\t${cleaned_vcf}" >> "${MISSING_REPORT}"
        log_info "Removed ${missing} all-missing variants; kept ${kept}. Cleaned: ${cleaned_vcf}"
    else
        log_info "Cleaned panel already exists for ${chr}: ${cleaned_vcf}"
    fi

    local first_contig=""
    first_contig=$(zcat "${cleaned_vcf}" 2>/dev/null | awk '!/^#/ {print $1; exit}')
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
    local cmd=(Rscript "${QUILT2_PREP_SCRIPT}"
        "--genetic_map_file=${GENETIC_MAP_FILE}"
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
log_info "Genetic map: ${GENETIC_MAP_FILE}"
log_info "Reference FASTA: ${REFERENCE_FASTA}"
if [[ -n "${BAMLIST}" ]]; then
    log_info "BAM list: ${BAMLIST}"
fi
if [[ -n "${TRUTH_VCF}" ]]; then
    log_info "Truth VCF for evaluation: ${TRUTH_VCF}"
fi

declare -A PANEL_CACHE=()
declare -a OUTPUT_VCFS=()

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
