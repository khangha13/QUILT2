#!/bin/bash
# Shared functions for QUILT2 SLURM-array pipeline
# Extracted from the legacy monolithic quilt2_pipeline.sh

# Logging (stderr to avoid contaminating command substitutions)
log_info() { echo "[INFO] $*" >&2; }
log_warn() { echo "[WARN] $*" >&2; }
log_error() { echo "[ERROR] $*" >&2; }

# Dry-run aware runner
run_cmd() {
    if [[ "${DRY_RUN:-false}" == "true" ]]; then
        echo "+ $*"
    else
        "$@"
    fi
}

require_cmd() {
    local cmd="$1"
    if ! command -v "${cmd}" >/dev/null 2>&1; then
        log_error "Command not found: ${cmd}"
        return 1
    fi
}

resolve_task_scratch_dir() {
    local scratch_root="$1"
    local fallback_root="$2"
    local prefix="${3:-quilt2}"
    local job_id="${SLURM_JOB_ID:-manual}"
    local task_id="${SLURM_ARRAY_TASK_ID:-0}"
    local base=""

    if [[ -n "${scratch_root}" ]]; then
        base="${scratch_root%/}"
    elif [[ -n "${TMPDIR:-}" ]]; then
        base="${TMPDIR%/}"
    else
        base="${fallback_root%/}/scratch"
    fi

    local task_dir="${base}/${prefix}_${job_id}_${task_id}"
    mkdir -p "${task_dir}" || {
        log_error "Failed to create scratch directory: ${task_dir}"
        return 1
    }
    echo "${task_dir}"
}

if command -v module >/dev/null 2>&1; then
    module purge
fi

ensure_bcftools() {
    # Load bcftools if absent; prefer BCFTOOLS_MODULE but fall back to the
    # cluster-safe default used elsewhere in the pipeline.
    local bcftools_module="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"

    if command -v bcftools >/dev/null 2>&1; then
        return 0
    fi
    if command -v module >/dev/null 2>&1; then
        if [[ -n "${bcftools_module}" ]]; then
            if module load "${bcftools_module}"; then
                log_info "Loaded ${bcftools_module} module for bcftools"
            else
                log_warn "Failed to load ${bcftools_module}; bcftools still unavailable"
            fi
        else
            log_warn "BCFTOOLS_MODULE not set; cannot auto-load bcftools."
        fi
    else
        log_warn "module command not found; cannot auto-load bcftools."
    fi

    if ! command -v bcftools >/dev/null 2>&1; then
        log_error "bcftools not found in PATH. Install it or load module ${bcftools_module:-<unset>}."
        return 1
    fi
}

load_quilt_env() {
    local conda_env="${QUILT2_CONDA_ENV:-quilt2}"
    local conda_sourced="false"

    if command -v module >/dev/null 2>&1; then
        if module load miniforge/25.3.0-3 >/dev/null 2>&1; then
            log_info "Loaded miniforge/25.3.0-3 module"
        else
            log_error "Failed to load miniforge/25.3.0-3 module"
            return 1
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
        return 0
    fi

    if conda activate "${conda_env}" >/dev/null 2>&1; then
        log_info "Activated conda environment: ${conda_env}"
    else
        log_error "Failed to activate conda environment: ${conda_env}"
        return 1
    fi
}

# Resolve genetic map for a chromosome
# args: chr, is_dir(bool), map_file, map_dir
resolve_genetic_map() {
    local chr="$1"
    local genetic_map_is_dir="$2"
    local genetic_map_file="$3"
    local genetic_map_dir="$4"

    if [[ "${genetic_map_is_dir}" == "true" ]]; then
        local -a patterns=(
            "${genetic_map_dir}/${chr}.map"
            "${genetic_map_dir}/${chr}.txt"
            "${genetic_map_dir}/${chr}.txt.gz"
            "${genetic_map_dir}/genetic_map.${chr}.txt"
            "${genetic_map_dir}/genetic_map.${chr}.txt.gz"
            "${genetic_map_dir}/genetic_map_${chr}.txt"
            "${genetic_map_dir}/genetic_map_${chr}.txt.gz"
            "${genetic_map_dir}/${chr}_genetic_map.txt"
            "${genetic_map_dir}/${chr}_genetic_map.txt.gz"
        )
        local -a globbed=()
        while IFS= read -r path; do
            globbed+=( "${path}" )
        done < <(find "${genetic_map_dir}" -maxdepth 1 -type f \( -name "*${chr}*.txt" -o -name "*${chr}*.txt.gz" -o -name "*${chr}*.map" \) 2>/dev/null | head -1)
        patterns+=( "${globbed[@]}" )

        for candidate in "${patterns[@]}"; do
            if [[ -f "${candidate}" ]]; then
                echo "${candidate}"
                return 0
            fi
        done
        log_error "Genetic map for ${chr} not found in ${genetic_map_dir}"
        return 1
    else
        echo "${genetic_map_file}"
    fi
}

# Generate a constant-rate (1.0 cM/Mb) dummy genetic map for one chromosome
# args: chr, fai_path, outfile
generate_dummy_genetic_map() {
    local chr="$1"
    local fai="$2"
    local outfile="$3"

    local len
    len=$(awk -v c="${chr}" '$1==c{print $2}' "${fai}")
    if [[ -z "${len}" ]]; then
        log_error "Chromosome ${chr} not found in reference FASTA index: ${fai}"
        return 1
    fi

    awk -v len="${len}" -v rate="1.0" -v step="1000" '
        BEGIN {
            print "position COMBINED_rate.cM.Mb. Genetic_Map.cM."
            first_pos = step
            for (p = step; p <= len; p += step) {
                cM = (p - first_pos) / 1000000.0
                printf "%d %.6f %.6f\n", p, rate, cM
            }
        }
    ' > "${outfile}"
}

# Panel VCF selection with safe globbing (avoids chr1 vs chr10 collisions)
pick_panel_vcf() {
    local reference_panel_dir="$1"
    local chr="$2"

    local best_phased=""
    local best_beagle=""
    local best_any=""

    local -a candidates=(
        "${reference_panel_dir}/nomiss/quilt.nomiss.${chr}.vcf.gz"
        "${reference_panel_dir}/standardised/${chr}_chr.vcf.gz"
        "${reference_panel_dir}/quilt.nomiss.${chr}.vcf.gz"
        "${reference_panel_dir}/${chr}_chr.vcf.gz"
        "${reference_panel_dir}/apple_panel.refpol.${chr}.vcf.gz"
        "${reference_panel_dir}/panel.snps.clean__${chr}.vcf.gz"
        "${reference_panel_dir}/${chr}.vcf.gz"
    )
    local -a globbed=()
    while IFS= read -r path; do
        globbed+=( "${path}" )
    done < <(find "${reference_panel_dir}" -maxdepth 1 \( -type f -o -type l \) \( -name "${chr}_*.vcf.gz" -o -name "${chr}.*.vcf.gz" -o -name "${chr}-*.vcf.gz" \) 2>/dev/null | sort)
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
    else
        echo "${best_any}"
    fi
}

# Normalize/filter panel VCF; returns path to usable VCF
normalize_panel_vcf() {
    local chr="$1"
    local reference_panel_dir="$2"

    local remove_missing="${REMOVE_MISSING:-false}"
    local min_valid_gt_rate="${MIN_VALID_GT_RATE:-0.95}"
    local chunk_file="${CHUNK_FILE:-}"
    local panel_out_dir="${PANEL_OUT_DIR:-}"
    local missing_report="${MISSING_REPORT:-}"

    local invcf
    invcf="$(pick_panel_vcf "${reference_panel_dir}" "${chr}")"
    if [[ -z "${invcf}" ]]; then
        log_error "Reference panel VCF for ${chr} not found in ${reference_panel_dir}."
        return 1
    fi

    log_info "Using Step1C/QUILT2 panel VCF for ${chr}: ${invcf}"
    if [[ "${DRY_RUN:-false}" != "true" ]]; then
        if [[ "${invcf}" =~ \.vcf\.gz$ ]]; then
            if ! bcftools index -f -c "${invcf}"; then
                log_error "Failed to index panel VCF: ${invcf}"
                return 1
            fi
        fi
    else
        echo "+ bcftools index -f -c \"${invcf}\""
    fi

    if [[ "${remove_missing}" != "true" ]]; then
        local first_contig=""
        first_contig=$(gzip -dc "${invcf}" 2>/dev/null | awk '!/^#/ {print $1; exit}')
        if [[ -n "${first_contig}" && -z "${chunk_file}" ]]; then
            log_warn "Detected contig name '${first_contig}' in ${invcf}. Ensure --chr/--region values match this naming."
        fi
        echo "${invcf}"
        return 0
    fi

    if [[ -z "${panel_out_dir}" ]]; then
        log_error "PANEL_OUT_DIR is not set but --remove-missing was requested."
        return 1
    fi
    mkdir -p "${panel_out_dir}"

    local cleaned_vcf="${panel_out_dir}/quilt.nomiss.${chr}.vcf.gz"
    # Filter expression: keep variants with high phased rate AND no missing genotypes (./. or .|.)
    local filter_expr='COUNT(GT~"[|]") >= '"${min_valid_gt_rate}"' * N_SAMPLES && COUNT(GT="mis") = 0'

    if [[ "${DRY_RUN:-false}" == "true" ]]; then
        cat <<EOF
+ bcftools view -H "${invcf}" | wc -l    # total variants
+ bcftools view -e '${filter_expr}' -H "${invcf}" | wc -l    # removed (below ${min_valid_gt_rate} phased rate or has missing GT)
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
            log_error "bcftools filter (min phased rate ${min_valid_gt_rate}, no missing GT) failed for ${invcf}"
            return 1
        fi
        run_cmd bcftools index -f -c "${cleaned_vcf}"
        kept=$((total - removed))
        if [[ -n "${missing_report}" ]]; then
            if [[ ! -f "${missing_report}" ]]; then
                echo -e "chromosome\tinput_vcf\ttotal\tremoved\tkept\tmin_valid_gt_rate\toutput_vcf" > "${missing_report}"
            fi
            echo -e "${chr}\t${invcf}\t${total}\t${removed}\t${kept}\t${min_valid_gt_rate}\t${cleaned_vcf}" >> "${missing_report}"
        fi
        log_info "Filtered ${chr}: ${removed}/${total} variants removed (below ${min_valid_gt_rate} phased rate or missing GT); ${kept} kept. Output: ${cleaned_vcf}"
    else
        log_info "Cleaned panel already exists for ${chr}: ${cleaned_vcf}"
    fi

    local first_contig=""
    first_contig=$(gzip -dc "${cleaned_vcf}" 2>/dev/null | awk '!/^#/ {print $1; exit}')
    if [[ -n "${first_contig}" && -z "${chunk_file}" ]]; then
        log_warn "Detected contig name '${first_contig}' in ${cleaned_vcf}. Ensure --chr/--region values match this naming."
    fi

    echo "${cleaned_vcf}"
}

# Validate core inputs (paths exist)
validate_quilt2_inputs() {
    local work_dir="$1"
    local reference_panel_dir="$2"
    local genetic_map="$3"
    local bamlist="$4"
    local truth_vcf="$5"

    if [[ ! -d "${work_dir}" ]]; then
        log_error "Working directory not found: ${work_dir}"
        return 1
    fi
    if [[ ! -d "${reference_panel_dir}" ]]; then
        log_error "Reference panel directory not found: ${reference_panel_dir}"
        return 1
    fi
    if [[ -z "${genetic_map}" ]]; then
        log_error "Genetic map is required (--genetic-map or QUILT2_GENETIC_MAP)."
        return 1
    fi
    if [[ ! -e "${genetic_map}" ]]; then
        log_error "Genetic map not found: ${genetic_map}"
        return 1
    fi
    if [[ -n "${bamlist}" ]]; then
        if [[ ! -f "${bamlist}" ]]; then
            log_error "BAM list file does not exist: ${bamlist}"
            return 1
        fi
    fi
    if [[ -n "${truth_vcf}" && ! -f "${truth_vcf}" ]]; then
        log_error "Truth VCF not found: ${truth_vcf}"
        return 1
    fi
}

# Persist chunk definitions to manifest (one per line)
create_chunk_manifest() {
    local manifest_path="$1"
    shift
    local -a chunks=( "$@" )

    : > "${manifest_path}"
    for chunk in "${chunks[@]}"; do
        echo "${chunk}" >> "${manifest_path}"
    done

    if [[ ! -s "${manifest_path}" ]]; then
        log_error "Chunk manifest ${manifest_path} is empty."
        return 1
    fi
}
