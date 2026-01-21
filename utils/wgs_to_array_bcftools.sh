#!/bin/bash
#SBATCH --job-name=wgs2array_bcf
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
# Submit with: sbatch utils/wgs_to_array_bcftools.sh /path/input.vcf.gz [output.vcf.gz]

# Recodes WGS GT to array-style using bcftools +setGT.
# Assumes biallelic SNPs and diploid samples; other sites are left unchanged.
# REF/ALT are preserved unless FORCE_REF_ALT=true.

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  bash utils/wgs_to_array_bcftools.sh /path/input.vcf.gz [output.vcf.gz]

Notes:
  - Requires bcftools with the +setGT plugin.
  - Input should be biallelic SNPs and diploid samples.
  - Output defaults to <input_basename>_array.vcf.gz in the same directory.
  - REF/ALT are preserved unless FORCE_REF_ALT=true (uses ARRAY_REF/ARRAY_ALT, default A/G).
USAGE
}

INPUT="${1:-}"
OUTPUT="${2:-}"

if [[ -z "${INPUT}" ]]; then
    usage
    exit 1
fi

if [[ ! -f "${INPUT}" ]]; then
    echo "[ERROR] Input VCF not found: ${INPUT}" >&2
    exit 1
fi

# Derive default output path if not supplied.
if [[ -z "${OUTPUT}" ]]; then
    base="$(basename "${INPUT}")"
    dir="$(cd "$(dirname "${INPUT}")" && pwd)"
    case "${base}" in
        *.vcf.gz) OUTPUT="${dir}/${base%.vcf.gz}_array.vcf.gz" ;;
        *.vcf) OUTPUT="${dir}/${base%.vcf}_array.vcf.gz" ;;
        *.bcf) OUTPUT="${dir}/${base%.bcf}_array.vcf.gz" ;;
        *) OUTPUT="${dir}/${base}.array.vcf.gz" ;;
    esac
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "[ERROR] bcftools is required but not found in PATH." >&2
    exit 1
fi

if ! bcftools +setGT --help >/dev/null 2>&1; then
    echo "[ERROR] bcftools +setGT plugin is not available." >&2
    exit 1
fi

THREADS="${THREADS:-1}"
FORCE_REF_ALT="${FORCE_REF_ALT:-false}"
ARRAY_REF="${ARRAY_REF:-A}"
ARRAY_ALT="${ARRAY_ALT:-G}"

if [[ "${FORCE_REF_ALT}" == "true" ]]; then
    if [[ ! "${ARRAY_REF}" =~ ^[ACGT]$ || ! "${ARRAY_ALT}" =~ ^[ACGT]$ ]]; then
        echo "[ERROR] ARRAY_REF/ARRAY_ALT must be single A/C/G/T bases." >&2
        exit 1
    fi
    if [[ "${ARRAY_REF}" == "${ARRAY_ALT}" ]]; then
        echo "[ERROR] ARRAY_REF and ARRAY_ALT must be different bases." >&2
        exit 1
    fi
fi

case "${OUTPUT}" in
    *.vcf.gz) OUTTYPE="-Oz" ;;
    *.vcf) OUTTYPE="-Ov" ;;
    *.bcf) OUTTYPE="-Ob" ;;
    *)
        echo "[WARN] Unrecognized output extension; writing bgzipped VCF." >&2
        OUTPUT="${OUTPUT}.vcf.gz"
        OUTTYPE="-Oz"
        ;;
 esac

AT_EXPR='REF~"^[AT]$" && ALT~"^[AT]$"'
CG_EXPR='REF~"^[CG]$" && ALT~"^[CG]$"'
CG_AT_EXPR='REF~"^[CG]$" && ALT~"^[AT]$"'
DIPLOID='GT~"^[0-9]+[/|][0-9]+$"'
HOMREF='GT~"^0[/|]0$"'
HOMALT='GT~"^1[/|]1$"'

recode_stream() {
    bcftools view -Ou --threads "${THREADS}" "${INPUT}" \
        | bcftools +setGT -Ou -- -n 0/0 -i "${AT_EXPR} && ${DIPLOID}" \
        | bcftools +setGT -Ou -- -n 1/1 -i "${CG_EXPR} && ${DIPLOID}" \
        | bcftools +setGT -Ou -- -n 1/1 -i "${CG_AT_EXPR} && ${HOMREF}" \
        | bcftools +setGT -Ou -- -n 0/0 -i "${CG_AT_EXPR} && ${HOMALT}"
}

if [[ "${FORCE_REF_ALT}" == "true" ]]; then
    if [[ "${OUTTYPE}" == "-Ov" ]]; then
        recode_stream \
            | bcftools view -Ov \
            | awk -v ref="${ARRAY_REF}" -v alt="${ARRAY_ALT}" 'BEGIN{OFS="\t"} /^#/ {print; next} {$4=ref; $5=alt; print}' \
            > "${OUTPUT}"
    else
        recode_stream \
            | bcftools view -Ov \
            | awk -v ref="${ARRAY_REF}" -v alt="${ARRAY_ALT}" 'BEGIN{OFS="\t"} /^#/ {print; next} {$4=ref; $5=alt; print}' \
            | bcftools view "${OUTTYPE}" -o "${OUTPUT}"
    fi
else
    recode_stream \
        | bcftools view "${OUTTYPE}" -o "${OUTPUT}"
fi

if [[ "${OUTPUT}" =~ \.vcf\.gz$ || "${OUTPUT}" =~ \.bcf$ ]]; then
    bcftools index -f -c "${OUTPUT}"
fi

echo "[INFO] Done. Output: ${OUTPUT}" >&2
