#!/bin/bash
# Ad-hoc WGS -> array recode without bcftools +setGT.
# Uses an embedded pysam-based recoder; REF/ALT forcing is optional.

set -euo pipefail

usage() {
    cat <<'USAGE'
Usage:
  bash utils/adhoc_wgs_to_array.sh /path/input.vcf.gz [output.vcf.gz]

Notes:
  - Requires python + pysam. bcftools is optional (used only for indexing if present).
  - Assumes biallelic SNPs and diploid samples.
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

PYTHON_BIN="${PYTHON_BIN:-python3}"
if ! command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
    echo "[ERROR] Python not found (tried ${PYTHON_BIN})." >&2
    exit 1
fi

if ! "${PYTHON_BIN}" - <<'PY' >/dev/null 2>&1
import pysam  # noqa: F401
PY
then
    echo "[ERROR] pysam is required in the active python environment." >&2
    exit 1
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

if [[ "${FORCE_REF_ALT}" == "true" ]]; then
    echo "[INFO] Recoding GT and forcing REF/ALT to ${ARRAY_REF}/${ARRAY_ALT}" >&2
else
    echo "[INFO] Recoding GT (REF/ALT preserved)" >&2
fi
"${PYTHON_BIN}" - "${INPUT}" "${OUTPUT}" "${ARRAY_REF}" "${ARRAY_ALT}" "${FORCE_REF_ALT}" <<'PY'
import sys

import pysam


def allele_group(base: str):
    if base in ("A", "T"):
        return "AT"
    if base in ("C", "G"):
        return "CG"
    return None


def recode_genotype(gt, allele_lookup):
    if gt is None or len(gt) != 2:
        return (None, None)

    bases = []
    for idx in gt:
        if idx is None or idx < 0:
            return (None, None)
        base = allele_lookup.get(idx)
        if base is None or len(base) != 1:
            return (None, None)
        group = allele_group(base.upper())
        if group is None:
            return (None, None)
        bases.append(group)

    if bases[0] == "AT" and bases[1] == "AT":
        return (0, 0)
    if bases[0] == "CG" and bases[1] == "CG":
        return (1, 1)
    return (0, 1)


def output_mode(path: str) -> str:
    if path.endswith(".vcf.gz"):
        return "wz"
    if path.endswith(".bcf"):
        return "wb"
    return "w"


def main():
    if len(sys.argv) != 6:
        sys.stderr.write("Usage: <input> <output> <array_ref> <array_alt> <force_ref_alt>\\n")
        sys.exit(2)
    input_vcf, output_vcf, array_ref, array_alt, force_flag = sys.argv[1:6]
    force_ref_alt = force_flag.lower() == "true"
    mode = output_mode(output_vcf)

    with pysam.VariantFile(input_vcf, "r") as reader:
        header = reader.header.copy()
        meta = "adhoc_wgs_to_array.sh recoded GT (AT->0/0, CG->1/1, mix->0/1)"
        if force_ref_alt:
            meta += f"; REF/ALT forced to {array_ref}/{array_alt}"
        else:
            meta += "; REF/ALT preserved"
        header.add_meta("source", meta)

        with pysam.VariantFile(output_vcf, mode, header=header) as writer:
            for record in reader:
                allele_lookup = {0: record.ref}
                for i, alt in enumerate(record.alts or [], start=1):
                    allele_lookup[i] = alt

                for sample in record.samples:
                    call = record.samples[sample]
                    phased = call.phased
                    new_gt = recode_genotype(call.get("GT"), allele_lookup)
                    call["GT"] = new_gt
                    call.phased = phased

                if force_ref_alt:
                    record.ref = array_ref
                    record.alts = (array_alt,)
                writer.write(record)


if __name__ == "__main__":
    main()
PY

if [[ "${OUTPUT}" =~ \.vcf\.gz$ || "${OUTPUT}" =~ \.bcf$ ]]; then
    if command -v bcftools >/dev/null 2>&1; then
        bcftools index -f -c "${OUTPUT}"
    elif [[ "${OUTPUT}" =~ \.vcf\.gz$ ]]; then
        "${PYTHON_BIN}" - "${OUTPUT}" <<'PY'
import sys
import pysam

path = sys.argv[1]
pysam.tabix_index(path, preset="vcf", force=True)
PY
    else
        echo "[WARN] Index not created (bcftools not available for BCF)." >&2
    fi
fi

echo "[INFO] Done. Output: ${OUTPUT}" >&2
