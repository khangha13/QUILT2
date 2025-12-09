#!/bin/bash
# Convert numeric contig names (1..17) in VCFs to Chr01..Chr17 and reindex.
# Keeps inputs untouched; writes renamed files with a suffix in the same or target directory.

set -euo pipefail

usage() {
    cat <<'EOF'
Usage:
  standardize_chr_names.sh -i <vcf_dir> [options]

Required:
  -i, --vcf-dir PATH        Directory containing *.vcf.gz to rename (numeric -> ChrNN).

Options:
  -o, --out-dir PATH        Output directory (default: same as --vcf-dir).
  -s, --suffix SUFFIX       Suffix inserted before .vcf.gz (default: _renamed).
  -f, --force               Overwrite existing outputs.
  --contigs LIST            Comma list of numeric contigs (default: 1..17).
  --prefix STR              Chromosome prefix (default: Chr).
  --pad N                   Zero-pad width (default: 2, so 1 -> Chr01).
  -h, --help                Show this message.

Notes:
  - Requires bcftools on PATH.
  - Creates a temporary rename map (e.g., "1<TAB>Chr01") and uses bcftools annotate --rename-chrs.
  - Input VCFs are not modified.
EOF
}

VCF_DIR=""
OUT_DIR=""
SUFFIX="_renamed"
CONTIG_LIST=""
PREFIX="Chr"
PAD=2
FORCE="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -i|--vcf-dir) VCF_DIR="$2"; shift 2 ;;
        -o|--out-dir) OUT_DIR="$2"; shift 2 ;;
        -s|--suffix) SUFFIX="$2"; shift 2 ;;
        --contigs) CONTIG_LIST="$2"; shift 2 ;;
        --prefix) PREFIX="$2"; shift 2 ;;
        --pad) PAD="$2"; shift 2 ;;
        -f|--force) FORCE="true"; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown option: $1" >&2; usage; exit 1 ;;
    esac
done

if [[ -z "${VCF_DIR}" ]]; then
    echo "ERROR: --vcf-dir is required" >&2
    usage
    exit 1
fi

if [[ ! -d "${VCF_DIR}" ]]; then
    echo "ERROR: VCF directory not found: ${VCF_DIR}" >&2
    exit 1
fi

if ! command -v bcftools >/dev/null 2>&1; then
    echo "ERROR: bcftools not found in PATH" >&2
    exit 1
fi

if [[ -z "${OUT_DIR}" ]]; then
    OUT_DIR="${VCF_DIR}"
fi
mkdir -p "${OUT_DIR}"

# Default contigs 1..17 for apple
if [[ -z "${CONTIG_LIST}" ]]; then
    CONTIG_LIST="$(seq 1 17 | paste -sd, -)"
fi

IFS=',' read -r -a CONTIGS <<< "${CONTIG_LIST}"

rename_map="$(mktemp)"
trap 'rm -f "${rename_map}"' EXIT

: > "${rename_map}"
for c in "${CONTIGS[@]}"; do
    # Normalize to integer then pad back out
    num=$((10#$c))
    printf "%d\t%s%0${PAD}d\n" "${num}" "${PREFIX}" "${num}" >> "${rename_map}"
done

shopt -s nullglob
vcfs=( "${VCF_DIR}"/*.vcf.gz )
if [[ "${#vcfs[@]}" -eq 0 ]]; then
    echo "ERROR: No .vcf.gz files found in ${VCF_DIR}" >&2
    exit 1
fi

for vcf in "${vcfs[@]}"; do
    base="$(basename "${vcf}" .vcf.gz)"
    # If the file name already suggests it was renamed (e.g., contains "_renamed"), skip unless forcing.
    if [[ "${base}" == *"_renamed"* && "${FORCE}" != "true" ]]; then
        echo "Skipping already-renamed file: ${vcf}" >&2
        continue
    fi
    # Skip files that already end with the suffix (to avoid suffix stacking)
    if [[ "${base}" == *"${SUFFIX}" ]]; then
        echo "Skipping already-standardized file: ${vcf}" >&2
        continue
    fi
    out="${OUT_DIR%/}/${base}${SUFFIX}.vcf.gz"
    if [[ -f "${out}" && "${FORCE}" != "true" ]]; then
        echo "Skipping existing: ${out} (use --force to overwrite)" >&2
        continue
    fi
    echo "Renaming ${vcf} -> ${out}"
    bcftools annotate --rename-chrs "${rename_map}" "${vcf}" -Oz -o "${out}"
    bcftools index -f -c "${out}"
done
