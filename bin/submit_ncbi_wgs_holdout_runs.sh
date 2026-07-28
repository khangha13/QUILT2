#!/usr/bin/env bash
set -euo pipefail

BASE="/scratch/project/bigdata_apple"
INPUT_DIR="${BASE}/downsampling/NCBI_downsampling_bam"
IMPUTED_BASE="${BASE}/downsampling/NCBI_downsampling_imputed"
RUN_QUILT2="${BASE}/QUILT2/bin/run_quilt2.sh"
BAMLIST="${BASE}/downsampling/NCBI_downsampling_imputed_biased/NCBI_bamlist.txt"
EXCLUDE="${INPUT_DIR}/accessions.txt"


## This script submits the 6 runs for set of WGS experiments

## This holds out the 7 downsampling samples from the combined, Liao, and NCBI reference panels using --exclude flag 
## (--reference-exclude-samples built-in from QUILT2)


# Filter
bash "${RUN_QUILT2}" \
  -i "${INPUT_DIR}" \
  --output-dir "${IMPUTED_BASE}/filter/NCBI_WGS_Combined_holdout" \
  --bamlist "${BAMLIST}" \
  --reference-panel-dir "${BASE}/phased_reference_panel/combined" \
  --exclude "${EXCLUDE}" \
  --genetic-map dummy \
  --auto-chunk-map

bash "${RUN_QUILT2}" \
  -i "${INPUT_DIR}" \
  --output-dir "${IMPUTED_BASE}/filter/NCBI_WGS_Liao_holdout" \
  --bamlist "${BAMLIST}" \
  --reference-panel-dir "${BASE}/phased_reference_panel/liao" \
  --exclude "${EXCLUDE}" \
  --genetic-map dummy \
  --auto-chunk-map

bash "${RUN_QUILT2}" \
  -i "${INPUT_DIR}" \
  --output-dir "${IMPUTED_BASE}/filter/NCBI_WGS_NCBI_holdout" \
  --bamlist "${BAMLIST}" \
  --reference-panel-dir "${BASE}/phased_reference_panel/NCBI" \
  --exclude "${EXCLUDE}" \
  --genetic-map dummy \
  --auto-chunk-map

# No filter
bash "${RUN_QUILT2}" \
  -i "${INPUT_DIR}" \
  --output-dir "${IMPUTED_BASE}/no_filter/NCBI_WGS_no_filter_Combined_holdout" \
  --bamlist "${BAMLIST}" \
  --reference-panel-dir "${BASE}/imputed_reference_panel/Combined" \
  --exclude "${EXCLUDE}" \
  --genetic-map dummy \
  --auto-chunk-map

bash "${RUN_QUILT2}" \
  -i "${INPUT_DIR}" \
  --output-dir "${IMPUTED_BASE}/no_filter/NCBI_WGS_no_filter_Liao_holdout" \
  --bamlist "${BAMLIST}" \
  --reference-panel-dir "${BASE}/imputed_reference_panel/Liao" \
  --exclude "${EXCLUDE}" \
  --genetic-map dummy \
  --auto-chunk-map

bash "${RUN_QUILT2}" \
  -i "${INPUT_DIR}" \
  --output-dir "${IMPUTED_BASE}/no_filter/NCBI_WGS_no_filter_NCBI_holdout" \
  --bamlist "${BAMLIST}" \
  --reference-panel-dir "${BASE}/imputed_reference_panel/NCBI" \
  --exclude "${EXCLUDE}" \
  --genetic-map dummy \
  --auto-chunk-map
