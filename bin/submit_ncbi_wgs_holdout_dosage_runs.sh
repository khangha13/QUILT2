#!/usr/bin/env bash
set -euo pipefail

TRUTH_DIR="/QRISdata/Q8367/WGS_Reference_Panel/NCBI_truth_set/7.Consolidated_VCF"
BASE="/scratch/project/bigdata_apple/downsampling/NCBI_downsampling_imputed"
DOSAGE_R2="/scratch/project/bigdata_apple/QUILT2/bin/dosage_r2_sbatch.sh"

# Filter
bash "${DOSAGE_R2}" \
  --truth-mode wgs \
  --chunks-dir "${BASE}/filter/NCBI_WGS_Combined_holdout/chunks/imputed" \
  --truth-dataset-dir "${TRUTH_DIR}" \
  -- --force

bash "${DOSAGE_R2}" \
  --truth-mode wgs \
  --chunks-dir "${BASE}/filter/NCBI_WGS_Liao_holdout/chunks/imputed" \
  --truth-dataset-dir "${TRUTH_DIR}" \
  -- --force

bash "${DOSAGE_R2}" \
  --truth-mode wgs \
  --chunks-dir "${BASE}/filter/NCBI_WGS_NCBI_holdout/chunks/imputed" \
  --truth-dataset-dir "${TRUTH_DIR}" \
  -- --force

# No filter
bash "${DOSAGE_R2}" \
  --truth-mode wgs \
  --chunks-dir "${BASE}/no_filter/NCBI_WGS_no_filter_Combined_holdout/chunks/imputed" \
  --truth-dataset-dir "${TRUTH_DIR}" \
  -- --force

bash "${DOSAGE_R2}" \
  --truth-mode wgs \
  --chunks-dir "${BASE}/no_filter/NCBI_WGS_no_filter_Liao_holdout/chunks/imputed" \
  --truth-dataset-dir "${TRUTH_DIR}" \
  -- --force

bash "${DOSAGE_R2}" \
  --truth-mode wgs \
  --chunks-dir "${BASE}/no_filter/NCBI_WGS_no_filter_NCBI_holdout/chunks/imputed" \
  --truth-dataset-dir "${TRUTH_DIR}" \
  -- --force
