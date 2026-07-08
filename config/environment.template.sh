#!/bin/bash
# =============================================================================
# QUILT2 PIPELINE ENVIRONMENT TEMPLATE
# =============================================================================
# Copy to config/environment.sh and adjust paths, tools, and behavior toggles
# for your cluster.
# The orchestrator will source config/environment.sh if present; otherwise it
# will fall back to this template.
# =============================================================================

# -----------------------------------------------------------------------------
# PIPELINE LOCATION & STORAGE
# -----------------------------------------------------------------------------

# Root of the QUILT2 pipeline checkout
QUILT2_ROOT="${QUILT2_ROOT:-$HOME/QUILT2_Pipeline_KH_v1}"

# Scratch base for optional task staging
QUILT2_SCRATCH_BASE="${QUILT2_SCRATCH_BASE:-/scratch/user/$(whoami)}"

# Persistent pipeline output directory. If empty, defaults to <input-dir>/quilt2_output.
QUILT2_OUTPUT_DIR="${QUILT2_OUTPUT_DIR:-}"

# Optional scratch/staging root. If empty, array tasks use $TMPDIR when present,
# otherwise <output-dir>/scratch. Do not store final results only in scratch.
QUILT2_SCRATCH_DIR="${QUILT2_SCRATCH_DIR:-}"

# -----------------------------------------------------------------------------
# REFERENCES
# -----------------------------------------------------------------------------

# Reference genome (FASTA with .fai/.dict)
QUILT2_REFERENCE_FASTA="${QUILT2_REFERENCE_FASTA:-/QRISdata/Q8367/Reference_Genome/GDDH13_1-1_formatted.fasta}"

# Genetic map (file or directory with per-chromosome maps); can be overridden via CLI
QUILT2_GENETIC_MAP="${QUILT2_GENETIC_MAP:-/QRISdata/Q8367/Genetic_Maps/PLEASE_SET}"

# Reference panel directory. Set this here or pass --reference-panel-dir.
QUILT2_REFERENCE_PANEL_DIR="${QUILT2_REFERENCE_PANEL_DIR:-}"

# SLURM resource defaults live in config/quilt2_config.sh.

# -----------------------------------------------------------------------------
# TOOLS / ENVIRONMENTS
# -----------------------------------------------------------------------------

# bcftools module name to load if bcftools not on PATH
# NOTE: Module names can be case-sensitive on HPC (see GATK pipeline incident log).
BCFTOOLS_MODULE="${BCFTOOLS_MODULE:-bcftools/1.18-gcc-12.3.0}"

# Conda environment for QUILT2 scripts
QUILT2_CONDA_ENV="${QUILT2_CONDA_ENV:-quilt2}"

# Optional overrides for QUILT2 R scripts or home directory
QUILT2_HOME="${QUILT2_HOME:-}"
QUILT2_PREP_SCRIPT="${QUILT2_PREP_SCRIPT:-}"
QUILT2_RUN_SCRIPT="${QUILT2_RUN_SCRIPT:-}"

# -----------------------------------------------------------------------------
# BEHAVIOR TOGGLES / DEFAULTS
# -----------------------------------------------------------------------------

# Chromosome list and chunking defaults
QUILT2_CHROMS="${QUILT2_CHROMS:-Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17}"
QUILT2_BUFFER="${QUILT2_BUFFER:-500000}"
QUILT2_NGEN="${QUILT2_NGEN:-100}"
QUILT2_REGION_START="${QUILT2_REGION_START:-1}"
QUILT2_REGION_END="${QUILT2_REGION_END:-}"       # empty means require chunk file/auto-chunk
QUILT2_CHUNK_FILE="${QUILT2_CHUNK_FILE:-}"
QUILT2_AUTO_CHUNK_MAP="${QUILT2_AUTO_CHUNK_MAP:-false}"

# Panel prep toggles
QUILT2_REMOVE_MISSING="${QUILT2_REMOVE_MISSING:-false}"
QUILT2_MIN_VALID_GT_RATE="${QUILT2_MIN_VALID_GT_RATE:-0.95}"
QUILT2_STANDARDISE_NAME="${QUILT2_STANDARDISE_NAME:-false}"
QUILT2_STANDARDISE_NAME_FORCE="${QUILT2_STANDARDISE_NAME_FORCE:-false}"

# Phase selection and dry-run
QUILT2_PREP_ONLY="${QUILT2_PREP_ONLY:-false}"
QUILT2_IMPUTE_ONLY="${QUILT2_IMPUTE_ONLY:-false}"
QUILT2_DRY_RUN="${QUILT2_DRY_RUN:-false}"

# Optional BAM list override (otherwise inferred from working directory)
QUILT2_BAMLIST="${QUILT2_BAMLIST:-}"

# -----------------------------------------------------------------------------
# LOGGING / NOTIFICATIONS
# -----------------------------------------------------------------------------

QUILT2_LOG_LEVEL="${QUILT2_LOG_LEVEL:-INFO}"
QUILT2_EMAIL_NOTIFICATIONS="${QUILT2_EMAIL_NOTIFICATIONS:-false}"
QUILT2_EMAIL_ADDRESS="${QUILT2_EMAIL_ADDRESS:-}"
