# QUILT2 Pipeline (KH v1)

SLURM-array wrapper around QUILT2 imputation for apple data. Mirrors the Step1C array pattern from the GATK workflow.

## Layout
- `bin/run_quilt2.sh` – orchestrator; builds chunk manifest, generates SLURM array script.
- `templates/quilt2_job.sh` – array worker; processes one chunk (prepare + impute).
- `lib/functions.sh` – shared helpers (env bootstrap, bcftools/QUILT checks, panel/map resolution).
- `config/quilt2_config.sh` – SLURM defaults (account/partition/qos/resources/array cap).
- `quilt2_pipeline.legacy.sh` – pre-array monolithic script (kept only for rollback).
- `quilt2_past_problem_and_solution.md` – troubleshooting history.

## Prerequisites
- Bash, SLURM.
- `Rscript` with QUILT2 scripts available (`QUILT2_prepare_reference.R`, `QUILT2.R`) via `--quilt2-home` or explicit `--quilt2-*` paths.
- `bcftools` on PATH or loadable via `BCFTOOLS_MODULE` (default `bcftools/1.18-gcc-12.3.0`).
- Optional: conda env name via `QUILT2_CONDA_ENV` (default `quilt2`); module `miniforge/25.3.0-3` loaded if present.

## Configuration (environment.sh)
- Copy `config/environment.template.sh` to `config/environment.sh` and set site defaults (scratch/log roots, reference FASTA, genetic map path, panel dir).
- SLURM defaults: `QUILT2_ACCOUNT`, `QUILT2_PARTITION`, `QUILT2_QOS`, `QUILT2_CPUS_PER_TASK`, `QUILT2_MEMORY`, `QUILT2_TIME_LIMIT`, `QUILT2_ARRAY_LIMIT` (0=no cap), `QUILT2_CONSTRAINT` (optional).
- Tooling: `BCFTOOLS_MODULE`, `QUILT2_CONDA_ENV`, optional `QUILT2_HOME`/`QUILT2_PREP_SCRIPT`/`QUILT2_RUN_SCRIPT`.
- Behavior toggles: `QUILT2_CHROMS`, `QUILT2_BUFFER`, `QUILT2_NGEN`, `QUILT2_AUTO_CHUNK_MAP`, `QUILT2_CHUNK_FILE`, `QUILT2_REGION_START/END`, `QUILT2_REMOVE_MISSING`, `QUILT2_MIN_VALID_GT_RATE`, `QUILT2_STANDARDISE_NAME`, `QUILT2_STANDARDISE_NAME_FORCE`, `QUILT2_PREP_ONLY`, `QUILT2_IMPUTE_ONLY`, `QUILT2_DRY_RUN`, `QUILT2_BAMLIST`.

## Inputs
- `--input-dir` (`WORK_DIR`) containing:
  - Panel VCFs under `8.Imputated_VCF_BEAGLE/` (preferred) or `7.Consolidated_VCF/`; otherwise `WORK_DIR` is searched.
  - `bamlist.txt` (or `bamlist.1.0.txt` / `bamlist.tsv`) unless running `--prepare-only`; **still required for `--impute-only`**.
- Genetic map: `--genetic-map` file or directory with per-chromosome maps. Names must match `--chr` values (`Chr01` vs `1`, etc.).
- Reference panel should be **phased**; use `--remove-missing` with `--min-valid-gt-rate` if needed.
- Ensure VCFs are indexed (`.tbi/.csi`).

## Quick Start
Full run (auto chunks from genetic map directory):
```bash
bash bin/run_quilt2.sh \
  -i /path/to/work_dir \
  --genetic-map /path/to/genetic_maps_dir \
  --auto-chunk-map \
  --remove-missing --min-valid-gt-rate 0.95
```

Fixed region for all chromosomes:
```bash
bash bin/run_quilt2.sh \
  -i /path/to/work_dir \
  --genetic-map /path/to/genetic_maps_dir \
  --region-start 1 --region-end 5000000
```

Impute-only (prepared references already exist):
```bash
bash bin/run_quilt2.sh \
  -i /path/to/work_dir \
  --genetic-map /path/to/genetic_maps_dir \
  --impute-only \
  --region-start 1 --region-end 5000000
```

Dry-run (generate SLURM script, no submission):
```bash
bash bin/run_quilt2.sh ... --dry-run
```

Remove-missing toggle:
```bash
  --remove-missing --min-valid-gt-rate 0.9   # keep variants with >=90% phased
```

Point to QUILT2 scripts:
```bash
  --quilt2-home /path/to/QUILT2_scripts
# or
  --quilt2-prepare-script /path/to/QUILT2_prepare_reference.R \
  --quilt2-run-script /path/to/QUILT2.R
```

SLURM overrides (env, matches Step1C style):
```bash
export QUILT2_ACCOUNT=youracct
export QUILT2_PARTITION=compute
export QUILT2_QOS=normal
export QUILT2_CPUS_PER_TASK=8
export QUILT2_MEMORY=48G
export QUILT2_TIME_LIMIT=12:00:00
export QUILT2_ARRAY_LIMIT=0   # 0 = no cap
export QUILT2_CONSTRAINT=epyc4
```

## Execution model (self-submit + phases)
- Default: `bin/run_quilt2.sh` self-submits via `sbatch` when not already in SLURM, then exits. Opt out with `--no-submit` or `--submit-self=false`. Master job ID: `quilt2_slurm/quilt2_master_job_id.txt`; logs: `quilt2_slurm/quilt2_master_%j.(output|error)`.
- Phase 1 (panel prep array): runs per chromosome when `--standardise-name` and/or `--remove-missing` are set. Outputs in `quilt2_output/panel/`:
  - Standardised: `<chr>_chr.vcf.gz` (+ index), skipped unless `--standardise-name-force`.
  - Filtered: `quilt.nomiss.<chr>.vcf.gz` (+ index) when `--remove-missing`.
  - Job ID/logs: `quilt2_slurm/quilt2_nomiss_job_id.txt`, `quilt2_slurm/quilt2_nomiss_%A_%a.(output|error)`.
- Phase 2 (chunk array): uses `quilt2_output/panel/` if Phase 1 ran; otherwise the original panel dir. Job ID/logs: `quilt2_slurm/quilt2_job_id.txt`, `quilt2_slurm/quilt2_%A_%a.(output|error)`.
- `--dry-run` prints the master sbatch command and exits; Phase 1/2 scripts are still generated for inspection.

Cache awareness (when changing inputs/settings):
- Remove cached chunks: `rm -f quilt2_output/tmp/quilt_auto_chunks.tsv`
- Remove filtered panels: `rm -f quilt2_output/panel/quilt.nomiss.*`
- Remove prepared references: `rm -f quilt2_output/RData/QUILT_prepared_reference.*`

## Troubleshooting
See `quilt2_past_problem_and_solution.md` for fixes on genetic map columns, chr naming, symlinks, chunk parsing, phased panel requirements, and cache invalidation. Use `--dry-run` first to ensure SLURM script generation succeeds before submitting.

## Concordance & dosage r2 (imputed vs truth)
- Utility: `utils/dosage_r2.sh` (bash) + `utils/dosage_r2.R` (R/data.table/ggplot2). Loads `miniforge/25.3.0-3` and bcftools; activates `QUILT2_CONDA_ENV` when available.
- Inputs: imputed VCF with `DS` (or `GP` fallback) and truth VCF with `GT`; both indexed. Samples default to the intersection; supply `--samples` to override.
- Behavior: intersects sites allele-aware (`bcftools isec`), restricts to biallelic SNPs by default, extracts DS/GP and GT, computes per-variant r2 and genotype concordance, plus MAF-binned summaries; optional plots via R.
- Example:
```bash
bash utils/dosage_r2.sh \
  --imputed /path/imputed.vcf.gz \
  --truth /path/truth.vcf.gz \
  --out-prefix results/dosage_eval \
  --region chr1 \
  --samples common_samples.txt
```
