# QUILT2 Pipeline (KH v1)

SLURM-array wrapper around QUILT2 imputation for apple data. Mirrors the Step1C array pattern from the GATK workflow.

## Layout
- `bin/run_quilt2.sh` – orchestrator; builds chunk manifest, generates SLURM array script, self-submits.
- `bin/dosage_r2_sbatch.sh` – SLURM submit wrapper for `modules/evaluate/dosage_r2.sh`.
- `templates/quilt2_job.sh` – Phase 2 array worker; processes one chunk (prepare + impute).
- `templates/quilt2_nomiss_job.sh` – Phase 1 array worker; standardises contig names and/or filters missing genotypes per chromosome.
- `lib/functions.sh` – shared helpers (env bootstrap, bcftools/QUILT checks, panel/map resolution).
- `config/quilt2_config.sh` – SLURM defaults (account/partition/qos/resources/array cap).
- `config/environment.template.sh` – site-specific environment template; copy to `config/environment.sh` and customise.
- `modules/evaluate/dosage_r2.sh` – post-imputation evaluation (concordance & dosage r/r²).
- `modules/evaluate/dosage_r2.R` – R companion script for metrics computation and plots.
- `quilt2_pipeline.legacy.sh` – pre-array monolithic script (kept only for rollback).
- `quilt2_past_problem_and_solution.md` – troubleshooting history.
- `dummy_map.md` – guide for creating dummy genetic maps when a species-specific map is unavailable.

## Prerequisites
- Bash, SLURM.
- `Rscript` with QUILT2 scripts available (`QUILT2_prepare_reference.R`, `QUILT2.R`) via `--quilt2-home` or explicit `--quilt2-*` paths.
- `bcftools` on PATH or loadable via `BCFTOOLS_MODULE` (default `bcftools/1.18-gcc-12.3.0`).
- Optional: conda env name via `QUILT2_CONDA_ENV` (default `quilt2`); module `miniforge/25.3.0-3` loaded if present.

## Configuration (environment.sh)
- Copy `config/environment.template.sh` to `config/environment.sh` and set site defaults (output/scratch roots, reference FASTA, genetic map path, panel dir).
- SLURM defaults: `QUILT2_ACCOUNT`, `QUILT2_PARTITION`, `QUILT2_QOS`, `QUILT2_NODES`, `QUILT2_NTASKS`, `QUILT2_CPUS_PER_TASK`, `QUILT2_MEMORY`, `QUILT2_TIME_LIMIT`, `QUILT2_ARRAY_MAX` (0=no cap; falls back to `QUILT2_ARRAY_LIMIT`), `QUILT2_CONSTRAINT` (optional).
- Tooling: `BCFTOOLS_MODULE`, `QUILT2_CONDA_ENV`, optional `QUILT2_HOME`/`QUILT2_PREP_SCRIPT`/`QUILT2_RUN_SCRIPT`.
- Paths and behavior toggles: `QUILT2_OUTPUT_DIR`, `QUILT2_SCRATCH_DIR`, `QUILT2_CHROMS`, `QUILT2_BUFFER`, `QUILT2_NGEN`, `QUILT2_AUTO_CHUNK_MAP`, `QUILT2_CHUNK_FILE`, `QUILT2_REGION_START/END`, `QUILT2_REMOVE_MISSING`, `QUILT2_MIN_VALID_GT_RATE`, `QUILT2_STANDARDISE_NAME`, `QUILT2_STANDARDISE_NAME_FORCE`, `QUILT2_PREP_ONLY`, `QUILT2_IMPUTE_ONLY`, `QUILT2_DRY_RUN`, `QUILT2_BAMLIST`.

## Inputs
- `--input-dir` (`WORK_DIR`): run directory for outputs, logs, temporary files, and default input discovery. This is **not necessarily** the reference panel directory.
- `--bamlist`: text file listing the low-pass BAMs to impute. It is required for normal runs and `--impute-only`. If omitted, the script searches `WORK_DIR` for `bamlist.txt`, `bamlist.1.0.txt`, then `bamlist.tsv`.
- `--reference-panel-dir`: directory containing the phased reference panel VCFs. This is required, unless `QUILT2_REFERENCE_PANEL_DIR` is set in `config/environment.sh`.
- `--output-dir`: persistent output directory. Defaults to `WORK_DIR/quilt2_output`.
- `--scratch-dir`: optional scratch/staging root. If omitted, SLURM tasks use `$TMPDIR` when available; otherwise they use `OUTPUT_DIR/scratch`. Scratch is for disposable task-local files only.
- `--genetic-map`: genetic map file or directory with per-chromosome maps. Names must match the chromosome names used for the run (`Chr01` vs `1`, etc.). See `dummy_map.md` for creating maps when none are available.
- Optional: `--reference-fasta` with `.fai` index, used during `--standardise-name` when VCF headers need contig repair.

Reference panel requirements:
- Use **VCF/BCF-style reference panel files**, normally compressed as `*.vcf.gz`. Do not pass gVCF files; QUILT2 expects called genotype records, not gVCF reference blocks.
- The panel must be **phased** and contain phased `GT` values such as `0|0`, `0|1`, or `1|0`.
- Panel VCFs should be split or named per chromosome. The script looks for names such as `apple_panel.refpol.Chr01.vcf.gz`, `panel.snps.clean__Chr01.vcf.gz`, `Chr01.vcf.gz`, or matching `Chr01_*.vcf.gz`, `Chr01.*.vcf.gz`, `Chr01-*.vcf.gz`.
- VCFs must be bgzip-compressed and indexed (`.tbi` or `.csi`). The script tries to index `*.vcf.gz`, but pre-indexing avoids cluster-time failures.
- Recommended chromosome naming is `Chr01`-`Chr17`. If the panel uses bare numeric contigs (`1`-`17`), add `--standardise-name` to create `ChrNN` panel VCFs in `OUTPUT_DIR/panel/standardised/`. If the panel uses another convention, pre-standardise it or make sure `--chr`, the panel VCFs, and genetic maps all use the same names.
- If panel variants contain missing or unphased genotypes, use `--remove-missing --min-valid-gt-rate <rate>` to create cleaned per-chromosome panel VCFs before imputation.

## Quick Start
Full run with an already phased, indexed, clean panel:
```bash
bash bin/run_quilt2.sh \
  -i /path/to/quilt2_run_dir \
  --output-dir /path/to/quilt2_output \
  --bamlist /path/to/bamlist.txt \
  --reference-panel-dir /path/to/phased_panel_vcfs \
  --genetic-map /path/to/genetic_maps_dir \
  --auto-chunk-map
```

Full run with panel cleanup first:
```bash
bash bin/run_quilt2.sh \
  -i /path/to/quilt2_run_dir \
  --output-dir /path/to/quilt2_output \
  --bamlist /path/to/bamlist.txt \
  --reference-panel-dir /path/to/phased_panel_vcfs \
  --genetic-map /path/to/genetic_maps_dir \
  --auto-chunk-map \
  --remove-missing \
  --min-valid-gt-rate 0.95
```

Fixed region for all chromosomes:
```bash
bash bin/run_quilt2.sh \
  -i /path/to/quilt2_run_dir \
  --bamlist /path/to/bamlist.txt \
  --reference-panel-dir /path/to/phased_panel_vcfs \
  --genetic-map /path/to/genetic_maps_dir \
  --region-start 1 --region-end 5000000
```

Impute-only (prepared references already exist):
```bash
bash bin/run_quilt2.sh \
  -i /path/to/quilt2_run_dir \
  --bamlist /path/to/bamlist.txt \
  --reference-panel-dir /path/to/phased_panel_vcfs \
  --genetic-map /path/to/genetic_maps_dir \
  --impute-only \
  --region-start 1 --region-end 5000000
```

Dry-run (generate SLURM script, no submission):
```bash
bash bin/run_quilt2.sh ... --dry-run
```

Panel cleanup (`--remove-missing`):
- Default: **off**. The pipeline uses the reference panel as provided unless this flag is set.
- What it does: runs a Phase 1 SLURM array over chromosomes and writes cleaned panel VCFs to `OUTPUT_DIR/panel/nomiss/quilt.nomiss.<chr>.vcf.gz`.
- Filter rule: keep sites with no missing `GT` calls and at least `--min-valid-gt-rate` of samples carrying phased genotypes (`|` in `GT`). Default rate is `0.95`.
- When to use it: use it when the reference panel may contain missing genotypes (`./.` or `.|.`), unphased genotypes (`0/1`), or mixed-quality phased output.
- When to skip it: skip it if your panel is already phased, bgzip-compressed, indexed, and known to be missing-free. This avoids an extra per-chromosome `bcftools` filtering job.
- Reuse behavior: if `OUTPUT_DIR/panel/nomiss/quilt.nomiss.<chr>.vcf.gz` already exists with an index, the cleanup step reuses it instead of re-filtering. Delete `OUTPUT_DIR/panel/nomiss/quilt.nomiss.*` if you changed the source panel or filtering threshold.

Point to QUILT2 scripts:
```bash
  --quilt2-home /path/to/QUILT2_scripts
# or
  --quilt2-prepare-script /path/to/QUILT2_prepare_reference.R \
  --quilt2-run-script /path/to/QUILT2.R
```

Per-chunk evaluation against a truth VCF (uses vcfppR):
```bash
bash bin/run_quilt2.sh \
  -i /path/to/quilt2_run_dir \
  --bamlist /path/to/bamlist.txt \
  --reference-panel-dir /path/to/phased_panel_vcfs \
  --genetic-map /path/to/genetic_maps_dir \
  --region-start 1 --region-end 5000000 \
  --truth-vcf /path/to/truth.vcf.gz \
  --eval-output /path/to/eval_dir
```

SLURM overrides (env, matches Step1C style):
```bash
export QUILT2_ACCOUNT=youracct
export QUILT2_PARTITION=compute
export QUILT2_QOS=normal
export QUILT2_CPUS_PER_TASK=8
export QUILT2_MEMORY=48G
export QUILT2_TIME_LIMIT=12:00:00
export QUILT2_ARRAY_MAX=0   # 0 = no cap
export QUILT2_CONSTRAINT=epyc4
```

## Execution model (self-submit + phases)
- Default: `bin/run_quilt2.sh` self-submits via `sbatch` when not already in SLURM, then exits. Opt out with `--no-submit` or `--submit-self=false`. Master job ID and logs are written under `OUTPUT_DIR/logs/`.
- Phase 1 (panel prep array): runs per chromosome when `--standardise-name` and/or `--remove-missing` are set. Uses `templates/quilt2_nomiss_job.sh`. Outputs:
  - Standardised: `OUTPUT_DIR/panel/standardised/<chr>_chr.vcf.gz` (+ index), skipped unless `--standardise-name-force`.
  - Filtered: `OUTPUT_DIR/panel/nomiss/quilt.nomiss.<chr>.vcf.gz` (+ index) when `--remove-missing`.
  - Job ID/logs: `OUTPUT_DIR/logs/quilt2_nomiss_job_id.txt`, `OUTPUT_DIR/logs/phase1_panel/quilt2_nomiss_%A_%a.(output|error)`.
- Phase 2 (chunk array): uses `panel/nomiss/` if `--remove-missing` ran, `panel/standardised/` if only `--standardise-name` ran, otherwise the original panel dir. Uses `templates/quilt2_job.sh`. Job ID/logs: `OUTPUT_DIR/logs/quilt2_job_id.txt`, `OUTPUT_DIR/logs/phase2_chunks/quilt2_%A_%a.(output|error)`.
- `--dry-run` prints the master sbatch command and exits; Phase 1/2 scripts are still generated for inspection.

Cache awareness (when changing inputs/settings):
- Remove cached chunk manifests: `rm -f OUTPUT_DIR/chunks/manifests/quilt_auto_chunks.tsv`
- Remove filtered panels: `rm -f OUTPUT_DIR/panel/nomiss/quilt.nomiss.*`
- Remove prepared references: `rm -f OUTPUT_DIR/prepared_reference/QUILT_prepared_reference.*`

## Output Layout
Persistent outputs are organised under `OUTPUT_DIR`:

```text
OUTPUT_DIR/
├── panel/
│   ├── standardised/
│   └── nomiss/
├── prepared_reference/
├── chunks/
│   ├── manifests/
│   └── imputed/
│       └── <chr>/
├── eval/
├── logs/
│   ├── scripts/
│   ├── master/
│   ├── phase1_panel/
│   └── phase2_chunks/
└── run_manifest.tsv
```

`--scratch-dir` is only used for task-local staging. Do not point downstream analysis at scratch files; use the persistent files under `OUTPUT_DIR`.

## Troubleshooting
See `quilt2_past_problem_and_solution.md` for fixes on genetic map columns, chr naming, symlinks, chunk parsing, phased panel requirements, and cache invalidation. Use `--dry-run` first to ensure SLURM script generation succeeds before submitting.

## Concordance & dosage r² (imputed vs truth)

Post-imputation evaluation comparing imputed genotypes against a truth (array) VCF.

### Scripts
- `modules/evaluate/dosage_r2.sh` (bash) + `modules/evaluate/dosage_r2.R` (R/data.table/ggplot2).
- `bin/dosage_r2_sbatch.sh` – SLURM submit wrapper (recommended for cluster runs).

### Environment
Loads `miniforge/25.3.0-3` and bcftools modules. Activates conda env `CONDA_ENV` (default `myenv_py310`); override with `CONDA_ENV`, `MINIFORGE_MODULE`, or `BCFTOOLS_MODULE` environment variables.

### How it works
Dosages are derived from GT fields only; DS/GP tags are not used. The pipeline:
1. Normalizes contig names to canonical ChrNN format.
2. Finds overlapping positions (CHROM + POS) between imputed and truth VCFs via position-only intersection.
3. Deduplicates multi-allelic positions (`bcftools norm -d snps`).
4. Removes strand-ambiguous loci (A/T, T/A, C/G, G/C) that cannot be reliably assigned to allele classes.
5. Translates both VCFs to A/B format using AT/CG nucleotide grouping (A,T → group A; C,G → group B).
6. Feeds A/B genotype TSVs to R for per-variant r (signed Pearson) and r², genotype concordance, MAF-binned summaries, per-sample metrics, and diagnostic plots.

### Inputs
- Imputed VCF and truth VCF (both with GT fields; both indexed).
- Samples default to the intersection; supply `--samples` to override.

### Options
- `--region STR` – limit evaluation to a region (e.g., `Chr01:1-1e6`).
- `--use-vcfpp` – additionally run vcfppR comparison on unambiguous VCFs.
- `--no-parquet` – skip writing the concordance Parquet file.
- `--no-biallelic-only` – do not restrict to biallelic SNPs (default: restrict).
- `--no-plots` – skip plotting; only produce metric TSVs.
- `--force` – re-run all steps even if output files already exist.
- `--keep-temp` – do not delete the temporary working directory.

### Outputs (PREFIX.*)
| File | Description |
|------|-------------|
| `metrics.tsv` | Per-variant r (signed Pearson), r², concordance, MAF |
| `per_sample_metrics.tsv` | Per-sample r and r², overall and per 0.1 MAF bin |
| `summary.tsv` | Overall r_mean, r_median, r2_mean, r2_median, concordance_mean |
| `maf_bins.tsv` | r_mean, r2_mean, concordance aggregated by MAF bins |
| `concordance.parquet` | Per-site per-sample concordance (0/1/NA) |
| `IMPUTED_overlapped_only.vcf.gz` | Imputed VCF at common positions (deduped) |
| `TRUTH_overlapped_only.vcf.gz` | Truth VCF at common positions (deduped) |
| `IMPUTED_overlapped_unambiguous_only.vcf.gz` | Imputed VCF after removing ambiguous loci |
| `TRUTH_overlapped_unambiguous_only.vcf.gz` | Truth VCF after removing ambiguous loci |
| `ambiguous_loci_removed.tsv` | Strand-ambiguous positions that were removed |
| `duplicates_removed.{imputed,truth}.tsv` | Duplicate positions removed |
| `imputed.AB_format.tsv` | Imputed genotypes in A/B format |
| `truth.AB_format.tsv` | Truth genotypes in A/B format |
| `translation_exceptions.tsv` | Unexpected GTs found during translation |
| `r2_hist.png` | Distribution of r² |
| `concordance_hist.png` | Distribution of concordance |
| `r2_vs_maf.png` | r² vs MAF heatmap |
| `r2_vs_maf_line.{png,tsv}` | Mean r² vs MAF line plot (0.01 bins) |
| `r2_per_chr_1Mb.{png,tsv}` | Mean r² per 1 Mb window, faceted by chromosome |
| `r2_per_sample.png` | Per-sample r² bar chart (sorted ascending by r²) |

### Examples

Direct execution:
```bash
bash modules/evaluate/dosage_r2.sh \
  --imputed /path/imputed.vcf.gz \
  --truth /path/truth.vcf.gz \
  --out-prefix results/dosage_eval \
  --region Chr01:1-1e6 \
  --samples common_samples.txt
```

Via SLURM (recommended):
```bash
bash bin/dosage_r2_sbatch.sh \
  --imputed /path/imputed.vcf.gz \
  --truth /path/truth.vcf.gz \
  --out-prefix results/dosage_eval \
  -- --region Chr01:1-1e6 --samples common_samples.txt
```
