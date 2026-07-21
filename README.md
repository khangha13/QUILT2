# QUILT2 Pipeline (KH v1)

SLURM-array wrapper around QUILT2 imputation for apple data. Mirrors the Step1C array pattern from the GATK workflow.

## Layout
- `bin/run_quilt2.sh` – orchestrator; builds chunk manifest, generates SLURM array script, self-submits.
- `bin/dosage_r2_sbatch.sh` – SLURM submit wrapper for array or WGS truth evaluation.
- `templates/quilt2_job.sh` – Phase 2 array worker; processes one chunk (prepare + impute).
- `templates/quilt2_nomiss_job.sh` – Phase 1 array worker; standardises contig names and/or filters missing genotypes per chromosome.
- `lib/functions.sh` – shared helpers (env bootstrap, bcftools/QUILT checks, panel/map resolution).
- `config/quilt2_config.sh` – SLURM defaults (account/partition/qos/resources/array cap).
- `config/environment.template.sh` – site-specific environment template; copy to `config/environment.sh` and customise.
- `modules/evaluate/concat_imputed.sh` – stitches per-chunk imputed VCFs into per-chromosome/genome-wide VCFs.
- `modules/evaluate/dosage_r2.sh` + `dosage_r2.R` – array-truth GT-to-GT evaluation.
- `modules/evaluate/dosage_r2_wgs.sh` + `dosage_r2_wgs.R` – WGS-truth DS-to-GT-dosage evaluation.
- `quilt2_pipeline.legacy.sh` – pre-array monolithic script (kept only for rollback).
- `quilt2_past_problem_and_solution.md` – troubleshooting history.
- `dummy_map.md` – guide for creating dummy genetic maps when a species-specific map is unavailable.

## Prerequisites
- Bash, SLURM.
- `Rscript` with QUILT2 scripts available (`QUILT2_prepare_reference.R`, `QUILT2.R`) via `--quilt2-home` or explicit `--quilt2-*` paths.
- `bcftools` on PATH or loadable via `BCFTOOLS_MODULE` (default `bcftools/1.18-gcc-12.3.0`).
- Optional: conda env name via `QUILT2_CONDA_ENV` (default `quilt2`); module `miniforge/25.3.0-3` loaded if present.

## Configuration
- Copy `config/environment.template.sh` to `config/environment.sh` and set site defaults for paths, tools, and behavior toggles (output/scratch roots, reference FASTA, genetic map path, panel dir).
- Edit `config/quilt2_config.sh` for SLURM resource defaults: `QUILT2_ACCOUNT`, `QUILT2_PARTITION`, `QUILT2_QOS`, `QUILT2_NODES`, `QUILT2_NTASKS`, `QUILT2_CPUS_PER_TASK`, `QUILT2_PHASE2_CPUS_PER_TASK`, `QUILT2_MEMORY`, `QUILT2_TIME_LIMIT`, `QUILT2_MASTER_TIME_LIMIT`, `QUILT2_ARRAY_MAX`, `QUILT2_CONSTRAINT`.
- Tooling: `BCFTOOLS_MODULE`, `QUILT2_CONDA_ENV`, optional `QUILT2_HOME`/`QUILT2_PREP_SCRIPT`/`QUILT2_RUN_SCRIPT`.
- Paths and behavior toggles: `QUILT2_OUTPUT_DIR`, `QUILT2_SCRATCH_DIR`, `QUILT2_CHROMS`, `QUILT2_BUFFER`, `QUILT2_NGEN`, `QUILT2_AUTO_CHUNK_MAP`, `QUILT2_CHUNK_FILE`, `QUILT2_REGION_START/END`, `QUILT2_REMOVE_MISSING`, `QUILT2_MIN_VALID_GT_RATE`, `QUILT2_STANDARDISE_NAME`, `QUILT2_STANDARDISE_NAME_FORCE`, `QUILT2_PREP_ONLY`, `QUILT2_IMPUTE_ONLY`, `QUILT2_DRY_RUN`, `QUILT2_BAMLIST`, and the `QUILT2_WGS_TRUTH_*` filters.

## Inputs
- `--input-dir` (`WORK_DIR`): run directory for outputs, logs, temporary files, and default input discovery. This is **not necessarily** the reference panel directory.
- `--bamlist`: text file listing the low-pass BAMs to impute. It is required for normal runs and `--impute-only`. If omitted, the script searches `WORK_DIR` for `bamlist.txt`, `bamlist.1.0.txt`, then `bamlist.tsv`.
- `--reference-panel-dir`: directory containing the phased reference panel VCFs. This is required, unless `QUILT2_REFERENCE_PANEL_DIR` is set in `config/environment.sh`.
- `--output-dir`: persistent output directory. Defaults to `WORK_DIR/quilt2_output`.
- `--scratch-dir`: optional scratch/staging root. If omitted, SLURM tasks use `$TMPDIR` when available; otherwise they use `OUTPUT_DIR/scratch`. Scratch is for disposable task-local files only.
- `--genetic-map`: genetic map file or directory with per-chromosome maps. Names must match the chromosome names used for the run (`Chr01` vs `1`, etc.). Pass `dummy` to auto-generate constant 1.0 cM/Mb maps into `OUTPUT_DIR/genetic_map/dummy/` (requires `--reference-fasta` with a `.fai` index). See `dummy_map.md` for details and for creating maps manually when none are available.
- Optional: `--reference-fasta` with `.fai` index, used during `--standardise-name` when VCF headers need contig repair.

Reference panel requirements:
- Use **VCF/BCF-style reference panel files**, normally compressed as `*.vcf.gz`. Do not pass gVCF files; QUILT2 expects called genotype records, not gVCF reference blocks.
- The panel must be **phased** and contain phased `GT` values such as `0|0`, `0|1`, or `1|0`.
- Panel VCFs should be split or named per chromosome. The script looks for names such as `apple_panel.refpol.Chr01.vcf.gz`, `panel.snps.clean__Chr01.vcf.gz`, `Chr01.vcf.gz`, or matching `Chr01_*.vcf.gz`, `Chr01.*.vcf.gz`, `Chr01-*.vcf.gz`.
- VCFs must be bgzip-compressed and indexed (`.tbi` or `.csi`). The script tries to index `*.vcf.gz`, but pre-indexing avoids cluster-time failures.
- Recommended chromosome naming is `Chr01`-`Chr17`. If the panel uses bare numeric contigs (`1`-`17`), the pipeline auto-detects this (by peeking at the first contig of each chromosome's panel VCF) and automatically renames them into `ChrNN` panel VCFs in `OUTPUT_DIR/panel/standardised/`. Pass `--standardise-name` to force this on regardless of detection, or `--no-standardise-name` to disable detection and always use the panel as-is. If the panel uses another convention, pre-standardise it or make sure `--chr`, the panel VCFs, and genetic maps all use the same names.
- If panel variants contain missing or unphased genotypes, use `--remove-missing --min-valid-gt-rate <rate>` to create cleaned per-chromosome panel VCFs before imputation.

## Worked Example
This example uses a fictional low-pass dataset named `Apple_LowPass_2026`.
Assume the inputs are organised as follows:

```text
/QRISdata/Q8367/WGS_Reference_Panel/Apple_LowPass_2026/
├── quilt2_run/
│   └── bamlist.txt
└── quilt2_output/                  # created by the pipeline

/QRISdata/Q8367/Reference_Panels/apple_phased/
├── Chr01.vcf.gz                    # plus Chr01.vcf.gz.tbi
├── Chr02.vcf.gz                    # plus Chr02.vcf.gz.tbi
└── ...                             # through Chr17

/QRISdata/Q8367/Genetic_Maps/apple/
├── Chr01.txt
├── Chr02.txt
└── ...                             # through Chr17
```

`bamlist.txt` contains one absolute BAM path per line, for example:

```text
/QRISdata/Q8367/WGS_Reference_Panel/Apple_LowPass_2026/4.BAM/Gala_01/Gala_01.bam
/QRISdata/Q8367/WGS_Reference_Panel/Apple_LowPass_2026/4.BAM/Fuji_02/Fuji_02.bam
```

Run the following commands from the `QUILT2_Pipeline_KH_v1` directory. Replace
the fictional paths with your own paths.

```bash
RUN_DIR=/QRISdata/Q8367/WGS_Reference_Panel/Apple_LowPass_2026/quilt2_run
OUTPUT_DIR=/QRISdata/Q8367/WGS_Reference_Panel/Apple_LowPass_2026/quilt2_output
PANEL_DIR=/QRISdata/Q8367/Reference_Panels/apple_phased
MAP_DIR=/QRISdata/Q8367/Genetic_Maps/apple

# 1. Validate the inputs and generate the Phase 1/2 SLURM scripts without
#    submitting jobs. --no-submit keeps the dry-run in the current shell.
bash bin/run_quilt2.sh \
  --input-dir "${RUN_DIR}" \
  --output-dir "${OUTPUT_DIR}" \
  --bamlist "${RUN_DIR}/bamlist.txt" \
  --reference-panel-dir "${PANEL_DIR}" \
  --genetic-map "${MAP_DIR}" \
  --auto-chunk-map \
  --no-submit --dry-run

# 2. Submit the full Chr01-Chr17 run using map-based automatic chunks.
bash bin/run_quilt2.sh \
  --input-dir "${RUN_DIR}" \
  --output-dir "${OUTPUT_DIR}" \
  --bamlist "${RUN_DIR}/bamlist.txt" \
  --reference-panel-dir "${PANEL_DIR}" \
  --genetic-map "${MAP_DIR}" \
  --auto-chunk-map

# 3. Alternatively, run a small Chr01 pilot over bases 1-5,000,000.
bash bin/run_quilt2.sh \
  --input-dir "${RUN_DIR}" \
  --output-dir "${OUTPUT_DIR}/Chr01_pilot" \
  --bamlist "${RUN_DIR}/bamlist.txt" \
  --reference-panel-dir "${PANEL_DIR}" \
  --genetic-map "${MAP_DIR}" \
  --chr Chr01 \
  --region-start 1 --region-end 5000000
```

If the phased panel contains missing or unphased genotypes, add
`--remove-missing --min-valid-gt-rate 0.95`. If no species-specific genetic map
is available, replace `--genetic-map "${MAP_DIR}"` with
`--genetic-map dummy --reference-fasta /path/to/reference.fasta`; the FASTA must
have a matching `.fai` index.

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

SLURM resources (edit `config/quilt2_config.sh`):
```bash
QUILT2_ACCOUNT="youracct"
QUILT2_PARTITION="compute"
QUILT2_QOS="normal"
QUILT2_CPUS_PER_TASK="8"
QUILT2_PHASE2_CPUS_PER_TASK="3"
QUILT2_MEMORY="48G"
QUILT2_TIME_LIMIT="12:00:00"
QUILT2_ARRAY_MAX="0"   # 0 = no cap
QUILT2_CONSTRAINT="epyc4"
```

## Execution model (self-submit + phases)
- Default: `bin/run_quilt2.sh` self-submits via `sbatch` when not already in SLURM, then exits. Opt out with `--no-submit` or `--submit-self=false`. Master job ID and logs are written under `OUTPUT_DIR/logs/`.
- Phase 1 (panel prep array): runs per chromosome when `--standardise-name` and/or `--remove-missing` are set. Uses `templates/quilt2_nomiss_job.sh`. Outputs:
  - Standardised: `OUTPUT_DIR/panel/standardised/<chr>_chr.vcf.gz` (+ index), skipped unless `--standardise-name-force`.
  - Filtered: `OUTPUT_DIR/panel/nomiss/quilt.nomiss.<chr>.vcf.gz` (+ index) when `--remove-missing`.
  - Job ID/logs: `OUTPUT_DIR/logs/quilt2_nomiss_job_id.txt`, `OUTPUT_DIR/logs/phase1_panel/quilt2_nomiss_%A_%a.(output|error)`.
- Phase 2 (chunk array): uses `panel/nomiss/` if `--remove-missing` ran, `panel/standardised/` if only `--standardise-name` ran, otherwise the original panel dir. Uses `templates/quilt2_job.sh`. Job ID/logs: `OUTPUT_DIR/logs/quilt2_job_id.txt`, `OUTPUT_DIR/logs/phase2_chunks/quilt2_%A_%a.(output|error)`.
- Phase 1 short-circuit: on a rerun, if every chromosome already has its expected standardised/filtered output (+ index), Phase 1 submission is skipped entirely (no `sbatch`, no wait). Use `--standardise-name-force` to force Phase 1 to resubmit and rebuild outputs regardless.
- `--dry-run` prints the master sbatch command and exits; Phase 1/2 scripts are still generated for inspection.

Cache awareness (when changing inputs/settings):
- Remove cached chunk manifests: `rm -f OUTPUT_DIR/chunks/manifests/quilt_auto_chunks.tsv`
- Remove filtered panels: `rm -f OUTPUT_DIR/panel/nomiss/quilt.nomiss.*`
- Remove prepared references: `rm -f OUTPUT_DIR/prepared_reference/RData/QUILT_prepared_reference.*`

## Output Layout
Persistent outputs are organised under `OUTPUT_DIR`:

```text
OUTPUT_DIR/
├── panel/
│   ├── standardised/
│   └── nomiss/
├── prepared_reference/
│   └── RData/
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

Post-imputation evaluation has two explicit truth modes:

- `array` (default) preserves the existing GT-to-GT A/B workflow.
- `wgs` compares QUILT2 `FORMAT/DS` with ALT dosage derived from filtered GATK `FORMAT/GT`.

### Scripts
- `modules/evaluate/dosage_r2.sh` (bash) + `modules/evaluate/dosage_r2.R` (R/data.table) for array mode.
- `modules/evaluate/dosage_r2_wgs.sh` (bash) + `modules/evaluate/dosage_r2_wgs.R` (base R) for WGS mode.
- `modules/evaluate/concat_imputed.sh` – stitches per-chunk imputed VCFs (`OUTPUT_DIR/chunks/imputed/<chr>/quilt2.diploid.<chr>.<start>-<end>.vcf.gz`) into per-chromosome and genome-wide VCFs, ordered via the run manifest (or numeric filename sort as a fallback). Adjacent chunks may overlap (e.g. `--auto-chunk-map`); each chunk is trimmed to end just before the next chunk's start before concatenating, so `bcftools index` on the result doesn't fail with "unsorted positions".
- `bin/dosage_r2_sbatch.sh` – SLURM submit wrapper (recommended for cluster runs); accepts `--chunks-dir` as an alternative to `--imputed` to chain concatenation directly into evaluation.

### Environment
Loads `miniforge/25.3.0-3` and bcftools modules. Activates conda env `CONDA_ENV` (default `myenv_py310`); override with `CONDA_ENV`, `MINIFORGE_MODULE`, or `BCFTOOLS_MODULE` environment variables.

### Array mode: how it works
Dosages are derived from GT fields only; DS/GP tags are not used in array mode. The pipeline:
1. Normalizes contig names to canonical ChrNN format.
2. Finds overlapping positions (CHROM + POS) between imputed and truth VCFs via position-only intersection.
3. Deduplicates multi-allelic positions (`bcftools norm -d snps`).
4. Removes strand-ambiguous loci (A/T, T/A, C/G, G/C) that cannot be reliably assigned to allele classes.
5. Translates both VCFs to A/B format using AT/CG nucleotide grouping (A,T → group A; C,G → group B).
6. Feeds A/B genotype TSVs to R for per-sample r² (overall and per 0.1 MAF bin).

### Array mode inputs
- Imputed VCF and truth VCF (both with GT fields; both indexed).
- Samples default to the intersection; supply `--samples` to override.

### Array mode options
- `--region STR` – limit evaluation to a region (e.g., `Chr01:1-1e6`).
- `--use-vcfpp` – additionally run vcfppR comparison on unambiguous VCFs.
- `--no-parquet` – skip writing the concordance Parquet file.
- `--no-biallelic-only` – do not restrict to biallelic SNPs (default: restrict).
- `--force` – re-run all steps even if output files already exist.
- `--keep-temp` – do not delete the temporary working directory.

### Outputs (`EVAL_DIR/`, e.g. `eval/dosage_eval/`)

Deliverables at the eval run root:

| File | Description |
|------|-------------|
| `per_sample_metrics.tsv` | Per-sample r and r², overall and per 0.1 MAF bin |
| `concordance.parquet` | Per-site per-sample concordance (0/1/NA); skip with `--no-parquet` |

Intermediates under `intermediate/`:

| Path | Description |
|------|-------------|
| `intermediate/vcfs/IMPUTED_overlapped_only.vcf.gz` | Imputed VCF at common positions (deduped) |
| `intermediate/vcfs/TRUTH_overlapped_only.vcf.gz` | Truth VCF at common positions (deduped) |
| `intermediate/vcfs/IMPUTED_overlapped_unambiguous_only.vcf.gz` | Imputed after removing ambiguous loci |
| `intermediate/vcfs/TRUTH_overlapped_unambiguous_only.vcf.gz` | Truth after removing ambiguous loci |
| `intermediate/ab/imputed.AB_format.tsv` | Imputed genotypes in A/B format |
| `intermediate/ab/truth.AB_format.tsv` | Truth genotypes in A/B format |
| `intermediate/qc/common_samples.txt` | Sample IDs evaluated |
| `intermediate/qc/duplicates_removed.{imputed,truth}.tsv` | Duplicate positions removed |
| `intermediate/qc/ambiguous_loci_removed.tsv` | Strand-ambiguous positions removed |
| `intermediate/qc/translation_exceptions.tsv` | Unexpected GTs during A/B decoding |
| `intermediate/qc/pipeline_audit.tsv` | Step log (timestamp, counts) |

> **Breaking change:** older runs wrote flat `{prefix}.*` files next to the eval directory name. After upgrading, re-run with `--force` or delete legacy flat files; skip checks use the new `intermediate/` layout.

### Examples

Array mode, direct execution:
```bash
bash modules/evaluate/dosage_r2.sh \
  --imputed /path/imputed.vcf.gz \
  --truth /path/truth.vcf.gz \
  --out-prefix results/dosage_eval \
  --region Chr01:1-1e6 \
  --samples common_samples.txt
```

Array mode via SLURM (recommended; `--truth-mode array` is optional because it is the default):
```bash
bash bin/dosage_r2_sbatch.sh \
  --truth-mode array \
  --imputed /path/imputed.vcf.gz \
  --truth /path/truth.vcf.gz \
  --out-prefix results/dosage_eval \
  -- --region Chr01:1-1e6 --samples common_samples.txt
```

Straight from raw chunk output (concatenates via `modules/evaluate/concat_imputed.sh` first, then evaluates):
```bash
bash bin/dosage_r2_sbatch.sh \
  --truth-mode array \
  --chunks-dir OUTPUT_DIR/chunks/imputed \
  --truth /path/truth.vcf.gz
```
`--out-prefix` defaults to `OUTPUT_DIR/eval/dosage_eval` in this mode; add `--chr LIST` to restrict chromosomes or `--concat-force` to re-concatenate existing outputs.

### WGS truth mode

WGS mode discovers numerically ordered, indexed `Chr*_consolidated.vcf.gz` files in the GATK pipeline's `7.Consolidated_VCF` directory. It normalizes both sides against the same reference, keeps original biallelic SNP records, and matches exact `CHROM:POS:REF:ALT` alleles. `Chr00` is excluded unless it is explicitly requested and present on both sides.

For each retained exact match, QUILT2 `DS` is compared with truth ALT dosage (`0/0 = 0`, heterozygous = 1, `1/1 = 2`). Site failures remove the site for all samples. A missing/invalid truth GT or a failing GQ/DP masks only that sample at that site; no call-rate filter is applied. Disabling WGS filtering skips the configurable site, GQ, and DP thresholds but still requires biallelic SNPs, exact allele matches, valid diploid truth GTs, and valid imputed DS values in `[0,2]`.

Configure WGS filtering in `config/environment.sh` (copied from `environment.template.sh`). Defaults are `QUAL >= 30`, `QD >= 2`, `SOR <= 3`, `FS <= 60`, `MQ >= 40`, `MQRankSum >= -12.5`, `ReadPosRankSum >= -8`, `GQ >= 60`, and `DP >= 10`. Missing QUAL, QD, or MQ fails a site; missing rank-sum annotations are allowed. Set `QUILT2_WGS_TRUTH_FILTER_ENABLED=false` to skip the configurable thresholds. These settings are ignored in array mode and are not duplicated as command-line options.

WGS mode from concatenated output:

```bash
bash bin/dosage_r2_sbatch.sh \
  --truth-mode wgs \
  --imputed OUTPUT_DIR/chunks/imputed/imputed.all_chroms.vcf.gz \
  --truth-dataset-dir ../../7.Consolidated_VCF \
  --out-prefix OUTPUT_DIR/eval/dosage_eval_wgs \
  -- --samples truth_samples.txt --region Chr01:1-10000000
```

WGS mode straight from chunk output:

```bash
bash bin/dosage_r2_sbatch.sh \
  --truth-mode wgs \
  --chunks-dir OUTPUT_DIR/chunks/imputed \
  --truth-dataset-dir ../../7.Consolidated_VCF
```

Both examples use `QUILT2_REFERENCE_FASTA` from `config/environment.sh`. Pass `--reference-fasta` only when a run needs to override that configured reference.

The default WGS output is `OUTPUT_DIR/eval/dosage_eval_wgs`. It contains `per_variant_metrics.tsv`, `per_sample_metrics.tsv`, `filter_summary.tsv`, `genotype_masking_summary.tsv`, variant matching reports, a settings `run_manifest.tsv`, and the two matrices `intermediate/imputed_ds.tsv` and `intermediate/truth_gt_dosage.tsv`. Both metric tables retain signed Pearson `r` and `r²`; per-variant values are `NA` below three usable sample pairs or when either dosage vector has zero variance.

On Bunya, run the synthetic acceptance test after loading the same modules used for evaluation:

```bash
bash tests/test_dosage_r2_wgs.sh
```
