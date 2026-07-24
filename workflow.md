# QUILT2 Imputation Pipeline — Detailed Workflow

> **Visual overview:** `workflow_diagram.html` in the repository root renders an interactive metromap diagram.
> **Evaluation tools:** `modules/evaluate/concat_imputed.sh` (stitch chunk VCFs) → the array evaluator (`dosage_r2.sh` + `dosage_r2.R`) or WGS evaluator (`dosage_r2_wgs.sh` + `dosage_r2_wgs.R`), plus `utils/test_concordance_check_with_array_validation.sh` for sample-identity matching against array truth.

---

## Table of Contents

1. [Overview](#1-overview)
2. [Repository Structure](#2-repository-structure)
3. [Stage 1 — Panel Preparation (Optional)](#3-stage-1--panel-preparation-optional)
4. [Stage 2 — Chunk Definition](#4-stage-2--chunk-definition)
5. [Stage 3 — Reference Preparation (SLURM array)](#5-stage-3--reference-preparation-slurm-array)
6. [Stage 4 — Imputation (SLURM array)](#6-stage-4--imputation-slurm-array)
7. [Stage 5 — Evaluation (Optional)](#7-stage-5--evaluation-optional)
   - [5a. Contig Name Normalisation](#5a-contig-name-normalisation)
   - [5b. Position-Only Overlap](#5b-position-only-overlap)
   - [5c. Biallelic Filter and Deduplication](#5c-biallelic-filter-and-deduplication)
   - [5d. Strand-Ambiguous Loci Removal](#5d-strand-ambiguous-loci-removal)
   - [5e. A/B Format Translation](#5e-ab-format-translation)
   - [5f. R Metrics and Plots](#5f-r-metrics-and-plots)
   - [5g. Standalone Concordance Check for Array Validation](#5g-standalone-concordance-check-for-array-validation)
8. [Configuration](#8-configuration)
9. [How to Run](#9-how-to-run)
10. [Output Files](#10-output-files)

---

## 1. Overview

This pipeline performs low-pass whole-genome imputation using **QUILT2** on a SLURM cluster. It is designed around the apple genome (17 chromosomes, `Chr01`–`Chr17`) 

The pipeline has five stages:

| Stage | Name | Script | SLURM | Mandatory |

|---|---|---|---|---|
| 1 | Panel Preparation | `bin/run_quilt2.sh` (Phase 1) | Array per chromosome | Optional |
| 2 | Chunk Definition | `bin/run_quilt2.sh` | Local / inline | Mandatory |
| 3 | Reference Preparation | `templates/quilt2_nomiss_job.sh` | Array per chunk | Mandatory |
| 4 | Imputation + Concat | `templates/quilt2_job.sh` | Array per chunk | Mandatory |
| 5 | Evaluation | `bin/dosage_r2_sbatch.sh` → array or WGS evaluator, or `utils/test_concordance_check_with_array_validation.sh` | Single job, chromosome array + finalizer, or local utility | Optional |

Stages 3 and 4 run as a **single SLURM array** where each task processes one chunk of one chromosome. Stage 1 runs as a separate SLURM array over chromosomes. Stage 5 runs as a single SLURM job.

---

## 2. Repository Structure

```
QUILT2_Pipeline_KH_v1/
├── bin/
│   ├── run_quilt2.sh          # Main orchestrator (submit master + array jobs)
│   └── dosage_r2_sbatch.sh    # Evaluation SLURM submit wrapper
├── utils/
│   └── test_concordance_check_with_array_validation.sh  # Standalone all-vs-all sample concordance utility
├── config/
│   ├── environment.template.sh  # Copy to environment.sh and fill in paths
│   ├── environment.sh           # (user-created, not committed)
│   └── quilt2_config.sh         # SLURM resource defaults
├── lib/
│   └── functions.sh             # Shared bash helper functions
├── modules/
│   └── evaluate/
│       ├── concat_imputed.sh    # Stitch per-chunk imputed VCFs into per-chr/genome-wide VCFs
│       ├── dosage_r2.sh         # Evaluation pipeline (bash) — called by bin/dosage_r2_sbatch.sh
│       ├── dosage_r2.R          # Array-mode metrics, r/r², concordance, and plots
│       ├── dosage_r2_wgs.sh     # WGS-truth normalization, filtering, and extraction
│       └── dosage_r2_wgs.R      # Exact-allele GT-to-GT metrics and reports
├── templates/
│   ├── quilt2_job.sh            # SLURM array job template (ref prep + imputation)
│   └── quilt2_nomiss_job.sh     # Variant: filters missing panel variants first
├── workflow.md                  # This file
└── workflow_diagram.html        # Interactive pipeline diagram
```

---

## 3. Stage 1 — Panel Preparation (Optional)

**Script:** `bin/run_quilt2.sh` (Phase 1, triggered by `--standardise-name` and/or `--remove-missing`)
**Execution:** SLURM array, one task per chromosome
**Purpose:** Prepare a clean reference panel before imputation.

### Sub-step 1a — Standardise Chromosome Names

If the reference panel uses numeric chromosome names (e.g., `1`, `2`, ..., `17`) or UCSC-style names (`chr1`), they are renamed to the canonical `ChrNN` format (e.g., `Chr01`, `Chr02`) expected by the downstream pipeline.

```bash
bcftools annotate --rename-chrs <chrmap.txt> panel.vcf.gz | bcftools view -Oz -o panel.ChrNN.vcf.gz
tabix -p vcf panel.ChrNN.vcf.gz
```

The contig map file (`chrmap.txt`) translates between naming conventions. Without this step, QUILT2 may silently fail to find variants on misnamed chromosomes.

**Trigger:** `--standardise-name` (or `QUILT2_STANDARDISE_NAME=true` in environment)

### Sub-step 1b — Filter Low-Quality Sites

Sites where fewer than `--min-valid-gt-rate` (default 0.95, i.e. 95%) of samples have a phased genotype are removed. This prevents rare, sparsely-phased variants from degrading the QUILT2 reference model.

```bash
bcftools view --min-af <rate> panel.vcf.gz | bcftools view -Oz -o panel.nomiss.vcf.gz
tabix -p vcf panel.nomiss.vcf.gz
```

A report of removed sites is written to `quilt2_output/panel/missing_sites_removed.tsv`.

**Trigger:** `--remove-missing` (or `QUILT2_REMOVE_MISSING=true`)

### Outputs

| File | Description |
|---|---|
| `quilt2_output/panel/<chr>.filtered.vcf.gz` | Filtered/renamed panel per chromosome |
| `quilt2_output/panel/missing_sites_removed.tsv` | Sites removed by quality filter |

> If Stage 1 is skipped, the raw panel VCF is used directly in Stage 3.

---

## 4. Stage 2 — Chunk Definition

**Script:** `bin/run_quilt2.sh` (inline, before array submission)
**Execution:** Local (no SLURM job)
**Purpose:** Divide each chromosome into non-overlapping chunks for parallel imputation.

QUILT2 requires genomic coordinates to be chunked because the HMM model becomes memory-intensive for very long regions. Three chunking strategies are available:

### Strategy A — Automatic Chunking (recommended)

Uses `QUILT::quilt_chunk_map()` (an R function from the QUILT package) to select chunk boundaries based on variant density and target chunk size.

```bash
bash bin/run_quilt2.sh --auto-chunk-map ...
```

### Strategy B — User-Supplied Chunk File

A TSV file specifying chunk boundaries explicitly:

```
chr    start    end    [buffer]
Chr01  1        5000000  500000
Chr01  5000001  10000000 500000
...
```

```bash
bash bin/run_quilt2.sh --chunk-file chunks.tsv ...
```

### Strategy C — Single Region

Apply a single start/end coordinate to all chromosomes (useful for testing):

```bash
bash bin/run_quilt2.sh --region-start 1 --region-end 5000000 ...
```

### Output

The chunk manifest is written to `OUTPUT_DIR/chunks/manifests/` and drives the SLURM array size in Stages 3 and 4.

---

## 5. Stage 3 — Reference Preparation (SLURM array)

**Script:** `templates/quilt2_nomiss_job.sh` / `templates/quilt2_job.sh` (Phase 1 of each array task)
**Execution:** SLURM array — one task per (chromosome × chunk)
**Tool:** `QUILT2_prepare_reference.R`

For each chunk, QUILT2 builds a compressed reference haplotype object (`.RData`) from:

- The phased reference panel VCF (from Stage 1 or raw panel)
- The genetic map for the chromosome
- The chunk start/end coordinates and buffer size

```bash
Rscript QUILT2_prepare_reference.R \
  --outputdir OUTPUT_DIR/prepared_reference \
  --chr Chr01 \
  --regionStart 1 \
  --regionEnd 5000000 \
  --buffer 500000 \
  --nGen 100 \
  --reference_vcf_file panel.Chr01.vcf.gz \
  --genetic_map_file genetic_map_Chr01.txt
```

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `--nGen` | 100 | Effective number of generations since divergence from panel |
| `--buffer` | 500000 bp | Overlap region on each side of a chunk |
| `--outputdir` | `OUTPUT_DIR/prepared_reference/` | Passed to QUILT2, which nests `.RData` objects under an `RData/` subdirectory of this path |

**Output:** `OUTPUT_DIR/prepared_reference/RData/QUILT_prepared_reference.Chr01.1.5000000.RData` (one per chunk)

> This phase can be re-used across multiple imputation runs without re-running if the panel is unchanged. Use `--impute-only` to skip it.

---

## 6. Stage 4 — Imputation (SLURM array)

**Script:** `templates/quilt2_nomiss_job.sh` / `templates/quilt2_job.sh` (Phase 2 of each array task)
**Execution:** SLURM array — one task per (chromosome × chunk)
**Tool:** `QUILT2.R` followed by `bcftools concat`

### 4a. Per-Chunk Imputation

Each array task imputes genotypes for all samples listed in `bamlist.txt` for one chunk, using the pre-built `.RData` reference:

```bash
Rscript QUILT2.R \
  --output_filename OUTPUT_DIR/chunks/imputed/Chr01/quilt2.diploid.Chr01.1-5000000.vcf.gz \
  --chr Chr01 \
  --regionStart 1 \
  --regionEnd 5000000 \
  --buffer 500000 \
  --nGen 100 \
  --bamlist bamlist.txt \
  --reference_haplotype_file OUTPUT_DIR/prepared_reference/RData/QUILT_prepared_reference.Chr01.1.5000000.RData
```

Output: a per-chunk VCF with imputed `GT`, `DS` (dosage), and `GP` (genotype probability) fields.

### 4b. Concatenation

Once all chunk jobs for a chromosome complete, per-chunk VCFs are available under `OUTPUT_DIR/chunks/imputed/<chr>/`. Chunk start coordinates are **not** zero-padded, so a naive filename glob/sort does not match genomic order (e.g. `"20000001-..."` sorts before `"5000001-..."`). Use `modules/evaluate/concat_imputed.sh` instead, which orders chunks from the run manifest (falling back to numeric-sorted filename parsing if no manifest is available) and checks for missing chunks before concatenating:

```bash
bash modules/evaluate/concat_imputed.sh \
  --chunks-dir OUTPUT_DIR/chunks/imputed \
  --chr Chr01
```

This prints the resulting VCF path (`OUTPUT_DIR/chunks/imputed/Chr01/imputed.Chr01.vcf.gz`) to stdout; omit `--chr` to concatenate every chromosome and also produce a genome-wide `OUTPUT_DIR/chunks/imputed/imputed.all_chroms.vcf.gz`. See [Stage 5](#7-stage-5--evaluation-optional) for how this is chained directly into evaluation via `bin/dosage_r2_sbatch.sh --chunks-dir`.

**Per-chunk output:** `quilt2.diploid.<chr>.<start>-<end>.vcf.gz` — one VCF per chunk.

---

## 7. Stage 5 — Evaluation (Optional)

**Script:** `bin/dosage_r2_sbatch.sh` → array or WGS evaluator
**Execution:** Single SLURM job for array mode; chromosome array plus finalizer for WGS mode
**Purpose:** Compare the imputed VCF against array truth through the existing GT-to-GT A/B route, or against consolidated WGS truth through exact-allele GT-to-GT ALT counts.

Stage 5 has these entry points:

- `bin/dosage_r2_sbatch.sh --truth-mode array` → `modules/evaluate/dosage_r2.sh` for the existing array A/B comparison.
- `bin/dosage_r2_sbatch.sh --truth-mode wgs` → `modules/evaluate/dosage_r2_wgs.sh` for exact-allele GT-to-GT metrics.
- `utils/test_concordance_check_with_array_validation.sh` for Quarto-ready all-vs-all concordance matching between a nucleotide query VCF and an array truth VCF whose GT indices encode A/B genotype classes.

### How to Submit

```bash
bash bin/dosage_r2_sbatch.sh \
  --truth-mode array \
  --imputed  imputed.Chr01.vcf.gz \
  --truth    array_truth.vcf.gz \
  --out-prefix results/eval/Chr01
```

Additional flags can be passed after `--`:

```bash
bash bin/dosage_r2_sbatch.sh \
  --truth-mode array \
  --imputed  imputed.Chr01.vcf.gz \
  --truth    array_truth.vcf.gz \
  --out-prefix results/eval/Chr01 \
  -- --samples sample_list.txt --use-vcfpp --no-parquet
```

Alternatively, pass `--chunks-dir` (instead of `--imputed`) to go straight from raw Phase 2 chunk output to evaluation in a single submission — `dosage_r2_sbatch.sh` runs `modules/evaluate/concat_imputed.sh` first (see [4b](#4b-concatenation)) and feeds its output into `dosage_r2.sh`:

```bash
bash bin/dosage_r2_sbatch.sh \
  --truth-mode array \
  --chunks-dir OUTPUT_DIR/chunks/imputed \
  --truth      array_truth.vcf.gz
```

`--out-prefix` defaults to `OUTPUT_DIR/eval/dosage_eval` in this mode (`OUTPUT_DIR` resolved from `run_manifest.tsv`, two levels up from `--chunks-dir`). Add `--chr LIST` to restrict to specific chromosomes, or `--concat-force` to re-concatenate even if the concat outputs already exist.

For WGS truth, point directly to the GATK pipeline's `7.Consolidated_VCF` directory:

```bash
bash bin/dosage_r2_sbatch.sh \
  --truth-mode wgs \
  --chunks-dir OUTPUT_DIR/chunks/imputed \
  --truth-dataset-dir ../../7.Consolidated_VCF \
  -- --samples truth_samples.txt
```

WGS mode uses `QUILT2_REFERENCE_FASTA` from `config/environment.sh` by default; `--reference-fasta` is an optional override. It defaults to `OUTPUT_DIR/eval/dosage_eval_wgs`, normalizes both VCF sources against the configured reference, retains biallelic SNPs, and runs `bcftools isec -c none` to require exact `CHROM:POS:REF:ALT` alleles. Both `FORMAT/GT` fields are converted to hard-call ALT counts; `FORMAT/DS` is not used. Site thresholds apply globally; invalid GT or truth GQ/DP failures mask only the affected sample. There is no locus call-rate filter. Filter defaults and overrides live only in `config/environment.sh`, and the effective settings are written to root `run_manifest.tsv` and `qc/filter_summary.tsv`.

The WGS wrapper submits one chromosome task per shared chromosome, capped by `QUILT2_WGS_EVAL_MAX_CONCURRENT_CHROMS` (default 4), followed by an `afterok` finalizer. With `--chunks-dir`, each task concatenates and evaluates only its chromosome. Exact common, imputed-only, and truth-only VCFs remain in task-local `$TMPDIR`; R reads GT only from the exact common pair. The finalizer combines small sufficient-statistic files rather than rereading genome-wide genotype tables. Root `per_sample_metrics.tsv` has the same overall and 0.1-MAF-bin columns as array mode. Variant metrics, per-variant concordance, and audit tables are Snappy-compressed Parquet datasets partitioned under `metrics/`. Retained matches are in `per_variant_metrics`, while rejected exact matches are in `site_filtered_variants`. No wide genotype or dosage matrix is persisted. See the README for the complete layout and Arrow extraction examples. The detailed steps below describe array mode.

In array mode, the truth VCF is expected to use array-style genotype coding where GT index `0` means allele `A` and GT index `1` means allele `B`. The truth-side decoder does not derive A/B labels from truth REF/ALT nucleotides.

The dosage evaluation pipeline (`modules/evaluate/dosage_r2.sh`) runs six sequential steps, each with a **skip-check**: if the output files for a step already exist, the step is skipped automatically. Use `--force` to override all skip-checks.

---

### 5a. Contig Name Normalisation

Both VCFs are inspected for chromosome naming style. If either uses a non-canonical style (e.g., `1`, `chr1`), its contigs are renamed to apple-style `ChrNN` using a dynamically-generated renaming map and `bcftools annotate --rename-chrs`. This ensures CHROM values are comparable across both files before any intersection step.

Supported input styles:

| Style | Example | Normalised to |
|---|---|---|
| ChrNN (canonical) | `Chr01` | `Chr01` (no change) |
| chrN | `chr1` | `Chr01` |
| N (bare numeric) | `1` | `Chr01` |

Only the canonical apple chromosomes `Chr01`–`Chr17` are carried forward into Stage 5. Non-apple contigs are excluded before overlap, concordance, and metric calculation.

---

### 5b. Position-Only Overlap

Before deriving common positions, each normalised VCF is first restricted to the evaluation sample set, optional region, biallelic SNPs, and canonical apple chromosomes. The two filtered VCFs are then reduced to their **intersection by (CHROM, POS)** — deliberately ignoring REF and ALT alleles at this stage.

```bash
# Restrict to the evaluation subset before deriving common positions
bcftools view -S samples.txt -m2 -M2 -v snps imputed.vcf.gz -Oz -o imputed.filtered.vcf.gz
bcftools view -S samples.txt -m2 -M2 -v snps truth.vcf.gz   -Oz -o truth.filtered.vcf.gz

# Extract sorted CHROM:POS lists from both filtered VCFs
bcftools query -f '%CHROM\t%POS\n' imputed.filtered.vcf.gz | sort -k1,1 -k2,2n > imputed.pos
bcftools query -f '%CHROM\t%POS\n' truth.filtered.vcf.gz   | sort -k1,1 -k2,2n > truth.pos

# Find common positions
comm -12 imputed.pos truth.pos > common.pos

# Filter each already-filtered VCF to common positions
bcftools view -T common.pos imputed.filtered.vcf.gz -Oz -o imputed.overlap.vcf.gz
bcftools view -T common.pos truth.filtered.vcf.gz   -Oz -o truth.overlap.vcf.gz
```

> **Why position-only?** The REF and ALT alleles may legitimately differ between a WGS-based imputed VCF and an array-based truth VCF (e.g., REF=A in WGS vs REF=G in the array due to different strand or reference conventions). Comparing only by position avoids discarding valid overlapping variants due to apparent REF/ALT mismatches.

**Checkpoint outputs** (`{eval}/intermediate/vcfs/`):

- `IMPUTED_overlapped_only.vcf.gz`
- `TRUTH_overlapped_only.vcf.gz`

---

### 5c. Deduplication

Each overlapping filtered VCF is then **deduplicated** to retain a single record per (CHROM, POS):

```bash
bcftools norm -d snps imputed.overlap.vcf.gz -Oz -o imputed.dedup.vcf.gz
```

Multi-allelic positions that still appear as duplicate rows after the side-specific filtering are collapsed to one representative record. Removed duplicates are recorded, and the script asserts that the two post-dedup position sets still match before continuing.

| Report file | Contents |
|---|---|
| `{eval}/intermediate/qc/duplicates_removed.imputed.tsv` | CHROM, POS, REF, ALT of duplicates removed from imputed VCF |
| `{eval}/intermediate/qc/duplicates_removed.truth.tsv` | Same for truth VCF |

---

### 5d. Strand-Ambiguous Loci Removal

Variants where REF and ALT are **complementary base pairs** (A/T, T/A, C/G, G/C) are removed from both VCFs:

```bash
# Remove positions where REF/ALT are a complementary pair
bcftools view --exclude 'REF="A" && ALT="T" || REF="T" && ALT="A" || REF="C" && ALT="G" || REF="G" && ALT="C"' ...
```

These "strand-ambiguous" loci cannot be reliably classified as A-group or B-group nucleotides (see Step 5e), so they must be excluded to avoid dosage assignment errors. After this step, every retained variant has **one allele in {A, T} and one allele in {C, G}**.

**Report:** `{eval}/intermediate/qc/ambiguous_loci_removed.tsv` — lists CHROM, POS, REF, ALT of all removed loci.

**Checkpoint outputs** (`{eval}/intermediate/vcfs/`):

- `IMPUTED_overlapped_unambiguous_only.vcf.gz`
- `TRUTH_overlapped_unambiguous_only.vcf.gz`

---

### 5e. A/B Format Translation

The imputed VCF is translated into a unified A/B genotype format. The truth VCF is decoded using array genotype indices, so truth GTs are mapped directly from index codes to `A/A`, `A/B`, `B/A`, and `B/B`.

#### The Rule

For the imputed VCF, each GT index is first decoded to its actual nucleotide using the variant's REF and ALT columns:

```
index 0 → REF nucleotide
index 1 → ALT nucleotide
```

The nucleotide is then grouped:

```
{A, T} → group A
{C, G} → group B
```

This produces:

| Nucleotide pair | A/B genotype | Dosage |
|---|---|---|
| Both in {A, T} | A/A | 0 |
| One {A,T} + one {C,G} | A/B | 1 |
| Both in {C, G} | B/B | 2 |

#### Truth VCF Handling

For the truth VCF:

- GT indices are decoded directly as array genotype classes:
  - `0/0 -> A/A`
  - `0/1 -> A/B`
  - `1/0 -> B/A`
  - `1/1 -> B/B`
- Truth `REF` and `ALT` are retained for filtering, position matching, and strand-ambiguity QC, but they are not used to derive truth A/B genotype labels.

This means Stage 5 assumes the array truth GT field already follows the A/B genotype-class convention used by the original pipeline.

#### Exception Handling

- Missing genotypes (`./. or .|.`) are recorded as `./.` in the output and are **not** treated as exceptions.
- Any imputed GT that cannot be parsed, or whose decoded nucleotide is not in `{A, T, C, G}`, is set to `./.` and recorded in `{eval}/intermediate/qc/translation_exceptions.tsv`.
- Any truth GT containing unsupported allele indices is recorded in `{eval}/intermediate/qc/translation_exceptions.tsv` and emitted with missing genotypes.

**Checkpoint outputs:**

- `{eval}/intermediate/ab/imputed.AB_format.tsv` — tab-separated: CHROM, POS, REF, ALT, ID, then one column per sample with A/A, A/B, or B/B genotypes.
- `{eval}/intermediate/ab/truth.AB_format.tsv` — same tabular layout, decoded from the array truth GT index classes.
- `{eval}/intermediate/qc/translation_exceptions.tsv` — imputed translation exceptions plus truth GT decode exceptions.

---

### 5f. R Metrics

`modules/evaluate/dosage_r2.R` reads the two A/B TSV files and computes per-sample imputation quality metrics.

#### Dosage Conversion

A/B genotypes are converted to numeric dosage:

| Genotype | Dosage |
|---|---|
| A/A | 0 |
| A/B or B/A | 1 |
| B/B | 2 |
| ./. | NA |

#### Per-Variant Metrics

For each variant *i*:

- **Dosage r²**: Pearson correlation squared between imputed and truth dosage vectors across all non-missing samples: `r²ᵢ = cor(imputed_dosageᵢ, truth_dosageᵢ)²`
- **Genotype concordance**: Fraction of samples where `round(imputed_dosage) == round(truth_dosage)`
- **MAF**: Minor allele frequency derived from truth dosages: `maf = min(mean(truth_dosage/2), 1 − mean(truth_dosage/2))`

#### Per-Sample Metrics

For each sample *j*, r² is computed across all variants (treating variant positions as observations): `r²ⱼ = cor(imputed_dosage[all variants, j], truth_dosage[all variants, j])²`. This identifies individual samples with systematically poor imputation quality. A per-0.1-MAF-bin breakdown is also provided in `{eval}/per_sample_metrics.tsv`.

#### Output layout (`EVAL_DIR`, e.g. `OUTPUT_DIR/eval/dosage_eval/`)

```
eval/dosage_eval/
  per_sample_metrics.tsv          # deliverable
  concordance.parquet             # optional (--no-parquet to skip)
  intermediate/
    vcfs/                         # overlapped + unambiguous checkpoints
    ab/                           # A/B genotype TSVs
    qc/                           # audit reports, common_samples.txt, pipeline_audit.tsv
```

| File | Description |
|---|---|
| `{eval}/per_sample_metrics.tsv` | Per-sample r (signed) and r² overall and per 0.1 MAF bin |
| `{eval}/concordance.parquet` | Per-site per-sample concordance matrix (0/1/NA); Arrow Parquet |

> **Breaking change:** pre-folder-layout runs wrote flat `{eval}.*` files. After upgrading, use `--force` once or remove legacy flat files before relying on skip/resume.

---

### 5g. Standalone Concordance Check for Array Validation

When the sequencing/imputed samples and array samples are the same biological individuals but use different sample IDs, use the standalone concordance matcher:

```bash
bash utils/test_concordance_check_with_array_validation.sh \
  --vcf1      imputed.Chr01.vcf.gz \
  --truth     array_truth_ab.vcf.gz \
  --out-prefix results/concordance/Chr01
```

This utility expects:

- `--vcf1`: a standard nucleotide VCF/BCF.
- `--truth` or `--vcf`: an array truth VCF/BCF whose GT indices encode A/B genotype classes.

It reuses the same overlap, deduplication, and strand-ambiguity filtering model as `dosage_r2.sh`, then:

1. translates only the query-side VCF into A/B genotype space,
2. decodes the truth-side array GT indices into A/B genotype classes without using truth REF/ALT for genotype translation,
3. computes all-vs-all pairwise genotype concordance across every `vcf1_sample × truth_sample` pair,
4. ranks the best and second-best truth match for each query sample.

Quarto-oriented outputs:

- `{prefix}.best_matches.tsv`
- `{prefix}.pairwise_concordance.tsv`
- `{prefix}.pipeline_audit.tsv`
- `{prefix}.output_manifest.tsv`

Downloadable intermediate outputs include:

- overlapped and unambiguous VCF checkpoints,
- query and truth A/B TSVs,
- duplicate-removal reports,
- ambiguous-loci report,
- query translation exceptions and truth GT decode exceptions.

---

## 8. Configuration

Configuration is split across two files:

### `config/environment.sh`

Paths, tools, and behavior settings. Copy from `config/environment.template.sh` and fill in:

```bash
# Conda/module setup
export MINIFORGE_MODULE="miniforge/25.3.0-3"
export CONDA_ENV="myenv_py310"
export BCFTOOLS_MODULE="bcftools/1.18-gcc-12.3.0"

# Paths
export QUILT2_HOME="/path/to/QUILT"                # Directory containing QUILT2.R
export QUILT2_GENETIC_MAP="/path/to/genetic_map"   # File or directory of per-chr maps
export QUILT2_REFERENCE_FASTA="/path/to/ref.fa"    # Optional; needed for contig header fixes

# WGS-truth evaluation (ignored by array mode)
export QUILT2_WGS_TRUTH_FILTER_ENABLED="true"
export QUILT2_WGS_TRUTH_MIN_GQ="60"
export QUILT2_WGS_TRUTH_MIN_DP="10"
# The remaining QUAL/QD/SOR/FS/MQ/rank-sum defaults are documented in
# config/environment.template.sh and can be overridden here in the same way.
```

### `config/quilt2_config.sh`

Canonical SLURM resource defaults. Edit this file to change resource requests;
`config/environment.sh` and exported shell variables do not override these values.

| Variable | Default | Description |
|---|---|---|
| `QUILT2_ACCOUNT` | `a_qaafi_cas` | SLURM account |
| `QUILT2_PARTITION` | `general` | SLURM partition |
| `QUILT2_QOS` | empty | Optional SLURM QoS |
| `QUILT2_NODES` | `1` | Nodes per job |
| `QUILT2_NTASKS` | `1` | Tasks per job |
| `QUILT2_CPUS_PER_TASK` | `2` | CPUs per Phase 1 panel-prep task |
| `QUILT2_PHASE2_CPUS_PER_TASK` | `3` | CPUs per Phase 2 imputation chunk task |
| `QUILT2_MEMORY` | `8G` | Memory per array task |
| `QUILT2_TIME_LIMIT` | `72:00:00` | Wall time per Phase 1/2 array task |
| `QUILT2_MASTER_TIME_LIMIT` | `336:00:00` | Wall time for the master polling job |
| `QUILT2_ARRAY_MAX` | `0` (no cap) | Max concurrent array tasks |
| `QUILT2_CONSTRAINT` | `epyc4` | Node constraint (avoids ISA issues) |

---

## 9. How to Run

### Prerequisites

1. Copy and fill in `config/environment.sh`.
2. Ensure the conda environment `myenv_py310` (or your chosen `CONDA_ENV`) contains:
   - `r-base`, `r-data.table`
   - Optional: `r-arrow` (Parquet output), `r-vcfppr`, `r-hexbin`
3. Ensure `bcftools` module is available via Lmod.
4. Ensure any array truth VCF used in Stage 5 has already been harmonised into A/B format before submission.

### Full Run (Panel Prep + Imputation)

```bash
bash bin/run_quilt2.sh \
  --input-dir /data/apple_lowpass \
  --output-dir /data/apple_lowpass/quilt2_output \
  --reference-panel-dir /data/panels/apple_phased \
  --genetic-map /data/maps/apple \
  --auto-chunk-map \
  --standardise-name \
  --remove-missing \
  --min-valid-gt-rate 0.95 \
  --n-gen 100 \
  --buffer 500000
```

### Imputation Only (Skip Panel Prep)

```bash
bash bin/run_quilt2.sh \
  --input-dir /data/apple_lowpass \
  --output-dir /data/apple_lowpass/quilt2_output \
  --reference-panel-dir /data/panels/apple_phased \
  --genetic-map /data/maps/apple \
  --chunk-file chunks.tsv \
  --impute-only
```

### Evaluation Only

```bash
bash bin/dosage_r2_sbatch.sh \
  --imputed  results/imputed.Chr01.vcf.gz \
  --truth    data/truth_array.vcf.gz \
  --out-prefix results/eval/Chr01 \
  -- --samples sample_ids.txt --use-vcfpp
```

### Standalone Concordance Matching Only

```bash
bash utils/test_concordance_check_with_array_validation.sh \
  --vcf1      results/imputed.Chr01.vcf.gz \
  --truth     data/truth_array_ab.vcf.gz \
  --out-prefix results/concordance/Chr01
```

### Dry Run (Preview SLURM Commands)

```bash
bash bin/run_quilt2.sh \
  --input-dir /data/apple_lowpass \
  --output-dir /data/apple_lowpass/quilt2_output \
  --reference-panel-dir /data/panels/apple_phased \
  --genetic-map /data/maps/apple \
  --auto-chunk-map \
  --dry-run
```

---

## 10. Output Files

### Imputation Outputs (`OUTPUT_DIR/`)

```
OUTPUT_DIR/
├── panel/
│   ├── standardised/
│   │   └── <chr>_chr.vcf.gz             # Standardised panel, if requested
│   └── nomiss/
│       └── quilt.nomiss.<chr>.vcf.gz    # Missing-filtered panel, if requested
├── prepared_reference/
│   └── RData/
│       └── QUILT_prepared_reference.<chr>.<start>.<end>.RData
├── chunks/
│   ├── manifests/
│   │   └── quilt2_chunks_<timestamp>.txt
│   └── imputed/
│       └── <chr>/
│           └── quilt2.diploid.<chr>.<start>-<end>.vcf.gz
├── eval/
│   └── dosage_eval/
│       ├── per_sample_metrics.tsv
│       ├── concordance.parquet
│       ├── intermediate/
│       │   ├── vcfs/
│       │   ├── ab/
│       │   └── qc/
│       └── slurm/                # when submitted via dosage_r2_sbatch.sh
├── logs/
│   ├── scripts/
│   ├── master/
│   ├── phase1_panel/
│   └── phase2_chunks/
└── run_manifest.tsv
```

`--scratch-dir` is optional and only used for disposable task staging. If omitted, SLURM tasks use `$TMPDIR` when available; otherwise they use `OUTPUT_DIR/scratch`.

### Evaluation Outputs (`EVAL_DIR/`)

See [Stage 5f](#5f-r-metrics). `--out-prefix` is the eval run directory (e.g. `OUTPUT_DIR/eval/dosage_eval`). Deliverables (`per_sample_metrics.tsv`, `concordance.parquet`) sit at that root; checkpoints and QC files live under `intermediate/{vcfs,ab,qc}/`.

### Standalone Concordance Outputs (`{prefix}.*`)

For `utils/test_concordance_check_with_array_validation.sh`, the main Quarto-facing outputs are:

- `{prefix}.best_matches.tsv`
- `{prefix}.pairwise_concordance.tsv`
- `{prefix}.pipeline_audit.tsv`
- `{prefix}.output_manifest.tsv`

The utility also keeps overlapped/unambiguous VCF checkpoints and intermediate TSVs for download and audit.

### SLURM Logs

SLURM logs and generated job scripts are written under `OUTPUT_DIR/logs/`.

For evaluation jobs:

```
<out_dir>/slurm/
├── dosage_r2_<jobid>.out
└── dosage_r2_<jobid>.err
```
