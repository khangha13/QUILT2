# QUILT2 Imputation Pipeline — Detailed Workflow

> **Visual overview:** `workflow_diagram.html` in the repository root renders an interactive metromap diagram.
> **Evaluation tools:** `modules/evaluate/dosage_r2.sh` → `modules/evaluate/dosage_r2.R`, plus `utils/test_concordance_check_with_array_validation.sh` for sample-identity matching against array truth.

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

This pipeline performs low-pass whole-genome imputation using **QUILT2** on a SLURM cluster. It is designed around the apple genome (17 chromosomes, `Chr01`–`Chr17`) but is generalisable to any organism.

The pipeline has five stages:

| Stage | Name | Script | SLURM | Mandatory |
|---|---|---|---|---|
| 1 | Panel Preparation | `bin/run_quilt2.sh` (Phase 1) | Array per chromosome | Optional |
| 2 | Chunk Definition | `bin/run_quilt2.sh` | Local / inline | Mandatory |
| 3 | Reference Preparation | `templates/quilt2_nomiss_job.sh` | Array per chunk | Mandatory |
| 4 | Imputation + Concat | `templates/quilt2_job.sh` | Array per chunk | Mandatory |
| 5 | Evaluation | `bin/dosage_r2_sbatch.sh` → `modules/evaluate/dosage_r2.sh`, or `utils/test_concordance_check_with_array_validation.sh` | Single job / local utility | Optional |

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
│       ├── dosage_r2.sh         # Evaluation pipeline (bash) — called by bin/dosage_r2_sbatch.sh
│       └── dosage_r2.R          # R metrics, r/r², concordance, and plots
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

The chunk manifest is written to `quilt2_output/tmp/chunks.tsv` and drives the SLURM array size in Stages 3 and 4.

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
  --outputdir quilt2_output/RData \
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
| `--outputdir` | `quilt2_output/RData/` | Where `.RData` objects are saved |

**Output:** `quilt2_output/RData/RData.Chr01.1.5000000.RData` (one per chunk)

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
  --outputdir quilt2_output/Chr01/1_5000000 \
  --chr Chr01 \
  --regionStart 1 \
  --regionEnd 5000000 \
  --buffer 500000 \
  --nGen 100 \
  --bamlist bamlist.txt \
  --reference_haplotype_file quilt2_output/RData/RData.Chr01.1.5000000.RData
```

Output: a per-chunk VCF with imputed `GT`, `DS` (dosage), and `GP` (genotype probability) fields.

### 4b. Concatenation

Once all chunk jobs for a chromosome complete, the orchestrator concatenates chunk VCFs into a final per-chromosome imputed VCF:

```bash
bcftools concat --naive -Oz -o imputed.Chr01.vcf.gz quilt2_output/Chr01/*/*.vcf.gz
tabix -p vcf imputed.Chr01.vcf.gz
```

**Final output:** `imputed.<chr>.vcf.gz` — one VCF per chromosome, containing imputed genotypes for all samples across all chunks.

---

## 7. Stage 5 — Evaluation (Optional)

**Script:** `bin/dosage_r2_sbatch.sh` → `modules/evaluate/dosage_r2.sh` → `modules/evaluate/dosage_r2.R`
**Execution:** Single SLURM job
**Purpose:** Compare the imputed VCF against a harmonised truth VCF to quantify imputation accuracy, and optionally perform all-vs-all sample identity matching when the array and sequencing sample names differ.

Stage 5 currently has two entry points:

- `bin/dosage_r2_sbatch.sh` → `modules/evaluate/dosage_r2.sh` for dosage r², concordance, and diagnostic plots.
- `utils/test_concordance_check_with_array_validation.sh` for Quarto-ready all-vs-all concordance matching between a nucleotide query VCF and an already A/B-encoded array truth VCF.

### How to Submit

```bash
bash bin/dosage_r2_sbatch.sh \
  --imputed  imputed.Chr01.vcf.gz \
  --truth    array_truth.vcf.gz \
  --out-prefix results/eval/Chr01
```

Additional flags can be passed after `--`:

```bash
bash bin/dosage_r2_sbatch.sh \
  --imputed  imputed.Chr01.vcf.gz \
  --truth    array_truth.vcf.gz \
  --out-prefix results/eval/Chr01 \
  -- --samples sample_list.txt --use-vcfpp --no-parquet
```

The evaluation truth VCF is expected to already be encoded in harmonised A/B format before it enters Stage 5. Raw vendor A/B exports are not translated here.

The dosage evaluation pipeline (`modules/evaluate/dosage_r2.sh`) runs six sequential steps, each with a **skip-check**: if the output files for a step already exist, the step is skipped automatically. Use `--force` to override all skip-checks.

---

### 5a. Contig Name Normalisation

Both VCFs are inspected for chromosome naming style. If either uses a non-canonical style (e.g., `1`, `chr1`), its contigs are renamed to `ChrNN` using a dynamically-generated renaming map and `bcftools annotate --rename-chrs`. This ensures CHROM values are comparable across both files before any intersection step.

Supported input styles:

| Style | Example | Normalised to |
|---|---|---|
| ChrNN (canonical) | `Chr01` | `Chr01` (no change) |
| chrN | `chr1` | `Chr01` |
| N (bare numeric) | `1` | `Chr01` |

---

### 5b. Position-Only Overlap

The two normalised VCFs are filtered to their **intersection by (CHROM, POS)** — deliberately ignoring REF and ALT alleles at this stage.

```bash
# Extract sorted CHROM:POS lists from both VCFs
bcftools query -f '%CHROM\t%POS\n' imputed.vcf.gz | sort -k1,1 -k2,2n > imputed.pos
bcftools query -f '%CHROM\t%POS\n' truth.vcf.gz   | sort -k1,1 -k2,2n > truth.pos

# Find common positions
comm -12 imputed.pos truth.pos > common.pos

# Filter each VCF to common positions (targets mode avoids pulling in multi-allelic neighbours)
bcftools view -T common.pos imputed.vcf.gz -Oz -o imputed.overlap.vcf.gz
bcftools view -T common.pos truth.vcf.gz   -Oz -o truth.overlap.vcf.gz
```

> **Why position-only?** The REF and ALT alleles may legitimately differ between a WGS-based imputed VCF and an array-based truth VCF (e.g., REF=A in WGS vs REF=G in the array due to different strand or reference conventions). Comparing only by position avoids discarding valid overlapping variants due to apparent REF/ALT mismatches.

**Checkpoint outputs:**

- `{prefix}.IMPUTED_overlapped_only.vcf.gz`
- `{prefix}.TRUTH_overlapped_only.vcf.gz`

---

### 5c. Biallelic Filter and Deduplication

Each overlapping VCF is filtered to **biallelic SNPs** only, then **deduplicated** to retain a single record per (CHROM, POS):

```bash
bcftools view -m2 -M2 -v snps imputed.overlap.vcf.gz | bcftools norm -d snps -Oz -o imputed.biallelic.vcf.gz
```

Multi-allelic positions (where two different ALT alleles are split into separate rows after decomposition) are collapsed to one representative record. Removed duplicates are recorded:

| Report file | Contents |
|---|---|
| `{prefix}.duplicates_removed.imputed.tsv` | CHROM, POS, REF, ALT of duplicates removed from imputed VCF |
| `{prefix}.duplicates_removed.truth.tsv` | Same for truth VCF |

---

### 5d. Strand-Ambiguous Loci Removal

Variants where REF and ALT are **complementary base pairs** (A/T, T/A, C/G, G/C) are removed from both VCFs:

```bash
# Remove positions where REF/ALT are a complementary pair
bcftools view --exclude 'REF="A" && ALT="T" || REF="T" && ALT="A" || REF="C" && ALT="G" || REF="G" && ALT="C"' ...
```

These "strand-ambiguous" loci cannot be reliably classified as A-group or B-group nucleotides (see Step 5e), so they must be excluded to avoid dosage assignment errors. After this step, every retained variant has **one allele in {A, T} and one allele in {C, G}**.

**Report:** `{prefix}.ambiguous_loci_removed.tsv` — lists CHROM, POS, REF, ALT of all removed loci.

**Checkpoint outputs:**

- `{prefix}.IMPUTED_overlapped_unambiguous_only.vcf.gz`
- `{prefix}.TRUTH_overlapped_unambiguous_only.vcf.gz`

---

### 5e. A/B Format Translation

The imputed VCF is translated into a unified A/B genotype format. The truth VCF is already expected to be in that harmonised A/B space, so this step only performs QC validation on the truth alleles before decoding GT indices.

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

- `REF` and `ALT` must already be distinct values in `{A, B}`.
- GT indices are decoded directly to `A/A`, `A/B`, `B/A`, or `B/B`.
- If a truth row fails the A/B allele QC, it is written to the exception report and emitted with missing genotypes rather than being retranslated from nucleotide alleles.

This means Stage 5 assumes the array truth has already been harmonised upstream into the same A/B space used for downstream dosage comparison.

#### Exception Handling

- Missing genotypes (`./. or .|.`) are recorded as `./.` in the output and are **not** treated as exceptions.
- Any imputed GT that cannot be parsed, or whose decoded nucleotide is not in `{A, T, C, G}`, is set to `./.` and recorded in `{prefix}.translation_exceptions.tsv`.
- Any truth row whose alleles are not distinct `{A, B}` values is recorded as a truth QC exception in `{prefix}.translation_exceptions.tsv` and emitted with missing genotypes.

**Checkpoint outputs:**

- `{prefix}.imputed.AB_format.tsv` — tab-separated: CHROM, POS, REF, ALT, ID, then one column per sample with A/A, A/B, or B/B genotypes.
- `{prefix}.truth.AB_format.tsv` — same tabular layout, decoded directly from the already A/B-encoded truth VCF.
- `{prefix}.translation_exceptions.tsv` — imputed translation exceptions plus truth A/B QC exceptions.

---

### 5f. R Metrics and Plots

`modules/evaluate/dosage_r2.R` reads the two A/B TSV files and computes imputation quality metrics.

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

#### MAF-Binned Summary

Variants are binned by MAF into standard bins: `[0, 0.01)`, `[0.01, 0.05)`, `[0.05, 0.1)`, `[0.1, 0.2)`, `[0.2, 0.3)`, `[0.3, 0.5]`. Mean r² and concordance are reported per bin. Rare variants (low MAF) are typically harder to impute and tend to have lower r².

#### Per-Sample Metrics

For each sample *j*, r² is computed across all variants (treating variant positions as observations): `r²ⱼ = cor(imputed_dosage[all variants, j], truth_dosage[all variants, j])²`. This identifies individual samples with systematically poor imputation quality. A per-0.1-MAF-bin breakdown is also provided.

#### Fine-Scale MAF Line Plot

Variants are binned in 0.01 MAF increments (50 bins from 0 to 0.5) and mean r² is plotted as a line, providing a smooth view of how imputation accuracy varies with allele frequency.

#### Per-Chromosome 1 Mb Plot

Variants are binned into 1 Mb windows along each chromosome and mean r² is plotted as a faceted line chart, revealing regional patterns of imputation accuracy (e.g., centromere effects, panel coverage gaps).

#### Output Files

| File | Description |
|---|---|
| `{prefix}.metrics.tsv` | Per-variant r (signed), r², concordance, MAF, n_non_missing |
| `{prefix}.summary.tsv` | Overall mean/median r², mean concordance, total variant counts |
| `{prefix}.maf_bins.tsv` | r (signed), r², and concordance aggregated by standard MAF bins |
| `{prefix}.per_sample_metrics.tsv` | Per-sample r (signed) and r² overall and per 0.1 MAF bin |
| `{prefix}.concordance.parquet` | Per-site per-sample concordance matrix (0/1/NA); Arrow Parquet |
| `{prefix}.r2_hist.png` | Histogram of per-variant r² distribution |
| `{prefix}.concordance_hist.png` | Histogram of per-variant concordance |
| `{prefix}.r2_vs_maf.png` | 2D density (hex or bin2d) of r² vs MAF |
| `{prefix}.r2_vs_maf_line.png` | Mean r² per 0.01 MAF bin (line plot) |
| `{prefix}.r2_vs_maf_line.tsv` | Data underlying the MAF line plot |
| `{prefix}.r2_per_chr_1Mb.png` | Mean r² per 1 Mb window, faceted by chromosome |
| `{prefix}.r2_per_chr_1Mb.tsv` | Data underlying the per-chromosome plot |
| `{prefix}.r2_per_sample.png` | Per-sample r² bar chart (sorted ascending) |
| `{prefix}.r_per_sample.png` | Per-sample signed r bar chart (blue = positive, red = negative) |

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
- `--truth` or `--vcf`: an already harmonised A/B truth VCF/BCF.

It reuses the same overlap, deduplication, and strand-ambiguity filtering model as `dosage_r2.sh`, then:

1. translates only the query-side VCF into A/B genotype space,
2. validates and decodes the truth-side A/B VCF without attempting nucleotide translation,
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
- query translation exceptions and truth A/B QC exceptions.

---

## 8. Configuration

Configuration is split across two files:

### `config/environment.sh`

Paths and environment settings. Copy from `config/environment.template.sh` and fill in:

```bash
# Conda/module setup
export MINIFORGE_MODULE="miniforge/25.3.0-3"
export CONDA_ENV="myenv_py310"
export BCFTOOLS_MODULE="bcftools/1.18-gcc-12.3.0"

# Paths
export QUILT2_HOME="/path/to/QUILT"                # Directory containing QUILT2.R
export QUILT2_GENETIC_MAP="/path/to/genetic_map"   # File or directory of per-chr maps
export QUILT2_REFERENCE_FASTA="/path/to/ref.fa"    # Optional; needed for contig header fixes
```

### `config/quilt2_config.sh`

SLURM resource defaults (can be overridden via environment variables):

| Variable | Default | Description |
|---|---|---|
| `QUILT2_ACCOUNT` | `a_qaafi_cas` | SLURM account |
| `QUILT2_PARTITION` | `general` | SLURM partition |
| `QUILT2_CPUS_PER_TASK` | `2` | CPUs per chunk task |
| `QUILT2_MEMORY` | `8G` | Memory per chunk task |
| `QUILT2_TIME_LIMIT` | `72:00:00` | Wall time per job |
| `QUILT2_ARRAY_MAX` | `0` (no cap) | Max concurrent array tasks |
| `QUILT2_CONSTRAINT` | `epyc4` | Node constraint (avoids ISA issues) |

---

## 9. How to Run

### Prerequisites

1. Copy and fill in `config/environment.sh`.
2. Ensure the conda environment `myenv_py310` (or your chosen `CONDA_ENV`) contains:
   - `r-base`, `r-data.table`, `r-ggplot2`
   - Optional: `r-arrow` (Parquet output), `r-vcfppr`, `r-hexbin`
3. Ensure `bcftools` module is available via Lmod.
4. Ensure any array truth VCF used in Stage 5 has already been harmonised into A/B format before submission.

### Full Run (Panel Prep + Imputation)

```bash
bash bin/run_quilt2.sh \
  --input-dir /data/apple_lowpass \
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
  --genetic-map /data/maps/apple \
  --auto-chunk-map \
  --dry-run
```

---

## 10. Output Files

### Imputation Outputs (`quilt2_output/`)

```
quilt2_output/
├── panel/
│   ├── <chr>.filtered.vcf.gz       # Panel after Stage 1 (if run)
│   └── missing_sites_removed.tsv   # Sites removed by quality filter
├── RData/
│   └── RData.<chr>.<start>.<end>.RData  # QUILT2 reference objects (one per chunk)
├── <chr>/
│   └── <start>_<end>/
│       └── quilt.output.vcf.gz     # Per-chunk imputed VCF
└── tmp/
    └── chunks.tsv                  # Chunk manifest used for array submission
```

Final per-chromosome imputed VCFs are concatenated to the `--input-dir` or a path specified by `--eval-output`.

### Evaluation Outputs (`{prefix}.*`)

See the table in [Stage 5f](#5f-r-metrics-and-plots). All files share a common prefix specified by `--out-prefix`. Intermediate files (overlapped, unambiguous, A/B TSVs) are retained alongside final metrics so that specific steps can be re-run or inspected.

### Standalone Concordance Outputs (`{prefix}.*`)

For `utils/test_concordance_check_with_array_validation.sh`, the main Quarto-facing outputs are:

- `{prefix}.best_matches.tsv`
- `{prefix}.pairwise_concordance.tsv`
- `{prefix}.pipeline_audit.tsv`
- `{prefix}.output_manifest.tsv`

The utility also keeps overlapped/unambiguous VCF checkpoints and intermediate TSVs for download and audit.

### SLURM Logs

```
quilt2_slurm/
├── quilt2_master_<jobid>.output     # Master job stdout
├── quilt2_master_<jobid>.error      # Master job stderr
├── quilt2_array_<jobid>_<task>.output  # Per-chunk stdout
└── quilt2_array_<jobid>_<task>.error   # Per-chunk stderr
```

For evaluation jobs:

```
<out_dir>/slurm/
├── dosage_r2_<jobid>.out
└── dosage_r2_<jobid>.err
```
