# Creating Dummy Genetic Maps for QUILT2

When a real genetic/linkage map is not available for your species, you can create
a dummy map with a constant recombination rate. This allows QUILT2 to run, though
imputation accuracy may be reduced compared to using a species-specific map.

## Automatic Generation (Recommended)

`bin/run_quilt2.sh` can generate the dummy maps for you. Pass `dummy` as the
`--genetic-map` value and provide `--reference-fasta` with an existing `.fai`
index:

```bash
bash bin/run_quilt2.sh \
    -i /path/to/work_dir \
    --reference-panel-dir /path/to/panel \
    --genetic-map dummy \
    --reference-fasta /path/to/reference.fasta \
    --auto-chunk-map \
    --bamlist bamlist.txt
```

This creates one file per chromosome in `CHR_LIST` at
`OUTPUT_DIR/genetic_map/dummy/<chr>.txt`, using the same constant 1.0 cM/Mb,
1000bp-step format described in Method 1 below (chromosome lengths come from
the FASTA's `.fai` index). Files are only generated once; re-running the
pipeline reuses existing files. If `--reference-fasta`/`.fai` is not
available, or you want a different rate/step or the variant-anchored Method 2
below, create maps manually with the scripts in this document and pass the
resulting directory to `--genetic-map` instead.

## QUILT2 Genetic Map Format

QUILT2 expects a **3-column space-delimited file**.

Based on testing and the QUILT documentation, the function `QUILT::quilt_chunk_map` expects:

| Column | Name | Description |
|--------|------|-------------|
| 1 | `position` | Physical position in base pairs (1-based) |
| 2 | `COMBINED_rate.cM.Mb.` | Recombination rate in cM/Mb |
| 3 | `Genetic_Map.cM.` | Cumulative genetic distance in cM |

### Required format

QUILT's internal functions (`QUILT::quilt_chunk_map`, `QUILT_prepare_reference`) expect **exact column names** with periods (not parentheses):

```
position COMBINED_rate.cM.Mb. Genetic_Map.cM.
1000 1.000000 0.000000
2000 1.000000 0.001000
3000 1.000000 0.002000
4000 1.000000 0.003000
```

> ⚠️ **Important**: Using simplified column names like `position rate cM` will cause the error:
> `Error in [.data.frame(genetic_map, , "Genetic_Map.cM.") : undefined columns selected`

**Important notes:**
- Header line is **required**
- Columns are **space-separated** (not tab)
- The `cM` column starts at 0 for the first position
- See: [QUILT paper (Davies et al., 2021)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7611184/)

## Creating Dummy Maps (Constant Rate)

A common assumption is **1 cM per 1 Mb** (rate = 1.0 cM/Mb). This means:
- `cM = (position - first_position) / 1,000,000`
- `rate = 1.0` (constant)

### Method 1: Per-chromosome maps from reference FASTA index

This creates one map file per chromosome, using chromosome lengths from the `.fai` index.

```bash
#!/bin/bash
# Generate dummy genetic maps with constant 1 cM/Mb rate
# Uses QUILT's required column names: position, COMBINED_rate.cM.Mb., Genetic_Map.cM.

rate_cM_per_Mb=1.0
ref="/QRISdata/Q8367/Reference_Genome/GDDH13_1-1_formatted.fasta"
outdir="dummy_map"
step=1000  # Position step in bp (1kb intervals)

mkdir -p "$outdir"

for chr in Chr01 Chr02 Chr03 Chr04 Chr05 Chr06 Chr07 Chr08 Chr09 Chr10 Chr11 Chr12 Chr13 Chr14 Chr15 Chr16 Chr17; do
    len=$(awk -v c="$chr" '$1==c{print $2}' "${ref}.fai")
    if [ -z "$len" ]; then
        echo "ERROR: $chr not found in ${ref}.fai"
        exit 1
    fi
    
    echo "Generating map for $chr (length: $len bp)"
    
    awk -v len="$len" -v rate="$rate_cM_per_Mb" -v step="$step" '
        BEGIN {
            # QUILT required column names (space-separated)
            print "position COMBINED_rate.cM.Mb. Genetic_Map.cM."
            first_pos = step
            for (p = step; p <= len; p += step) {
                # cM starts at 0 for first position
                cM = (p - first_pos) / 1000000.0
                printf "%d %.6f %.6f\n", p, rate, cM
            }
        }
    ' > "${outdir}/${chr}.txt"
done

echo "Done. Maps saved to ${outdir}/"
```

### Method 2: Maps anchored to variant positions from VCF

If you want map entries only at variant sites (slightly cleaner):

```bash
#!/bin/bash
# Generate dummy map from VCF variant positions
# Uses QUILT's required column names: position, COMBINED_rate.cM.Mb., Genetic_Map.cM.

rate_cM_per_Mb=1.0
vcf="/path/to/panel.vcf.gz"
outdir="dummy_map"

mkdir -p "$outdir"

# Extract chromosomes from VCF
chroms=$(zcat "$vcf" | grep -v "^#" | cut -f1 | sort -u)

for chr in $chroms; do
    echo "Generating map for $chr"
    
    zcat "$vcf" | awk -v chr="$chr" -v rate="$rate_cM_per_Mb" '
        BEGIN {
            print "position COMBINED_rate.cM.Mb. Genetic_Map.cM."
            first_pos = 0
        }
        !/^#/ && $1 == chr {
            if (first_pos == 0) first_pos = $2
            cM = ($2 - first_pos) / 1000000.0
            printf "%d %.6f %.6f\n", $2, rate, cM
        }
    ' > "${outdir}/${chr}.txt"
done

echo "Done. Maps saved to ${outdir}/"
```

## File Naming Conventions

The QUILT2 pipeline (`bin/run_quilt2.sh`) searches for per-chromosome map files
with these naming patterns (in order):

1. `Chr01.map`
2. `Chr01.txt`
3. `Chr01.txt.gz`
4. `genetic_map.Chr01.txt`
5. `genetic_map.Chr01.txt.gz`
6. `genetic_map_Chr01.txt`
7. `genetic_map_Chr01.txt.gz`
8. `Chr01_genetic_map.txt`
9. `Chr01_genetic_map.txt.gz`
10. Any file matching `*Chr01*.txt`, `*Chr01*.txt.gz`, or `*Chr01*.map`

**Recommended:** Use `{chr}.txt` format (e.g., `Chr01.txt`, `Chr02.txt`, ...).

## Chromosome Naming Consistency

> ⚠️ **Critical**: Chromosome names must match across all files:
> - Genetic map filenames (`Chr01.txt` or `1.txt`)
> - VCF chromosome IDs (in header and data)
> - BAM file chromosome names (`@SQ SN:Chr01` or `@SQ SN:1`)
> - The `--chr` argument passed to the pipeline

If your VCF uses numeric names (`1`, `2`, ...) but BAMs use `Chr01`, `Chr02`, ..., `bin/run_quilt2.sh` auto-detects this by peeking at the first contig of each chromosome's panel VCF and automatically renames them (equivalent to `--standardise-name`) before imputation — see `README.md`. Pass `--no-standardise-name` if you'd rather manage naming yourself. To do it manually instead:

1. **Rename VCF chromosomes** to match BAMs:
   ```bash
   for i in 01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17; do
       num=$((10#$i))
       zcat ${num}_phased.vcf.gz | \
           sed "s/^${num}\t/Chr${i}\t/; s/ID=${num}/ID=Chr${i}/g" | \
           bgzip > Chr${i}_phased.vcf.gz
       tabix -p vcf Chr${i}_phased.vcf.gz
   done
   ```

2. **Create genetic maps for both naming schemes**:
   ```bash
   # If you have Chr01.txt, create 1.txt as a copy
   for i in 01 02 ... 17; do
       num=$((10#$i))
       cp Chr${i}.txt ${num}.txt
   done
   ```

## Usage with quilt2_pipeline.sh

Pass the directory containing per-chromosome maps:

```bash
bash bin/run_quilt2.sh \
    -i /path/to/work_dir \
    --genetic-map /path/to/dummy_map/ \
    --auto-chunk-map \
    --bamlist bamlist.txt
```

Or pass a single concatenated map file (less common):

```bash
bash bin/run_quilt2.sh \
    -i /path/to/work_dir \
    --genetic-map /path/to/single_map.txt \
    --chunk-file chunks.tsv \
    --bamlist bamlist.txt
```

## Verifying Map Format

Check that your map file has the correct format:

```bash
# View first 5 lines
head -5 dummy_map/Chr01.txt

# Should show something like:
# position COMBINED_rate.cM.Mb. Genetic_Map.cM.
# 1000 1.000000 0.000000
# 2000 1.000000 0.001000
# 3000 1.000000 0.002000
# 4000 1.000000 0.003000

# Check column count (should be 3)
awk '{print NF}' dummy_map/Chr01.txt | sort -u
# Should output: 3

# Check for correct header (QUILT required column names)
head -1 dummy_map/Chr01.txt
# Should be: position COMBINED_rate.cM.Mb. Genetic_Map.cM.

# Verify it's space-separated (not tab)
head -2 dummy_map/Chr01.txt | cat -A
# Should show spaces, not ^I (tab characters)
```

## Limitations of Dummy Maps

Using a constant-rate dummy map has limitations:

1. **Reduced accuracy**: Real recombination rates vary across the genome
   (hotspots vs coldspots). A constant rate ignores this variation.

2. **Suboptimal chunk boundaries**: `--auto-chunk-map` uses the genetic map to
   place chunk boundaries at recombination hotspots. With a constant rate, chunks
   are placed uniformly, which may not be optimal.

3. **Phase accuracy**: Haplotype phasing benefits from accurate recombination
   rate information.

**Recommendation**: If a species-specific genetic map becomes available, use it
instead of the dummy map for better imputation accuracy.

## References

- [QUILT2 Tutorial](https://github.com/rwdavies/QUILT/blob/master/README_QUILT2.org)
- [QUILT GitHub Repository](https://github.com/rwdavies/QUILT)
