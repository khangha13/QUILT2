# QUILT2 Pipeline: Past Problems and Solutions

This document records issues encountered while developing and testing `quilt2_pipeline.sh` and their resolutions.

---

## 1. Genetic Map Column Names

### Problem
```
Error in `[.data.frame`(genetic_map, , "Genetic_Map.cM.") : 
  undefined columns selected
```
QUILT expects specific column names in the genetic map file.

### Cause
The genetic map files used simple column names (`position rate cM`) instead of QUILT's expected format.

### Solution
Use the exact column names QUILT expects:
```
position COMBINED_rate.cM.Mb. Genetic_Map.cM.
1000 1.000000 0.000000
2000 1.000000 0.001000
```

**Script to regenerate dummy maps:**
```bash
for chr in Chr01 Chr02 ... Chr17; do
    len=$(awk -v c="$chr" '$1==c{print $2}' reference.fasta.fai)
    awk -v len="$len" '
        BEGIN {
            print "position COMBINED_rate.cM.Mb. Genetic_Map.cM."
            for (p = 1000; p <= len; p += 1000) {
                cM = (p - 1000) / 1000000.0
                printf "%d 1.000000 %.6f\n", p, cM
            }
        }
    ' > dummy_map/${chr}.txt
done
```

---

## 2. Log Messages Contaminating VCF Path

### Problem
```
[E::hts_open_format] Failed to open file "[INFO] Using Step1C panel VCF for Chr01: /path/to/vcf.gz
/path/to/vcf.gz"
```
The VCF file path variable contained log messages.

### Cause
`log_info()` and `log_warn()` functions were printing to **stdout** instead of **stderr**. When functions were called in command substitutions like `vcf=$(normalize_panel_vcf ...)`, log messages were captured along with the return value.

### Solution
Changed logging functions to write to stderr:
```bash
# Before (BUG):
log_info() { echo "[INFO] $*"; }
log_warn() { echo "[WARN] $*"; }

# After (FIXED):
log_info() { echo "[INFO] $*" >&2; }
log_warn() { echo "[WARN] $*" >&2; }
```

---

## 3. Scientific Notation in Chunk Coordinates

### Problem
```
[ERROR] Chunk 0 has non-numeric coordinates: start=1 end=5e+06
```
QUILT's `quilt_chunk_map` outputs coordinates in scientific notation.

### Cause
R outputs large numbers in scientific notation (e.g., `5e+06`), but the bash validation regex `^[0-9]+$` only accepts pure integers.

### Solution
Modified the R script to use `sprintf` with `%d` format to force integer output:
```r
start_int <- as.integer(as.numeric(reg_parts[2]))
end_int <- as.integer(as.numeric(reg_parts[3]))
lines <- c(lines, sprintf("%s\t%s\t%d\t%d", chunk, chr, start_int, end_int))
```

---

## 4. Symbolic Links Not Found by `find`

### Problem
```
[ERROR] Reference panel VCF for 1 not found
```
VCF symlinks were not being discovered.

### Cause
The `find` command used `-type f` which only finds regular files, not symbolic links.

### Solution
Changed to find both files and symlinks:
```bash
# Before:
find "${REFERENCE_PANEL_DIR}" -maxdepth 1 -type f -name "${chr}*.vcf.gz"

# After:
find "${REFERENCE_PANEL_DIR}" -maxdepth 1 \( -type f -o -type l \) -name "${chr}*.vcf.gz"
```

---

## 5. Chromosome Name Matching Too Broad

### Problem
```
[INFO] Using Step1C panel VCF for 1: .../10_phased.vcf.gz
```
When searching for chromosome `1`, it incorrectly matched `10_phased.vcf.gz`.

### Cause
The glob pattern `${chr}*.vcf.gz` matched any file starting with the chromosome name, so `1*.vcf.gz` matched both `1_phased.vcf.gz` and `10_phased.vcf.gz`. Alphabetical sorting put `10` before `1_`.

### Solution
Changed glob patterns to require a separator character after the chromosome name:
```bash
# Before:
-name "${chr}*.vcf.gz"

# After:
\( -name "${chr}_*.vcf.gz" -o -name "${chr}.*.vcf.gz" -o -name "${chr}-*.vcf.gz" \)
```

---

## 6. Chunk File Parsing with Numeric Chromosomes

### Problem
```
[INFO] ---- Chunk 0_1_1 (0:1-1, buffer=5000000) ----
```
Chunks were being parsed incorrectly when chromosome names were numeric.

### Cause
The chunk file format detection checked if column 2 started with `chr*`. With numeric chromosome names like `1`, `2`, etc., this check failed and the parser fell back to wrong column assignments.

### Solution
Added heuristic to detect 4-column format by checking if column 4 looks like a coordinate:
```bash
elif [[ "${c4}" =~ ^[0-9]+$ && "${c4}" -gt 1000 ]]; then
    # Format: chunk chr start end [buffer]
    add_chunk "${c1}" "${c2}" "${c3}" "${c4}" "${c5:-${BUFFER}}"
```

---

## 7. VCF and BAM Chromosome Name Mismatch

### Problem
```
Error in quilt_get_chromosome_length(...) : 
  Could not find chromosome length for file:.../BAM/sample.bam
```

### Cause
- BAM files were aligned to reference using `Chr01`, `Chr02`, etc.
- VCF files (from NCBI) used `1`, `2`, etc.
- QUILT couldn't find matching chromosomes.

### Solution
Rename chromosomes in VCF files to match BAM files:
```bash
# Create rename file
cat > chr_rename.txt << 'EOF'
1 Chr01
2 Chr02
...
17 Chr17
EOF

# Rename using sed (handles header too)
for i in 01 02 ... 17; do
    num=$((10#$i))
    zcat Chr${i}.vcf.gz_phased.vcf.gz | \
        sed "s/^${num}\t/Chr${i}\t/; s/ID=${num}/ID=Chr${i}/g" | \
        bgzip > Chr${i}_phased_renamed.vcf.gz
    tabix -p vcf Chr${i}_phased_renamed.vcf.gz
done
```

---

## 8. Unphased VCF Reference Panel

### Problem
```
Error in QUILT_prepare_reference(...) : 
  There are no variants in the region you are imputing.
```

### Cause
The VCF used as reference panel contained mostly **unphased** genotypes (`0/1`, `1/1` with `/`) and many missing genotypes (`./.|.`). QUILT requires **phased** haplotypes (using `|`).

### Solution
1. Use properly phased VCF files (e.g., from Beagle phasing)
2. Or use `--remove-missing` with `--min-phased-rate` to filter:
```bash
bash quilt2.sh ... --remove-missing --min-phased-rate 0.5
```

The `--min-phased-rate` filter keeps only variants where ≥X% of samples have phased genotypes.

---

## 9. Cached Files Not Regenerated

### Problem
Pipeline reuses old cached files even after configuration changes.

### Cause
The pipeline caches:
- Auto-generated chunk files: `quilt2_output/tmp/quilt_auto_chunks.tsv`
- Filtered panel VCFs: `quilt2_output/panel/quilt.nomiss.*.vcf.gz`
- Prepared references: `quilt2_output/RData/QUILT_prepared_reference.*.RData`

### Solution
Delete cached files before re-running with new settings:
```bash
rm -f <work_dir>/quilt2_output/tmp/quilt_auto_chunks.tsv
rm -f <work_dir>/quilt2_output/panel/quilt.nomiss.*.vcf.gz*
rm -f <work_dir>/quilt2_output/RData/QUILT_prepared_reference.*.RData
```

---

## 10. Genetic Map Files for Different Chromosome Naming

### Problem
Genetic map files named `Chr01.txt` but pipeline searching for `1.txt` (or vice versa).

### Cause
Genetic map files must match the chromosome names used in `--chr` argument.

### Solution
Create copies of genetic map files with matching names:
```bash
# If using numeric chromosomes (1, 2, ...)
for i in 01 02 ... 17; do
    num=$((10#$i))
    cp Chr${i}.txt ${num}.txt
done

# Or vice versa if using Chr format
```

---

## 11. `--remove-missing` now runs in two SLURM phases

### Problem
- The `--remove-missing` filter ran sequentially inside `run_quilt2.sh`, so per-chromosome filtering could not use SLURM and the master could launch Phase 2 without knowing Phase 1 status.
- Failures during filtering were silent (no job ID, no failure flag), making diagnosis hard when panel inputs were bad.

### Solution
- Added a dedicated Phase 1 SLURM array (`templates/quilt2_nomiss_job.sh`) that filters the panel per chromosome and writes cleaned VCFs to `quilt2_output/panel/` (`quilt.nomiss.<chr>.vcf.gz` + index).
- The master (`bin/run_quilt2.sh`) now:
  - Submits Phase 1 when `--remove-missing` is set, records the job ID in `quilt2_slurm/quilt2_nomiss_job_id.txt`, and waits on the job (`squeue` poll).
  - Uses a failure flag `quilt2_slurm/quilt2_nomiss_failed.flag` that Phase 1 workers touch on any error (trap on `ERR`).
  - Verifies every cleaned VCF and index exist before submitting Phase 2 (chunk array).
- Phase 2 still passes `--remove-missing`; workers detect existing `quilt.nomiss.*` files and skip re-filtering.

### Notes / How to debug
- Phase 1 SLURM logs: `quilt2_slurm/quilt2_nomiss_%A_%a.output|error`.
- If Phase 1 fails: check the failure flag, then inspect per-task logs and the panel directory for the missing chromosome.
- Job IDs:
  - Phase 1: `quilt2_slurm/quilt2_nomiss_job_id.txt`
  - Phase 2: `quilt2_slurm/quilt2_job_id.txt`

---

## 12. `--standardise-name` now runs in Phase 1 (with `--remove-missing`)

### Problem
- Contig standardisation ran sequentially from the master and could not be combined with `--remove-missing`; per-chrom jobs were not SLURM-parallelised, and the master could not validate per-chrom outputs before Phase 2.

### Solution
- Phase 1 array now handles both operations per chromosome when requested:
  - If `--standardise-name`: create `quilt2_output/panel/<chr>_chr.vcf.gz` (+ index), skipping unless `--standardise-name-force`.
  - If `--remove-missing`: create `quilt2_output/panel/quilt.nomiss.<chr>.vcf.gz` (+ index).
  - When both are set: each task runs standardise → remove-missing on the same chromosome.
- The master runs Phase 1 when either flag is set, waits on the job ID, and validates outputs (standardised and/or filtered) before submitting Phase 2. Phase 2 consumes `quilt2_output/panel/` to avoid re-running standardise/remove-missing.

### Notes / How to debug
- Phase 1 logs: `quilt2_slurm/quilt2_nomiss_%A_%a.output|error`; failure flag: `quilt2_slurm/quilt2_nomiss_failed.flag`.
- Outputs checked before Phase 2:
  - Standardised: `quilt2_output/panel/<chr>_chr.vcf.gz` (+ .csi/.tbi) when `--standardise-name`.
  - Filtered: `quilt2_output/panel/quilt.nomiss.<chr>.vcf.gz` (+ .csi/.tbi) when `--remove-missing`.
  - Phase 2 panel input directory: `quilt2_output/panel/` when Phase 1 ran; original directory otherwise.

---

## Summary Checklist

Before running QUILT2 pipeline, verify:

1. ✅ **Chromosome names match** across: VCF, BAM files, genetic maps, and `--chr` argument
2. ✅ **Genetic map format** uses correct column names: `position COMBINED_rate.cM.Mb. Genetic_Map.cM.`
3. ✅ **Reference panel is phased** (genotypes use `|` not `/`)
4. ✅ **VCF files are indexed** (`.csi` or `.tbi` index files exist)
5. ✅ **Delete cached files** when changing configuration
