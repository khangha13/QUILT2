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

## 13. Master self-submit (detaching run_quilt2)

### Problem
- Running `bin/run_quilt2.sh` directly would block on the login node while Phase 1/2 SLURM jobs ran, unlike the GATK master which self-submits and detaches.

### Solution
- `bin/run_quilt2.sh` now defaults to self-submit when not already in SLURM:
  - Adds `--submit-self[=bool]` (default true) and `--no-submit`/`--submit-self=false` to opt out.
  - Creates `quilt2_slurm/quilt2_master_*.sh`, submits via `sbatch`, records job ID in `quilt2_slurm/quilt2_master_job_id.txt`, prints job ID + SLURM script/log paths, then exits.
  - When inside SLURM or when self-submit is disabled, runs inline as before.

### Notes
- DRY-RUN prints the sbatch command and exits without submitting.
- SLURM logs: `quilt2_slurm/quilt2_master_%j.(output|error)`.

---

## 14. Phase 1 failures: env not loaded + standardise-only validation

### Symptoms
- Master: `Cannot locate quilt2_prepare (QUILT2_prepare_reference.R). Set --quilt2-prepare-script or --quilt2-home.`
- Phase 1 task (standardise-only): `Filtered VCF missing or empty for ChrXX: <empty>` even though `--remove-missing` was **not** set.

### Causes
- The master resolved `QUILT2_prepare_reference.R` before the conda env was activated, so scripts provided by the env were not found.
- Phase 1 validation always expected a “filtered” VCF, even when only `--standardise-name` was requested, so `cleaned_vcf` was empty and validation failed.

### Fixes
- Load the conda env at the start of `bin/run_quilt2.sh` so `QUILT2_prepare_reference.R`/`QUILT2.R` are on PATH before resolution (or set `--quilt2-home`/`--quilt2-prepare-script` explicitly).
- In `templates/quilt2_nomiss_job.sh`:
  - Treat the standardised (or original) VCF as the Phase 1 output when `--remove-missing` is **false**; only expect `quilt.nomiss.*` when the flag is **true**.
  - Sort and re-index standardised VCFs before validation; if an existing `_chr` VCF cannot be indexed, rebuild it.

### Notes for future runs
- Ensure the quilt2 conda env is activated (or `--quilt2-home`/`--quilt2-*-script` are set) so the prepare/impute scripts are locatable.
- For standardise-only runs, Phase 1 success is defined by the presence of `<chr>_chr.vcf.gz` (+ index), not `quilt.nomiss.*`.

---

## 15. `--standardise-name` fails: "Contig 'X' is not defined in the header"

### Symptom
When running with `--standardise-name`, Phase 1 fails with:
```
[W::vcf_parse] Contig '1' is not defined in the header. (Quick workaround: index the file with tabix.)
Encountered an error, cannot proceed. Please check the error output above.
[E::bcf_hdr_read] Input is not detected as bcf or vcf format
Could not read VCF/BCF headers from -
[ERROR] Failed to standardise+sort Chr01
```

### Cause
The input VCF has variant records with contig names (e.g., `1`) that are not defined in the VCF header (`##contig=<ID=1,...>`). `bcftools annotate --rename-chrs` fails when it encounters undefined contigs, producing invalid output that breaks the downstream `bcftools sort`.

### Fix
The pipeline now automatically creates a **tabix index** on the source VCF before running `bcftools annotate --rename-chrs`. This allows bcftools to handle VCFs with undefined contig headers.

### How it works
1. Before renaming contigs, the pipeline runs `tabix -f -p vcf <source.vcf.gz>` to create a `.tbi` index.
2. With a tabix index present, bcftools can process the VCF even if contig headers are missing.
3. If the VCF already has `Chr*`-prefixed contigs, the rename step is skipped entirely.

### Notes
- **Smart skip:** If the VCF already has `Chr*`-prefixed contigs (e.g., `Chr01`), the entire standardisation is skipped — the file is just copied and indexed.
- **No reference FASTA required:** The tabix workaround doesn't require a reference FASTA, unlike `bcftools reheader`.
- The tabix index is created on the source VCF (in place), not a temporary copy.

---

## 16. `concat_imputed.sh`: "Unsorted positions" when indexing concatenated per-chromosome VCF

### Problem
```
Concatenating .../quilt2.diploid.Chr01.4990000-8990000.vcf.gz	0.007837 seconds
Concatenating .../quilt2.diploid.Chr01.8980000-12980000.vcf.gz	0.006933 seconds
[E::hts_idx_push] Unsorted positions on sequence #1: 8989860 followed by 8980005
index: failed to create index for ".../imputed.Chr01.vcf.gz"
```
`bcftools concat --naive` succeeds (headers compatible, all files "concatenated"), but `bcftools index` on the result fails.

### Cause
Adjacent chunks are not just listed in the right order — their *content* overlaps. Chunk boundaries from `QUILT::quilt_chunk_map()` (`--auto-chunk-map`) share a small buffer between neighbors (e.g. chunk 1 = `1-5000000`, chunk 2 = `4990000-8990000`: a 10,000bp overlap). `templates/quilt2_job.sh` passes these directly as QUILT2.R's `--regionStart`/`--regionEnd`, and QUILT2.R writes out every variant across that full region — so both chunks independently call genotypes for the shared window. Naively concatenating full chunk files then places the tail of chunk *N* (higher positions, up to its end) immediately before the head of chunk *N+1* (lower positions, from its start), which is locally out-of-order — `bcftools index` correctly rejects this.

This is **not** the same issue as #3 (scientific notation) or the earlier "chunk ordering" fix in `concat_imputed.sh` (manifest-based ordering). Ordering was already correct; the problem is overlapping *content* between correctly-ordered neighbors.

### Fix
`concat_imputed.sh` now trims each chunk to end just before the next chunk's start before concatenating, using the manifest's (or filename's) own start/end values — no hardcoded overlap size is assumed:
```bash
# For chunk i with (start_i, end_i) and next chunk's start_{i+1}:
trimmed_end_i = min(end_i, start_{i+1} - 1)
bcftools view -t "chr:start_i-trimmed_end_i" -Oz -o trimmed.vcf.gz chunk_i.vcf.gz
```
- Uses `bcftools view -t` (targets, streaming) rather than `-r` (regions), so no index is required on the raw per-chunk files.
- If chunks don't actually overlap, `trimmed_end_i == end_i` and the original file is used unchanged (no-op, no extra `bcftools view` call).
- The last chunk of a chromosome is never trimmed (no next chunk to overlap with).

### Notes
- This only affects concatenation of **overlapping-chunk** pipelines (e.g. `--auto-chunk-map`). Pipelines using non-overlapping/contiguous `--chunk-file` regions are unaffected (trim is a no-op).
- If you ever see `bcftools index` complain about unsorted positions again after concatenation, check whether the chunk manifest's regions overlap before assuming it's a different bug.

---

## Summary Checklist

Before running QUILT2 pipeline, verify:

1. ✅ **Chromosome names match** across: VCF, BAM files, genetic maps, and `--chr` argument
2. ✅ **Genetic map format** uses correct column names: `position COMBINED_rate.cM.Mb. Genetic_Map.cM.`
3. ✅ **Reference panel is phased** (genotypes use `|` not `/`)
4. ✅ **VCF files are indexed** (`.csi` or `.tbi` index files exist)
5. ✅ **Delete cached files** when changing configuration
6. ✅ **`--standardise-name`** auto-creates tabix index to handle VCFs with undefined contig headers
7. ✅ **Overlapping chunks** (e.g. `--auto-chunk-map`) are safe to concatenate via `modules/evaluate/concat_imputed.sh`, which trims each chunk's overlap before concatenating
