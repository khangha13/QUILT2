#!/bin/bash
# End-to-end synthetic test for exact-isec WGS GT-to-GT evaluation.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WRAPPER="${ROOT_DIR}/bin/dosage_r2_sbatch.sh"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/test_dosage_r2_wgs.XXXXXX")"
trap 'rm -rf "${WORK_DIR}"' EXIT

for cmd in bcftools bgzip Rscript; do
    command -v "${cmd}" >/dev/null 2>&1 || { echo "SKIP: ${cmd} is unavailable"; exit 0; }
done
Rscript -e 'quit(status=if (requireNamespace("data.table", quietly=TRUE) && requireNamespace("arrow", quietly=TRUE)) 0 else 1)' \
    || { echo "SKIP: R packages data.table and arrow are unavailable"; exit 0; }

fail() { echo "FAIL: $*" >&2; exit 1; }
assert_file() { [[ -s "$1" ]] || fail "missing or empty file: $1"; }
assert_tsv() {
    local file="$1" expression="$2" message="$3"
    awk -F '\t' "${expression}" "${file}" || fail "${message}"
}
assert_parquet() {
    local dataset="$1"
    [[ -n "$(find "${dataset}" -type f -name 'part-*.parquet' -print -quit 2>/dev/null)" ]] \
        || fail "missing Parquet parts: ${dataset}"
}
assert_parquet_r() {
    local dataset="$1" expression="$2" message="$3"
    Rscript - "${dataset}" "${expression}" <<'RS' || fail "${message}"
args <- commandArgs(trailingOnly = TRUE)
files <- list.files(args[[1]], pattern = "[.]parquet$", recursive = TRUE, full.names = TRUE)
if (!length(files)) quit(status = 1)
x <- data.table::rbindlist(lapply(files, function(path) {
  part <- data.table::as.data.table(arrow::read_parquet(path))
  part[, CHROM := sub(".*CHROM=([^/]+).*", "\\1", path)]
  part
}), fill = TRUE)
if (!isTRUE(eval(parse(text = args[[2]])))) quit(status = 1)
RS
}

mkdir -p "${WORK_DIR}/truth/7.Consolidated_VCF"
{
    printf '>Chr01\n'
    awk 'BEGIN { for (i=0; i<1000; i++) printf "A"; printf "\n" }'
    printf '>Chr02\n'
    awk 'BEGIN { for (i=0; i<1000; i++) printf "A"; printf "\n" }'
} > "${WORK_DIR}/reference.fa"
if command -v samtools >/dev/null 2>&1; then
    samtools faidx "${WORK_DIR}/reference.fa"
else
    {
        printf 'Chr01\t1000\t7\t1000\t1001\n'
        printf 'Chr02\t1000\t1015\t1000\t1001\n'
    } > "${WORK_DIR}/reference.fa.fai"
fi

cat > "${WORK_DIR}/imputed.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=1000>
##contig=<ID=Chr02,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
Chr01	100	v100	A	C	.	PASS	.	GT	0|0	0|1	1|1	0|0	0|1	1|1
Chr01	110	v110	A	C	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	120	v120	A	C	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	130	v130	A	C	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	140	v140	A	C	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	150	v150	A	C	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	160	v160	A	C	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	170	v170i	A	G	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	175	v175_shared	A	C	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	175	v175_extra	A	G	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	180	v180	A	C	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	200	v200	A	C,G	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr01	210	v210	A	C	.	PASS	.	GT	./.	0	1/1	0/0	0/1	1/1
Chr01	220	v220	A	C	.	PASS	.	GT	0/1	0/1	0/1	0/1	0/1	0/1
Chr01	230	v230	A	C	.	PASS	.	GT	0/0	1/1	0/1	0/0	0/1	1/1
Chr01	240	v240	A	C	.	PASS	.	GT	0/0	0/1	1/1	0/0	0/1	1/1
Chr02	100	v2_100	A	C	.	PASS	.	GT	0|0	0|1	1|1	0|0	0|1	1|1
Chr02	110	v2_110	A	C	.	PASS	.	GT	1/1	0/1	0/0	1/1	0/1	0/0
Chr02	120	v2_120	A	C	.	PASS	.	GT	0/1	0/0	0/1	0/1	1/1	1/1
EOF

awk '{ gsub(/\\t/, "\t"); print }' "${WORK_DIR}/imputed.vcf" > "${WORK_DIR}/imputed.tabs.vcf"
mv "${WORK_DIR}/imputed.tabs.vcf" "${WORK_DIR}/imputed.vcf"

cat > "${WORK_DIR}/truth/Chr01_consolidated.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=1000>
##INFO=<ID=QD,Number=1,Type=Float,Description="Quality by depth">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Strand odds ratio">
##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="MQ rank sum">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Read position rank sum">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
Chr01	100	t100	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	110	t110	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:59:20	0/1:60:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	120	t120	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:9	0/1:99:10	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	130	t130	A	C	29	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	140	t140	A	C	30	PASS	QD=2;SOR=3;FS=60;MQ=40;MQRankSum=-12.5;ReadPosRankSum=-8	GT:GQ:DP	0/0:60:10	0/1:60:10	1/1:60:10	0/0:60:10	0/1:60:10	1/1:60:10
Chr01	150	t150	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	160	t160	A	C	100	PASS	SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	170	t170	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	175	t175_shared	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	190	t190	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	200	t200	A	C,G	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	210	t210	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	220	t220	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/1:99:20	0/1:99:20	0/1:99:20	0/1:99:20	0/1:99:20	0/1:99:20
Chr01	230	t230	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	1/1:99:20	0/1:59:20	0/0:59:20	0/1:59:20	1/1:59:20
Chr01	240	t240	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0:99:20	./.:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
EOF

cat > "${WORK_DIR}/truth/Chr02_consolidated.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=Chr02,length=1000>
##INFO=<ID=QD,Number=1,Type=Float,Description="Quality by depth">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Strand odds ratio">
##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand">
##INFO=<ID=MQ,Number=1,Type=Float,Description="Mapping quality">
##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="MQ rank sum">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Read position rank sum">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
Chr02	100	t2_100	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr02	110	t2_110	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	1/1:99:20	0/1:99:20	0/0:99:20	1/1:99:20	0/1:99:20	0/0:99:20
Chr02	120	t2_120	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/0:99:20	0/1:99:20	0/1:99:20	1/1:99:20	1/1:99:20
EOF

awk '{ gsub(/\\t/, "\t"); print }' "${WORK_DIR}/truth/Chr01_consolidated.vcf" > "${WORK_DIR}/truth/Chr01_consolidated.tabs.vcf"
mv "${WORK_DIR}/truth/Chr01_consolidated.tabs.vcf" "${WORK_DIR}/truth/Chr01_consolidated.vcf"

awk '{ gsub(/\\t/, "\t"); print }' "${WORK_DIR}/truth/Chr02_consolidated.vcf" > "${WORK_DIR}/truth/Chr02_consolidated.tabs.vcf"
mv "${WORK_DIR}/truth/Chr02_consolidated.tabs.vcf" "${WORK_DIR}/truth/Chr02_consolidated.vcf"

bgzip -c "${WORK_DIR}/imputed.vcf" > "${WORK_DIR}/imputed.vcf.gz"
bcftools index -t "${WORK_DIR}/imputed.vcf.gz"
for chromosome in Chr01 Chr02; do
    mkdir -p "${WORK_DIR}/chunks/imputed/${chromosome}"
    chunk_vcf="${WORK_DIR}/chunks/imputed/${chromosome}/quilt2.diploid.${chromosome}.1-1000.vcf.gz"
    bcftools view -r "${chromosome}" -Oz -o "${chunk_vcf}" "${WORK_DIR}/imputed.vcf.gz"
    bcftools index -t "${chunk_vcf}"
done
cat > "${WORK_DIR}/imputed.no_gt.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=1000>
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
Chr01	100	no_gt	A	C	.	PASS	.	DP	10	10	10	10	10	10
EOF
awk '{ gsub(/\\t/, "\t"); print }' "${WORK_DIR}/imputed.no_gt.vcf" > "${WORK_DIR}/imputed.no_gt.tabs.vcf"
mv "${WORK_DIR}/imputed.no_gt.tabs.vcf" "${WORK_DIR}/imputed.no_gt.vcf"
bgzip -c "${WORK_DIR}/imputed.no_gt.vcf" > "${WORK_DIR}/imputed.no_gt.vcf.gz"
bcftools index -t "${WORK_DIR}/imputed.no_gt.vcf.gz"
bgzip -c "${WORK_DIR}/truth/Chr01_consolidated.vcf" > "${WORK_DIR}/truth/7.Consolidated_VCF/Chr01_consolidated.vcf.gz"
bcftools index -t "${WORK_DIR}/truth/7.Consolidated_VCF/Chr01_consolidated.vcf.gz"
bgzip -c "${WORK_DIR}/truth/Chr02_consolidated.vcf" > "${WORK_DIR}/truth/7.Consolidated_VCF/Chr02_consolidated.vcf.gz"
bcftools index -t "${WORK_DIR}/truth/7.Consolidated_VCF/Chr02_consolidated.vcf.gz"

OUT_ENABLED="${WORK_DIR}/output/dosage_eval_wgs"
QUILT2_REFERENCE_FASTA="${WORK_DIR}/reference.fa" \
QUILT2_WGS_KEEP_DOSAGE_MATRICES=true \
SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs \
    --imputed "${WORK_DIR}/imputed.vcf.gz" \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --out-prefix "${OUT_ENABLED}" 2> "${WORK_DIR}/obsolete-setting.log"
grep -Fq "QUILT2_WGS_KEEP_DOSAGE_MATRICES is obsolete and will be ignored" "${WORK_DIR}/obsolete-setting.log" \
    || fail "obsolete dosage-matrix setting was not reported"

assert_file "${OUT_ENABLED}/per_sample_metrics.tsv"
assert_file "${OUT_ENABLED}/run_manifest.tsv"
assert_file "${OUT_ENABLED}/qc/filter_summary.tsv"
assert_file "${OUT_ENABLED}/qc/genotype_masking_summary.tsv"
for dataset in \
    metrics/per_variant_metrics metrics/site_filtered_variants metrics/imputed_only_variants \
    metrics/truth_only_variants metrics/allele_mismatches; do
    assert_parquet "${OUT_ENABLED}/${dataset}"
done
[[ ! -e "${OUT_ENABLED}/per_variant_metrics.tsv" ]] || fail "large WGS TSV outputs should not be written"
[[ ! -e "${OUT_ENABLED}/intermediate/imputed_ds" ]] || fail "obsolete imputed DS dataset was written"
[[ ! -e "${OUT_ENABLED}/intermediate/truth_gt_dosage" ]] || fail "obsolete truth dosage dataset was written"

expected_header=$'sample\tr_overall\tr2_overall\tn_variants\tr_maf_[0.0,0.1)\tr2_maf_[0.0,0.1)\tn_maf_[0.0,0.1)\tr_maf_[0.1,0.2)\tr2_maf_[0.1,0.2)\tn_maf_[0.1,0.2)\tr_maf_[0.2,0.3)\tr2_maf_[0.2,0.3)\tn_maf_[0.2,0.3)\tr_maf_[0.3,0.4)\tr2_maf_[0.3,0.4)\tn_maf_[0.3,0.4)\tr_maf_[0.4,0.5]\tr2_maf_[0.4,0.5]\tn_maf_[0.4,0.5]'
[[ "$(head -n 1 "${OUT_ENABLED}/per_sample_metrics.tsv")" == "${expected_header}" ]] \
    || fail "WGS per-sample header does not match array mode"
assert_tsv "${OUT_ENABLED}/per_sample_metrics.tsv" 'NR>1 && $4>=3 {ok=1} END {exit !ok}' "chromosome statistics were not aggregated"

# A chromosome-only reference case checks vectorized variant/sample correlations against cor().
OUT_KNOWN="${WORK_DIR}/output/dosage_eval_wgs_known"
QUILT2_REFERENCE_FASTA="${WORK_DIR}/reference.fa" SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" --chr Chr02 \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --out-prefix "${OUT_KNOWN}"
Rscript - "${OUT_KNOWN}" <<'RS' || fail "vectorized GT correlations differ from cor()"
out <- commandArgs(trailingOnly = TRUE)[[1]]
read_dataset <- function(path) {
  files <- list.files(path, pattern = "[.]parquet$", recursive = TRUE, full.names = TRUE)
  data.table::rbindlist(lapply(files, function(x) {
    part <- data.table::as.data.table(arrow::read_parquet(x))
    part[, CHROM := sub(".*CHROM=([^/]+).*", "\\1", x)]
    part
  }), fill = TRUE)
}
variant <- read_dataset(file.path(out, "metrics", "per_variant_metrics"))
same <- function(a, b) (is.na(a) && is.na(b)) || isTRUE(all.equal(a, b, tolerance = 1e-12))
truth <- rbind(
  c(0, 1, 2, 0, 1, 2),
  c(2, 1, 0, 2, 1, 0),
  c(0, 0, 1, 1, 2, 2)
)
imputed <- rbind(
  c(0, 1, 2, 0, 1, 2),
  c(2, 1, 0, 2, 1, 0),
  c(1, 0, 1, 1, 2, 2)
)
data.table::setorder(variant, POS)
if (!identical(as.integer(variant$POS), c(100L, 110L, 120L))) quit(status = 1)
for (i in seq_len(nrow(variant))) {
  expected_r <- cor(imputed[i, ], truth[i, ])
  expected_concordance <- mean(imputed[i, ] == truth[i, ])
  if (!same(expected_r, variant$r[[i]])) quit(status = 1)
  if (!same(expected_r^2, variant$r2[[i]])) quit(status = 1)
  if (!same(expected_concordance, variant$concordance[[i]])) quit(status = 1)
}
sample_metrics <- data.table::fread(file.path(out, "per_sample_metrics.tsv"), check.names = FALSE)
for (j in seq_len(ncol(truth))) {
  row <- sample_metrics[sample == paste0("S", j)]
  expected_r <- cor(imputed[, j], truth[, j])
  if (!same(expected_r, row$r_overall)) quit(status = 1)
  if (!same(expected_r^2, row$r2_overall)) quit(status = 1)
  if (row$n_variants != 3L) quit(status = 1)
  if (!same(expected_r, row[["r_maf_[0.4,0.5]"]])) quit(status = 1)
  if (!same(expected_r^2, row[["r2_maf_[0.4,0.5]"]])) quit(status = 1)
  if (row[["n_maf_[0.4,0.5]"]] != 3L) quit(status = 1)
}
RS

# Chunk input concatenates and evaluates only the selected chromosome.
OUT_CHUNKS="${WORK_DIR}/output/dosage_eval_wgs_chunks"
printf 'S1\nS2\nS3\n' > "${WORK_DIR}/selected_samples.txt"
QUILT2_REFERENCE_FASTA="${WORK_DIR}/reference.fa" SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --chunks-dir "${WORK_DIR}/chunks/imputed" --chr Chr02 \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --out-prefix "${OUT_CHUNKS}" -- --samples "${WORK_DIR}/selected_samples.txt"
assert_file "${OUT_CHUNKS}/per_sample_metrics.tsv"
[[ "$(awk 'END {print NR-1}' "${OUT_CHUNKS}/per_sample_metrics.tsv")" == "3" ]] \
    || fail "sample restriction was not applied"
grep -Fq $'chunks_dir\t'"${WORK_DIR}/chunks/imputed" "${OUT_CHUNKS}/run_manifest.tsv" \
    || fail "chunk input was not recorded"
assert_parquet_r "${OUT_CHUNKS}/metrics/per_variant_metrics" \
    'all(x$CHROM=="Chr02") && nrow(x)==3' "chunk mode evaluated the wrong chromosome"

assert_parquet_r "${OUT_ENABLED}/metrics/per_variant_metrics" \
    'identical(names(x), c("POS","REF","ALT","ID","ID_truth","n_pairs","truth_alt_frequency","truth_maf","r","r2","concordance","CHROM"))' "per-variant schema is incorrect"
assert_parquet_r "${OUT_ENABLED}/metrics/per_variant_metrics" \
    'any(x$CHROM=="Chr01" & x$POS==100 & abs(x$r2-1)<1e-12 & x$concordance==1, na.rm=TRUE)' "known perfect GT result was incorrect"
assert_parquet_r "${OUT_ENABLED}/metrics/per_variant_metrics" \
    'any(x$CHROM=="Chr01" & x$POS==100 & x$ID=="v100" & x$ID_truth=="t100")' "retained match IDs were not preserved"
assert_parquet_r "${OUT_ENABLED}/metrics/per_variant_metrics" \
    'any(x$CHROM=="Chr01" & x$POS==110 & x$n_pairs==5, na.rm=TRUE)' "GQ 59 should mask one sample only"
assert_parquet_r "${OUT_ENABLED}/metrics/per_variant_metrics" \
    'any(x$CHROM=="Chr01" & x$POS==230 & x$n_pairs==2 & is.na(x$r2))' "per-variant r2 must be NA below three pairs"
assert_parquet_r "${OUT_ENABLED}/metrics/per_variant_metrics" \
    'any(x$CHROM=="Chr01" & x$POS==220 & is.na(x$r) & is.na(x$r2))' "zero-variance metric must be NA"
assert_parquet_r "${OUT_ENABLED}/metrics/site_filtered_variants" \
    'any(x$CHROM=="Chr01" & x$POS==130 & !x$site_filter_pass)' "QUAL failure did not remove the site globally"
assert_parquet_r "${OUT_ENABLED}/metrics/per_variant_metrics" \
    'any(x$CHROM=="Chr01" & x$POS==140)' "exact threshold boundary should pass"
assert_parquet_r "${OUT_ENABLED}/metrics/per_variant_metrics" \
    'any(x$CHROM=="Chr01" & x$POS==150)' "missing rank-sum annotations should not fail"
assert_parquet_r "${OUT_ENABLED}/metrics/site_filtered_variants" \
    'any(x$CHROM=="Chr01" & x$POS==160 & !x$site_filter_pass)' "missing QD should fail"
assert_tsv "${OUT_ENABLED}/qc/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_GQ_below_min" && $4==1 {ok=1} END {exit !ok}' "GQ boundary masking count is wrong"
assert_tsv "${OUT_ENABLED}/qc/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_DP_below_min" && $4==1 {ok=1} END {exit !ok}' "DP boundary masking count is wrong"
assert_tsv "${OUT_ENABLED}/qc/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="imputed_missing_GT" && $4==1 {ok=1} END {exit !ok}' "missing imputed GT was not masked"
assert_tsv "${OUT_ENABLED}/qc/genotype_masking_summary.tsv" '$1=="S2" && $2=="ALL" && $3=="imputed_invalid_GT" && $4==1 {ok=1} END {exit !ok}' "haploid imputed GT was not masked"
assert_parquet_r "${OUT_ENABLED}/metrics/allele_mismatches" 'any(x$CHROM=="Chr01" & x$POS==170)' "allele mismatch was not reported"
assert_parquet_r "${OUT_ENABLED}/metrics/per_variant_metrics" 'any(x$CHROM=="Chr01" & x$POS==175 & x$ALT=="C")' "shared allele at a multi-record position was not retained"
assert_parquet_r "${OUT_ENABLED}/metrics/allele_mismatches" 'any(x$CHROM=="Chr01" & x$POS==175 & x$ALT_imputed=="G" & x$ALT_truth=="C")' "common-plus-extra allele mismatch was not reported"
assert_parquet_r "${OUT_ENABLED}/metrics/imputed_only_variants" 'any(x$CHROM=="Chr01" & x$POS==180)' "imputed-only variant was not reported"
assert_parquet_r "${OUT_ENABLED}/metrics/truth_only_variants" 'any(x$CHROM=="Chr01" & x$POS==190)' "truth-only variant was not reported"
assert_parquet_r "${OUT_ENABLED}/metrics/per_variant_metrics" '!any(x$CHROM=="Chr01" & x$POS==200)' "multiallelic site should have been structurally excluded"
grep -Fq $'QUILT2_WGS_TRUTH_MIN_GQ\t60' "${OUT_ENABLED}/run_manifest.tsv" || fail "default GQ was not recorded"
grep -Fq $'QUILT2_WGS_TRUTH_MIN_DP\t10' "${OUT_ENABLED}/run_manifest.tsv" || fail "default DP was not recorded"
grep -Fq $'output_schema\twgs-gt-isec-v4' "${OUT_ENABLED}/run_manifest.tsv" || fail "GT/isec schema was not recorded"
grep -Fq $'comparison_field\tGT' "${OUT_ENABLED}/run_manifest.tsv" || fail "GT comparison field was not recorded"
grep -Fq $'genotype_encoding\tALT_COUNT_0_1_2' "${OUT_ENABLED}/run_manifest.tsv" || fail "GT encoding was not recorded"
grep -Fq $'intersection_key\tCHROM:POS:REF:ALT' "${OUT_ENABLED}/run_manifest.tsv" || fail "exact intersection key was not recorded"
grep -Fq $'intersection_tool\tbcftools_isec' "${OUT_ENABLED}/run_manifest.tsv" || fail "intersection tool was not recorded"
! grep -q 'KEEP_DOSAGE\\|imputed_ds\\|truth_gt_dosage' "${OUT_ENABLED}/run_manifest.tsv" || fail "obsolete dosage outputs remain in the manifest"
grep -Fq $'reference_fasta\t'"${WORK_DIR}/reference.fa" "${OUT_ENABLED}/run_manifest.tsv" || fail "configured reference FASTA was not used"

# Resume should reuse only the completed WGS-mode output.
resume_message="$(SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --reference-fasta "${WORK_DIR}/reference.fa" --out-prefix "${OUT_ENABLED}" 2>&1)"
[[ "${resume_message}" == *"Complete WGS evaluation already exists"* ]] || fail "resume did not reuse complete output"

# An incomplete run should skip valid chromosome checkpoints and rebuild only the missing chromosome.
rm -f "${OUT_ENABLED}/.complete" "${OUT_ENABLED}/intermediate/checkpoints/Chr02.done"
rm -f "${OUT_ENABLED}/metrics/per_variant_metrics/CHROM=Chr02/part-000.parquet"
partial_resume="$(SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --reference-fasta "${WORK_DIR}/reference.fa" --out-prefix "${OUT_ENABLED}" 2>&1)"
[[ "${partial_resume}" == *"skipping Chr01"* ]] || fail "resume did not skip the completed chromosome"
[[ "${partial_resume}" == *"Completed WGS chromosome task 2: Chr02"* ]] || fail "resume did not rebuild the incomplete chromosome"
assert_file "${OUT_ENABLED}/.complete"

# Disabling filtering retains site failures and skips GQ/DP masking, but invalid GT still masks.
OUT_DISABLED="${WORK_DIR}/output/dosage_eval_wgs_unfiltered"
QUILT2_WGS_TRUTH_FILTER_ENABLED=false SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --reference-fasta "${WORK_DIR}/reference.fa" --out-prefix "${OUT_DISABLED}"
assert_parquet_r "${OUT_DISABLED}/metrics/per_variant_metrics" \
    'any(x$CHROM=="Chr01" & x$POS==130)' "disabled filtering still applied QUAL"
assert_parquet_r "${OUT_DISABLED}/metrics/per_variant_metrics" \
    'any(x$CHROM=="Chr01" & x$POS==160)' "disabled filtering still applied missing QD"
assert_tsv "${OUT_DISABLED}/qc/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_GQ_below_min" && $4==0 {ok=1} END {exit !ok}' "disabled filtering still applied GQ"
assert_tsv "${OUT_DISABLED}/qc/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_invalid_GT" && $4==1 {ok=1} END {exit !ok}' "disabled filtering weakened GT validation"

# Custom thresholds are loaded from the environment, applied, and recorded.
OUT_CUSTOM="${WORK_DIR}/output/dosage_eval_wgs_custom"
QUILT2_WGS_TRUTH_MIN_GQ=59 QUILT2_WGS_TRUTH_MIN_DP=9 SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" --chr Chr01 \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --reference-fasta "${WORK_DIR}/reference.fa" --out-prefix "${OUT_CUSTOM}" \
    -- --region Chr01:100-240
assert_tsv "${OUT_CUSTOM}/qc/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_GQ_below_min" && $4==0 {ok=1} END {exit !ok}' "custom GQ was not applied"
assert_tsv "${OUT_CUSTOM}/qc/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_DP_below_min" && $4==0 {ok=1} END {exit !ok}' "custom DP was not applied"
grep -Fq $'QUILT2_WGS_TRUTH_MIN_GQ\t59' "${OUT_CUSTOM}/run_manifest.tsv" || fail "custom GQ was not recorded"
grep -Fq $'QUILT2_WGS_TRUTH_MIN_DP\t9' "${OUT_CUSTOM}/run_manifest.tsv" || fail "custom DP was not recorded"

if QUILT2_WGS_TRUTH_MIN_GQ=61 QUILT2_WGS_TRUTH_MIN_DP=11 SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" --chr Chr01 \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --reference-fasta "${WORK_DIR}/reference.fa" --out-prefix "${OUT_CUSTOM}" \
    -- --region Chr01:100-240 >/dev/null 2>&1; then
    fail "run-signature mismatch was accepted without --force"
fi

# --force recomputes the same WGS cache with new configured values.
QUILT2_WGS_TRUTH_MIN_GQ=61 QUILT2_WGS_TRUTH_MIN_DP=11 SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" --chr Chr01 \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --reference-fasta "${WORK_DIR}/reference.fa" --out-prefix "${OUT_CUSTOM}" \
    -- --region Chr01:100-240 --force
grep -Fq $'QUILT2_WGS_TRUTH_MIN_GQ\t61' "${OUT_CUSTOM}/run_manifest.tsv" || fail "--force did not refresh configured values"

# Invalid configuration must fail before VCF processing.
if QUILT2_WGS_TRUTH_MIN_GQ=100 bash "${ROOT_DIR}/modules/evaluate/dosage_r2_wgs.sh" \
    --imputed /does/not/exist --truth-dataset-dir /does/not/exist \
    --reference-fasta /does/not/exist --out-prefix "${WORK_DIR}/invalid" >/dev/null 2>&1; then
    fail "invalid GQ configuration was accepted"
fi

# GT is mandatory on both sides; DS or another FORMAT field is not a fallback.
if QUILT2_REFERENCE_FASTA="${WORK_DIR}/reference.fa" SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.no_gt.vcf.gz" --chr Chr01 \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --out-prefix "${WORK_DIR}/output/no_gt" >/dev/null 2>&1; then
    fail "imputed VCF without FORMAT/GT was accepted"
fi

# A non-WGS cache cannot be reused by WGS mode.
mkdir -p "${WORK_DIR}/output/collision"
printf 'array output\n' > "${WORK_DIR}/output/collision/per_sample_metrics.tsv"
if SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --reference-fasta "${WORK_DIR}/reference.fa" --out-prefix "${WORK_DIR}/output/collision" >/dev/null 2>&1; then
    fail "cross-mode cache reuse was accepted"
fi

# Array routing remains on the original evaluator and does not validate WGS settings.
array_route="$(QUILT2_WGS_TRUTH_MIN_GQ=not-a-number SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode array --imputed "${WORK_DIR}/imputed.vcf.gz" \
    --truth "${WORK_DIR}/truth/7.Consolidated_VCF/Chr01_consolidated.vcf.gz" \
    --out-prefix "${WORK_DIR}/output/array_route" -- --no-parquet 2>&1 || true)"
[[ "${array_route}" == *"modules/evaluate/dosage_r2.sh"* ]] || fail "array mode did not route to the original evaluator"
[[ "${array_route}" != *"MIN_GQ must"* ]] || fail "array mode evaluated WGS-only configuration"

# Submission mode should create a capped chromosome array and an afterok finalizer.
MOCK_BIN="${WORK_DIR}/mock_bin"
MOCK_SBATCH_LOG="${WORK_DIR}/mock_sbatch.log"
MOCK_SBATCH_COUNTER="${WORK_DIR}/mock_sbatch.counter"
mkdir -p "${MOCK_BIN}"
cat > "${MOCK_BIN}/sbatch" <<'EOF'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$*" >> "${MOCK_SBATCH_LOG}"
count=0
[[ ! -f "${MOCK_SBATCH_COUNTER}" ]] || count="$(<"${MOCK_SBATCH_COUNTER}")"
count=$((count + 1))
printf '%s\n' "${count}" > "${MOCK_SBATCH_COUNTER}"
printf '%s\n' "$((9000 + count))"
EOF
chmod +x "${MOCK_BIN}/sbatch"
export MOCK_SBATCH_LOG MOCK_SBATCH_COUNTER
OUT_SUBMIT="${WORK_DIR}/output/dosage_eval_wgs_submit"
PATH="${MOCK_BIN}:${PATH}" QUILT2_REFERENCE_FASTA="${WORK_DIR}/reference.fa" bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" --out-prefix "${OUT_SUBMIT}" >/dev/null
sed -n '1p' "${MOCK_SBATCH_LOG}" | grep -Fq -- '--array=1-2%4' || fail "WGS chromosome array was not capped at four"
sed -n '2p' "${MOCK_SBATCH_LOG}" | grep -Fq -- '--dependency=afterok:9001' || fail "WGS finalizer dependency is incorrect"

echo "PASS: exact-isec WGS GT-to-GT integration tests"
