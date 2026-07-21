#!/bin/bash
# End-to-end synthetic test for WGS DS-versus-GT dosage evaluation.

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WRAPPER="${ROOT_DIR}/bin/dosage_r2_sbatch.sh"
WORK_DIR="$(mktemp -d "${TMPDIR:-/tmp}/test_dosage_r2_wgs.XXXXXX")"
trap 'rm -rf "${WORK_DIR}"' EXIT

for cmd in bcftools bgzip Rscript; do
    command -v "${cmd}" >/dev/null 2>&1 || { echo "SKIP: ${cmd} is unavailable"; exit 0; }
done

fail() { echo "FAIL: $*" >&2; exit 1; }
assert_file() { [[ -s "$1" ]] || fail "missing or empty file: $1"; }
assert_tsv() {
    local file="$1" expression="$2" message="$3"
    awk -F '\t' "${expression}" "${file}" || fail "${message}"
}

mkdir -p "${WORK_DIR}/truth/7.Consolidated_VCF"
{
    printf '>Chr01\n'
    awk 'BEGIN { for (i=0; i<1000; i++) printf "A"; printf "\n" }'
} > "${WORK_DIR}/reference.fa"
if command -v samtools >/dev/null 2>&1; then
    samtools faidx "${WORK_DIR}/reference.fa"
else
    printf 'Chr01\t1000\t7\t1000\t1001\n' > "${WORK_DIR}/reference.fa.fai"
fi

cat > "${WORK_DIR}/imputed.vcf" <<'EOF'
##fileformat=VCFv4.2
##contig=<ID=Chr01,length=1000>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=1,Type=Float,Description="ALT dosage">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3	S4	S5	S6
Chr01	100	v100	A	C	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	110	v110	A	C	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	120	v120	A	C	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	130	v130	A	C	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	140	v140	A	C	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	150	v150	A	C	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	160	v160	A	C	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	170	v170i	A	G	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	180	v180	A	C	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	200	v200	A	C,G	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	210	v210	A	C	.	PASS	.	GT:DS	0/0:.	0/1:2.5	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	220	v220	A	C	.	PASS	.	GT:DS	0/1:1.0	0/1:1.0	0/1:1.0	0/1:1.0	0/1:1.0	0/1:1.0
Chr01	230	v230	A	C	.	PASS	.	GT:DS	0/0:0.0	1/1:2.0	0/1:1.0	0/0:0.0	0/1:1.0	1/1:2.0
Chr01	240	v240	A	C	.	PASS	.	GT:DS	0/0:0.0	0/1:1.0	1/1:2.0	0/0:0.0	0/1:1.0	1/1:2.0
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
Chr01	190	t190	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	200	t200	A	C,G	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	210	t210	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	0/1:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
Chr01	220	t220	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/1:99:20	0/1:99:20	0/1:99:20	0/1:99:20	0/1:99:20	0/1:99:20
Chr01	230	t230	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0/0:99:20	1/1:99:20	0/1:59:20	0/0:59:20	0/1:59:20	1/1:59:20
Chr01	240	t240	A	C	100	PASS	QD=10;SOR=1;FS=1;MQ=60;MQRankSum=0;ReadPosRankSum=0	GT:GQ:DP	0:99:20	./.:99:20	1/1:99:20	0/0:99:20	0/1:99:20	1/1:99:20
EOF

awk '{ gsub(/\\t/, "\t"); print }' "${WORK_DIR}/truth/Chr01_consolidated.vcf" > "${WORK_DIR}/truth/Chr01_consolidated.tabs.vcf"
mv "${WORK_DIR}/truth/Chr01_consolidated.tabs.vcf" "${WORK_DIR}/truth/Chr01_consolidated.vcf"

bgzip -c "${WORK_DIR}/imputed.vcf" > "${WORK_DIR}/imputed.vcf.gz"
bcftools index -t "${WORK_DIR}/imputed.vcf.gz"
bgzip -c "${WORK_DIR}/truth/Chr01_consolidated.vcf" > "${WORK_DIR}/truth/7.Consolidated_VCF/Chr01_consolidated.vcf.gz"
bcftools index -t "${WORK_DIR}/truth/7.Consolidated_VCF/Chr01_consolidated.vcf.gz"

OUT_ENABLED="${WORK_DIR}/output/dosage_eval_wgs"
QUILT2_REFERENCE_FASTA="${WORK_DIR}/reference.fa" SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs \
    --imputed "${WORK_DIR}/imputed.vcf.gz" \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --out-prefix "${OUT_ENABLED}"

for file in \
    per_variant_metrics.tsv per_sample_metrics.tsv filter_summary.tsv genotype_masking_summary.tsv \
    matched_variants.tsv imputed_only_variants.tsv truth_only_variants.tsv allele_mismatches.tsv \
    intermediate/imputed_ds.tsv intermediate/truth_gt_dosage.tsv run_manifest.tsv; do
    assert_file "${OUT_ENABLED}/${file}"
done

assert_tsv "${OUT_ENABLED}/per_variant_metrics.tsv" '$2==100 && $10==1 {ok=1} END {exit !ok}' "known perfect r2 was not 1"
assert_tsv "${OUT_ENABLED}/per_variant_metrics.tsv" '$2==110 && $6==5 {ok=1} END {exit !ok}' "GQ 59 should mask one sample only"
assert_tsv "${OUT_ENABLED}/per_variant_metrics.tsv" '$2==230 && $6==2 && $10=="NA" {ok=1} END {exit !ok}' "per-variant r2 must be NA below three pairs"
assert_tsv "${OUT_ENABLED}/per_variant_metrics.tsv" '$2==220 && $9=="NA" && $10=="NA" {ok=1} END {exit !ok}' "zero-variance metric must be NA"
assert_tsv "${OUT_ENABLED}/matched_variants.tsv" '$2==130 && $7=="FALSE" {ok=1} END {exit !ok}' "QUAL failure did not remove the site globally"
assert_tsv "${OUT_ENABLED}/matched_variants.tsv" '$2==140 && $7=="TRUE" {ok=1} END {exit !ok}' "exact threshold boundary should pass"
assert_tsv "${OUT_ENABLED}/matched_variants.tsv" '$2==150 && $7=="TRUE" {ok=1} END {exit !ok}' "missing rank-sum annotations should not fail"
assert_tsv "${OUT_ENABLED}/matched_variants.tsv" '$2==160 && $7=="FALSE" {ok=1} END {exit !ok}' "missing QD should fail"
assert_tsv "${OUT_ENABLED}/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_GQ_below_min" && $4==1 {ok=1} END {exit !ok}' "GQ boundary masking count is wrong"
assert_tsv "${OUT_ENABLED}/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_DP_below_min" && $4==1 {ok=1} END {exit !ok}' "DP boundary masking count is wrong"
assert_tsv "${OUT_ENABLED}/genotype_masking_summary.tsv" '$1=="S2" && $2=="ALL" && $3=="imputed_invalid_DS" && $4==1 {ok=1} END {exit !ok}' "invalid DS was not masked"
assert_tsv "${OUT_ENABLED}/allele_mismatches.tsv" '$2==170 {ok=1} END {exit !ok}' "allele mismatch was not reported"
assert_tsv "${OUT_ENABLED}/imputed_only_variants.tsv" '$2==180 {ok=1} END {exit !ok}' "imputed-only variant was not reported"
assert_tsv "${OUT_ENABLED}/truth_only_variants.tsv" '$2==190 {ok=1} END {exit !ok}' "truth-only variant was not reported"
if awk -F '\t' '$2==200 {found=1} END {exit !found}' "${OUT_ENABLED}/matched_variants.tsv"; then
    fail "multiallelic site should have been structurally excluded"
fi
grep -Fq $'QUILT2_WGS_TRUTH_MIN_GQ\t60' "${OUT_ENABLED}/run_manifest.tsv" || fail "default GQ was not recorded"
grep -Fq $'QUILT2_WGS_TRUTH_MIN_DP\t10' "${OUT_ENABLED}/run_manifest.tsv" || fail "default DP was not recorded"
grep -Fq $'reference_fasta\t'"${WORK_DIR}/reference.fa" "${OUT_ENABLED}/run_manifest.tsv" || fail "configured reference FASTA was not used"

# Resume should reuse only the completed WGS-mode output.
resume_message="$(SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --reference-fasta "${WORK_DIR}/reference.fa" --out-prefix "${OUT_ENABLED}" 2>&1)"
[[ "${resume_message}" == *"Complete WGS evaluation already exists"* ]] || fail "resume did not reuse complete output"

# Disabling filtering retains site failures and skips GQ/DP masking, but invalid GT still masks.
OUT_DISABLED="${WORK_DIR}/output/dosage_eval_wgs_unfiltered"
QUILT2_WGS_TRUTH_FILTER_ENABLED=false SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --reference-fasta "${WORK_DIR}/reference.fa" --out-prefix "${OUT_DISABLED}"
assert_tsv "${OUT_DISABLED}/matched_variants.tsv" '$2==130 && $7=="TRUE" {ok=1} END {exit !ok}' "disabled filtering still applied QUAL"
assert_tsv "${OUT_DISABLED}/matched_variants.tsv" '$2==160 && $7=="TRUE" {ok=1} END {exit !ok}' "disabled filtering still applied missing QD"
assert_tsv "${OUT_DISABLED}/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_GQ_below_min" && $4==0 {ok=1} END {exit !ok}' "disabled filtering still applied GQ"
assert_tsv "${OUT_DISABLED}/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_invalid_GT" && $4==1 {ok=1} END {exit !ok}' "disabled filtering weakened GT validation"

# Custom thresholds are loaded from the environment, applied, and recorded.
OUT_CUSTOM="${WORK_DIR}/output/dosage_eval_wgs_custom"
QUILT2_WGS_TRUTH_MIN_GQ=59 QUILT2_WGS_TRUTH_MIN_DP=9 SLURM_JOB_ID=test bash "${WRAPPER}" \
    --truth-mode wgs --imputed "${WORK_DIR}/imputed.vcf.gz" --chr Chr01 \
    --truth-dataset-dir "${WORK_DIR}/truth/7.Consolidated_VCF" \
    --reference-fasta "${WORK_DIR}/reference.fa" --out-prefix "${OUT_CUSTOM}" \
    -- --region Chr01:100-240
assert_tsv "${OUT_CUSTOM}/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_GQ_below_min" && $4==0 {ok=1} END {exit !ok}' "custom GQ was not applied"
assert_tsv "${OUT_CUSTOM}/genotype_masking_summary.tsv" '$1=="S1" && $2=="ALL" && $3=="truth_DP_below_min" && $4==0 {ok=1} END {exit !ok}' "custom DP was not applied"
grep -Fq $'QUILT2_WGS_TRUTH_MIN_GQ\t59' "${OUT_CUSTOM}/run_manifest.tsv" || fail "custom GQ was not recorded"
grep -Fq $'QUILT2_WGS_TRUTH_MIN_DP\t9' "${OUT_CUSTOM}/run_manifest.tsv" || fail "custom DP was not recorded"

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

echo "PASS: WGS dosage-r2 integration tests"
