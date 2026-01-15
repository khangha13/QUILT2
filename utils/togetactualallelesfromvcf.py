# vcf_to_alleles.py by Shashi Goonetilleke

from cyvcf2 import VCF
import sys

def convert_genotypes_to_alleles(vcf_path, output_path=None):
    vcf = VCF(vcf_path)
    samples = vcf.samples

    results = []

    for variant in vcf:
        alleles = [variant.REF] + variant.ALT
        line = {
            "CHROM": variant.CHROM,
            "POS": variant.POS,
            "REF": variant.REF,
            "ALT": ",".join(variant.ALT)
        }

        for i, genotype in enumerate(variant.genotypes):
            gt_indices = genotype[:2]
            sample_name = samples[i]
            gt_alleles = [alleles[idx] if idx is not None and idx >= 0 else '.' for idx in gt_indices]
            line[sample_name] = "/".join(gt_alleles)

        results.append(line)

    headers = ["CHROM", "POS", "REF", "ALT"] + samples
    if output_path:
        with open(output_path, "w") as f:
            f.write("\t".join(headers) + "\n")
            for row in results:
                f.write("\t".join(str(row[h]) for h in headers) + "\n")
    else:
        print("\t".join(headers))
        for row in results:
            print("\t".join(str(row[h]) for h in headers))

# Run if script is called from command line
if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python vcf_to_alleles.py input.vcf [output.tsv]")
        sys.exit(1)
    input_vcf = sys.argv[1]
    output_tsv = sys.argv[2] if len(sys.argv) > 2 else None
    convert_genotypes_to_alleles(input_vcf, output_tsv)