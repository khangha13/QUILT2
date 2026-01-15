#!/usr/bin/env python3
"""
Translate a WGS-style VCF (GT coded as allele indexes) into SNP array-style
genotype codes:

  - A/A, A/T, T/A, T/T -> 0/0
  - C/C, C/G, G/C, G/G -> 1/1
  - any mix of (A|T) with (C|G) -> 0/1
  - anything else (non-ACGT, non-diploid, missing) -> ./.

Assumes the input VCF is already biallelic SNP-only (as produced by the pipeline).
All INFO/FORMAT fields are preserved, only GT is recoded.
"""

import argparse
import sys

try:
    import pysam
except ImportError:  # pragma: no cover
    # Keep dependency failure obvious and friendly.
    sys.stderr.write("pysam is required to run this script: pip install pysam\n")
    raise


def allele_group(base: str):
    """Return group tag for allele or None if not A/C/G/T."""
    if base in ("A", "T"):
        return "AT"
    if base in ("C", "G"):
        return "CG"
    return None


def recode_genotype(gt, allele_lookup):
    """
    Recode a diploid genotype tuple of allele indexes to array-style GT.
    Returns a tuple of ints/None suitable for pysam (None -> '.').
    """
    if gt is None or len(gt) != 2:
        return (None, None)

    bases = []
    for idx in gt:
        if idx is None or idx < 0:
            return (None, None)
        base = allele_lookup.get(idx)
        if base is None or len(base) != 1:
            return (None, None)
        group = allele_group(base.upper())
        if group is None:
            return (None, None)
        bases.append(group)

    if bases[0] == "AT" and bases[1] == "AT":
        return (0, 0)
    if bases[0] == "CG" and bases[1] == "CG":
        return (1, 1)
    return (0, 1)


def translate_vcf(input_vcf: str, output_vcf: str):
    mode = "wz" if output_vcf.endswith(".gz") else "w"

    with pysam.VariantFile(input_vcf, "r") as reader:
        header = reader.header.copy()
        header.add_meta(
            "source",
            "wgs_to_array_vcf.py recoded GT (AT->0/0, CG->1/1, mix->0/1)",
        )

        with pysam.VariantFile(output_vcf, mode, header=header) as writer:
            for record in reader:
                allele_lookup = {0: record.ref}
                for i, alt in enumerate(record.alts or [], start=1):
                    allele_lookup[i] = alt

                for sample in record.samples:
                    call = record.samples[sample]
                    phased = call.phased
                    new_gt = recode_genotype(call.get("GT"), allele_lookup)
                    call["GT"] = new_gt
                    call.phased = phased  # preserve original phase indicator

                writer.write(record)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Recodes WGS VCF GT to SNP array-style VCF GT."
    )
    parser.add_argument("-i", "--input", required=True, help="Input VCF/BCF path.")
    parser.add_argument("-o", "--output", required=True, help="Output VCF/BCF path.")
    return parser.parse_args()


def main():
    args = parse_args()
    translate_vcf(args.input, args.output)


if __name__ == "__main__":
    main()
