"""
Module that contains the command line app, so we can still import __main__
without executing side effects
"""

import argparse
import vcf

try:
    from utils import all_combinations
except ModuleNotFoundError:
    from haploblock_shuffler.utils import all_combinations


def main():
    parser = argparse.ArgumentParser(
        description="Generate all possible phaseing combination"
    )
    parser.add_argument("vcf")
    args = parser.parse_args()

    vcf_in = vcf.Reader(filename=args.vcf)

    all_variants = list(vcf_in)

    for i, variants in enumerate(all_combinations(all_variants)):
        with open(f"out_{i}.vcf", "w") as fout:
            vcf_out = vcf.Writer(fout, template=vcf_in)
            flatten = [record for group in variants for record in group]
            for variant in sorted(flatten, key=lambda x: (x.CHROM, x.POS)):
                vcf_out.write_record(variant)


if __name__ == "__main__":
    main()
