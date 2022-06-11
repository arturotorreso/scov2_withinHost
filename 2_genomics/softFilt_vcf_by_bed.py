#!/usr/bin/env python3

import pysam
import argparse
from pysam import VariantFile


def main(args):

    # VCF file reading
    vcf = pysam.VariantFile(args.vcf, "r")

    bed = []
    with open(args.bed, "r") as f:
        for line in f:
            bed.append(line.strip().split())

    if args.tag not in vcf.header.filters:
        vcf.header.filters.add(args.tag, None, None, args.filt)

    if args.out_vcf is False:
        vcf_out = VariantFile('-', 'w', header=vcf.header)
    else:
        vcf_out = VariantFile(args.out_vcf, 'w', header=vcf.header)

    for variant in vcf:
        pos = variant.pos
        pos_filt = any([pos >= int(reg[1]) and pos <= int(reg[2]) for reg in bed])

        if (pos_filt):
            variant.filter.add(args.tag)

        vcf_out.write(variant)


parser = argparse.ArgumentParser(description='Soft-filter VCF file using a BED file.')
parser.add_argument("-i", "--vcf", dest="vcf", required=True, help="Input VCF file", type=str)
parser.add_argument("-b", "--bed", dest="bed", required=True, help="Input BED file", type=str)
parser.add_argument("-t", "--tag", dest="tag", help="Tag to add", default='LowComplex', type=str)
parser.add_argument("-f", "--filter", dest="filt", help="Tag to add", default="Position in a low complexity region", type=str)
parser.add_argument("-o", "--output", dest="out_vcf", help="Output VCF file", default=False)
args = parser.parse_args()

main(args)
