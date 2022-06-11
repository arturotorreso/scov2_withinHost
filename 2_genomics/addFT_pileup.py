#!/usr/bin/env python

import pysam
import argparse
from pysam import VariantFile


def main(args):

    # VCF file reading
    vcf = pysam.VariantFile(args.vcf, "r")
    vcf = format_header(vcf)

    if args.out_vcf is False:
        vcf_out = VariantFile('-', 'w', header=vcf.header)
    else:
        vcf_out = VariantFile(args.out_vcf, 'w', header=vcf.header)

    for variant in vcf:
        for sample in variant.samples:

            if variant.alts is not None:
                if len(variant.alts) > 1:
                    variant.alts = variant.alts[0]

            if len(variant.samples[sample]['AD']) > 2:
                variant.samples[sample]['AD'] = variant.samples[sample]['AD'][0:2]
                variant.samples[sample]['ADR'] = variant.samples[sample]['ADR'][0:2]
                variant.samples[sample]['ADF'] = variant.samples[sample]['ADF'][0:2]

            dp = variant.samples[sample]['DP']
            ad = variant.samples[sample]['AD']
            gt = variant.samples[sample]['GT'][0]
            adr = variant.samples[sample]['ADR']
            adf = variant.samples[sample]['ADF']

            # Missing genotype
            if gt is None or gt == '.':
                variant.samples[sample]['FT'] = 0
                variant.samples[sample]['GT'] = ('.', '.')
                continue

            # Known genotype: HAPLOID
            ad0 = ad[0]
            ad1 = 0 if len(ad) == 1 else ad[1]

            af1 = ad1/dp
            af0 = ad0/dp

            Maf = max(af0, af1)
            maf = min(af0, af1)

            # One of the alleles passes Max Allele Frequency for Homozygous
            if Maf >= float(args.freq):
                gtmax = [af0, af1].index(max([af0, af1]))
                gt = '{}/{}'.format(gtmax, gtmax)
                variant.samples[sample]['AF'] = Maf
            # Both alleles pass minDP
            elif all([i >= args.het_dp for i in [ad0, ad1]]):
                gt = '0/1'
                variant.samples[sample]['AF'] = maf
            # One of them passess filter
            elif any([i >= args.het_dp for i in [ad0, ad1]]):
                gtmax = [af0, af1].index(max([af0, af1]))
                gt = '{}/{}'.format(gtmax, gtmax)
                variant.samples[sample]['AF'] = Maf
            else:
                gt = '0/0'
                variant.samples[sample]['AF'] = maf
            variant.samples[sample]['GT'] = (int(gt.split('/')[0]), int(gt.split('/')[1]))
            if 'PASS' not in variant.filter:
                variant.samples[sample]['FT'] = 0
                continue

            # Homozygous genotype
            if gt.split('/')[0] == gt.split('/')[1]:
                gt = int(gt.split('/')[0])
                if dp >= 20 and max([ad0, ad1]) >= 20 and adr[gt] >= 5 and adf[gt] >= 5:
                    variant.samples[sample]['FT'] = 1
                    continue
                else:
                    variant.samples[sample]['FT'] = 0
                    continue
            # Heterozygous genotype
            elif (gt.split('/')[0] != gt.split('/')[1]):
                gt0 = int(gt.split('/')[0])
                gt1 = int(gt.split('/')[1])
                if dp >= 100 and ad0 >= 20 and ad1 >= 20 and adr[gt0] >= 5 and adr[gt1] >= 5 and adf[gt0] >= 5 and adf[gt1] >= 5:
                    variant.samples[sample]['FT'] = 1
                    continue
                else:
                    variant.samples[sample]['FT'] = 0
                    continue

        vcf_out.write(variant)


def format_header(vcf):
    # Add the FT field to header.
    if 'FT' not in vcf.header.formats:
        vcf.header.formats.add("FT", "1", "Integer", "Whether a sample was a Pass(1) or fail (0) based on FORMAT values")
    if 'AD' not in vcf.header.formats:
        vcf.header.formats.add("AD", "R", "Integer", "Allelic depths (high-quality bases)")
    if 'DP' not in vcf.header.formats:
        vcf.header.formats.add("DP", "1", "Integer", "Number of high-quality bases")
    if 'AF' not in vcf.header.formats:
        vcf.header.formats.add("AF", "1", "Float", "Frequency of minor allele")
    return(vcf)


parser = argparse.ArgumentParser(description='Soft-filter VCF file using a BED file.')
parser.add_argument("-i", "--vcf", dest="vcf", required=True, help="Input VCF file", type=str)
parser.add_argument("-o", "--output", dest="out_vcf", help="Output VCF file", default=False)
parser.add_argument("-f", "--freq", dest="freq", help="Frequency for consensus call (1/1)", default=0.8, type=float)
parser.add_argument("-d", "--hetDP", dest="het_dp", help="Min DP to consider allele", default=30, type=int)
args = parser.parse_args()

main(args)
