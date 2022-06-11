#!/usr/bin/env python3

import sys
import argparse
from Bio import AlignIO
from Bio.SeqIO import FastaIO
import pysam


def main(args):
    if args.input_fa == '-':
        with sys.stdin as reader:
            alignment = AlignIO.read(reader, "fasta")
    else:
        alignment = AlignIO.read(args.input_fa, "fasta")

    dictvcf = load_vcf(args.vcf_file)

    for p in range(alignment.get_alignment_length()):
        # Position is ambiguous base, change to -
        base = alignment[0].seq[p].upper()
        base_ext = ambiguous_dna_values[base].upper()

        if base not in ['A', 'C', 'G', 'T', 'N', '-']:
            # Get number of reads for each base
            base1 = dictvcf[p+1][base_ext[0]]
            base2 = dictvcf[p+1][base_ext[1]]

            # Get index of order based on number of reads
            base_index = sorted(range(len([base1, base2])), key=[base1, base2].__getitem__, reverse=True)
            # Reorder the bases in MajMin form
            new_char = list(base_ext)[base_index[0]].upper() + list(base_ext)[base_index[1]].lower()
        else:
            new_char = base_ext
        base_ext = maj_min_dna_values[new_char]
        if base not in ['A', 'C', 'G', 'T']:
            # Change it in position
            alignment[0].seq = alignment[0].seq.tomutable()
            alignment[0].seq[p] = base_ext
            alignment[0].seq = alignment[0].seq.toseq()

    if args.output_fa:
        output_handle = open(args.output_fa, "w")
        fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_out.write_file(alignment)
        # AlignIO.write(alignment, output_handle, "fasta")
        output_handle.close()
    else:
        print(alignment.format("fasta"))


def load_vcf(vcf_file):
    outvcf = {}
    invcf = pysam.VariantFile(vcf_file, "r")
    for variant in invcf:
        pos = variant.pos
        outvcf[pos] = {}
        outvcf[pos][variant.ref] = variant.samples[0]['AD'][0]
        outvcf[pos][variant.alts[0]] = variant.samples[0]['AD'][1]
    return(outvcf)


ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "N": "N",
    "-": "-"
}


maj_min_dna_values = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'N': '-',
    '-': '-',
    'Ac': 'M',
    'Ag': 'R',
    'At': 'W',
    'Cg': 'S',
    'Ct': 'Y',
    'Gt': 'K',
    'Ca': 'N',
    'Ga': 'D',
    'Ta': 'Q',
    'Gc': 'E',
    'Tc': 'H',
    'Tg': 'I'
}

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_fa', required=True,
                        help='Input pseudosequence file (- for stdin)')
    parser.add_argument('-v', '--vcf', dest='vcf_file', required=True,
                        help='''File with sites to convert to maj/min hets.''')
    parser.add_argument('-o', '--output', dest='output_fa', default=False,
                        help='Output pseudosequence file (Default: stdout)')

    args = parser.parse_args()

    if not len(sys.argv) > 1:
        parser.print_help()
        exit(1)

    main(args)
