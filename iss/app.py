#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import generator
from Bio import SeqIO

import argparse


def main():
    parser = argparse.ArgumentParser(
        prog='InSilicoSeq',
        usage='iss [options]',
        description='InSilicoSeq: A sequencing simulator'
    )

    parser.add_argument(
        '--genome',
        '-g',
        metavar='<fasta>',
        help='Input genome from where the reads will originate (Required)',
        required=True
    )
    parser.add_argument(
        '--length',
        '-l',
        metavar='<int>',
        type=int,
        default=150,
        help='Read length (default: %(default)s)'
        )
    parser.add_argument(
        '--coverage',
        '-c',
        metavar='<int>',
        type=int,
        default=10,
        help='Coverage (default: %(default)s)'
        )
    parser.add_argument(
        '--output',
        '-o',
        metavar='<fastq>',
        help='Output file (Required)',
        required=True
        )
    parser._optionals.title = 'arguments'
    args = parser.parse_args()

    with open(args.genome, 'r') as f:
        fasta_file = SeqIO.parse(f, 'fasta')
        for record in fasta_file:
            read_gen = generator.reads(
                record,
                args.length,
                args.coverage
                )
            break  # breaks after the 1st record. We only want a single-fasta!
        generator.to_fastq(read_gen, args.output)
