#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import generator

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

    read_gen = generator.reads(args.genome, args.length, args.coverage)
    generator.to_fastq(read_gen, args.output)
