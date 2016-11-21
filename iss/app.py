#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import abundance
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
        help='Input genome(s) from where the reads will originate (Required)',
        required=True
    )
    parser.add_argument(
        '--abundance',
        '-a',
        metavar='<txt>',
        help='abundance file for coverage calculations (default: %(default)s)'
        )
    parser.add_argument(
        '--read_length',
        '-l',
        metavar='<int>',
        type=int,
        default=150,
        help='Read length (default: %(default)s)'
        )
    parser.add_argument(
        '--n_reads',
        '-n',
        metavar='<int>',
        type=int,
        default=1000000,
        help='Number of reads to generate (default: %(default)s)'
        )
    parser.add_argument(
        '--insert_size',
        '-i',
        metavar='<int>',
        type=int,
        default=200,
        help='Insert size for paired-end data (default: %(default)s)'
        )
    parser.add_argument(
        '--output',
        '-o',
        metavar='<fastq>',
        help='Output file prefix (Required)',
        required=True
        )
    parser._optionals.title = 'arguments'
    args = parser.parse_args()

    abundance_dic = abundance.parse_abundance_file(args.abundance)
    with open(args.genome, 'r') as f:
        fasta_file = SeqIO.parse(f, 'fasta')
        for record in fasta_file:
            species_abundance = abundance_dic[record.id]
            genome_size = len(record.seq)
            coverage = abundance.to_coverage(
                args.n_reads,
                species_abundance,
                args.read_length,
                genome_size
                )

            read_gen = generator.reads(
                record,
                args.read_length,
                coverage,
                args.insert_size
                )
            # for read_pairs in read_gen:
            #     print(read_pairs)
            #     break
            generator.to_fastq(read_gen, args.output)
