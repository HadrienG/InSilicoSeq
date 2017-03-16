#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import bam
from iss import abundance
from iss import generator
from iss import error_model
from Bio import SeqIO

import argparse


def generate_reads(args):
    abundance_dic = abundance.parse_abundance_file(args.abundance)
    with open(args.genome, 'r') as f:
        fasta_file = SeqIO.parse(f, 'fasta')
        for record in fasta_file:
            species_abundance = abundance_dic[record.id]
            genome_size = len(record.seq)
            coverage = abundance.to_coverage(
                args.n_reads,
                species_abundance,
                125,
                genome_size
                )

            read_gen = generator.reads(
                record,
                125,
                coverage,
                args.insert_size,
                30
                    )

            generator.to_fastq(read_gen, args.output)


def model_from_bam(args):
    i_size = bam.get_insert_size(args.bam)
    hist_forward, hist_reverse = bam.quality_distribution(args.bam)
    read_length = len(hist_forward)
    sub_forward, sub__reverse = bam.substitutions(args.bam, read_length)
    bam.write_to_file(
        read_length,
        hist_forward,
        hist_reverse,
        sub_forward,
        sub__reverse,
        i_size,
        args.output + '.npz')


def main():
    parser = argparse.ArgumentParser(
        prog='InSilicoSeq',
        usage='iss [options]',
        description='InSilicoSeq: A sequencing simulator'
    )
    subparsers = parser.add_subparsers(
            title='available commands',
            metavar=''
    )

    parser_mod = subparsers.add_parser(
        'mod',
        prog='iss mod',
        description='generate an error model from a bam file',
        help='generate an error model from a bam file'
    )
    parser_gen = subparsers.add_parser(
        'gen',
        prog='iss gen',
        description='simulate reads from an error model',
        help='simulate reads from an error model'
    )

    # arguments form the read generator module
    parser_gen.add_argument(
        '--genome',
        '-g',
        metavar='<fasta>',
        help='Input genome(s) from where the reads will originate (Required)',
        required=True
    )
    parser_gen.add_argument(
        '--abundance',
        '-a',
        metavar='<txt>',
        help='abundance file for coverage calculations (default: %(default)s)'
    )
    parser_gen.add_argument(
        '--n_reads',
        '-n',
        metavar='<int>',
        type=int,
        default=1000000,
        help='Number of reads to generate (default: %(default)s)'
    )
    parser_gen.add_argument(
        '--insert_size',
        '-i',
        metavar='<int>',
        type=int,
        default=200,
        help='Insert size for paired-end data (default: %(default)s)'
    )
    parser_gen.add_argument(
        '--output',
        '-o',
        metavar='<fastq>',
        help='Output file prefix (Required)',
        required=True
    )
    parser_gen._optionals.title = 'arguments'
    parser_gen.set_defaults(func=generate_reads)

    # arguments for the error model module
    parser_mod.add_argument(
        '--bam',
        '-b',
        metavar='<bam>',
        help='aligned reads from which the model will be inferred (Required)',
        required=True
    )
    parser_mod.add_argument(
        '--output',
        '-o',
        metavar='<npy>',
        help='Output file prefix (Required)',
        required=True
    )
    parser_mod._optionals.title = 'arguments'
    parser_mod.set_defaults(func=model_from_bam)

    args = parser.parse_args()
    args.func(args)
