#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import bam
from iss import abundance
from iss import generator
from Bio import SeqIO

import argparse
import sys


def generate_reads(args):
    try:  # try to load the correct error model
        if args.model == 'cdf':
            from iss.error_models import cdf
            npz = args.model_file
            err_mod = cdf.CDFErrorModel(npz)
        elif args.model == 'kde':
            from iss.error_models import kde
            npz = args.model_file
            err_mod = kde.KDErrorModel(npz)
        elif args.model == 'basic':
            from iss.error_models import basic
            err_mod = basic.BasicErrorModel()
    except ImportError as e:
        print('Error:', e)
        sys.exit(1)
    except FileNotFoundError as e:
        print('Error:', e)
        sys.exit(1)

    # read the abundance file
    abundance_dic = abundance.parse_abundance_file(args.abundance)

    try:  # read genomes and generate reads
        f = open(args.genomes, 'r')
    except IOError as e:
        print('Error', e)
        sys.exit(1)
    else:
        with f:
            fasta_file = SeqIO.parse(f, 'fasta')
            for record in fasta_file:
                try:
                    species_abundance = abundance_dic[record.id]
                except KeyError as e:
                    print('Error:', e)
                    sys.exit(1)
                else:
                    genome_size = len(record.seq)
                    coverage = abundance.to_coverage(
                        args.n_reads,
                        species_abundance,
                        err_mod.read_length,
                        genome_size
                        )

                    read_gen = generator.reads(
                        record,
                        coverage,
                        err_mod
                        )

                    generator.to_fastq(read_gen, args.output)


def model_from_bam(args):
    i_size = bam.get_insert_size(args.bam)
    hist_forward, hist_reverse = bam.quality_distribution(args.model, args.bam)
    read_length = len(hist_forward)
    sub_forward, sub_reverse, ins_forward, \
        ins_reverse, del_forward, del_reverse = \
        bam.get_mismatches(args.bam, read_length)
    bam.write_to_file(
        read_length,
        hist_forward,
        hist_reverse,
        sub_forward,
        sub_reverse,
        ins_forward,
        ins_reverse,
        del_forward,
        del_reverse,
        i_size,
        args.output + '.npz')


def main():
    parser = argparse.ArgumentParser(
        prog='InSilicoSeq',
        usage='iss [subcammand] [options]',
        description='InSilicoSeq: A sequencing simulator'
    )
    subparsers = parser.add_subparsers(
            title='available subcommands',
            metavar=''
    )

    parser_mod = subparsers.add_parser(
        'model',
        prog='iss model',
        description='generate an error model from a bam file',
        help='generate an error model from a bam file'
    )
    parser_gen = subparsers.add_parser(
        'generate',
        prog='iss generate',
        description='simulate reads from an error model',
        help='simulate reads from an error model'
    )

    # arguments form the read generator module
    parser_gen.add_argument(
        '--genomes',
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
        '--model',
        '-m',
        metavar='[\'cdf\', \'kde\', \'basic\']',
        choices=['cdf', 'kde', 'basic'],
        default='kde',
        help='Error model. If not specified, using kernel density estimation \
        (default: %(default)s). Can be \'kde\', \'cdf\' or \'basic\''
    )
    parser_gen.add_argument(
        '--model_file',
        '-f',
        metavar='<npz>',
        default=None,
        help='Error model file. If not specified, using a basic \
        error model instead (default: %(default)s)'
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
        '--model',
        '-m',
        metavar='[\'cdf\', \'kde\']',
        choices=['cdf', 'kde'],
        default='kde',
        help='Error model to generate. (default: %(default)s). \
        Can be \'kde\' or \'cdf\''
    )
    parser_mod.add_argument(
        '--output',
        '-o',
        metavar='<npz>',
        help='Output file prefix (Required)',
        required=True
    )
    parser_mod._optionals.title = 'arguments'
    parser_mod.set_defaults(func=model_from_bam)

    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError:
        parser.print_help()
