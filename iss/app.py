#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import bam
from iss import abundance
from iss import generator
from Bio import SeqIO

import os
import sys
import logging
import argparse


def generate_reads(args):
    """Main function for the `iss generate` submodule

    This submodule generates reads from an ErrorModel and write them to
        args.output + _R(1|2).fastq

    Args:
        args (object): the command-line arguments from argparse
    """
    logger = logging.getLogger(__name__)
    logger.debug('Using verbose logger')

    try:  # try to load the correct error model
        logger.info('Starting iss generate')
        logger.info('Using %s ErrorModel' % args.model)
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
        logger.error('Failed to import ErrorModel module: %s' % e)
        sys.exit(1)

    # read the abundance file
    abundance_dic = abundance.parse_abundance_file(args.abundance)

    try:  # read genomes and generate reads
        assert os.stat(args.genomes).st_size != 0
        f = open(args.genomes, 'r')
    except IOError as e:
        logger.error('Failed to open genome(s) file:%s' % e)
        sys.exit(1)
    except AssertionError as e:
        logger.error('Genome(s) file seems empty: %s' % args.genomes)
        sys.exit(1)
    else:
        with f:
            fasta_file = SeqIO.parse(f, 'fasta')
            n_records = 0
            for record in fasta_file:
                try:
                    n_records += 1
                    species_abundance = abundance_dic[record.id]
                except KeyError as e:
                    logger.error(
                        'Fasta record not found in abundance file: %s' % e)
                    sys.exit(1)
                else:
                    logger.info('Generating reads for record: %s' % record.id)
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

            # check if at least one record was in fasta file
            try:
                assert n_records != 0
            except AssertionError as e:
                logger.error(
                    'Failed to find records in genome(s) file:%s'
                    % args.genomes)
                sys.exit(1)
            else:
                logger.info('Read generation complete')


def model_from_bam(args):
    """Main function for the `iss model` submodule

    This submodule write all variables necessary for building an ErrorModel
    to args.output + .npz

    Args:
        args (object): the command-line arguments from argparse
    """
    logger = logging.getLogger(__name__)
    logger.debug('Using verbose logger')

    try:
        logger.info('Starting iss model')
        from iss import bam
    except ImportError as e:
        logger.error('Failed to import bam module: %s' % e)
        sys.exit(1)
    else:
        logger.info('Using %s ErrorModel' % args.model)
        bam.to_model(args.bam, args.model, args.output)
        logger.info('Model generation complete')


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
        '--quiet',
        '-q',
        action='store_true',
        default=False,
        help='Disable info logging. (default: %(default)s).'
    )
    parser_gen.add_argument(
        '--debug',
        '-d',
        action='store_true',
        default=False,
        help='Enable debug logging. (default: %(default)s).'
    )
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
        '--quiet',
        '-q',
        action='store_true',
        default=False,
        help='Disable info logging. (default: %(default)s).'
    )
    parser_mod.add_argument(
        '--debug',
        '-d',
        action='store_true',
        default=False,
        help='Enable debug logging. (default: %(default)s).'
    )
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

    # set logger
    try:
        if args.quiet:
            logging.basicConfig(level=logging.ERROR)
        elif args.debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        args.func(args)
    except AttributeError as e:
        logger = logging.getLogger(__name__)
        logger.debug(e)
        parser.print_help()
