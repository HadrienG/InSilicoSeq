#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import bam
from iss import util
from iss import download
from iss import abundance
from iss import generator
from iss.version import __version__

from Bio import SeqIO
from joblib import Parallel, delayed, load, dump

import gc
import os
import sys
import pickle
import random
import logging
import argparse
import numpy as np


def generate_reads(args):
    """Main function for the `iss generate` submodule

    This submodule generates reads from an ErrorModel and write them to
        args.output + _R(1|2).fastq

    Args:
        args (object): the command-line arguments from argparse
    """
    logger = logging.getLogger(__name__)
    logger.debug('iss version %s' % __version__)
    logger.debug('Using verbose logger')

    try:  # try to import and load the correct error model
        logger.info('Starting iss generate')
        logger.info('Using %s ErrorModel' % args.mode)
        if args.seed:
            logger.info('Setting random seed to %i' % args.seed)
            random.seed(args.seed)
            np.random.seed(args.seed)
        if args.mode == 'kde':
            from iss.error_models import kde
            if args.model is None:
                logger.error('--model is required in --mode kde')
                sys.exit(1)
            elif args.model.lower() == 'hiseq':
                npz = os.path.join(
                    os.path.dirname(__file__),
                    'profiles/HiSeq')
            elif args.model.lower() == 'novaseq':
                npz = os.path.join(
                    os.path.dirname(__file__),
                    'profiles/NovaSeq')
            elif args.model.lower() == 'miseq':
                npz = os.path.join(
                    os.path.dirname(__file__),
                    'profiles/MiSeq')
            else:
                npz = args.model
            err_mod = kde.KDErrorModel(npz)
        elif args.mode == 'basic':
            if args.model is not None:
                logger.warning('--model %s will be ignored in --mode %s' % (
                    args.model, args.mode))
            from iss.error_models import basic
            err_mod = basic.BasicErrorModel()
        elif args.mode == 'perfect':
            if args.model is not None:
                logger.warning('--model %s will be ignored in --mode %s' % (
                    args.model, args.mode))
            from iss.error_models import perfect
            err_mod = perfect.PerfectErrorModel()
    except ImportError as e:
        logger.error('Failed to import ErrorModel module: %s' % e)
        sys.exit(1)

    try:  # try to read genomes and concatenate --genomes and --ncbi genomes
        if args.genomes or args.draft or args.ncbi:
            genome_files = []
            if args.genomes:
                genome_files.extend(args.genomes)
            if args.draft:
                genome_files.extend(args.draft)
            if args.ncbi and args.n_genomes_ncbi:
                util.genome_file_exists(args.output + '_ncbi_genomes.fasta')
                total_genomes_ncbi = []
                try:
                    assert len(*args.ncbi) == len(*args.n_genomes_ncbi)
                except AssertionError as e:
                    logger.error(
                        '--ncbi and --n_genomes_ncbi of unequal lengths. \
                        Aborting')
                    sys.exit(1)
                for g, n in zip(*args.ncbi, *args.n_genomes_ncbi):
                    genomes_ncbi = download.ncbi(
                        g, n, args.output + '_ncbi_genomes.fasta')
                genome_files.append(genomes_ncbi)
            if args.ncbi and not args.n_genomes_ncbi:
                logger.error(
                    '--ncbi/-k requires --n_genomes_ncbi/-U. Aborting.')
                sys.exit(1)

        else:
            logger.error("One of --genomes/-g, --draft, --ncbi/-k is required")
            sys.exit(1)

        genome_file = args.output + '.iss.tmp.genomes.fasta'
        util.concatenate(
            genome_files,
            output=genome_file)

        # for n_genomes we use reservoir sampling to draw random genomes
        # from the concatenated genome file. We then override the file.
        if args.n_genomes and not args.draft and not args.ncbi:
            genome_count = util.count_records(genome_file)
            genome_files = [genome for genome in util.reservoir(
                SeqIO.parse(genome_file, 'fasta'),
                genome_count,
                args.n_genomes)]
            SeqIO.write(genome_files, genome_file, 'fasta')

        assert os.stat(genome_file).st_size != 0
        f = open(genome_file, 'r')
        with f:  # count the number of records
            genome_list = util.count_records(f)
    except IOError as e:
        logger.error('Failed to open genome(s) file:%s' % e)
        sys.exit(1)
    except AssertionError as e:
        logger.error('Genome(s) file seems empty: %s' % genome_file)
        sys.exit(1)
    except KeyboardInterrupt as e:
        logger.error('iss generate interrupted: %s' % e)
        sys.exit(1)
    else:
        abundance_dispatch = {
            'uniform': abundance.uniform,
            'halfnormal': abundance.halfnormal,
            'exponential': abundance.exponential,
            'lognormal': abundance.lognormal,
            'zero_inflated_lognormal': abundance.zero_inflated_lognormal
        }
        # read the abundance file
        if args.abundance_file:
            logger.info('Using abundance file:%s' % args.abundance_file)
            if args.draft:
                abundance_dic_short = abundance.parse_abundance_file(
                    args.abundance_file)
                complete_genomes_dic = {k: v for
                                        k, v in abundance_dic_short.items()
                                        if k not in args.draft}
                draft_dic = abundance.expand_draft_abundance(
                    abundance_dic_short,
                    args.draft)
                abundance_dic = {**complete_genomes_dic,
                                 **draft_dic}
            else:
                abundance_dic = abundance.parse_abundance_file(
                    args.abundance_file)
        elif args.coverage_file:
            logger.warning('--coverage_file is an experimental feature')
            logger.warning('--coverage_file disables --n_reads')
            logger.info('Using coverage file:%s' % args.coverage_file)
            if args.draft:
                coverage_dic = abundance.parse_abundance_file(
                    args.coverage_file)
                complete_genomes_dic = {k: v for
                                        k, v in coverage_dic.items()
                                        if k not in args.draft}
                draft_dic = abundance.expand_draft_abundance(
                    abundance_dic_short,
                    args.draft,
                    mode="coverage")
                abundance_dic = {**complete_genomes_dic,
                                 **draft_dic}
            else:
                abundance_dic = abundance.parse_abundance_file(
                    args.coverage_file)
        elif args.coverage in abundance_dispatch:
            # todo coverage distribution with --draft
            logger.warning('--coverage is an experimental feature')
            logger.info('Using %s coverage distribution' % args.coverage)
            if args.draft:
                abundance_dic = abundance.draft(
                    genome_list,
                    args.draft,
                    abundance_dispatch[args.abundance],
                    args.output,
                    mode="coverage")
            else:
                abundance_dic = abundance_dispatch[
                    args.coverage](genome_list)
            if args.n_reads:
                n_reads = util.convert_n_reads(args.n_reads)
                logger.info('scaling coverage to %s reads' % n_reads)
                abundance_dic = abundance.coverage_scaling(n_reads,
                                                           abundance_dic,
                                                           genome_file,
                                                           err_mod.read_length)
            abundance.to_file(abundance_dic, args.output, mode="coverage")
        elif args.abundance in abundance_dispatch:
            logger.info('Using %s abundance distribution' % args.abundance)
            if args.draft:
                abundance_dic = abundance.draft(
                    genome_list,
                    args.draft,
                    abundance_dispatch[args.abundance],
                    args.output)
            else:
                abundance_dic = abundance_dispatch[
                    args.abundance](genome_list)
                abundance.to_file(abundance_dic, args.output)
        else:
            logger.error('Could not get abundance')
            sys.exit(1)

        cpus = args.cpus
        logger.info('Using %s cpus for read generation' % cpus)

        if not (args.coverage or args.coverage_file):
            n_reads = util.convert_n_reads(args.n_reads)
            logger.info('Generating %s reads' % n_reads)

        try:
            temp_file_list = []  # list holding the prefix of all temp files
            f = open(genome_file, 'r')  # re-opens the file
            with f:
                fasta_file = SeqIO.parse(f, 'fasta')

                for record in fasta_file:
                    # generate reads for records
                    try:
                        species_abundance = abundance_dic[record.id]
                    except KeyError as e:
                        logger.error(
                            'Fasta record not found in abundance file: %s' % e)
                        sys.exit(1)
                    else:
                        logger.info('Generating reads for record: %s'
                                    % record.id)
                        genome_size = len(record.seq)

                        if args.coverage or args.coverage_file:
                            coverage = species_abundance
                        else:
                            coverage = abundance.to_coverage(
                                n_reads,
                                species_abundance,
                                err_mod.read_length,
                                genome_size
                            )
                        n_pairs = int(round(
                            (coverage *
                                len(record.seq)) / err_mod.read_length) / 2)
                        # skip record if n_reads == 0
                        if n_pairs == 0:
                            continue

                        # exact n_reads for each cpus
                        if n_pairs % cpus == 0:
                            n_pairs_per_cpu = [(n_pairs // cpus)
                                               for _ in range(cpus)]
                        else:
                            n_pairs_per_cpu = [(n_pairs // cpus)
                                               for _ in range(cpus)]
                            n_pairs_per_cpu[-1] += n_pairs % cpus

                        # due to a bug in multiprocessing
                        # https://bugs.python.org/issue17560
                        # we can't send records taking more than 2**31 bytes
                        # through serialisation.
                        # In those cases we use memmapping
                        if sys.getsizeof(str(record.seq)) >= 2**31 - 1:
                            logger.warning(
                                "record %s unusually big." % record.id)
                            logger.warning("Using a memory map.")
                            mode = "memmap"

                            record_mmap = "%s.memmap" % args.output
                            if os.path.exists(record_mmap):
                                os.unlink(record_mmap)
                            util.dump(record, record_mmap)
                            del record
                            record = record_mmap
                            gc.collect()
                        else:
                            mode = "default"

                        record_file_name_list = Parallel(
                            n_jobs=cpus)(
                                delayed(generator.reads)(
                                    record, err_mod,
                                    n_pairs_per_cpu[i], i, args.output,
                                    args.seed,
                                    args.gc_bias, mode) for i in range(cpus))
                        temp_file_list.extend(record_file_name_list)
        except KeyboardInterrupt as e:
            logger.error('iss generate interrupted: %s' % e)
            temp_file_unique = list(set(temp_file_list))
            temp_R1 = [temp_file + '_R1.fastq' for temp_file in temp_file_list]
            temp_R2 = [temp_file + '_R2.fastq' for temp_file in temp_file_list]
            full_tmp_list = temp_R1 + temp_R2
            full_tmp_list.append(genome_file)
            if os.path.exists("%s.memmap" % args.output):
                full_tmp_list.append("%s.memmap" % args.output)
            util.cleanup(full_tmp_list)
            sys.exit(1)
        else:
            # remove the duplicates in file list and cleanup
            # we remove the duplicates in case two records had the same header
            # and reads were appended to the same temp file.
            temp_file_unique = list(set(temp_file_list))
            temp_R1 = [temp_file + '_R1.fastq' for temp_file in temp_file_list]
            temp_R2 = [temp_file + '_R2.fastq' for temp_file in temp_file_list]
            util.concatenate(temp_R1, args.output + '_R1.fastq')
            util.concatenate(temp_R2, args.output + '_R2.fastq')
            full_tmp_list = temp_R1 + temp_R2
            full_tmp_list.append(genome_file)
            if os.path.exists("%s.memmap" % args.output):
                full_tmp_list.append("%s.memmap" % args.output)
            util.cleanup(full_tmp_list)
            if args.compress:
                util.compress(args.output + '_R1.fastq')
                util.compress(args.output + '_R2.fastq')
            logger.info('Read generation complete')


def model_from_bam(args):
    """Main function for the `iss model` submodule

    This submodule write all variables necessary for building an ErrorModel
    to args.output + .npz

    Args:
        args (object): the command-line arguments from argparse
    """
    logger = logging.getLogger(__name__)
    logger.debug('iss version %s' % __version__)
    logger.debug('Using verbose logger')

    try:  # try to import bam module and write model data to file
        logger.info('Starting iss model')
        from iss import bam
    except ImportError as e:
        logger.error('Failed to import bam module: %s' % e)
        sys.exit(1)
    else:
        logger.info('Using KDE ErrorModel')
        bam.to_model(args.bam, args.output)
        logger.info('Model generation complete')


def main():
    parser = argparse.ArgumentParser(
        prog='InSilicoSeq',
        usage='iss [subcommand] [options]',
        description='InSilicoSeq: A sequencing simulator'
    )
    parser.add_argument(
        '-v',
        '--version',
        action='store_true',
        default=False,
        help='print software version and exit'
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
    param_logging = parser_gen.add_mutually_exclusive_group()
    # --genomes and --ncbi should not be exclusive anymore
    # input_genomes = parser_gen.add_mutually_exclusive_group()
    input_abundance = parser_gen.add_mutually_exclusive_group()
    param_logging.add_argument(
        '--quiet',
        '-q',
        action='store_true',
        default=False,
        help='Disable info logging. (default: %(default)s).'
    )
    param_logging.add_argument(
        '--debug',
        '-d',
        action='store_true',
        default=False,
        help='Enable debug logging. (default: %(default)s).'
    )
    parser_gen.add_argument(
        '--seed',
        type=int,
        metavar='<int>',
        help='Seed all the random number generators',
        default=None
    )
    parser_gen.add_argument(
        '--cpus',
        '-p',
        default=2,
        type=int,
        metavar='<int>',
        help='number of cpus to use. (default: %(default)s).'
    )
    parser_gen.add_argument(
        '--genomes',
        '-g',
        metavar='<genomes.fasta>',
        nargs="+",
        help='Input genome(s) from where the reads will originate'
    )
    parser_gen.add_argument(
        '--draft',
        metavar='<draft.fasta>',
        nargs="+",
        help='Input draft genome(s) from where the reads will originate'
    )
    parser_gen.add_argument(
        '--n_genomes',
        '-u',
        type=int,
        metavar='<int>',
        help='How many genomes will be used for the simulation. is set with \
            --genomes/-g or/and --draft to take random genomes from the \
            input multifasta'
    )
    parser_gen.add_argument(
        '--ncbi',
        '-k',
        choices=['bacteria', 'viruses', 'archaea'],
        action='append',
        nargs='*',
        metavar='<str>',
        help='Download input genomes from NCBI. Requires --n_genomes/-u\
            option. Can be bacteria, viruses, archaea or a combination of the\
            three (space-separated)'
    )
    parser_gen.add_argument(
        '--n_genomes_ncbi',
        '-U',
        type=int,
        action='append',
        metavar='<int>',
        nargs='*',
        help='How many genomes will be downloaded from NCBI. Required if\
            --ncbi/-k is set. If more than one kingdom is set with --ncbi,\
            multiple values are necessary (space-separated).'
    )
    input_abundance.add_argument(
        '--abundance',
        '-a',
        choices=['uniform', 'halfnormal',
                 'exponential', 'lognormal', 'zero_inflated_lognormal'],
        metavar='<str>',
        default='lognormal',
        help='abundance distribution (default: %(default)s). Can be uniform,\
            halfnormal, exponential, lognormal or zero-inflated-lognormal.'
    )
    input_abundance.add_argument(
        '--abundance_file',
        '-b',
        metavar='<abundance.txt>',
        help='abundance file for coverage calculations (default: %(default)s).'
    )
    input_abundance.add_argument(
        '--coverage',
        '-C',
        choices=['uniform', 'halfnormal',
                 'exponential', 'lognormal', 'zero_inflated_lognormal'],
        metavar='<str>',
        help='coverage distribution. Can be uniform,\
            halfnormal, exponential, lognormal or zero-inflated-lognormal.'
    )
    input_abundance.add_argument(
        '--coverage_file',
        '-D',
        metavar='<coverage.txt>',
        help='file containing coverage information (default: %(default)s).'
    )
    parser_gen.add_argument(
        '--n_reads',
        '-n',
        metavar='<int>',
        default='1000000',
        help='Number of reads to generate (default: %(default)s). Allows \
        suffixes k, K, m, M, g and G (ex 0.5M for 500000).'
    )
    parser_gen.add_argument(
        '--mode',
        '-e',
        metavar='<str>',
        choices=['kde', 'basic', 'perfect'],
        default='kde',
        help='Error model. If not specified, using kernel density estimation \
        (default: %(default)s). Can be kde, basic or perfect'
    )
    parser_gen.add_argument(
        '--model',
        '-m',
        metavar='<npz>',
        default=None,
        help='Error model file. (default: %(default)s). Use HiSeq, NovaSeq or \
        MiSeq for a pre-computed error model provided with the software, or a \
        file generated with iss model. If you do not wish to use a model, use \
        --mode basic or --mode perfect. The name of the built-in models are  \
        case insensitive.'
    )
    parser_gen.add_argument(
        '--gc_bias',
        '-c',
        action='store_true',
        default=False,
        help='If set, may fail to sequence reads with abnormal GC content. \
        (default: %(default)s)'
    )
    parser_gen.add_argument(
        '--compress',
        '-z',
        action='store_true',
        default=False,
        help='Compress the output in gzip format (default: %(default)s).'
    )
    parser_gen.add_argument(
        '--output',
        '-o',
        metavar='<fastq>',
        help='Output file path and prefix (Required)',
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
        '--output',
        '-o',
        metavar='<npz>',
        help='Output file path and prefix (Required)',
        required=True
    )
    parser_mod._optionals.title = 'arguments'
    parser_mod.set_defaults(func=model_from_bam)
    args = parser.parse_args()

    # set logger and display version if args.version
    try:
        if args.version:
            print('iss version %s' % __version__)
            sys.exit(0)
        elif args.quiet:
            logging.basicConfig(level=logging.ERROR)
        elif args.debug:
            logging.basicConfig(level=logging.DEBUG)
        else:
            logging.basicConfig(level=logging.INFO)

        args.func(args)
        logging.shutdown()
    except AttributeError as e:
        logger = logging.getLogger(__name__)
        logger.debug(e)
        parser.print_help()
        # raise  # extra traceback to uncomment if all hell breaks lose
