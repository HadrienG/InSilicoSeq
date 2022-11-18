#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util
from iss import download
from iss import abundance
from iss import generator
from iss.version import __version__

from Bio import SeqIO
from joblib import Parallel, delayed

import gc
import os
import sys
import random
import logging
import argparse
import numpy as np
from typing import List, Optional

LOGGER = logging.getLogger(__name__)


class ReadGenerator:
    def __init__(self, args):
        LOGGER.debug('iss version %s' % __version__)
        LOGGER.debug('Using verbose logger')

        if args.seed:
            self.seed = args.seed
            LOGGER.info('Setting random seed to %i' % self.seed)
            random.seed(self.seed)
            np.random.seed(self.seed)
        else: 
            self.seed = None

        # initialize error model
        self.mode: str = args.mode
        self.model: str = args.model
        self.store_mutations: bool = args.store_mutations != "none"
        self.err_mod = self._load_model()

        # initialize genomes
        self.genomes: List[str] = args.genomes
        self.draft: List[str] = args.draft
        self.ncbi: str = args.ncbi
        self.output: str = args.output
        self.n_genomes_ncbi: List[int] = args.n_genomes_ncbi
        self.n_genomes: int = args.n_genomes
        self.genome_file = self._load_genome()

        # initialize abundance / coverage
        self.abundance_dispatch = {
            'uniform': abundance.uniform,
            'halfnormal': abundance.halfnormal,
            'exponential': abundance.exponential,
            'lognormal': abundance.lognormal,
            'zero_inflated_lognormal': abundance.zero_inflated_lognormal
        }
        self.abundance_file: str = args.abundance_file
        self.coverage_file: str = args.coverage_file
        self.readcount_file: str = args.readcount_file
        self.coverage: str = args.coverage
        self.abundance: str = args.abundance
        self.n_reads: int = args.n_reads

        self.abundance_dic = self._load_abundance_file()
        self.readcount_dic = self._load_readcount_file()
        # print(*self.readcount_dic.items(), sep="\n")

        self.cpus: int = args.cpus

        # initialize generation variables
        self.sequence_type: str = args.sequence_type
        self.gc_bias: bool = args.gc_bias
        self.compress: bool = args.compress
        self.parallel_genomes = self.sequence_type == "amplicon"
        self.total_reads_generated = 0
        self.total_reads_generated_unrounded = 0
        self.mutations_format: str = args.store_mutations  # VCF
        if self.store_mutations:
            LOGGER.info(f'Storing mutations to {self.mutations_format} format')
        self.mutations = []


    def _load_model(self):
        try:  # try to import and load the correct error model
            LOGGER.info('Starting iss generate')
            LOGGER.info('Using %s ErrorModel' % self.mode)
            
            if self.mode == 'kde':
                from iss.error_models import kde
                if self.model is None:
                    LOGGER.error('--model is required in --mode kde')
                    sys.exit(1)
                elif self.model.lower() == 'hiseq':
                    npz = os.path.join(
                        os.path.dirname(__file__),
                        'profiles/HiSeq')
                elif self.model.lower() == 'novaseq':
                    npz = os.path.join(
                        os.path.dirname(__file__),
                        'profiles/NovaSeq')
                elif self.model.lower() == 'miseq':
                    npz = os.path.join(
                        os.path.dirname(__file__),
                        'profiles/MiSeq')
                else:
                    npz = self.model
                err_mod = kde.KDErrorModel(npz, self.store_mutations)
            elif self.mode == 'basic':
                if self.model is not None:
                    LOGGER.warning('--model %s will be ignored in --mode %s' % (
                        self.model, self.mode))
                from iss.error_models import basic
                err_mod = basic.BasicErrorModel(self.store_mutations)
            elif self.mode == 'perfect':
                if self.model is not None:
                    LOGGER.warning('--model %s will be ignored in --mode %s' % (
                        self.model, self.mode))
                from iss.error_models import perfect
                err_mod = perfect.PerfectErrorModel()
            else:
                raise RuntimeError(f"mode '{self.mode}' is not supported")
        except ImportError as e:
            LOGGER.error('Failed to import ErrorModel module: %s' % e)
            sys.exit(1)
        return err_mod

    def _load_genome(self):
        try:  # try to read genomes and concatenate --genomes and --ncbi genomes
            if self.genomes or self.draft or self.ncbi:
                genome_files = []
                if self.genomes:
                    genome_files.extend(self.genomes)
                if self.draft:
                    genome_files.extend(self.draft)
                if self.ncbi and self.n_genomes_ncbi:
                    util.genome_file_exists(self.output + '_ncbi_genomes.fasta')
                    try:
                        assert len(*self.ncbi) == len(*self.n_genomes_ncbi)
                    except AssertionError as e:
                        LOGGER.error(
                            '--ncbi and --n_genomes_ncbi of unequal lengths. \
                            Aborting')
                        sys.exit(1)
                    for g, n in zip(*self.ncbi, *self.n_genomes_ncbi):
                        self.genomes_ncbi = download.ncbi(
                            g, n, self.output + '_ncbi_genomes.fasta')
                    genome_files.append(self.genomes_ncbi)
                if self.ncbi and not self.n_genomes_ncbi:
                    LOGGER.error(
                        '--ncbi/-k requires --n_genomes_ncbi/-U. Aborting.')
                    sys.exit(1)

            else:
                LOGGER.error("One of --genomes/-g, --draft, --ncbi/-k is required")
                sys.exit(1)

            genome_file = self.output + '.iss.tmp.genomes.fasta'
            util.concatenate(
                genome_files,
                output=genome_file)

            # for n_genomes we use reservoir sampling to draw random genomes
            # from the concatenated genome file. We then override the file.
            if self.n_genomes and not self.draft and not self.ncbi:
                genome_count = util.count_records(genome_file)
                genome_files = [genome for genome in util.reservoir(
                    SeqIO.parse(genome_file, 'fasta'),
                    genome_count,
                    self.n_genomes)]
                SeqIO.write(genome_files, genome_file, 'fasta')

            assert os.stat(genome_file).st_size != 0
            f = open(genome_file, 'r')
            with f:  # count the number of records
                self.genome_list = util.count_records(f)
        except IOError as e:
            LOGGER.error('Failed to open genome(s) file:%s' % e)
            sys.exit(1)
        except AssertionError as e:
            LOGGER.error('Genome(s) file seems empty: %s' % genome_file)
            sys.exit(1)
        except KeyboardInterrupt as e:
            LOGGER.error('iss generate interrupted: %s' % e)
            sys.exit(1)
        return genome_file

    def _load_abundance_file(self):
        # read the abundance, coverage file
        if self.abundance_file:
            LOGGER.info('Using abundance file:%s' % self.abundance_file)
            if self.draft:
                abundance_dic_short = abundance.parse_abundance_file(
                    self.abundance_file)
                complete_genomes_dic = {k: v for
                                        k, v in abundance_dic_short.items()
                                        if k not in self.draft}
                draft_dic = abundance.expand_draft_abundance(
                    abundance_dic_short,
                    self.draft)
                abundance_dic = {**complete_genomes_dic,
                                    **draft_dic}
            else:
                abundance_dic = abundance.parse_abundance_file(
                    self.abundance_file)
        elif self.coverage_file:
            LOGGER.warning('--coverage_file is an experimental feature')
            LOGGER.warning('--coverage_file disables --n_reads')
            LOGGER.info('Using coverage file:%s' % self.coverage_file)
            if self.draft:
                coverage_dic = abundance.parse_abundance_file(
                    self.coverage_file)
                complete_genomes_dic = {k: v for
                                        k, v in coverage_dic.items()
                                        if k not in self.draft}
                draft_dic = abundance.expand_draft_abundance(
                    abundance_dic_short,
                    self.draft,
                    mode="coverage")
                abundance_dic = {**complete_genomes_dic,
                                    **draft_dic}
            else:
                abundance_dic = abundance.parse_abundance_file(
                    self.coverage_file)
        elif self.coverage in self.abundance_dispatch:
            # todo coverage distribution with --draft
            LOGGER.warning('--coverage is an experimental feature')
            LOGGER.info('Using %s coverage distribution' % self.coverage)
            if self.draft:
                abundance_dic = abundance.draft(
                    self.genome_list,
                    self.draft,
                    self.abundance_dispatch[self.abundance],
                    self.output,
                    mode="coverage")
            else:
                abundance_dic = self.abundance_dispatch[
                    self.coverage](self.genome_list)
            if self.n_reads:
                n_reads = util.convert_n_reads(self.n_reads)
                LOGGER.info('scaling coverage to %s reads' % n_reads)
                abundance_dic = abundance.coverage_scaling(n_reads,
                                                            abundance_dic,
                                                            self.genome_file,
                                                            self.err_mod.read_length)
            abundance.to_file(abundance_dic, self.output, mode="coverage")
        elif self.abundance in self.abundance_dispatch:
            LOGGER.info('Using %s abundance distribution' % self.abundance)
            if self.draft:
                abundance_dic = abundance.draft(
                    self.genome_list,
                    self.draft,
                    self.abundance_dispatch[self.abundance],
                    self.output)
            else:
                abundance_dic = self.abundance_dispatch[
                    self.abundance](self.genome_list)
                abundance.to_file(abundance_dic, self.output)
        else:
            LOGGER.error('Could not get abundance')
            sys.exit(1)
        return abundance_dic

    def _load_readcount_file(self):
        if self.readcount_file:
            LOGGER.warning('--readcount_file is an experimental feature')
            LOGGER.warning('--readcount_file disables --n_reads')
            if self.draft:
                raise RuntimeError("readcount_file is only supported using --genomes, not --draft")
            readcount_dic =  abundance.parse_readcount_file(self.readcount_file)
            self.n_reads = sum(readcount_dic.values())
            LOGGER.info(f'Will generate a total of {self.n_reads} read pairs based on the readcount_file')
            return readcount_dic
        return None

    def _get_n_pairs_for_record(self, record: SeqIO.SeqRecord):
        if self.readcount_dic:
            n_pairs = self.readcount_dic.get(record.id, 0)
            self.total_reads_generated += n_pairs
        else:
            try:
                species_abundance = self.abundance_dic[record.id]
            except KeyError as e:
                LOGGER.error(
                    'Fasta record not found in abundance file: %s' % e)
                sys.exit(1)

            genome_size = len(record.seq)
            if self.coverage or self.coverage_file:
                coverage = species_abundance
            else:
                coverage = abundance.to_coverage(
                    self.n_reads,
                    species_abundance,
                    self.err_mod.read_length,
                    genome_size
                )
            n_pairs_unrounded = (
                (coverage * len(record.seq)) /
                self.err_mod.read_length) / 2
            n_pairs = round(n_pairs_unrounded)

            # check that the rounding does not cause to drop
            self.total_reads_generated_unrounded += n_pairs_unrounded
            self.total_reads_generated += n_pairs
            if round(self.total_reads_generated_unrounded) > self.total_reads_generated:
                LOGGER.debug(
                    "Adding a pair to correct rounding error")
                n_pairs += 1
                self.total_reads_generated += 1
        return n_pairs


    def generate_for_record(
            self,
            record: SeqIO.SeqRecord,
            n_pairs: int
        ) -> Optional[str]:
        if n_pairs == 0:
            LOGGER.debug(
                f"Skipping record {record.id} because n_pairs == 0")
            return None

        # due to a bug in multiprocessing
        # https://bugs.python.org/issue17560
        # we can't send records taking more than 2**31 bytes
        # through serialisation.
        # In those cases we use memmapping
        if sys.getsizeof(str(record.seq)) >= 2**31 - 1:
                LOGGER.warning(
                    "record %s unusually big." % record.id)
                LOGGER.warning("Using a memory map.")
                mode = "memmap"

                record_mmap = "%s.memmap" % self.output
                if os.path.exists(record_mmap):
                    os.unlink(record_mmap)
                util.dump(record, record_mmap)
                del record
                record = record_mmap
                gc.collect()
        else:
            mode = "default"

        LOGGER.debug(f'Generating {n_pairs} read pair(s) for record: {record.id}')

        read_file = generator.reads(
            record, self.err_mod,
            n_pairs, 0, self.output,
            self.seed, self.sequence_type,
            self.gc_bias, mode, self.store_mutations
        )
        return read_file

    def generate_for_record_multi(
            self,
            record: SeqIO.SeqRecord,
            n_pairs: int
        ) -> List[Optional[str]]:
        # n_pairs = self._get_n_pairs_for_record(record)
        # skip record if n_pairs == 0
        if n_pairs == 0:
            LOGGER.debug(
                f"Skipping record {record.id} because n_pairs == 0")
            return [None]

        # due to a bug in multiprocessing
        # https://bugs.python.org/issue17560
        # we can't send records taking more than 2**31 bytes
        # through serialisation.
        # In those cases we use memmapping
        if sys.getsizeof(str(record.seq)) >= 2**31 - 1:
                LOGGER.warning(
                    "record %s unusually big." % record.id)
                LOGGER.warning("Using a memory map.")
                mode = "memmap"

                record_mmap = "%s.memmap" % self.output
                if os.path.exists(record_mmap):
                    os.unlink(record_mmap)
                util.dump(record, record_mmap)
                del record
                record = record_mmap
                gc.collect()
        else:
            mode = "default"

        LOGGER.debug(f'Generating {n_pairs} read pair(s) for record: {record.id}')

        # exact self.n_reads for each cpus
        if n_pairs % self.cpus == 0:
            n_pairs_per_cpu = [(n_pairs // self.cpus)
                                for _ in range(self.cpus)]
        else:
            n_pairs_per_cpu = [(n_pairs // self.cpus)
                                for _ in range(self.cpus)]
            n_pairs_per_cpu[-1] += n_pairs % self.cpus

        record_file_name_list = Parallel(
            n_jobs=self.cpus)(
                delayed(generator.reads)(
                    record, self.err_mod,
                    n_pairs_per_cpu[i], i, self.output,
                    self.seed, self.sequence_type,
                    self.gc_bias, mode) for i in range(self.cpus)
                )
        assert record_file_name_list is not None
        return record_file_name_list


    def generate(self):
        if not (self.coverage or self.coverage_file or self.readcount_file):
            self.n_reads = util.convert_n_reads(self.n_reads)
            LOGGER.info('Generating %s reads' % self.n_reads)
        gen_func = self.generate_for_record if self.parallel_genomes else self.generate_for_record_multi

        try:
            f = open(self.genome_file, 'r')  # re-opens the file
            with f:
                fasta_file = SeqIO.parse(f, 'fasta')
                tasks = []
                for record in fasta_file:
                    n_pairs = self._get_n_pairs_for_record(record)
                    if n_pairs == 0:
                        LOGGER.debug(
                            f"Skipping record {record.id} because n_pairs == 0")
                        continue
                    tasks.append((gen_func, record, n_pairs))
                temp_file_list = Parallel(
                    n_jobs=self.cpus, verbose=2, batch_size=1000
                    )(delayed(func)(record, n_pairs)
                        for (func, record, n_pairs) in tasks
                )

        except KeyboardInterrupt as e:
            LOGGER.error('iss generate interrupted: %s' % e)
            temp_file_list = list(filter(None, temp_file_list))
            temp_R1 = [temp_file + '_R1.fastq' for temp_file in temp_file_list]
            temp_R2 = [temp_file + '_R2.fastq' for temp_file in temp_file_list]
            temp_mut = [temp_file + '.' + self.mutations_format for temp_file in temp_file_list]
            full_tmp_list = temp_R1 + temp_R2 + temp_mut
            full_tmp_list.append(self.genome_file)
            if os.path.exists("%s.memmap" % self.output):
                full_tmp_list.append("%s.memmap" % self.output)
            util.cleanup(full_tmp_list)
            sys.exit(1)
        else:
            # remove the duplicates in file list and cleanup
            # we remove the duplicates in case two records had the same header
            # and reads were appended to the same temp file.
            temp_file_list = list(filter(None, temp_file_list))

            temp_R1 = [temp_file + '_R1.fastq' for temp_file in temp_file_list]
            temp_R2 = [temp_file + '_R2.fastq' for temp_file in temp_file_list]
            temp_mut = [temp_file + '.' + self.mutations_format for temp_file in temp_file_list] if self.store_mutations else []
            util.concatenate(temp_R1, self.output + '_R1.fastq')
            util.concatenate(temp_R2, self.output + '_R2.fastq')
            if self.store_mutations:
                util.concatenate(
                    temp_mut,
                    self.output + '.' + self.mutations_format,
                    "##fileformat=VCFv4.1\n" + "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"])
                )
            full_tmp_list = temp_R1 + temp_R2  + temp_mut
            full_tmp_list.append(self.genome_file)
            if os.path.exists("%s.memmap" % self.output):
                full_tmp_list.append("%s.memmap" % self.output)
            util.cleanup(full_tmp_list)
            if self.compress:
                util.compress(self.output + '_R1.fastq')
                util.compress(self.output + '_R2.fastq')
                if self.store_mutations:
                    util.compress(self.output + '.' + self.mutations_format)
            LOGGER.info(f'Read generation complete, {self.total_reads_generated} reads generated.')



def generate_reads(args):
    read_generator = ReadGenerator(args)
    read_generator.generate()


def model_from_bam(args):
    """Main function for the `iss model` submodule

    This submodule write all variables necessary for building an ErrorModel
    to args.output + .npz

    Args:
        args (object): the command-line arguments from argparse
    """
    LOGGER.debug('iss version %s' % __version__)
    LOGGER.debug('Using verbose logger')

    try:  # try to import bam module and write model data to file
        LOGGER.info('Starting iss model')
        from iss import bam
    except ImportError as e:
        LOGGER.error('Failed to import bam module: %s' % e)
        sys.exit(1)
    else:
        LOGGER.info('Using KDE ErrorModel')
        bam.to_model(args.bam, args.output)
        LOGGER.info('Model generation complete')


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
    input_abundance.add_argument(
        "--readcount_file",
        "-R",
        metavar='<readcount.txt>',
        help='file containing read_count information (default: %(default)s).'
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
    parser_gen.add_argument(
        '--sequence_type',
        '-t',
        choices=['metagenomics', 'amplicon'],
        required=True,
        help='Type of sequencing. Can be metagenomics or amplicon.'
    )
    parser_gen.add_argument(
        '--store_mutations',
        '-M',
        choices=['none', 'vcf'],
        default="none",
        help='Enables output of inserted mutations in format (default: %(choises)s).'
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
