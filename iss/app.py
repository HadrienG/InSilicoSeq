#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import multiprocessing as mp
import os
import sys

from iss import util
from iss.generator import (
    generate_work_divider,
    load_error_model,
    load_genomes,
    load_readcount_or_abundance,
    worker_iterator,
)
from iss.version import __version__


def generate_reads(args):
    """Main function for the `iss generate` submodule

    This submodule generates reads from an ErrorModel and write them to
        args.output + _R(1|2).fastq

    Args:
        args (object): the command-line arguments from argparse
    """
    logger = logging.getLogger(__name__)
    logger.debug("iss version %s" % __version__)
    logger.debug("Using verbose logger")
    logger.info("Starting iss generate")

    error_model = load_error_model(
        args.mode, args.seed, args.model, args.fragment_length, args.fragment_length_sd, args.store_mutations
    )

    genome_list, genome_file = load_genomes(
        args.genomes, args.draft, args.ncbi, args.n_genomes_ncbi, args.output, args.n_genomes
    )

    readcount_dic, abundance_dic = load_readcount_or_abundance(
        args.readcount_file,
        args.abundance_file,
        args.coverage_file,
        args.coverage,
        args.abundance,
        args.draft,
        genome_list,
        genome_file,
        args.n_reads,
        args.output,
        error_model,
    )

    if args.store_mutations:
        logger.info(f"Storing inserted sequence errors in {args.output}.vcf")

    logger.info("Using %s cpus for read generation" % args.cpus)

    if readcount_dic is not None:
        n_reads = sum(readcount_dic.values())
        logger.info("Generating %s reads" % n_reads)
    else:
        n_reads = util.convert_n_reads(args.n_reads)
        logger.info("Generating %s reads" % n_reads)

    try:
        # list holding the prefix for each cpu's temp file
        temp_file_list = [f"{args.output}.iss.tmp.{i}" for i in range(args.cpus)]

        # Calculate how many reads we want each cpu to generate
        n_read_pairs = n_reads // 2
        chunk_size = -((n_read_pairs) // -args.cpus)  # this is ceildiv, see https://stackoverflow.com/a/17511341
        logger.debug("Chunk size: %s" % chunk_size)

        # Divide the work of generating n_reads for each record into chunks
        work_chunks = generate_work_divider(
            genome_file,
            readcount_dic,
            abundance_dic,
            n_reads,
            args.coverage,
            args.coverage_file,
            error_model,
            chunk_size,
        )

        # Generate reads for each chunk in parallel
        with mp.Pool(args.cpus) as pool:
            pool.starmap(
                worker_iterator,
                [
                    (
                        work,
                        error_model,
                        cpu_number,
                        worker_prefix,
                        args.seed,
                        args.sequence_type,
                        args.gc_bias,
                        genome_file,
                    )
                    for cpu_number, (work, worker_prefix) in enumerate(zip(work_chunks, temp_file_list))
                ],
            )

    except KeyboardInterrupt as e:
        logger.error("iss generate interrupted: %s" % e)
        temp_R1 = [temp_file + "_R1.fastq" for temp_file in temp_file_list]
        temp_R2 = [temp_file + "_R2.fastq" for temp_file in temp_file_list]
        temp_mut = [temp_file + ".vcf" for temp_file in temp_file_list]
        full_tmp_list = temp_R1 + temp_R2 + temp_mut
        full_tmp_list.append(genome_file)
        if os.path.exists("%s.memmap" % args.output):
            full_tmp_list.append("%s.memmap" % args.output)
        util.cleanup(full_tmp_list)
        sys.exit(1)
    else:
        # remove the duplicates in file list and cleanup
        # we remove the duplicates in case two records had the same header
        # and reads were appended to the same temp file.
        temp_R1 = [temp_file + "_R1.fastq" for temp_file in temp_file_list]
        temp_R2 = [temp_file + "_R2.fastq" for temp_file in temp_file_list]
        temp_mut = [temp_file + ".vcf" for temp_file in temp_file_list]
        util.concatenate(temp_R1, args.output + "_R1.fastq")
        util.concatenate(temp_R2, args.output + "_R2.fastq")
        if args.store_mutations:
            util.concatenate(
                temp_mut,
                args.output + ".vcf",
                "##fileformat=VCFv4.1\n" + "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]),
            )
        full_tmp_list = temp_R1 + temp_R2 + temp_mut
        full_tmp_list.append(genome_file)
        if os.path.exists("%s.memmap" % args.output):
            full_tmp_list.append("%s.memmap" % args.output)
        util.cleanup(full_tmp_list)
        if args.compress:
            util.compress(args.output + "_R1.fastq")
            util.compress(args.output + "_R2.fastq")
            if args.store_mutations:
                util.compress(args.output + ".vcf")
        logger.info("Read generation complete")


def model_from_bam(args):
    """Main function for the `iss model` submodule

    This submodule write all variables necessary for building an ErrorModel
    to args.output + .npz

    Args:
        args (object): the command-line arguments from argparse
    """
    logger = logging.getLogger(__name__)
    logger.debug("iss version %s" % __version__)
    logger.debug("Using verbose logger")

    try:  # try to import bam module and write model data to file
        logger.info("Starting iss model")
        from iss import bam
    except ImportError as e:
        logger.error("Failed to import bam module: %s" % e)
        sys.exit(1)
    else:
        logger.info("Using KDE ErrorModel")
        bam.to_model(args.bam, args.output)
        logger.info("Model generation complete")


def main():
    parser = argparse.ArgumentParser(
        prog="InSilicoSeq",
        usage="iss [subcommand] [options]",
        description="InSilicoSeq: A sequencing simulator",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="store_true",
        default=False,
        help="print software version and exit",
    )
    subparsers = parser.add_subparsers(title="available subcommands", metavar="")

    parser_mod = subparsers.add_parser(
        "model",
        prog="iss model",
        description="generate an error model from a bam file",
        help="generate an error model from a bam file",
    )
    parser_gen = subparsers.add_parser(
        "generate",
        prog="iss generate",
        description="simulate reads from an error model",
        help="simulate reads from an error model",
    )

    # arguments form the read generator module
    param_logging = parser_gen.add_mutually_exclusive_group()
    # --genomes and --ncbi should not be exclusive anymore
    # input_genomes = parser_gen.add_mutually_exclusive_group()
    input_abundance = parser_gen.add_mutually_exclusive_group()
    param_logging.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        default=False,
        help="Disable info logging. (default: %(default)s).",
    )
    param_logging.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        help="Enable debug logging. (default: %(default)s).",
    )
    parser_gen.add_argument(
        "--seed",
        type=int,
        metavar="<int>",
        help="Seed all the random number generators",
        default=None,
    )
    parser_gen.add_argument(
        "--cpus",
        "-p",
        default=2,
        type=int,
        metavar="<int>",
        help="number of cpus to use. (default: %(default)s).",
    )
    parser_gen.add_argument(
        "--genomes",
        "-g",
        metavar="<genomes.fasta>",
        nargs="+",
        help="Input genome(s) from where the reads will originate",
    )
    parser_gen.add_argument(
        "--draft",
        metavar="<draft.fasta>",
        nargs="+",
        help="Input draft genome(s) from where the reads will originate",
    )
    parser_gen.add_argument(
        "--n_genomes",
        "-u",
        type=int,
        metavar="<int>",
        help="How many genomes will be used for the simulation. is set with \
            --genomes/-g or/and --draft to take random genomes from the \
            input multifasta",
    )
    parser_gen.add_argument(
        "--ncbi",
        "-k",
        choices=["bacteria", "viruses", "archaea"],
        action="append",
        nargs="*",
        metavar="<str>",
        help="Download input genomes from NCBI. Requires --n_genomes/-u\
            option. Can be bacteria, viruses, archaea or a combination of the\
            three (space-separated)",
    )
    parser_gen.add_argument(
        "--n_genomes_ncbi",
        "-U",
        type=int,
        action="append",
        metavar="<int>",
        nargs="*",
        help="How many genomes will be downloaded from NCBI. Required if\
            --ncbi/-k is set. If more than one kingdom is set with --ncbi,\
            multiple values are necessary (space-separated).",
    )
    input_abundance.add_argument(
        "--abundance",
        "-a",
        choices=[
            "uniform",
            "halfnormal",
            "exponential",
            "lognormal",
            "zero_inflated_lognormal",
        ],
        metavar="<str>",
        default="lognormal",
        help="abundance distribution (default: %(default)s). Can be uniform,\
            halfnormal, exponential, lognormal or zero-inflated-lognormal.",
    )
    input_abundance.add_argument(
        "--abundance_file",
        "-b",
        metavar="<abundance.txt>",
        help="abundance file for coverage calculations (default: %(default)s).",
    )
    input_abundance.add_argument(
        "--coverage",
        "-C",
        choices=[
            "uniform",
            "halfnormal",
            "exponential",
            "lognormal",
            "zero_inflated_lognormal",
        ],
        metavar="<str>",
        help="coverage distribution. Can be uniform,\
            halfnormal, exponential, lognormal or zero-inflated-lognormal.",
    )
    input_abundance.add_argument(
        "--coverage_file",
        "-D",
        metavar="<coverage.txt>",
        help="file containing coverage information (default: %(default)s).",
    )
    input_abundance.add_argument(
        "--readcount_file",
        "-R",
        metavar="<readcount.txt>",
        help="file containing read_count information (default: %(default)s).",
    )
    parser_gen.add_argument(
        "--n_reads",
        "-n",
        metavar="<int>",
        default="1000000",
        help="Number of reads to generate (default: %(default)s). Allows \
        suffixes k, K, m, M, g and G (ex 0.5M for 500000).",
    )
    parser_gen.add_argument(
        "--mode",
        "-e",
        metavar="<str>",
        choices=["kde", "basic", "perfect"],
        default="kde",
        help="Error model. If not specified, using kernel density estimation \
        (default: %(default)s). Can be kde, basic or perfect",
    )
    parser_gen.add_argument(
        "--model",
        "-m",
        metavar="<npz>",
        default=None,
        help="Error model file. (default: %(default)s). Use HiSeq, NextSeq, NovaSeq, \
        MiSeq or Miseq-[20,24,28,32] for a pre-computed error model provided with the \
        software, or a file generated with iss model. If you do not wish to use a model, use \
        --mode basic or --mode perfect. The name of the built-in models are  \
        case insensitive.",
    )
    parser_gen.add_argument(
        "--gc_bias",
        "-c",
        action="store_true",
        default=False,
        help="If set, may fail to sequence reads with abnormal GC content. \
        (default: %(default)s)",
    )
    parser_gen.add_argument(
        "--compress",
        "-z",
        action="store_true",
        default=False,
        help="Compress the output in gzip format (default: %(default)s).",
    )
    parser_gen.add_argument(
        "--output",
        "-o",
        metavar="<fastq>",
        help="Output file path and prefix (Required)",
        required=True,
    )
    parser_gen.add_argument(
        "--sequence_type",
        "-t",
        choices=["metagenomics", "amplicon"],
        default="metagenomics",
        required=False,
        help="Type of sequencing. Can be metagenomics or amplicon (default: %(default)s).",
    )
    parser_gen.add_argument(
        "--fragment-length",
        "-l",
        metavar="<int>",
        required=False,
        type=int,
        help="Fragment length for metagenomics sequencing (default: %(default)s).",
    )
    parser_gen.add_argument(
        "--fragment-length-sd",
        "-s",
        metavar="<int>",
        required=False,
        type=int,
        help="Fragment length standard deviation for metagenomics sequencing (default: %(default)s).",
    )
    parser_gen.add_argument(
        "--store_mutations",
        "-M",
        action="store_true",
        default=False,
        help="Generates an additional VCF file with the mutations introduced in the reads",
    )
    parser_gen._optionals.title = "arguments"
    parser_gen.set_defaults(func=generate_reads)

    # arguments for the error model module
    parser_mod.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        default=False,
        help="Disable info logging. (default: %(default)s).",
    )
    parser_mod.add_argument(
        "--debug",
        "-d",
        action="store_true",
        default=False,
        help="Enable debug logging. (default: %(default)s).",
    )
    parser_mod.add_argument(
        "--bam",
        "-b",
        metavar="<bam>",
        help="aligned reads from which the model will be inferred (Required)",
        required=True,
    )
    parser_mod.add_argument(
        "--output",
        "-o",
        metavar="<npz>",
        help="Output file path and prefix (Required)",
        required=True,
    )
    parser_mod._optionals.title = "arguments"
    parser_mod.set_defaults(func=model_from_bam)
    args = parser.parse_args()

    # set logger and display version if args.version
    try:
        if args.version:
            print("iss version %s" % __version__)
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
