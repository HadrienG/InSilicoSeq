#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, unicode_literals
from builtins import range

from iss.util import rev_comp

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from shutil import copyfileobj

import os
import sys
import random
import logging
import numpy as np


def reads(record, ErrorModel, n_pairs, cpu_number, output, seed,
          gc_bias=False):
    """Simulate reads from one genome (or sequence) according to an ErrorModel

    This function makes use of the `simulate_read` function to simulate reads
    and save them in a fastq file

    Args:
        record (SeqRecord): sequence or genome of reference
        ErrorModel (ErrorModel): an ErrorModel
        n_pairs (int): the number of reads to generate
        cpu_number (int): an int indentifying the cpu that is used by the
            function. Is used for naming the output file
        output (str): the output file prefix
        seed (int): random seed to use
        gc_bias (bool): if set, the function may skip a read due to abnormal
            GC content

    Returns:
        str: the name of the output file
    """
    logger = logging.getLogger(__name__)
    if seed is not None:
        random.seed(seed + cpu_number)
        np.random.seed(seed + cpu_number)
    logger.debug(
        'Cpu #%s: Generating %s read pairs'
        % (cpu_number, n_pairs))
    read_tuple_list = []
    i = 0
    while i < n_pairs:
        # try:
        #     forward, reverse = simulate_read(record, ErrorModel, i)
        # except ValueError as e:
        #     logger.error('Skipping this record: %s' % record.id)
        #     return
        forward, reverse = simulate_read(record, ErrorModel, i, cpu_number)
        if gc_bias:
            stiched_seq = forward.seq + reverse.seq
            gc_content = GC(stiched_seq)
            if 40 < gc_content < 60:
                read_tuple_list.append((forward, reverse))
                i += 1
            elif np.random.rand() < 0.90:
                read_tuple_list.append((forward, reverse))
                i += 1
            else:
                continue
        else:
            read_tuple_list.append((forward, reverse))
            i += 1

    temp_file_name = output + '.iss.tmp.%s.%s' % (record.id, cpu_number)
    to_fastq(read_tuple_list, temp_file_name)

    return temp_file_name


def simulate_read(record, ErrorModel, i, cpu_number):
    """From a read pair from one genome (or sequence) according to an
    ErrorModel

    Each read is a SeqRecord object
    returns a tuple containing the forward and reverse read.

    Args:
        record (SeqRecord): sequence or genome of reference
        ErrorModel (ErrorModel): an ErrorModel class
        i (int): a number identifying the read
        cpu_number (int): cpu number. Is added to the read id.

    Returns:
        tuple: tuple containg a forward read and a reverse read
    """
    logger = logging.getLogger(__name__)
    sequence = record.seq
    header = record.id

    read_length = ErrorModel.read_length
    insert_size = ErrorModel.random_insert_size()
    # generate the forward read
    try:  # a ref sequence has to be longer than 2 * read_length + i_size
        assert read_length < len(record.seq)
        forward_start = random.randrange(
            0, len(record.seq) - (2 * read_length + insert_size))
    except AssertionError as e:
        logger.error(
            '%s shorter than read length for this ErrorModel:%s'
            % (e, record.id))
        sys.exit(1)
    except ValueError as e:
        logger.debug(
            '%s shorter than template length for this ErrorModel:%s'
            % (record.id, e))
        forward_start = max(0, random.randrange(
            0, len(record.seq) - read_length))
        # raise

    forward_end = forward_start + read_length
    bounds = (forward_start, forward_end)
    # create a perfect read
    forward = SeqRecord(
        Seq(str(sequence[forward_start:forward_end]),
            IUPAC.unambiguous_dna
            ),
        id='%s_%s_%s/1' % (header, i, cpu_number),
        description=''
    )
    # add the indels, the qual scores and modify the record accordingly
    forward.seq = ErrorModel.introduce_indels(
        forward, 'forward', sequence, bounds)
    forward = ErrorModel.introduce_error_scores(forward, 'forward')
    forward.seq = ErrorModel.mut_sequence(forward, 'forward')

    # generate the reverse read
    try:
        reverse_start = forward_end + insert_size
        reverse_end = reverse_start + read_length
        assert reverse_end < len(record.seq)
    except AssertionError as e:
        # we use random insert when the modelled template length distribution
        # is too large
        reverse_end = random.randrange(read_length, len(record.seq))
        reverse_start = reverse_end - read_length
    bounds = (reverse_start, reverse_end)
    # create a perfect read
    reverse = SeqRecord(
        Seq(rev_comp(str(sequence[reverse_start:reverse_end])),
            IUPAC.unambiguous_dna
            ),
        id='%s_%s_%s/2' % (header, i, cpu_number),
        description=''
    )
    # add the indels, the qual scores and modify the record accordingly
    reverse.seq = ErrorModel.introduce_indels(
        reverse, 'reverse', sequence, bounds)
    reverse = ErrorModel.introduce_error_scores(reverse, 'reverse')
    reverse.seq = ErrorModel.mut_sequence(reverse, 'reverse')

    return (forward, reverse)


def to_fastq(generator, output):
    """Write reads to a fastq file

    Take a generator or a list containing read pairs (tuples) and write them
        in two fastq files: output_R1.fastq and output_R2.fastq

    Args:
        generator (generator): a read generator (or list)
        output (string): the output files prefix
    """
    logger = logging.getLogger(__name__)
    # define name of output files
    output_forward = output + '_R1.fastq'
    output_reverse = output + '_R2.fastq'

    try:
        f = open(output_forward, 'a')
        r = open(output_reverse, 'a')
    except PermissionError as e:
        logger.error('Failed to open output file(s): %s' % e)
        sys.exit(1)
    else:
        with f, r:
            for read_tuple in generator:
                SeqIO.write(read_tuple[0], f, 'fastq-sanger')
                SeqIO.write(read_tuple[1], r, 'fastq-sanger')
