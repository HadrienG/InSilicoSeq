#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, unicode_literals
from builtins import range

from iss.util import rev_comp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqUtils import GC
from shutil import copyfileobj

import os
import sys
import random
import logging
import numpy as np


def reads(record, ErrorModel, n_pairs, cpu_number, gc_bias=False):
    """Simulate reads from one genome (or sequence) according to an ErrorModel

    Each read is a SeqRecord object
    Return a generator of tuples containing the forward and
    reverse read.

    Args:
        record (SeqRecord): sequence or genome of reference
        coverage (float): desired coverage of the genome
        ErrorModel (ErrorModel): an ErrorModel class
        gc_bias (bool): if set, the function may skip a read due to abnormal
            GC content

    Yields:
        tuple: tuple containg a forward read and a reverse
            read
    """
    logger = logging.getLogger(__name__)
    logger.debug(
        'Cpu #%s: Generating %s read pairs'
        % (cpu_number, n_pairs))
    read_tuple_list = []
    for i in range(n_pairs):
        try:
            forward, reverse = simulate_read(record, ErrorModel, i)
        except ValueError as e:
            logger.error('Skipping this record: %s' % record.id)
            return
        # forward, reverse = simulate_read(record, ErrorModel, i)
        if gc_bias:
            stiched_seq = forward.seq + reverse.seq
            gc_content = GC(stiched_seq)
            if 40 < gc_content < 60:
                read_tuple_list.append((forward, reverse))
            elif np.random.rand() < 0.90:
                read_tuple_list.append((forward, reverse))
            else:
                continue
        else:
            read_tuple_list.append((forward, reverse))

    temp_file_name = '.iss.tmp.%s.%s' % (record.id, cpu_number)
    to_fastq(read_tuple_list, temp_file_name)

    return temp_file_name


def simulate_read(record, ErrorModel, i):
    """From a sequence record and an ErrorModel, generate a read pair

    EXPERIMENTAL. SHOULD BE MULTI-THREADABLE
    """
    logger = logging.getLogger(__name__)
    sequence = record.seq
    header = record.id

    read_length = ErrorModel.read_length
    insert_size = ErrorModel.random_insert_size()
    # generate the forward read
    try:  # a ref sequence has to be longer than 2 * read_length + i_size
        forward_start = random.randrange(
            0, len(record.seq) - (2 * read_length + insert_size))
    except ValueError as e:
        logger.error(
            '%s shorter than template length for this ErrorModel:%s'
            % (record.id, e))
        raise

    forward_end = forward_start + read_length
    bounds = (forward_start, forward_end)
    # create a perfect read
    forward = SeqRecord(
        Seq(str(sequence[forward_start:forward_end]),
            IUPAC.unambiguous_dna
            ),
        id='%s_%s_1' % (header, i),
        description=''
    )
    # add the indels, the qual scores and modify the record accordingly
    forward.seq = ErrorModel.introduce_indels(
        forward, 'forward', sequence, bounds)
    forward = ErrorModel.introduce_error_scores(forward, 'forward')
    forward.seq = ErrorModel.mut_sequence(forward, 'forward')

    # generate the reverse read
    reverse_start = forward_end + insert_size
    reverse_end = reverse_start + read_length
    bounds = (reverse_start, reverse_end)
    # create a perfect read
    reverse = SeqRecord(
        Seq(rev_comp(str(sequence[reverse_start:reverse_end])),
            IUPAC.unambiguous_dna
            ),
        id='%s_%s_2' % (header, i),
        description=''
    )
    # add the indels, the qual scores and modify the record accordingly
    reverse.seq = ErrorModel.introduce_indels(
        reverse, 'reverse', sequence, bounds)
    reverse = ErrorModel.introduce_error_scores(reverse, 'reverse')
    reverse.seq = ErrorModel.mut_sequence(reverse, 'reverse')

    return (forward, reverse)


def to_fastq(generator, output):
    """Write reads to fastq

    Take the read generator and write read pairs in two fastq files:
    output_R1.fastq and output_R2.fastq

    Args:
        generator (generator): the read generator
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


def concatenate(file_list, output):
    """Concatenate fastq files together

    Outputs two files: output_R1.fastq and output_R2.fastq

    Args:
        file_list (list): the list of input files prefix
        output (string): the output files prefix
    """
    logger = logging.getLogger(__name__)
    logger.info('Stitching temporary files together')
    # define name of output files
    output_forward = output + '_R1.fastq'
    output_reverse = output + '_R2.fastq'
    try:
        out_f = open(output_forward, 'wb')
        out_r = open(output_reverse, 'wb')
    except PermissionError as e:
        logger.error('Failed to open output file(s): %s' % e)
        sys.exit(1)

    with out_f, out_r:
        for file_name in file_list:
            if file_name is not None:
                temp_f = file_name + '_R1.fastq'
                temp_r = file_name + '_R2.fastq'
                with open(temp_f, 'rb') as f, open(temp_r, 'rb') as r:
                    copyfileobj(f, out_f)
                    copyfileobj(r, out_r)


def cleanup(file_list):
    """remove temporary files

    Args:
        file_list (list): a list of files to be removed
    """
    logger = logging.getLogger(__name__)
    logger.info('Cleaning up')
    for temp_file in file_list:
        if temp_file is not None:
            try:
                os.remove(temp_file + '_R1.fastq')
                os.remove(temp_file + '_R2.fastq')
            except FileNotFoundError as e:
                logger.error('Temporary file not found: %s' % temp_file)
                logger.error('You may have to remove temporary files manually')
                sys.exit(1)
