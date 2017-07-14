#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division, unicode_literals
from builtins import range

from iss.util import rev_comp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import sys
import random
import logging


def reads(record, coverage, ErrorModel):
    """Simulate reads from one genome (or sequence) according to an ErrorModel

    Each read is a SeqRecord object
    Return a generator of tuples containing the forward and
    reverse read.

    Args:
        record (SeqRecord): sequence or genome of reference
        coverage (float): desired coverage of the genome
        ErrorModel (ErrorModel): an ErrorModel class

    Yields:
        tuple: tuple containg a forward read and a reverse
            read
    """
    logger = logging.getLogger(__name__)
    header = record.id
    sequence = record.seq

    read_length = ErrorModel.read_length
    insert_size = ErrorModel.insert_size

    n_pairs = int(round((coverage * len(sequence)) / read_length) / 2)

    for i in range(n_pairs):
        # generate the forward read
        try:  # a ref sequence has to be longer than 2 * read_length + i_size
            forward_start = random.randrange(
                0, len(sequence) - (2 * read_length + insert_size))
        except ValueError as e:
            logger.error(
                '%s too small for this ErrorModel:%s' % (record.id, e))
            sys.exit(1)

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

        yield(forward, reverse)


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
