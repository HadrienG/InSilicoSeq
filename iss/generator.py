#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import error_model
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import random


def reads(record, coverage, ErrorModel):
    """Simulate reads from one genome (or sequence). Each read is a SeqRecord
    object. Return a generator of tuples containing the forward and
    reverse read.

    Arguments:
    record -- sequence of reference (from where the reads will
    originate). Must be a SeqRecord object.
    coverage -- desired coverage of the genome
    ErrorModel -- an ErrorModel class
    """
    header = record.id
    sequence = record.seq

    read_length = ErrorModel.read_length
    insert_size = ErrorModel.insert_size

    n_pairs = int(round((coverage * len(sequence)) / read_length) / 2)

    for i in range(n_pairs):
        # generate the forward read
        forward_start = random.randrange(0, len(sequence) - read_length)
        forward_end = forward_start + read_length
        forward = SeqRecord(
            Seq(str(sequence[forward_start:forward_end]),
                IUPAC.unambiguous_dna
                ),
            id='%s_%s_1' % (header, i),
            description=''
        )
        # add the quality and modify the nucleotides accordingly
        forward = ErrorModel.introduce_error_scores(forward, 'forward')
        forward.seq = ErrorModel.mut_sequence(forward, 'forward')

        # generate the reverse read
        reverse_start = forward_start + insert_size
        reverse_end = reverse_start + read_length
        reverse = SeqRecord(
            Seq(str(sequence[forward_start:forward_end]),
                IUPAC.unambiguous_dna
                ),
            id='%s_%s_1' % (header, i),
            description=''
        )
        # add the quality and modify the nucleotides accordingly
        reverse = ErrorModel.introduce_error_scores(reverse, 'reverse')
        reverse.seq = ErrorModel.mut_sequence(reverse, 'reverse')

        yield(forward, reverse.reverse_complement(
            id='%s_%s_2' % (header, i),
            description=''
            ))


def to_fastq(generator, output):
    """Take a generator and write to a file in fastq format"""
    # define name of output files
    output_forward = output + '_R1.fastq'
    output_reverse = output + '_R2.fastq'

    with open(output_forward, 'a') as f, open(output_reverse, 'a') as r:
        for read_tuple in generator:
            SeqIO.write(read_tuple[0], f, 'fastq-sanger')
            SeqIO.write(read_tuple[1], r, 'fastq-sanger')
