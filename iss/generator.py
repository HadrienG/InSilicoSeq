#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import error_model
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import random


def reads(record, read_length, coverage, insert_size, mean_qual):
    """Simulate perfect reads. Each read is a SeqRecord object. Return a
    generator of tuples containing the forward and reverse read.

    Arguments:
    input_record -- sequence of reference (from where the reads will
    originate). Must be a SeqRecord object.
    read_length -- desired read length (int)
    coverage -- desired coverage of the genome
    insert_size -- insert size between the pairs
    mean_qual -- mean quality score
    """
    header = record.id
    sequence = record.seq

    n_pairs = int(round((coverage * len(sequence)) / read_length) / 2)
    for i in range(n_pairs):
        forward_start = random.randrange(0, len(sequence) - read_length)
        forward_end = forward_start + read_length
        forward = SeqRecord(
            Seq(str(sequence[forward_start:forward_end]),
                IUPAC.unambiguous_dna
                ),
            id='%s_%s_1' % (header, i),
            description=''
        )
        # add the quality
        forward = error_model.introduce_errors(forward, mean_qual)

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
        # add the quality
        reverse = error_model.introduce_errors(reverse, mean_qual)

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
