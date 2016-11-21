#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import random


def reads(record, read_length, coverage):
    """Simulates perfect reads. Each read is a
    SeqRecord object. Return a generator.

    Arguments:
    input_record -- sequence of reference (from where the reads will
    originate). Must be a SeqRecord object.
    read_length -- desired read length (int)
    coverage -- desired coverage of the genome
    """
    header = record.id
    sequence = record.seq

    n_reads = int(round((coverage * len(sequence)) / read_length))
    for _ in range(n_reads):
        read_start = random.randrange(0, len(sequence) - read_length)
        read_end = read_start + read_length
        read = SeqRecord(
                        Seq(str(sequence[read_start:read_end]),
                            IUPAC.unambiguous_dna
                            ),
                        id='%s_%s' % (header, _),
                        description=''
                        )
        # add the quality
        read.letter_annotations["phred_quality"] = [40] * len(read)
        yield(read)


def to_fastq(generator, output):
    """Take a generator and write to a file in fastq format"""
    with open(output, 'a') as o:
        for read in generator:
            SeqIO.write(read, o, 'fastq-sanger')
