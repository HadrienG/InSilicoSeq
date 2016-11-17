#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import random


def reads(input_genome, read_length, coverage):
    """Simulates perfect reads. Each read is a
    SeqRecord object. Return a generator.

    Arguments:
    input_genome -- genome of reference (from where the reads will originate).
    Must be in fasta format.
    read_length -- desired read length (int)
    coverage -- desired coverage of the genome (int)
    """
    with open(input_genome, 'r') as f:
        fasta_file = SeqIO.parse(f, 'fasta')
        for record in fasta_file:
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
            break  # we only want a single-fasta so we break after 1 record


def to_fastq(generator, output):
    """Take a generator and write to a file in fastq format"""
    with open(output, 'w') as o:
        for read in generator:
            SeqIO.write(read, o, 'fastq-sanger')
