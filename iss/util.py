#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from builtins import dict

from Bio import SeqIO

import os
import sys
import logging
import numpy as np


def phred_to_prob(q):
    """Convert a phred score (Sanger or modern Illumina) in probabilty

    Given a phred score q, return the probabilty p
    of the call being right

    Args:
        q (int): phred score

    Returns:
        float: probabilty of basecall being right
    """
    p = 10 ** (-q / 10)
    return 1 - p


def prob_to_phred(p):
    """Convert a probabilty into a phred score (Sanger or modern Illumina)

    Given a probabilty p of the basecall being right, return the
    phred score q

    Args:
        p (int): probabilty of basecall being right

    Returns:
        int: phred score
    """
    q = int(round(-10 * np.log10(1 - p)))
    return q


def rev_comp(s):
    """A simple reverse complement implementation working on strings

    Args:
        s (string): a DNA sequence (IUPAC, can be ambiguous)

    Returns:
        list: reverse complement of the input sequence
    """
    bases = {
        "a": "t", "c": "g", "g": "c", "t": "a", "y": "r", "r": "y", "w": "w",
        "s": "s", "k": "m", "m": "k", "n": "n", "A": "T", "C": "G", "G": "C",
        "T": "A", "Y": "R", "R": "Y", "W": "W", "S": "S", "K": "M", "M": "K",
        "N": "N"}
    sequence = list(s)
    complement = "".join([bases[b] for b in sequence])
    reverse_complement = complement[::-1]
    return reverse_complement


def count_records(fasta_file):
    """Count the number of records in a fasta file and return a list of
    recods id

    Args:
        fasta_file (string): the path to a fasta file

    Returns:
        list: a list of record ids
    """
    logger = logging.getLogger(__name__)
    record_list = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        record_list.append(record.id)
    try:
        assert len(record_list) != 0
    except AssertionError as e:
        logger.error(
            'Failed to find records in genome(s) file:%s' % fasta_file)
        sys.exit(1)
    else:
        return record_list


def split_list(l, n_parts=1):
    """Split a list in a number of parts

    Args:
        l (list): a list
        n_parts (in): the number of parts to split the list in

    Returns:
        list: a list of n_parts lists
    """
    length = len(l)
    return [l[i * length // n_parts: (i + 1) * length // n_parts]
            for i in range(n_parts)]


def nplog(type, flag):
    logger = logging.getLogger(__name__)
    logger.debug("FloatingPointError (%s), with flag %s" % (type, flag))


def convert_n_reads(unit):
    """For strings representing a number of bases and ending with k, K, m, M,
    g, and G converts to a plain old number

    Args:
        n (str): a string representing a number ending with a suffix
    Returns:
        float: a number of reads
    """
    logger = logging.getLogger(__name__)
    suffixes = {'k': 3, 'm': 6, 'g': 9}
    if unit[-1].isdigit():
        try:
            unit_int = int(unit)
        except ValueError as e:
            logger.error('%s is not a valid number of reads' % unit)
            sys.exit(1)
    elif unit[-1].lower() in suffixes:
        number = unit[:-1]
        exponent = suffixes[unit[-1].lower()]
        unit_int = int(float(number) * 10**exponent)
    else:
        logger.error('%s is not a valid number of reads' % unit)
        sys.exit(1)
    return unit_int


def genome_file_exists(filename):
    logger = logging.getLogger(__name__)
    try:
        assert os.path.exists(filename) == False
    except AssertionError as e:
        logger.error('%s already exists. Aborting.' % filename)
        logger.error('Maybe use --genomes %s' % filename)
        logger.error('or use --ncbi with another output prefix')
        sys.exit(1)
