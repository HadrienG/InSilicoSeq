#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO

import os
import sys
import gzip
import pickle
import random
import logging
import numpy as np

from shutil import copyfileobj


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
        "s": "s", "k": "m", "m": "k", "n": "n", "b": "v", "v": "b", "d": "h",
        "h": "d", "A": "T", "C": "G", "G": "C", "T": "A", "Y": "R", "R": "Y",
        "W": "W", "S": "S", "K": "M", "M": "K", "N": "N", "B": "V", "V": "B",
        "D": "H", "H": "D"}
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
    """Checks if the output file from the --ncbi option already exists

    Args:
        filename (str): a file name
    """
    logger = logging.getLogger(__name__)
    try:
        assert os.path.exists(filename) is False
    except AssertionError as e:
        logger.error('%s already exists. Aborting.' % filename)
        logger.error('Maybe use another --output prefix')
        sys.exit(1)


def reservoir(records, record_list, n=None):
    """yield a number of records from a fasta file using reservoir sampling

    Args:
        records (obj): fasta records from SeqIO.parse

    Yields:
        record (obj): a fasta record
    """
    logger = logging.getLogger(__name__)
    if n is not None:
        try:
            total = len(record_list)
            assert n < total
        except AssertionError as e:
            logger.error(
                '-u should be strictly smaller than total number of records.')
            sys.exit(1)
        else:
            random.seed()
            x = 0
            samples = sorted(random.sample(range(0, total - 1), n))
            for sample in samples:
                while x < sample:
                    x += 1
                    _ = records.__next__()

                record = records.__next__()
                x += 1
                yield record
    else:
        for record in records:
            yield record


def concatenate(file_list, output):
    """Concatenate files together

    Args:
        file_list (list): the list of input files (can be a generator)
        output (string): the output file name
    """
    logger = logging.getLogger(__name__)
    logger.info('Stitching input files together')
    try:
        out_file = open(output, 'wb')
    except (IOError, OSError) as e:
        logger.error('Failed to open output file: %s' % e)
        sys.exit(1)

    with out_file:
        for file_name in file_list:
            if file_name is not None:
                with open(file_name, 'rb') as f:
                    copyfileobj(f, out_file)


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
                os.remove(temp_file)
            except (IOError, OSError) as e:
                logger.error('Could not read temporary file: %s' % temp_file)
                logger.error('You may have to remove temporary files manually')
                sys.exit(1)


def compress(filename, remove=True):
    """gzip a file

    Args:
        filename (string): name of file to be compressed
    """
    logger = logging.getLogger(__name__)
    logger.info('Compressing %s' % filename)
    outfile = filename + '.gz'
    with open(filename, 'rb') as i, gzip.open(outfile, 'wb') as o:
        copyfileobj(i, o)
    if remove:
        cleanup([filename])
    return outfile


def dump(object, output):
    """dump an object, like pickle.dump.
    This function uses pickle.dumps to dump large objects

    Args:
        object (object): a python object
    """
    MAX_BYTES = 2**31 - 1
    pickled_object = pickle.dumps(object, protocol=pickle.HIGHEST_PROTOCOL)
    size = sys.getsizeof(pickled_object)

    with open(output, 'wb') as out_file:
        for i in range(0, size, MAX_BYTES):
            out_file.write(pickled_object[i:i + MAX_BYTES])


def load(filename):
    """load a pickle from disk
    This function uses pickle.loads to load large objects

    Args:
        filename (string): the path of the pickle to load
    """
    MAX_BYTES = 2**31 - 1

    size = os.path.getsize(filename)
    bytes = bytearray(0)

    with open(filename, 'rb') as f:
        for _ in range(0, size, MAX_BYTES):
            bytes += f.read(MAX_BYTES)
    object = pickle.loads(bytes)

    return object
