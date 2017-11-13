#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

import os
import sys
import logging
import numpy as np

from scipy import stats


def parse_abundance_file(abundance_file):
    """Parse an abundance file

    The abundance file is a flat file of the format "genome_id<TAB>abundance"

    Args:
        abundance_file (string): the path to the abundance file

    Returns:
        dict: genome_id as keys, abundance as
            values
    """
    logger = logging.getLogger(__name__)
    abundance_dic = {}
    try:
        assert os.stat(abundance_file).st_size != 0
        f = open(abundance_file, 'r')
    except IOError as e:
        logger.error('Failed to open abundance file:%s' % e)
        sys.exit(1)
    except AssertionError as e:
        logger.error('Abundance file seems empty: %s' % abundance_file)
        sys.exit(1)
    else:
        with f:
            for line in f:
                try:
                    genome_id = line.split()[0]
                    abundance = float(line.split()[1])
                except IndexError as e:
                    logger.error('Failed to read abundance file: %s' % e)
                else:
                    abundance_dic[genome_id] = abundance
    logger.info('Loaded abundance file: %s' % abundance_file)
    return abundance_dic


def uniform(record_list):
    """Generate uniform abundance distribution from a number of records

    Args:
        record_list (list): a list of record.id

    Returns:
        list: a list of floats
    """
    abundance_dic = {}
    n_records = len(record_list)
    for record in record_list:
        abundance_dic[record] = 1 / n_records

    return abundance_dic


def halfnormal(record_list):
    """Generate scaled halfnormal abundance distribution from a number of
        records

    Args:
        record_list (list): a list of record.id

    Returns:
        list: a list of floats
    """
    abundance_dic = {}
    n_records = len(record_list)
    dist = stats.halfnorm.rvs(loc=0.00, scale=1.00, size=n_records)
    dist_scaled = dist / sum(dist)
    for record, abundance in zip(record_list, dist_scaled):
        abundance_dic[record] = abundance

    return abundance_dic


def exponential(record_list):
    """Generate scaled exponential abundance distribution from a number of
        records

    Args:
        record_list (list): a list of record.id

    Returns:
        list: a list of floats
    """
    abundance_dic = {}
    n_records = len(record_list)
    dist = np.random.exponential(size=n_records)
    dist_scaled = dist / sum(dist)
    for record, abundance in zip(record_list, dist_scaled):
        abundance_dic[record] = abundance

    return abundance_dic


def lognormal(record_list):
    """Generate scaled lognormal abundance distribution from a number of
        records

    Args:
        record_list (list): a list of record.id

    Returns:
        list: a list of floats
    """
    abundance_dic = {}
    n_records = len(record_list)
    dist = np.random.lognormal(size=n_records)
    dist_scaled = dist / sum(dist)
    for record, abundance in zip(record_list, dist_scaled):
        abundance_dic[record] = abundance

    return abundance_dic


def zero_inflated_lognormal(record_list):
    """Generate scaled zero-inflated lognormal abundance distribution from a
        number of records

    Args:
        record_list (list): a list of record.id

    Returns:
        list: a list of floats
    """
    abundance_dic = {}
    n_records = len(record_list)
    zero_inflated = stats.bernoulli.rvs(p=0.2, size=n_records)
    dist = (1 - zero_inflated) * np.random.lognormal(size=n_records)
    dist_scaled = dist / sum(dist)
    for record, abundance in zip(record_list, dist_scaled):
        abundance_dic[record] = abundance

    return abundance_dic


def to_coverage(total_n_reads, species_abundance, read_length, genome_size):
    """Calculate the coverage of a genome in a metagenome given its size and
    abundance

    Args:
        total_n_reads (int): total amount of reads in the dataset
        species_abundance (float): abundance of the species, between 0 and 1
        read_length (int): length of the reads in the dataset
        genome_size (int): size of the genome

    Returns:
        float: genome coverage
    """
    n_reads = total_n_reads * species_abundance
    coverage = (n_reads * read_length) / genome_size
    return coverage


def to_file(abundance_dic, output):
    """write the abundance dictionary to a file

    Args:
        abundance_dic (dict): the abundance dictionary
        output (str): the output file name
    """
    logger = logging.getLogger(__name__)
    output_abundance = output + '_abundance.txt'
    try:
        f = open(output_abundance, 'w')
    except PermissionError as e:
        logger.error('Failed to open output file: %s' % e)
        sys.exit(1)
    else:
        with f:
            for record, abundance in abundance_dic.items():
                f.write('%s\t%s\n' % (record, abundance))
