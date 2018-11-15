#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

from scipy import stats

import os
import sys
import logging
import numpy as np

from Bio import SeqIO


def parse_abundance_file(abundance_file):
    """Parse an abundance or coverage file

    The abundance/coverage file is a flat file of the format
    "genome_id<TAB>abundance"
    or
    "genome_id<TAB>coverage"

    Args:
        abundance_file (string): the path to the abundance file

    Returns:
        dict: genome_id as keys, abundance as values
    """
    logger = logging.getLogger(__name__)
    abundance_dic = {}
    try:
        assert os.stat(abundance_file).st_size != 0
        f = open(abundance_file, 'r')
    except (IOError, OSError) as e:
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
                    sys.exit(1)
                else:
                    abundance_dic[genome_id] = abundance
    logger.info('Loaded abundance file: %s' % abundance_file)
    return abundance_dic


def uniform(record_list):
    """Generate uniform abundance distribution from a number of records

    Args:
        record_list (list): a list of record.id

    Returns:
        dict: a dictionary with records as keys, abundance as values
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
        dict: a dictionary with records as keys, abundance as values
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
        dict: a dictionary with records as keys, abundance as values
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
        dict: a dictionary with records as keys, abundance as values
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
        dict: a dictionary with records as keys, abundance as values
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
        float: coverage of the genome
    """
    n_reads = total_n_reads * species_abundance
    coverage = (n_reads * read_length) / genome_size
    return coverage


def to_file(abundance_dic, output):
    """Write the abundance dictionary to a file

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


def draft(genomes, draft, distribution, output):
    """Computes the abundance dictionary for a mix of complete and
    draft genomes

    Args:
        genomes (list): list of all input records
        draft (list): draft genome files
        distribution (function): distribution function to use
        output (str): output file

    Returns:
        dict: the abundance dictionary
    """
    # first we get a list of contig names in draft genomes
    draft_records = []
    for d in draft:
        draft_records.extend(
            [record.name for record in SeqIO.parse(d, 'fasta')])
    genomes = list(set(genomes) - set(draft_records))
    abundance_dic = distribution(genomes + draft)
    complete_genomes_abundance = {k: v for
                                  k, v in abundance_dic.items()
                                  if k not in draft}
    to_file(abundance_dic, output)
    for key, abundance in abundance_dic.items():
        draft_dic = {}
        if key in draft:
            try:
                f = open(key, 'r')
                with f:
                    fasta_file = SeqIO.parse(f, 'fasta')
                    draft_records = [record for record in fasta_file]
                    total_length = sum(len(record) for record in draft_records)
            except Exception as e:
                raise
            else:
                total_length = sum(len(record) for record in draft_records)
                for record in draft_records:
                    length = len(record)
                    contig_abundance = abundance * (length / total_length)
                    # print(key, record.id, contig_abundance)
                    draft_dic[record.id] = contig_abundance
    # python 3.5 +
    # full_abundance_dic = {**complete_genomes_abundance, **draft_dic}
    full_abundance_dic = complete_genomes_abundance.copy()
    full_abundance_dic.update(draft_dic)
    return full_abundance_dic
