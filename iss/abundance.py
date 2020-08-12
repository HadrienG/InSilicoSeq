#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from scipy import stats

import os
import sys
import logging
import numpy as np


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
        logger.error('Failed to open file:%s' % e)
        sys.exit(1)
    except AssertionError as e:
        logger.error('File seems empty: %s' % abundance_file)
        sys.exit(1)
    else:
        with f:
            for line in f:
                try:
                    genome_id = line.split()[0]
                    abundance = float(line.split()[1])
                except IndexError as e:
                    logger.error('Failed to read file: %s' % e)
                    sys.exit(1)
                else:
                    abundance_dic[genome_id] = abundance
    logger.debug('Loaded abundance/coverage file: %s' % abundance_file)
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


def coverage_scaling(total_n_reads, abundance_dic, genome_file, read_length):
    """Scale coverage distribution according to the n_reads parameter

    Args:
        total_n_reads (int): total amount of reads to simulate
        abundance_dic (dict): a dictionary with records as keys, coverage as
            values
        genome_file (str): path to input fasta file containing genomes
        read_length (int): length of the reads in the dataset

    Returns:
        dict: scaled coverage dictionary
    """
    total_reads = 0
    f = open(genome_file, 'r')  # re-opens the file
    with f:
        fasta_file = SeqIO.parse(f, 'fasta')
        for record in fasta_file:
            try:
                species_coverage = abundance_dic[record.id]
            except KeyError as e:
                logger.error(
                    'Fasta record not found in abundance file: %s' % e)
                sys.exit(1)

            genome_size = len(record.seq)
            reads_g = species_coverage * genome_size / read_length / 2
            total_reads += reads_g

    scale_factor = total_n_reads / total_reads
    for key in abundance_dic:
        abundance_dic[key] *= scale_factor
    return abundance_dic


def to_file(abundance_dic, output, mode="abundance"):
    """Write the abundance dictionary to a file

    Args:
        abundance_dic (dict): the abundance dictionary
        output (str): the output file name
    """
    logger = logging.getLogger(__name__)
    if mode == "abundance":
        output_abundance = output + '_abundance.txt'
    else:
        output_abundance = output + '_coverage.txt'
    try:
        f = open(output_abundance, 'w')
    except PermissionError as e:
        logger.error('Failed to open output file: %s' % e)
        sys.exit(1)
    else:
        with f:
            for record, abundance in abundance_dic.items():
                f.write('%s\t%s\n' % (record, abundance))


def draft(genomes, draft, distribution, output, mode="abundance"):
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
    draft_dic = expand_draft_abundance(abundance_dic, draft, mode)
    full_abundance_dic = {**complete_genomes_abundance, **draft_dic}
    return full_abundance_dic


def expand_draft_abundance(abundance_dic, draft, mode="abundance"):
    """Calculate abundance for each contig of a draft genome
    The function takes the abundance dictionary and automatically
    detects draft genomes. In coverage mode the function simply assign
    the coverage value to each contig

    Args:
        abundance_dic (dict): dict with genome (paths or id) as key and
            abundance as value
        draft (list): draft genome files
        mode (str): abundance or coverage


    Returns:
        dict: abundance dictionary with abundance value for each contig
    """
    draft_dic = {}
    for key, abundance in abundance_dic.items():
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
                    if mode == "abundance":
                        length = len(record)
                        contig_abundance = abundance * (length / total_length)
                        # print(key, record.id, contig_abundance)
                        draft_dic[record.id] = contig_abundance
                    elif mode == "coverage":
                        draft_dic[record.id] = abundance
    return draft_dic
