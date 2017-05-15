#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys


def parse_abundance_file(abundance_file):
    """Parse an abundance file and return a dictionnary. The abundance file
    is a flat file of the format "genome_id\tabundance\n"

    Arguments:
    abundance_file -- the path to the abundance file
    """
    abundance_dic = {}
    try:
        f = open(abundance_file, 'r')
    except IOError as e:
        print('Error:', e)
        sys.exit(1)
    else:
        with f:
            for line in f:
                genome_id = line.split()[0]
                abundance = float(line.split()[1])
                abundance_dic[genome_id] = abundance
    return abundance_dic


def to_coverage(total_n_reads, species_abundance, read_length, genome_size):
    """Calculate the coverage of a genome in a metagenome given its size and
    abundance

    Arguments:
    total_n_reads -- total amount of reads in the dataset
    species_abundance -- abundance of the species
    read_length -- length of the reads in the dataset
    genome_size -- size of the genome
    """
    n_reads = total_n_reads * species_abundance
    coverage = (n_reads * read_length) / genome_size
    return coverage
