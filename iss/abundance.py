#!/usr/bin/env python
# -*- coding: utf-8 -*-


def parse_abundance_file(abundance_file):
    abundance_dic = {}
    with open(abundance_file, 'r') as f:
        for line in f:
            genome_id = line.split()[0]
            abundance = float(line.split()[1])
            abundance_dic[genome_id] = abundance
    return abundance_dic


def to_coverage(total_n_reads, species_abundance, read_length, genome_size):
    n_reads = total_n_reads * species_abundance
    coverage = (n_reads * read_length) / genome_size
    return coverage
