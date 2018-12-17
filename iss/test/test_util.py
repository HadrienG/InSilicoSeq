#!/usr/bin/env python
# -*- coding: utf-8 -*-

from nose.tools import raises, with_setup

from iss import util

from Bio import SeqIO

import random


def setup_function():
    output_file_prefix = 'data/.test'


def teardown_function():
    util.cleanup(['data/test_concat.iss.tmp.genomes.fasta'])


def test_phred_conversions():
    assert util.phred_to_prob(40) == 0.9999
    assert util.phred_to_prob(30) == 0.999
    assert util.phred_to_prob(20) == 0.99
    assert util.phred_to_prob(10) == 0.9

    assert util.prob_to_phred(0.9999) == 40
    assert util.prob_to_phred(0.999) == 30
    assert util.prob_to_phred(0.99) == 20
    assert util.prob_to_phred(0.9) == 10


def test_rev_comp():
    lowercase_seq = 'attgctat'
    uppercase_seq = 'CCGATTAC'
    assert util.rev_comp(lowercase_seq) == 'atagcaat'
    assert util.rev_comp(uppercase_seq) == 'GTAATCGG'


def test_count_records():
    f = open('data/genomes.fasta', 'r')
    with f:  # count the number of records
        record_list = util.count_records(f)
    assert record_list == [
        'genome_A', 'genome_T', 'genome_GC', 'genome_ATCG', 'genome_TA']


@raises(SystemExit)
def test_count_records_empty():
    util.count_records('data/empty_file')


def test_split_list():
    test_range = range(10)
    test_list = list(range(10))
    # tests on range
    assert util.split_list(test_range) == [range(10)]
    assert util.split_list(test_range, n_parts=2) == [
        range(0, 5),
        range(5, 10)
    ]
    assert util.split_list(test_range, n_parts=3) == [
        range(0, 3),
        range(3, 6),
        range(6, 10)
    ]
    # tests on lists
    assert util.split_list(test_list) == [list(range(10))]
    assert util.split_list(test_list, n_parts=2) == [
        list(range(0, 5)),
        list(range(5, 10))]
    assert util.split_list(test_list, n_parts=3) == [
        list(range(0, 3)),
        list(range(3, 6)),
        list(range(6, 10))
    ]


def test_convert_n_reads():
    simple_number = '10000'
    kilo_lower = '100k'
    kilo_upper = '100K'
    mega_lower = '20m'
    mega_upper = '20M'
    giga_lower = '1.5g'
    giga_upper = '0.5G'

    assert util.convert_n_reads(simple_number) == 10000
    assert util.convert_n_reads(kilo_lower) == 100000
    assert util.convert_n_reads(kilo_upper) == 100000
    assert util.convert_n_reads(mega_lower) == 20000000
    assert util.convert_n_reads(mega_upper) == 20000000
    assert util.convert_n_reads(giga_lower) == 1500000000
    assert util.convert_n_reads(giga_upper) == 500000000


@raises(SystemExit)
def test_convert_n_reads_bad_not_a_number():
    not_valid = 'rocket0'
    util.convert_n_reads(not_valid)


@raises(SystemExit)
def test_convert_n_reads_bad_suffix():
    not_valid = '0.12Q'
    util.convert_n_reads(not_valid)


@raises(SystemExit)
def test_genome_file_exists():
    my_important_file = 'data/ecoli.fasta'
    util.genome_file_exists(my_important_file)


def test_reservoir():
    samples = []
    genome_file = 'data/genomes.fasta'
    with open(genome_file, 'r') as f:
        record_list = util.count_records(f)
    n = 2
    with open(genome_file, 'r') as f:
        fasta_file = SeqIO.parse(f, 'fasta')
        for record in util.reservoir(fasta_file, record_list, n):
            samples.append(record.id)
    assert len(samples) == 2


@raises(SystemExit)
def test_reservoir_invalid_input():
    genome_file = 'data/ecoli.fasta'
    record_list = ['NC_002695.1']
    n = 4
    with open(genome_file, 'r') as f:
        fasta_file = SeqIO.parse(f, 'fasta')
        for record in util.reservoir(fasta_file, record_list, n):
            pass


@with_setup(setup_function, teardown_function)
def test_concatenate():
    genome_files = ['data/ecoli.fasta'] * 2
    util.concatenate(genome_files, 'data/test_concat.iss.tmp.genomes.fasta')
    with open('data/test_concat.iss.tmp.genomes.fasta', 'rb') as f:
        assert len(f.readlines()) == 40


@raises(SystemExit)
def test_concatenate_read_only():
    genome_files = ['data/ecoli.fasta'] * 2
    util.concatenate(genome_files, 'data/read_only.fasta')
