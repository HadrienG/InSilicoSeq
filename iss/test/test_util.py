#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util


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
