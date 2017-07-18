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
