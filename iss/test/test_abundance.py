#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util
from iss import abundance
from nose.tools import raises, with_setup

import numpy as np


def setup_function():
    output_file_prefix = 'data/.test'


def teardown_function():
    util.cleanup(['data/test_abundance.txt'])


def test_parsing():
    abundance_dic = abundance.parse_abundance_file('data/abundance.txt')
    assert abundance_dic == {
        'genome_ATCG': 0.1,
        'genome_TA': 0.1,
        'genome_A': 0.2,
        'genome_GC': 0.4,
        'genome_T': 0.2
    }


@raises(SystemExit)
def test_parsing_empty():
    abundance_dic = abundance.parse_abundance_file('data/empty_file')


@raises(SystemExit)
def test_parsing_no_exists():
    abundance_dic = abundance.parse_abundance_file('data/does_not_exist')


@raises(SystemExit)
def test_parsing_bad_abundance():
    abundance_dic = abundance.parse_abundance_file('data/bad_abundance.txt')


def test_cov_calc():
    coverage_ecoli = abundance.to_coverage(
        10000000,
        0.08,
        150,
        4639221
    )
    assert round(coverage_ecoli, 3) == 25.866


def test_distributions():
    np.random.seed(42)
    f = open('data/genomes.fasta', 'r')
    with f:  # count the number of records
        record_list = util.count_records(f)

    uniform_dic = abundance.uniform(record_list)
    halfnormal_dic = abundance.halfnormal(record_list)
    exponential_dic = abundance.exponential(record_list)
    lognormal_dic = abundance.lognormal(record_list)

    np.random.seed(42)  # reset the seed to get 0s in zero_inflated_lognormal
    zero_inflated_lognormal_dic = abundance.zero_inflated_lognormal(
        record_list)
    assert list(uniform_dic.values()) == [0.2] * 5
    assert round(halfnormal_dic['genome_A'], 2) == 0.16
    assert round(exponential_dic['genome_A'], 2) == 0.01
    assert round(lognormal_dic['genome_T'], 2) == 0.19
    assert zero_inflated_lognormal_dic['genome_T'] == 0.0
    assert round(zero_inflated_lognormal_dic['genome_A'], 2) == 0.44


@with_setup(setup_function, teardown_function)
def test_abunance_draft():
    abundance_dic = {'genome_A': 0.15511887441170918,
                     'genome_T': 0.08220476760848751,
                     'genome_GC': 0.18039811160555874,
                     'genome_ATCG': 0.4329003045949206,
                     'genome_TA': 0.07468835777633397,
                     'contig_1': 0.02776920430880394,
                     'contig_2': 0.011490705231229217,
                     'contig_3': 0.03542967446295675}
    np.random.seed(42)
    f = open('data/genomes.fasta', 'r')
    with f:  # count the number of records
        complete_genomes = util.count_records(f)
    draft_genomes = ['data/draft.fasta']
    ab = abundance.draft(
        complete_genomes,
        draft_genomes,
        abundance.lognormal,
        'data/test')
    for tv, v in zip(abundance_dic.values(), ab.values()):
        assert round(tv) == round(v)
    # assert_almost_equal(ab, abundance_dic)


def test_coverage_scaling():
    abundance_dic = abundance.parse_abundance_file('data/abundance.txt')
    coverage_dic = abundance.coverage_scaling(10000, abundance_dic,
                                              'data/genomes.fasta', 25)
    assert coverage_dic == {'genome_A': 136.6120218579235,
                            'genome_T': 136.6120218579235,
                            'genome_GC': 273.224043715847,
                            'genome_ATCG': 68.30601092896175,
                            'genome_TA': 68.30601092896175
                            }
