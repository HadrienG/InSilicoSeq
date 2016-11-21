#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import abundance


def test_parsing():
    abundance_dic = abundance.parse_abundance_file('data/abundance.txt')
    assert abundance_dic == {
        'genome_ATCG': 0.1,
        'genome_TA': 0.1,
        'genome_A': 0.2,
        'genome_GC': 0.4,
        'genome_T': 0.2
        }


def test_cov_calc():
    coverage_ecoli = abundance.to_coverage(
        10000000,
        0.08,
        150,
        4639221
        )
    assert round(coverage_ecoli, 3) == 25.866
