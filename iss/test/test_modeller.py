#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import bam
from iss import modeller
from nose.tools import assert_almost_equals

import numpy as np


def test_cdf_qualities():
    quality_distribution = [
        [40, 30],
        [40, 30],
        [20, 20],
        [40, 10],
        [10, 10]]

    weights_first_base, weights_sec_base = \
        modeller.raw_qualities_to_histogram(quality_distribution, 'cdf')
    assert weights_first_base[1][10] == 0.2
    assert weights_first_base[1][20] == 0.2
    assert weights_first_base[1][-1] == 0.6

    assert weights_sec_base[1][10] == 0.4
    assert weights_sec_base[1][20] == 0.2
    assert weights_sec_base[1][30] == 0.4


def test_kde_qualities():
    quality_distribution = [
        [40, 30],
        [40, 30],
        [20, 20],
        [40, 10],
        [10, 10]]
    cdf_first_base, cdf_sec_base = \
        modeller.raw_qualities_to_histogram(quality_distribution, 'kde')
    assert_almost_equals(cdf_first_base[10], 0.199999702)
    assert_almost_equals(cdf_first_base[20], 0.400000149)
    assert cdf_first_base[-1] == 1

    assert_almost_equals(cdf_sec_base[10], 0.399998509)
    assert_almost_equals(cdf_sec_base[20], 0.599999255)
    assert_almost_equals(cdf_sec_base[30], 0.999998509)
    assert cdf_sec_base[-1] == 1


def test_substitutions():
    subst_matrix = np.zeros([20, 16])
    bam_file = 'data/ecoli.bam'
    bam_reader = bam.read_bam(bam_file)
    for _ in range(2):
        bam_reader.__next__()
    read = bam_reader.__next__()  # read_1_2
    alignment = read.get_aligned_pairs(matches_only=True, with_seq=True)
    read_has_indels = False
    for base in alignment:
        pos, subst, read_has_indels = modeller.dispatch_subst(
            base, read, read_has_indels)
        subst_matrix[pos, subst] += 1
    choices = modeller.subst_matrix_to_choices(subst_matrix, 20)
    assert read_has_indels is False
    assert subst_matrix[0][1] == 1
    assert choices[0]['A'] == (['T', 'C', 'G'], [1.0, 0.0, 0.0])


def test_indels():
    indel_matrix = np.zeros([20, 9])
    bam_file = 'data/ecoli.bam'
    bam_reader = bam.read_bam(bam_file)
    for _ in range(8):
        bam_reader.__next__()
    read = bam_reader.__next__()  # read_4_1
    for pos, indel in modeller.dispatch_indels(read):
        indel_matrix[pos, indel] += 1
    for position in range(20):
        indel_matrix[position][0] = 5
    insertion, deletion = modeller.indel_matrix_to_choices(
        indel_matrix, 20)
    assert round(insertion[6]['T'], 2) == 0.2
    assert indel_matrix[6][2] == 1
