#!/usr/bin/env python
# -*- coding: utf-8 -*-

from builtins import next

from iss import bam
from iss import modeller
from nose.tools import assert_almost_equals

import sys
import numpy as np


def test_kde_qualities():
    quality_distribution = [
        [40, 30],
        [40, 30],
        [20, 20],
        [40, 10],
        [10, 10]]
    cdf_list = modeller.raw_qualities_to_histogram(quality_distribution)
    assert_almost_equals(cdf_list[0][-2], 0.500002794)
    assert cdf_list[-1][0] == 0.0
    assert cdf_list[-1][-1] == 1
    assert len(cdf_list) == 5


def test_substitutions():
    subst_matrix = np.zeros([20, 16])
    bam_file = 'data/substitutions_test.bam'
    bam_reader = bam.read_bam(bam_file)
    if sys.version_info > (3,):
        for _ in range(2):
            bam_reader.__next__()
        read = bam_reader.__next__()  # read_1_2
    else:
        for _ in range(2):
            bam_reader.next()
        read = bam_reader.next()  # read_1_2
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
    bam_file = 'data/substitutions_test.bam'
    bam_reader = bam.read_bam(bam_file)
    if sys.version_info > (3,):
        for _ in range(8):
            bam_reader.__next__()
        read = bam_reader.__next__()  # read_4_1
    else:
        for _ in range(8):
            bam_reader.next()
        read = bam_reader.next()  # read_4_1
    for pos, indel in modeller.dispatch_indels(read):
        indel_matrix[pos, indel] += 1
    for position in range(20):
        indel_matrix[position][0] = 5
    insertion, deletion = modeller.indel_matrix_to_choices(
        indel_matrix, 20)
    assert round(insertion[6]['T'], 2) == 0.2
    assert indel_matrix[6][2] == 1
