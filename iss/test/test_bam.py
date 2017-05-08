#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import bam
from nose.tools import assert_almost_equals


def test_mismatches():
    bam_file = 'data/SRR1660402_mapped.bam'
    sub_forward, sub_reverse = bam.get_mismatches(bam_file, 150)
    assert len(sub_forward) == 150
    assert sub_forward[-1]['G'] == (
        ['A', 'T', 'C'],
        [0.17763157894736842, 0.20394736842105263, 0.61842105263157898])


def test_substitutions_to_choices():
    dispatch = np.array(
        [2.77297000e+05, 2.58000000e+02, 3.22000000e+02, 7.37000000e+02,
            2.81467000e+05, 2.18000000e+02, 4.17000000e+02, 3.27000000e+02,
            2.87455000e+05, 2.55000000e+02, 1.96000000e+02, 2.65000000e+02,
            2.49238000e+05, 2.40000000e+02, 4.10000000e+02, 2.63000000e+02]
        )
    assert subst_matrix_to_choices(dispatch)['A'] == (
        ['T', 'C', 'G'],
        [0.1958997722095672, 0.24449506454062261, 0.55960516324981013])


def test_quality_distribution_cdf():
    bam_file = 'data/SRR1660402_mapped.bam'
    hist_forward, hist_reverse = bam.quality_distribution('cdf', bam_file)
    assert len(hist_forward) == 150
    assert len(hist_reverse) == 140
    assert hist_forward[0][0][2] == 2
    assert_almost_equals(hist_reverse[139][1][39], 0.03439153439153439)


def test_get_insert_size():
    bam_file = 'data/SRR1660402_mapped.bam'
    assert bam.get_insert_size(bam_file) == 227
