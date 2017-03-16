#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import bam


def test_substitutions():
    bam_file = 'data/SRR1660402_mapped.bam'
    sub_forward, sub_reverse = bam.substitutions(bam_file, 150)
    assert len(sub_forward) == 150
    assert sub_forward[0][0] == 701


def test_quality_distribution():
    bam_file = 'data/SRR1660402_mapped.bam'
    hist_forward, hist_reverse = bam.quality_distribution(bam_file)
    assert len(hist_forward) == 150
    assert len(hist_reverse) == 140
    assert hist_forward[0][0][2] == 719
    assert hist_reverse[139][0][39] == 156


def test_get_insert_size():
    bam_file = 'data/SRR1660402_mapped.bam'
    assert bam.get_insert_size(bam_file) == 227
