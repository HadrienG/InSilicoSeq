#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import bam

from nose.tools import raises

import os
import sys


@raises(SystemExit)
def test_read_fail():
    bam_file = 'data/empty_file'
    bam_reader = bam.read_bam(bam_file)
    for read in bam_reader:
        print(read)


def test_to_model():
    bam_file = 'data/ecoli.bam'
    output = 'data/test_bam'
    bam.to_model(bam_file, output)
    os.remove(output + '.npz')
