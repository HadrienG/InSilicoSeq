#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

import pytest

from iss import bam


def test_read_fail():
    with pytest.raises(SystemExit):
        bam_file = "data/empty_file"
        bam_reader = bam.read_bam(bam_file)
        for read in bam_reader:
            print(read)


def test_to_model():
    bam_file = "data/ecoli.bam"
    output = "data/test_bam"
    bam.to_model(bam_file, output)
    os.remove(output + ".npz")
