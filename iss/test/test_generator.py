#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import unicode_literals

from iss import generator
from iss.error_models import ErrorModel, basic, kde

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from nose.tools import with_setup, raises

import os
import sys
import random
import numpy as np

# due to inconsistent seeding between python 2 and 3, some of the following
# tests are disabled with python2


def setup_function():
    output_file_prefix = 'data/.test'


def teardown_function():
    generator.cleanup(['data/.test.iss.tmp.my_genome.0'])


@raises(SystemExit)
def test_cleanup_fail():
    generator.cleanup('data/does_not_exist')


@with_setup(setup_function, teardown_function)
def test_concatenate():
    output = 'data/.test.iss.tmp.my_genome.0'
    file_list = ['data/ecoli', 'data/ecoli']
    generator.concatenate(file_list, output)


@with_setup(setup_function, teardown_function)
def test_simulate_and_save():
    err_mod = basic.BasicErrorModel()
    ref_genome = SeqRecord(
        Seq(str('AAAAACCCCC' * 100),
            IUPAC.unambiguous_dna
            ),
        id='my_genome',
        description='test genome'
        )
    generator.reads(ref_genome, err_mod, 1000, 0, 'data/.test', True)


@raises(ValueError)
def test_small_input():
    err_mod = kde.KDErrorModel('data/ecoli.npz')
    ref_genome = SeqRecord(
        Seq(str('AAAAACCCCC'),
            IUPAC.unambiguous_dna
            ),
        id='my_genome',
        description='test genome'
        )
    generator.simulate_read(ref_genome, err_mod, 1)


def test_basic():
    if sys.version_info > (3,):
        random.seed(42)
        np.random.seed(42)
        err_mod = basic.BasicErrorModel()
        ref_genome = SeqRecord(
            Seq(str('AAAAACCCCC' * 100),
                IUPAC.unambiguous_dna
                ),
            id='my_genome',
            description='test genome'
            )
        read_tuple = generator.simulate_read(ref_genome, err_mod, 1)
        big_read = ''.join(str(read_tuple[0].seq) + str(read_tuple[1].seq))
        assert big_read[-15:] == 'TTTTGGGGGTTTTTG'


def test_kde():
    if sys.version_info > (3,):
        random.seed(42)
        np.random.seed(42)
        err_mod = kde.KDErrorModel('data/ecoli.npz')
        ref_genome = SeqRecord(
            Seq(str('CGTTTCAACC' * 400),
                IUPAC.unambiguous_dna
                ),
            id='my_genome',
            description='test genome'
            )
        read_tuple = generator.simulate_read(ref_genome, err_mod, 1)
        big_read = ''.join(str(read_tuple[0].seq) + str(read_tuple[1].seq))
        assert big_read[:15] == 'CCGTTTCAACCCGTT'
