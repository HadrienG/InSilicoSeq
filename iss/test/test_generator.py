#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pytest


from iss import generator
from iss.util import cleanup
from iss.error_models import basic, kde

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import sys
import os
import random
import numpy as np

# due to inconsistent seeding between python 2 and 3, some of the following
# tests are disabled with python2


def setup_function():
    output_file_prefix = 'data/.test'


def teardown_cleanup():
    cleanup(['data/.test.iss.tmp.my_genome.0_R1.fastq',
             'data/.test.iss.tmp.my_genome.0_R2.fastq'])


@pytest.fixture
def setup_and_teardown():
    setup_function()
    yield
    teardown_cleanup()


def test_cleanup_fail():
    with pytest.raises(SystemExit):
        cleanup('data/does_not_exist')


def test_simulate_and_save(setup_and_teardown):
    err_mod = basic.BasicErrorModel()
    ref_genome = SeqRecord(
        Seq(str('AAAAACCCCC' * 100)),
        id='my_genome',
        description='test genome'
    )
    generator.reads(ref_genome, err_mod, 1000, 0, 'data/.test', 0, 'metagenomics', gc_bias = True)


def test_simulate_and_save_short(setup_and_teardown):
    err_mod = basic.BasicErrorModel()
    ref_genome = SeqRecord(
        Seq(str('AACCC' * 100)),
        id='my_genome',
        description='test genome'
    )
    generator.reads(ref_genome, err_mod, 1000, 0, 'data/.test', 0, 'metagenomics', gc_bias = True)


def test_small_input():
    with pytest.raises(AssertionError):
        err_mod = kde.KDErrorModel('data/ecoli.npz')
        ref_genome = SeqRecord(
            Seq(str('AAAAACCCCC')),
            id='my_genome',
            description='test genome'
        )
        generator.simulate_read(ref_genome, err_mod, 1, 0)


def test_basic():
    if sys.version_info > (3,):
        random.seed(42)
        np.random.seed(42)
        err_mod = basic.BasicErrorModel()
        ref_genome = SeqRecord(
            Seq(str('AAAAACCCCC' * 100)),
            id='my_genome',
            description='test genome'
        )
        read_tuple = generator.simulate_read(ref_genome, err_mod, 1, 0, 'metagenomics')
        big_read = ''.join(str(read_tuple[0].seq) + str(read_tuple[1].seq))
        assert big_read[-15:] == 'TTTTGGGGGTTTTTG'


def test_kde():
    if sys.version_info > (3,):
        random.seed(42)
        np.random.seed(42)
        err_mod = kde.KDErrorModel('data/ecoli.npz')
        ref_genome = SeqRecord(
            Seq(str('CGTTTCAACC' * 400)),
            id='my_genome',
            description='test genome'
        )
        read_tuple = generator.simulate_read(ref_genome, err_mod, 1, 0, 'metagenomics')
        big_read = ''.join(str(read_tuple[0].seq) + str(read_tuple[1].seq))
        assert big_read[:15] == 'CCGTTTCAACCCGTT'


def test_kde_short():
    if sys.version_info > (3,):
        random.seed(42)
        np.random.seed(42)
        err_mod = kde.KDErrorModel('data/ecoli.npz')
        ref_genome = SeqRecord(
            Seq(str('AAACC' * 100)),
            id='my_genome',
            description='test genome'
        )
        read_tuple = generator.simulate_read(ref_genome, err_mod, 1, 0, 'metagenomics')
        big_read = ''.join(str(read_tuple[0].seq) + str(read_tuple[1].seq))
        assert big_read == 'ACCAAACCAAACCAAACCAAGGTTTGGTTTGGTTTGGTGT'


def test_amp():
    if sys.version_info > (3,):
        random.seed(42)
        np.random.seed(42)
        err_mod = kde.KDErrorModel(
            os.path.join(
                os.path.dirname(__file__),
                '../profiles/MiSeq'
            )
        )
        print(err_mod.read_length)
        ref_genome = SeqRecord(
            Seq((
                "TTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGG"
                "CCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTT"
                )),
            id='my_amplicon',
            description='test amplicon'
        )
        read_tuple = generator.simulate_read(ref_genome, err_mod, 1, 0, 'amplicon')
        assert len(read_tuple[0].seq) == 301
        assert read_tuple[0].seq.startswith("TTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
        assert len(read_tuple[1].seq) == 301
        assert read_tuple[1].seq.startswith("AAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
