#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import generator
from iss.error_models import ErrorModel, basic, cdf, kde

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import random
import numpy as np


def test_basic():
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
    read_gen = generator.reads(ref_genome, 2, err_mod)
    big_read = ''.join(
        str(read_tuple[0].seq) + str(read_tuple[1].seq)
        for read_tuple in read_gen)
    assert big_read[1890:1910] == 'GGGGGTGTTTGGGGGTTTTT'


def test_cdf():
    random.seed(42)
    np.random.seed(42)
    err_mod = cdf.CDFErrorModel('data/ecoli_cdf.npz')
    ref_genome = SeqRecord(
        Seq(str('CGTTTCAACC' * 40),
            IUPAC.unambiguous_dna
            ),
        id='my_genome',
        description='test genome'
        )
    read_gen = generator.reads(ref_genome, 2, err_mod)
    big_read = ''.join(
        str(read_tuple[0].seq) + str(read_tuple[1].seq)
        for read_tuple in read_gen)
    assert big_read[140:170] == 'TTGAAACGGGTTGAAACGGGGTTTCAACCC'


def test_kde():
    random.seed(42)
    np.random.seed(42)
    err_mod = kde.KDErrorModel('data/ecoli_kde.npz')
    ref_genome = SeqRecord(
        Seq(str('CGTTTCAACC' * 40),
            IUPAC.unambiguous_dna
            ),
        id='my_genome',
        description='test genome'
        )
    read_gen = generator.reads(ref_genome, 2, err_mod)
    big_read = ''.join(
        str(read_tuple[0].seq) + str(read_tuple[1].seq)
        for read_tuple in read_gen)
    assert big_read[140:170] == 'TTGAAACGGGTTGAAACGGGGTTTCAACCC'
