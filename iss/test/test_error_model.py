#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss.error_models import ErrorModel, basic, kde

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from nose.tools import raises

import random
import numpy as np


def test_basic_phred():
    np.random.seed(42)
    err_mod = basic.BasicErrorModel()

    distribution = err_mod.gen_phred_scores(20, 'forward')[:10]
    assert distribution == [23, 19, 25, 40, 19, 19, 40, 26, 18, 23]


def test_kde_phred():
    np.random.seed(42)
    err_mod = kde.KDErrorModel('data/ecoli.npz')
    distribution = err_mod.gen_phred_scores(err_mod.quality_reverse,
                                            'reverse')[10:]
    assert distribution == [40, 40, 40, 40, 40, 40, 40, 40, 10, 10]


def test_introduce_errors():
    np.random.seed(42)
    err_mod = basic.BasicErrorModel()

    read = SeqRecord(
        Seq(str('AATGC' * 25),
            IUPAC.unambiguous_dna
            ),
        id='read_1',
        description='test read'
        )
    read = err_mod.introduce_error_scores(read, 'forward')
    qualities = read.letter_annotations["phred_quality"][:10]
    assert qualities == [40, 26, 40, 40, 25, 25, 40, 40, 22, 40]


def test_mut_sequence():
    random.seed(42)
    np.random.seed(42)

    err_mod = basic.BasicErrorModel()

    read = SeqRecord(
        Seq(str('AAAAA' * 25),
            IUPAC.unambiguous_dna
            ),
        id='read_1',
        description='test read'
        )
    read.letter_annotations["phred_quality"] = [5] * 125
    read.seq = err_mod.mut_sequence(read, 'forward')
    assert str(read.seq[:10]) == 'AAAACAGAAA'


def test_introduce_indels():
        random.seed(42)
        np.random.seed(42)

        err_mod = basic.BasicErrorModel()
        err_mod.ins_for[1]['G'] = 1.0
        err_mod.del_for[0]['A'] = 1.0
        bounds = (5, 130)
        read = SeqRecord(
            Seq(str('ATATA' * 25),
                IUPAC.unambiguous_dna
                ),
            id='read_1',
            description='test read'
            )
        ref_genome = SeqRecord(
            Seq(str('ATATA' * 100),
                IUPAC.unambiguous_dna
                ),
            id='ref_genome',
            description='test reference'
            )
        read.seq = err_mod.introduce_indels(
            read, 'forward', ref_genome, bounds)
        assert len(read.seq) == 125
        assert read.seq[:10] == 'ATGATAATAT'


@raises(SystemExit)
def test_bad_err_mod():
    err_mod = kde.KDErrorModel('data/empty_file')
