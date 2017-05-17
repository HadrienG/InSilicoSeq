#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss.error_models import ErrorModel, basic, cdf, kde

import numpy as np


def test_basic_phred():
    np.random.seed(42)
    err_mod = basic.BasicErrorModel()

    distribution = err_mod.gen_phred_scores(20)[:10]
    assert distribution == [23, 19, 25, 40, 19, 19, 40, 26, 18, 23]


def test_cdf_phred():
    np.random.seed(42)
    err_mod = cdf.CDFErrorModel('data/ecoli_cdf.npz')
    distribution = err_mod.gen_phred_scores(err_mod.quality_reverse)[:10]
    assert distribution == [11, 21, 40, 40, 31, 31, 31, 40, 40, 40]


def test_kde_phred():
    np.random.seed(42)
    err_mod = kde.KDErrorModel('data/ecoli_kde.npz')
    distribution = err_mod.gen_phred_scores(err_mod.quality_reverse)[:10]
    assert distribution == [10, 20, 40, 40, 30, 30, 30, 40, 40, 40]
