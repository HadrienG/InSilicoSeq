#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss.error_models import ErrorModel, basic, cdf
from iss import util

import random
import numpy as np


def test_phred_conversions():
    assert util.phred_to_prob(40) == 0.9999
    assert util.phred_to_prob(30) == 0.999
    assert util.phred_to_prob(20) == 0.99
    assert util.phred_to_prob(10) == 0.9

    assert util.prob_to_phred(0.9999) == 40
    assert util.prob_to_phred(0.999) == 30
    assert util.prob_to_phred(0.99) == 20
    assert util.prob_to_phred(0.9) == 10


def test_basic_error_model():
    np.random.seed(42)
    err_mod = basic.BasicErrorModel()

    distribution = err_mod.gen_phred_scores(0.99)[:10]
    assert distribution == [23, 19, 25, 40, 19, 19, 40, 26, 18, 23]


def test_cdf_error_model_substitutions():
    err_mod = cdf.CDFErrorModel(
        npz_path='profiles/SRR5166376.npz'
        )
    dispatch = np.array(
        [2.77297000e+05, 2.58000000e+02, 3.22000000e+02, 7.37000000e+02,
            2.81467000e+05, 2.18000000e+02, 4.17000000e+02, 3.27000000e+02,
            2.87455000e+05, 2.55000000e+02, 1.96000000e+02, 2.65000000e+02,
            2.49238000e+05, 2.40000000e+02, 4.10000000e+02, 2.63000000e+02]
        )
    assert err_mod.subst_matrix_to_choices(dispatch)['A'] == (
        ['T', 'C', 'G'],
        [0.1958997722095672, 0.24449506454062261, 0.55960516324981013])
