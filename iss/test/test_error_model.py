#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import error_model

import random
import numpy as np


def test_phred_conversions():
    assert error_model.phred_to_prob(40) == 0.9999
    assert error_model.phred_to_prob(30) == 0.999
    assert error_model.phred_to_prob(20) == 0.99
    assert error_model.phred_to_prob(10) == 0.9

    assert error_model.prob_to_phred(0.9999) == 40
    assert error_model.prob_to_phred(0.999) == 30
    assert error_model.prob_to_phred(0.99) == 20
    assert error_model.prob_to_phred(0.9) == 10


def test_basic_error_model():
    np.random.seed(42)
    err_mod = error_model.BasicErrorModel()

    distribution = err_mod.gen_phred_scores(0.99)[:10]
    assert distribution == [23, 19, 25, 40, 19, 19, 40, 26, 18, 23]
