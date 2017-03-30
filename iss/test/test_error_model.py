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
