#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss.error_models import ErrorModel, basic, cdf

import numpy as np


def test_basic_error_model():
    np.random.seed(42)
    err_mod = basic.BasicErrorModel()

    distribution = err_mod.gen_phred_scores(0.99)[:10]
    assert distribution == [23, 19, 25, 40, 19, 19, 40, 26, 18, 23]
