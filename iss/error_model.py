#!/usr/bin/env python
# -*- coding: utf-8 -*-


import numpy as np


def phred_to_prob(q):
    """Given a phred score q, return the probabilty p
    of the call being RIGHT"""
    p = 10 ** (-q / 10)
    return 1 - p


def prob_to_phred(p):
    """Given the probablity p of a basecall being RIGHT
    return the phred score"""
    q = int(round(-10 * np.log10(1 - p)))
    return q


def introduce_errors(seq, mean_qual):
    seq.letter_annotations["phred_quality"] = basic(
        phred_to_prob(mean_qual), 0.01, len(seq))
    return seq


def basic(mean, stdev, length):
    """Generate a normal distribution, transform to phred scores"""
    # rate = 1
    norm = [min(q, 0.9999) for q in np.random.normal(mean, stdev, length)]
    # inverse transform sampling ? or lognormal distribution ?
    # uni = [min(q, 0.9999) for q in np.random.uniform()]
    # exp = [np.log(1 - q) / (- rate) for q in norm]
    phred = [prob_to_phred(p) for p in norm]
    return phred
