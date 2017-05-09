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


def rev_comp(s):
    """A simple reverse complement implementation working on strings"""
    bases = {
        "a": "t", "c": "g", "g": "c", "t": "a", "y": "r", "r": "y", "w": "w",
        "s": "s", "k": "m", "m": "k", "n": "n", "A": "T", "C": "G", "G": "C",
        "T": "A", "Y": "R", "R": "Y", "W": "W", "S": "S", "K": "M", "M": "K",
        "N": "N"}
    sequence = list(s)
    complement = "".join([bases[b] for b in sequence])
    reverse_complement = complement[::-1]
    return reverse_complement
