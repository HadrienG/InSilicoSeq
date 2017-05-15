#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util
from iss.error_models import ErrorModel
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

import random
import numpy as np


class BasicErrorModel(ErrorModel):
    """Basic Error Model class

    Basic error model. The phred scores are based on a normal distribution.
    Only substitutions errors occur. The substitution rate is assumed
    equal between all nucleotides."""
    def __init__(self):
        super().__init__()
        self.read_length = 125
        self.insert_size = 200
        self.quality_forward, self.quality_reverse = 30
        self.subst_choices_for, self.subst_choices_rev = [{
            'A': (['T', 'C', 'G'], [1/3, 1/3, 1/3]),
            'T': (['A', 'C', 'G'], [1/3, 1/3, 1/3]),
            'C': (['A', 'T', 'G'], [1/3, 1/3, 1/3]),
            'G': (['A', 'T', 'C'], [1/3, 1/3, 1/3])
        } for _ in range(125)]

    def gen_phred_scores(self, mean_quality):
        """Generate a normal distribution, transform to phred scores"""
        norm = [min(q, 0.9999) for q in np.random.normal(
            util.phred_to_prob(mean_quality), 0.01, self.read_length)]
        phred = [util.prob_to_phred(p) for p in norm]
        return phred
