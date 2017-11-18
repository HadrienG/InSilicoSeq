#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from builtins import dict, range, super

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
    equal between all nucleotides.
    """
    def __init__(self):
        super().__init__()
        self.read_length = 125
        self.insert_size = 200
        self.quality_forward = self.quality_reverse = 30
        self.subst_choices_for = self.subst_choices_rev = [{
            'A': (['T', 'C', 'G'], [1/3, 1/3, 1/3]),
            'T': (['A', 'C', 'G'], [1/3, 1/3, 1/3]),
            'C': (['A', 'T', 'G'], [1/3, 1/3, 1/3]),
            'G': (['A', 'T', 'C'], [1/3, 1/3, 1/3])
        } for _ in range(self.read_length)]

        self.ins_for = self.ins_rev = self.del_for = self.del_rev = [{
            'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0
        } for _ in range(125)]

    def gen_phred_scores(self, mean_quality, orientation):
        """Generate a normal distribution, transform to phred scores

        Generate a list of phred score according to a normal distribution
        centered around the ErrorModel quality

        Args:
            mean_quality (int): mean phred score

        Returns:
            list: list of phred scores following a normal distribution
        """
        norm = [min(q, 0.9999) for q in np.random.normal(
            util.phred_to_prob(mean_quality), 0.01, self.read_length)]
        phred = [util.prob_to_phred(p) for p in norm]
        return phred

    def random_insert_size(self):
        """Fake random function returning the default insert size of the
        basic arror model

        Returns:
            int: insert size
        """
        return self.insert_size
