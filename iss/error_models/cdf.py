#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division
from builtins import super

from iss import util
from iss.error_models import ErrorModel
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

import random
import numpy as np


class CDFErrorModel(ErrorModel):
    """CDFErrorModel class.

    Error model based on an .npz files derived from read alignments.
    the npz file must contain:

    - the length of the reads
    - the mean insert size
    - the distribution of qualities for each position (for R1 and R2) in form
        of weights and indices from a histogram
    - the substitution for each nucleotide at each position (for R1 and R2)
    - the insertion and deletion rates for each position (for R1 and R2)
    """
    def __init__(self, npz_path):
        super().__init__()
        self.npz_path = npz_path
        self.error_profile = self.load_npz(npz_path, 'cdf')

        self.read_length = self.error_profile['read_length']
        self.insert_size = self.error_profile['insert_size']

        self.quality_forward = self.error_profile['quality_hist_forward']
        self.quality_reverse = self.error_profile['quality_hist_reverse']

        self.subst_choices_for = self.error_profile['subst_choices_forward']
        self.subst_choices_rev = self.error_profile['subst_choices_forward']

        self.ins_for = self.error_profile['ins_forward']
        self.ins_rev = self.error_profile['ins_reverse']
        self.del_for = self.error_profile['del_forward']
        self.del_rev = self.error_profile['del_reverse']

    def gen_phred_scores(self, histograms):
        """Generate a list of phred scores based on a histogram from real data

        Args:
            histograms (ndarray): a list of weight, indices for each position
                of the read

        Returns:
            list: a list of phred scores
        """
        phred_list = []
        for w in histograms:
            random_quality = np.random.choice(
                w[0][1:], p=w[1]
            )
            phred_list.append(round(random_quality))
        return phred_list
