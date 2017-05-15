#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util
from iss.error_models import ErrorModel
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from scipy import stats

import random
import numpy as np


class KDErrorModel(ErrorModel):
    """KDErrorModel class

    Error model based on .npz files derived from alignment with bowtie2.
    the npz file must contain:

    - the length of the reads
    - the mean insert size
    - the raw counts of qualities (for R1 and R2)
    - the substitution for each nucleotide at each position (for R1 and R2)"""
    def __init__(self, npz_path):
        super().__init__()
        self.npz_path = npz_path
        self.error_profile = self.load_npz(npz_path)

        self.read_length = self.error_profile['read_length']
        self.insert_size = self.error_profile['insert_size']

        self.quality_forward = self.error_profile['quality_hist_forward']
        self.quality_reverse = self.error_profile['quality_hist_reverse']

        self.subst_choices_for = self.error_profile['subst_choices_forward']
        self.subst_choices_rev = self.error_profile['subst_choices_reverse']

        self.ins_for = self.error_profile['ins_forward']
        self.ins_rev = self.error_profile['ins_reverse']
        self.del_for = self.error_profile['del_forward']
        self.del_rev = self.error_profile['del_reverse']

    def gen_phred_scores(self, cdfs):
        """Generate a list of phred scores based on real datasets"""
        phred_list = []
        for cdf in cdfs:
            random_quality = np.searchsorted(cdf, np.random.rand())
            phred_list.append(random_quality)
        return phred_list
