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

    def load_npz(self, npz_path):
        """load the error profile npz file"""
        error_profile = np.load(npz_path)
        return error_profile

    def gen_phred_scores(self, qual_data):
        """Generate a list of phred scores based on real datasets"""
        phred_list = []
        for x in qual_data:
            kde = stats.gaussian_kde(x, bw_method=0.2 / x.std(ddof=1))
            kde = kde.evaluate(range(41))
            cdf = np.cumsum(kde)
            cdf = cdf / cdf[-1]
            random_quality = np.searchsorted(cdf, np.random.rand())
            phred_list.append(random_quality)
        return phred_list

    def introduce_error_scores(self, record, orientation):
        """Add phred scores to a SeqRecord according to the error_model"""
        if orientation == 'forward':
            record.letter_annotations["phred_quality"] = self.gen_phred_scores(
                self.quality_hist_forward)
        elif orientation == 'reverse':
            record.letter_annotations["phred_quality"] = self.gen_phred_scores(
                self.quality_hist_reverse)
        else:
            print('bad orientation. Fatal')  # add an exit here

        return record

    def subst_matrix_to_choices(self, subst_dispatch_dict):
        """from the raw substitutions at one position, returns nucleotides
        and probabilties of state change"""
        sums = {
            'A': sum(subst_dispatch_dict[1:4]),
            'T': sum(subst_dispatch_dict[5:8]),
            'C': sum(subst_dispatch_dict[9:12]),
            'G': sum(subst_dispatch_dict[13:])
        }

        nucl_choices = {
            'A': (
                ['T', 'C', 'G'],
                [count / sums['A'] for count in subst_dispatch_dict[1:4]]
                ),
            'T': (
                ['A', 'C', 'G'],
                [count / sums['T'] for count in subst_dispatch_dict[5:8]]
                ),
            'C': (
                ['A', 'T', 'G'],
                [count / sums['C'] for count in subst_dispatch_dict[9:12]]
                ),
            'G': (
                ['A', 'T', 'C'],
                [count / sums['G'] for count in subst_dispatch_dict[13:]]
                )
        }
        return nucl_choices

    def mut_sequence(self, record, orientation):
        # TODO
        """modify the nucleotides of a SeqRecord according to the phred scores.
        Return a sequence"""

        # get the right subst_matrix
        if orientation == 'forward':
            subst_matrix = self.subst_matrix_forward
        elif orientation == 'reverse':
            subst_matrix = self.subst_matrix_reverse
        else:
            print('this is bad')  # TODO error message and proper logging

        mutable_seq = record.seq.tomutable()
        quality_list = record.letter_annotations["phred_quality"]
        position = 0
        for nucl, qual in zip(mutable_seq, quality_list):
            nucl_choices = self.subst_matrix_to_choices(subst_matrix[position])
            if random.random() > util.phred_to_prob(qual):
                mutable_seq[position] = np.random.choice(
                    nucl_choices[nucl][0],
                    p=nucl_choices[nucl][1])
            position += 1
        return mutable_seq.toseq()
