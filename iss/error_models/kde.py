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

        self.quality_hist_for = self.error_profile['quality_hist_forward']
        self.quality_hist_rev = self.error_profile['quality_hist_reverse']

        self.subst_choices_for = self.error_profile['subst_choices_forward']
        self.subst_choices_rev = self.error_profile['subst_choices_reverse']

        self.indels_for = self.error_profile['indels_forward']
        self.indels_rev = self.error_profile['indels_reverse']

    def load_npz(self, npz_path):
        """load the error profile npz file"""
        error_profile = np.load(npz_path)
        return error_profile

    def gen_phred_scores(self, cdfs):
        """Generate a list of phred scores based on real datasets"""
        phred_list = []
        for cdf in cdfs:
            random_quality = np.searchsorted(cdf, np.random.rand())
            phred_list.append(random_quality)
        return phred_list

    def introduce_indels(self, record, orientation, full_sequence):
        """Introduce insertions or deletions in a sequence.
        Return a sequence"""
        # get the right indel array
        if orientation == 'forward':
            indels = self.indels_for
        elif orientation == 'reverse':
            indels = self.indels_rev
        else:
            print('this is bad')  # TODO error message and proper logging

        mutable_seq = record.seq.tomutable()
        position = 0
        for nucl in mutable_seq:
            if random.random() < indels[position][nucl][1][0]:  # insertions
                mutable_seq.insert(position, 'N')  # TODO shouldn't be Ns
            if random.random() < indels[position][nucl][1][1]:  # deletions
                mutable_seq.pop(position)
            position += 1

        if len(mutable_seq) == self.read_length:
            return mutable_seq.toseq()
        elif len(mutable_seq) > self.read_length:
            while len(mutable_seq) > self.read_length:
                mutable_seq.pop()
            return mutable_seq.toseq()
        else:  # len smaller
            while len(mutable_seq) < self.read_length:
                mutable_seq.append('N')  # TODO: shouldn't be Ns.
            return mutable_seq.toseq()

    def introduce_error_scores(self, record, orientation):
        """Add phred scores to a SeqRecord according to the error_model"""
        if orientation == 'forward':
            record.letter_annotations["phred_quality"] = self.gen_phred_scores(
                self.quality_hist_for)
        elif orientation == 'reverse':
            record.letter_annotations["phred_quality"] = self.gen_phred_scores(
                self.quality_hist_rev)
        else:
            print('bad orientation. Fatal')  # add an exit here

        return record

    def mut_sequence(self, record, orientation):
        """modify the nucleotides of a SeqRecord according to the phred scores.
        Return a sequence"""

        # get the right subst_matrix
        if orientation == 'forward':
            nucl_choices = self.subst_choices_for
        elif orientation == 'reverse':
            nucl_choices = self.subst_choices_rev
        else:
            print('this is bad')  # TODO error message and proper logging

        mutable_seq = record.seq.tomutable()
        quality_list = record.letter_annotations["phred_quality"]
        position = 0
        for nucl, qual in zip(mutable_seq, quality_list):
            if random.random() > util.phred_to_prob(qual):
                mutable_seq[position] = np.random.choice(
                    nucl_choices[position][nucl][0],
                    p=nucl_choices[position][nucl][1])
            position += 1
        return mutable_seq.toseq()
