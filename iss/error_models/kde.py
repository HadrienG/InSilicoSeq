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

    def introduce_indels(self, record, orientation, full_sequence, bounds):
        """Introduce insertions or deletions in a sequence.
        Return a sequence"""
        full_sequence_start, full_sequence_end = bounds

        # get the right indel arrays
        if orientation == 'forward':
            insertions = self.ins_for
            deletions = self.del_for
        elif orientation == 'reverse':
            insertions = self.ins_rev
            deletions = self.del_rev
        else:
            print('this is bad')  # TODO error message and proper logging

        mutable_seq = record.seq.tomutable()
        position = 0
        for nucl in mutable_seq:
            for nucl_to_insert, prob in insertions[position].items():  # ins
                if random.random() < prob:
                    # we want to insert after the base read, hence position + 1
                    mutable_seq.insert(position + 1, nucl_to_insert)
            if random.random() < deletions[position][nucl]:  # del
                mutable_seq.pop(position)
            position += 1

        if len(mutable_seq) == self.read_length:
            return mutable_seq.toseq()
        elif len(mutable_seq) > self.read_length:
            while len(mutable_seq) > self.read_length:
                mutable_seq.pop()
            return mutable_seq.toseq()
        else:  # len smaller
            to_add = self.read_length - len(mutable_seq)
            if orientation == 'forward':
                for i in range(to_add):
                    nucl_to_add = full_sequence[full_sequence_end + i]
                    mutable_seq.append(nucl_to_add)
            elif orientation == 'reverse':
                for i in range(to_add):
                    nucl_to_add = util.rev_comp(
                        full_sequence[full_sequence_start - i]
                    )
                    mutable_seq.append(nucl_to_add)
            return mutable_seq.toseq()

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
