#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord

import random
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


class ErrorModel(object):
    """Main ErrorModel Class

    This class is used to create inheriting classes
    """
    def __init__(self):
        self.read_length = int
        self.insert_size = int


class BasicErrorModel(ErrorModel):
    """Basic Error Model class

    Basic error model. The phred scores are based on a normal distribution.
    Only substitutions errors occur. The substitution rate is assumed
    equal between all nucleotides."""
    def __init__(self):
        super().__init__()
        self.read_length = 125
        self.insert_size = 200
        self.quality = 30

    def gen_phred_scores(self, mean_quality):
        """Generate a normal distribution, transform to phred scores"""
        norm = [min(q, 0.9999) for q in np.random.normal(
            mean_quality, 0.01, self.read_length)]
        phred = [prob_to_phred(p) for p in norm]
        return phred

    def introduce_error_scores(self, record):
        """Add phred scores to a SeqRecord"""
        record.letter_annotations["phred_quality"] = self.gen_phred_scores(
            phred_to_prob(self.quality))
        return record

    def mut_sequence(self, record):
        """modify the nucleotides of a SeqRecord according to the phred scores.
        Return a sequence"""
        nucl_choices = {
            'A': ['T', 'C', 'G'],
            'T': ['A', 'C', 'G'],
            'C': ['A', 'T', 'G'],
            'G': ['A', 'T', 'C']
            }
        mutable_seq = record.seq.tomutable()
        quality_list = record.letter_annotations["phred_quality"]
        position = 0
        for nucl, qual in zip(mutable_seq, quality_list):
            if random.random() > phred_to_prob(qual):
                mutable_seq[position] = random.choice(nucl_choices[nucl])
            position += 1
        return mutable_seq.toseq()


class KernelDensityErrorModel(ErrorModel):
    """KernelDensityErrorModel class.

    Error model based on .npz files derived from alignment with bowtie2.
    the npz file must contain:

    - the length of the reads
    - the mean insert size
    - the distribution of qualities for each position (for R1 and R2)
    - the substitution for each nucleotide at each position (for R1 and R2)"""
    def __init__(self, npz_path):
        super().__init__()
        self.npz_path = npz_path
        self.error_profile = self.load_npz(npz_path)

        self.read_length = error_profile['read_length']
        self.insert_size = error_profile['insert_size']

        self.quality_hist_forward = error_profile['quality_hist_forward']
        self.quality_hist_reverse = error_profile['quality_hist_reverse']

        self.subst_matrix_forward = error_profile['subst_matrix_forward']
        self.subst_matrix_reverse = error_profile['subst_matrix_reverse']

    def load_npz(self, npz_path):
        """load the error profile npz file"""
        error_profile = np.load(numpy_file)
        return error_profile

    def gen_phred_scores(self, histograms):
        """Generate a list of phred scores based on real datasets"""
        phred_list = []
        for hist in histograms:
            values, indices = hist
            weights = values / np.sum(values)
            random_quality = np.random.choice(
                indices[1:], p=weights
            )
            phred_list.append(round(random_quality))
        return phred_list

    def introduce_error_scores(self, record, histograms, orientation):
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

    def mut_sequence(self, record, subst_matrix, orientation):
        # TODO
        """modify the nucleotides of a SeqRecord according to the phred scores.
        Return a sequence"""
        nucl_choices = {
            'A': ['T', 'C', 'G'],
            'T': ['A', 'C', 'G'],
            'C': ['A', 'T', 'G'],
            'G': ['A', 'T', 'C']
            }
        mutable_seq = record.seq.tomutable()
        quality_list = record.letter_annotations["phred_quality"]
        position = 0
        for nucl, qual in zip(mutable_seq, quality_list):
            if random.random() > phred_to_prob(qual):
                mutable_seq[position] = random.choice(nucl_choices[nucl])
            position += 1
        return mutable_seq.toseq()
