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
        self.quality_forward = 30
        self.quality_reverse = 30

    def gen_phred_scores(self, mean_quality):
        """Generate a normal distribution, transform to phred scores"""
        norm = [min(q, 0.9999) for q in np.random.normal(
            util.phred_to_prob(mean_quality), 0.01, self.read_length)]
        phred = [util.prob_to_phred(p) for p in norm]
        return phred

    def mut_sequence(self, record, orientation):
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
            if random.random() > util.phred_to_prob(qual):
                mutable_seq[position] = random.choice(nucl_choices[nucl])
            position += 1
        return mutable_seq.toseq()
