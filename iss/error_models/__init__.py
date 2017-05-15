#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import sys


class ErrorModel(object):
    """Main ErrorModel Class

    This class is used to create inheriting classes
    """
    def load_npz(self, npz_path):
        """load the error profile npz file"""
        try:
            error_profile = np.load(npz_path)
        except OSError as e:
            print('Error:', e)
            sys.exit(1)
        return error_profile

    def introduce_error_scores(self, record, orientation):
        """Add phred scores to a SeqRecord according to the error_model"""
        if orientation == 'forward':
            record.letter_annotations["phred_quality"] = self.gen_phred_scores(
                self.quality_forward)
        elif orientation == 'reverse':
            record.letter_annotations["phred_quality"] = self.gen_phred_scores(
                self.quality_reverse)
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
