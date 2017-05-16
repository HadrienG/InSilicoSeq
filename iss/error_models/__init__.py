#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util

import sys
import logging
import random
import numpy as np


class ErrorModel(object):
    """Main ErrorModel Class

    This class is used to create inheriting classes
    """
    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def load_npz(self, npz_path, model):
        """load the error profile npz file"""
        try:
            error_profile = np.load(npz_path)
            assert error_profile['model'] == model
        except OSError as e:
            self.logger.error('Failed to read ErrorModel file: %s' % e)
            sys.exit(1)
        except AssertionError as e:
            self.logger.error(
                'Trying to load a %s ErrorModel in %s mode' % (
                    error_profile['model'], model))
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
        return record

    def mut_sequence(self, record, orientation):
        """modify the nucleotides of a SeqRecord according to the phred scores.
        Return a sequence"""

        # get the right subst_matrix
        if orientation == 'forward':
            nucl_choices = self.subst_choices_for
        elif orientation == 'reverse':
            nucl_choices = self.subst_choices_rev

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

    def adjust_seq_length(self, mut_seq, orientation, full_sequence, bounds):
        """take a mutable sequence after indel treatment. Return a seq of
        the correct read length"""
        read_start, read_end = bounds
        if len(mut_seq) == self.read_length:
            return mut_seq.toseq()
        elif len(mut_seq) > self.read_length:
            while len(mut_seq) > self.read_length:
                mut_seq.pop()
            return mut_seq.toseq()
        else:  # len smaller
            to_add = self.read_length - len(mut_seq)
            if orientation == 'forward':
                for i in range(to_add):
                    nucl_to_add = str(full_sequence[read_end + i])
                    mut_seq.append(nucl_to_add)
            elif orientation == 'reverse':
                for i in range(to_add):
                    nucl_to_add = util.rev_comp(
                        full_sequence[read_end + i]
                    )
                    mut_seq.append(nucl_to_add)
            return mut_seq.toseq()

    def introduce_indels(self, record, orientation, full_seq, bounds):
        """Introduce insertions or deletions in a sequence.
        Return a sequence"""

        # get the right indel arrays
        if orientation == 'forward':
            insertions = self.ins_for
            deletions = self.del_for
        elif orientation == 'reverse':
            insertions = self.ins_rev
            deletions = self.del_rev

        mutable_seq = record.seq.tomutable()
        position = 0
        for nucl in range(self.read_length - 1):
            for nucl_to_insert, prob in insertions[position].items():  # ins
                if random.random() < prob:
                    # we want to insert after the base read, hence position + 1
                    mutable_seq.insert(position + 1, nucl_to_insert)
            if random.random() < deletions[position][mutable_seq[nucl]]:  # del
                mutable_seq.pop(position)
            position += 1

        seq = self.adjust_seq_length(
            mutable_seq, orientation, full_seq, bounds)
        return seq
