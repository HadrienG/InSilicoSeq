#!/usr/bin/env python
# -*- coding: utf-8 -*-


class ErrorModel(object):
    """Main ErrorModel Class

    This class is used to create inheriting classes
    """
    def __init__(self):
        self.read_length = int
        self.insert_size = int

    def load_npz(self, npz_path):
        """load the error profile npz file"""
        error_profile = np.load(npz_path)
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
