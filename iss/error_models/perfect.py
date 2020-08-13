#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util
from iss.error_models import ErrorModel


class PerfectErrorModel(ErrorModel):
    """Perfect error model class

    Perfect Error Model. This is a model without errors.
    All Phred score are 40. No errors are introduced at all.
    """

    def __init__(self):
        super().__init__()
        self.read_length = 125
        self.insert_size = 200
        self.quality_forward = self.quality_reverse = 40

        self.subst_choices_for = self.subst_choices_rev = [{
            'A': (['A', 'T', 'C', 'G'], [1, 0, 0, 0]),
            'T': (['A', 'T', 'C', 'G'], [0, 1, 0, 0]),
            'C': (['A', 'T', 'C', 'G'], [0, 0, 1, 0]),
            'G': (['A', 'T', 'C', 'G'], [0, 0, 0, 1])
        } for _ in range(self.read_length)]

        self.ins_for = self.ins_rev = self.del_for = self.del_rev = [{
            'A': 0.0, 'T': 0.0, 'C': 0.0, 'G': 0.0
        } for _ in range(self.read_length)]

    def gen_phred_scores(self, mean_quality, orientation):
        """Fake randorm function returning the distribution of Phred
        scores. Score will be 40 for all positions

        Returns:
            list: list of phred scores (40 along the whole read)
        """
        return [40 for _ in range(self.read_length)]

    def random_insert_size(self):
        """Fake random function returning the default insert size of the
        basic arror model

        Returns:
            int: insert size
        """
        return self.insert_size
