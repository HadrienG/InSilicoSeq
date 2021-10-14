#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util
from iss.error_models import ErrorModel

from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from scipy import stats

import sys
import random
import numpy as np


class KDErrorModel(ErrorModel):
    """KDErrorModel class.

    Error model based on an .npz files derived from read alignments.
    the npz file must contain:

    - the length of the reads
    - the mean insert size
    - the size of mean sequence quality bins (for R1 and R2)
    - a cumulative distribution function of quality scores for each position
        (for R1 and R2)
    - the substitution for each nucleotide at each position (for R1 and R2)
    - the insertion and deletion rates for each position (for R1 and R2)
    """

    def __init__(self, npz_path):
        super().__init__()
        self.npz_path = npz_path
        self.error_profile = self.load_npz(npz_path, 'kde')

        self.read_length = self.error_profile['read_length']
        self.i_size_cdf = self.error_profile['insert_size']

        self.mean_forward = self.error_profile['mean_count_forward']
        self.mean_reverse = self.error_profile['mean_count_reverse']

        self.quality_forward = self.error_profile['quality_hist_forward']
        self.quality_reverse = self.error_profile['quality_hist_reverse']

        self.subst_choices_for = self.error_profile['subst_choices_forward']
        self.subst_choices_rev = self.error_profile['subst_choices_reverse']

        self.ins_for = self.error_profile['ins_forward']
        self.ins_rev = self.error_profile['ins_reverse']
        self.del_for = self.error_profile['del_forward']
        self.del_rev = self.error_profile['del_reverse']

        self.error_profile = ''  # discard the error profile after reading

    def gen_phred_scores(self, cdfs, orientation, quality_bin):
        """Generate a list of phred scores based on cdfs and mean bins

        For each position, draw a phred score from the cdf and append to the
        phred score list

        Args:
            cdfs (ndarray): array containing the cdfs
            orientation (string): orientation of the read. Can be 'forward' or
                'reverse'
            quality_bin (string): level of quality for the quality scores. Can be 'auto',
                'low', 'middle_low', 'middle_high', 'high' (default: 'auto').

        Returns:
            list: a list of phred scores
        """
        # choose which mean bin to use
        if orientation == 'forward':
            mean = self.mean_forward
        elif orientation == 'reverse':
            mean = self.mean_reverse

        bin_types = ['low', 'middle_low', 'middle_high', 'high']
        # picks the quality bin based on the user preference
        # if quality_bin == auto, picks the bin based on the probability
        if quality_bin == 'auto':
            norm_mean = mean / sum(mean)
            # quality_bin = np.searchsorted(norm_mean, np.random.rand())
            quality_bin = np.random.choice(range(len(norm_mean)), p=norm_mean)
            # downgrades index out of bound (ex rand is 1, last bin in searchsorted
            # is 0.95) to best quality bin
            if quality_bin == 4:
                quality_bin = 3
        # checks if the quality bin is in one of the bin types
        elif quality_bin in bin_types:
            # lists the bins that are not empty
            full_bins = [ind for ind, bin in enumerate(bin_types) if len(cdfs[ind]) != 0]

            quality_bin = bin_types.index(quality_bin)

            #checks if the quality bin is not empty.
            if quality_bin not in full_bins:
                quality_bin = full_bins[quality_bin]

                if not full_bins:
                    raise ValueError('Error: the model doesn\'t contain any cdfs.')

                start_bin = quality_bin
                #if bin is empty, look for the first available bin of lower quality
                #in the case of the starting bin being 0, will look for the first available bin of higher quality
                while len(cdfs[quality_bin]) == 0:
                    if start_bin == 0:
                        quality_bin += 1
                    #if the lowest possible bin is reached, cycle back to the highest possible bin
                    elif quality_bin == 0:
                        quality_bin == 3
                    else:
                        quality_bin -= 1

#        print(f'The quality scores used are of {bin_types[quality_bin]} quality.')

        cdfs_bin = cdfs[quality_bin]

        phred_list = []
        for cdf in cdfs_bin:
            random_quality = np.searchsorted(cdf, np.random.rand())
            phred_list.append(random_quality)
        return phred_list[:self.read_length]

    def random_insert_size(self):
        """Draw a random insert size from the insert size cdf

        Args:
            i_size_cdf: cumulative distribution function of the insert size

        Returns:
            int: an insert size
        """
        insert_size = np.searchsorted(self.i_size_cdf, np.random.rand())
        return insert_size
