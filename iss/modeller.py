#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy import stats

import numpy as np


def raw_qualities_to_histogram(qualities, model):
    """Generate numpy histograms for each position of the input mapped reads.
    A histogram contains the distribution of the phred scores for one position
    in all the reads. Returns a numpy array of the histograms for each position

    Arguments:
        bam_file (:obj:`str`): input bam file

    Returns:
        tuple: (histograms_forward, histograms_reverse)
    """
    if model == 'cdf':
        histograms = [np.histogram(
            i, bins=range(0, 41)) for i in zip(*qualities)]
        weights_list = []
        for hist in histograms:
            values, indices = hist
            weights = values / np.sum(values)
            weights_list.append((indices, weights))
        return weights_list
    elif model == 'kde':
        quals = [i for i in zip(*qualities)]
        cdfs_list = []
        for q in quals:
            kde = stats.gaussian_kde(q, bw_method=0.2 / np.std(q, ddof=1))
            kde = kde.evaluate(range(41))
            cdf = np.cumsum(kde)
            cdf = cdf / cdf[-1]
            cdfs_list.append(cdf)
        return cdfs_list


def get_mismatches(alignment):
    """Get substitution and indel rate for reads mapped to a reference genome.

    Arguments:
        bam_file (:obj:`str`): input bam file
        read_length (:obj"`int`): length of the mapped reads
    Returns:
        dict: for each nucleotide (keys) the values are a tuple containing the
        choices and probabilties of transiting to another nucleotide.
    """
    subst_array_f = np.zeros([read_length, 16])
    subst_array_r = np.zeros([read_length, 16])
    indel_array_f = np.zeros([read_length, 9])
    indel_array_r = np.zeros([read_length, 9])

    dispatch_subst = {
        'AA': 0,
        'aT': 1,
        'aG': 2,
        'aC': 3,
        'TT': 4,
        'tA': 5,
        'tG': 6,
        'tC': 7,
        'CC': 8,
        'cA': 9,
        'cT': 10,
        'cG': 11,
        'GG': 12,
        'gA': 13,
        'gT': 14,
        'gC': 15
    }

    has_indels = False
    for base in alignment:
        if read.seq[base[0]] != 'N':  # let's not deal with Ns
            query_pos = base[0]
            query_base = read.seq[query_pos]
            ref_base = base[2]
            dispatch_key = ref_base + query_base
            if dispatch_key not in dispatch_subst:
                # flag reads that have one or more indels
                has_indels = True
            if read.is_read1 and has_indels is False:
                subst_array_f[
                    query_pos,
                    dispatch_subst[dispatch_key]] += 1
            elif read.is_read2 and has_indels is False:
                subst_array_r[
                    query_pos,
                    dispatch_subst[dispatch_key]] += 1
