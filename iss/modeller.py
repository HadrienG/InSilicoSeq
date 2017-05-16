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


def dispatch_subst(base, read, read_has_indels):
    """Return the x and y position of a substitution to be inserted in the
    subst matrix.

    Arguments:
        base (:obj"`tuple`): one base from an aligmnent. alignment is a list of
            tuples: aligned read (query) and reference positions. the parameter
            with_seq adds the ref sequence as the 3rd element of the tuples.
            substitutions are lower-case.
    Returns:
        tuple: x and y position for incrementing the substitution matrix. and a
        third element: True if an indel has been detected, False otherwise
    """
    dispatch_dict = {
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

    query_pos = base[0]
    query_base = read.seq[query_pos]
    ref_base = base[2]
    dispatch_key = ref_base + query_base
    if dispatch_key not in dispatch_dict:
        # flag reads that have one or more indels
        read_has_indels = True  # flag the read for later indel treatment
        substitution = None  # flag this base to skip substitution treatment
    else:
        substitution = dispatch_dict[dispatch_key]
    return (query_pos, substitution, read_has_indels)


def subst_matrix_to_choices(substitution_matrix, read_length):
    """from the raw mismatches at one position, returns nucleotides
    and probabilties of state change (substitutions)"""
    nucl_choices_list = []
    for pos in range(read_length):
        sums = {
            'A': np.sum(substitution_matrix[pos][1:4]),
            'T': np.sum(substitution_matrix[pos][5:8]),
            'C': np.sum(substitution_matrix[pos][9:12]),
            'G': np.sum(substitution_matrix[pos][13:])
        }
        nucl_choices = {
            'A': (
                ['T', 'C', 'G'],
                [count / sums['A'] for count in substitution_matrix[pos][1:4]]
                ),
            'T': (
                ['A', 'C', 'G'],
                [count / sums['T'] for count in substitution_matrix[pos][5:8]]
                ),
            'C': (
                ['A', 'T', 'G'],
                [count / sums['C'] for count in substitution_matrix[pos][9:12]]
                ),
            'G': (
                ['A', 'T', 'C'],
                [count / sums['G'] for count in substitution_matrix[pos][13:]]
                )
        }
        nucl_choices_list.append(nucl_choices)
    return nucl_choices_list


def dispatch_indels(read):
    dispatch_indels = {
        0: 0,
        'A1': 1,
        'T1': 2,
        'C1': 3,
        'G1': 4,
        'A2': 5,
        'T2': 6,
        'C2': 7,
        'G2': 8
    }

    position = 0
    for (cigar_type, cigar_length) in read.cigartuples:
        if cigar_type == 0:  # match
            position += cigar_length
            continue
        elif cigar_type == 1:  # insertion
            query_base = read.query_sequence[position]
            insertion = query_base.upper() + '1'
            indel = dispatch_indels[insertion]
            dispatch_tuple = (position, indel)
            position += cigar_length
        elif cigar_type == 2:  # deletion
            ref_base = read.query_alignment_sequence[position]
            deletion = ref_base.upper() + '2'
            indel = dispatch_indels[deletion]
            dispatch_tuple = (position, indel)
            position -= cigar_length
        yield dispatch_tuple


def indel_matrix_to_choices(indel_matrix, read_length):
    """from the raw deletion rates at one position, returns nucleotides
    and probabilties of deletion"""
    ins_choices = []
    del_choices = []
    for pos in range(read_length):
        insertions = {
            'A': indel_matrix[pos][1] / indel_matrix[pos][0],
            'T': indel_matrix[pos][2] / indel_matrix[pos][0],
            'C': indel_matrix[pos][3] / indel_matrix[pos][0],
            'G': indel_matrix[pos][4] / indel_matrix[pos][0]
        }
        deletions = {
            'A': indel_matrix[pos][5] / indel_matrix[pos][0],
            'T': indel_matrix[pos][6] / indel_matrix[pos][0],
            'C': indel_matrix[pos][7] / indel_matrix[pos][0],
            'G': indel_matrix[pos][8] / indel_matrix[pos][0]
        }
        ins_choices.append(insertions)
        del_choices.append(deletions)
    return (ins_choices, del_choices)
