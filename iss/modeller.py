#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util
from scipy import stats
from joblib import Parallel, delayed

import logging
import numpy as np


def insert_size(insert_size_distribution):
    """Calculate cumulative distribution function from the raw insert size
    distributin. Uses 1D kernel density estimation.

    Args:
        insert_size_distribution (list): list of insert sizes from aligned
        read pairs

    Returns:
        1darray: a cumulative density function
    """
    kde = stats.gaussian_kde(
        insert_size_distribution,
        bw_method=0.2 / np.std(insert_size_distribution, ddof=1))
    x_grid = np.linspace(
        min(insert_size_distribution),
        max(insert_size_distribution), 1000)
    kde = kde.evaluate(x_grid)
    cdf = np.cumsum(kde)
    cdf = cdf / cdf[-1]
    return cdf


def divide_qualities_into_bins(qualities, n_bins=4):
    """Divides the raw quality scores in bins according to the mean phred
    quality of the sequence they come from

    Args:
        qualities (list): raw count of all the phred scores and mean sequence
            quality
        n_bins (int): number of bins to create (default: 4)

    Returns:
        list: a list of lists containing the binned quality scores
    """
    logger = logging.getLogger(__name__)
    logger.debug('Dividing qualities into mean clusters')
    bin_lists = [[] for _ in range(n_bins)]  # create list of `n_bins` list
    ranges = np.split(np.array(range(40)), n_bins)
    for quality in qualities:
        mean = int(quality[0][1])  # mean is at 1 and same regardless of b pos
        which_array = 0
        for array in ranges:
            if mean in array:
                read = np.fromiter((q[0] for q in quality), dtype=np.float)
                bin_lists[which_array].append(read)
            which_array += 1
    return bin_lists


def quality_bins_to_histogram(bin_lists):
    """Wrapper function to generate cdfs for each quality bins

    Generate cumulative distribution functions for a number of mean quality
    bins

    Args:
        bins_lists (list): list of list containing raw count of all phred
        scores

    Returns:
        list: a list of lists containg cumulative density functions
    """
    logger = logging.getLogger(__name__)
    cdf_bins = []
    i = 0
    for qual_bin in bin_lists:
        if len(qual_bin) > 1:
            logger.debug('Transposing matrix for mean cluster #%s' % i)
            # quals = np.asarray(qual_bin).T  # seems to make clunkier models
            quals = [q for q in zip(*qual_bin)]
            logger.debug(
                'Modelling quality distribution for mean cluster #%s' % i)
            cdfs_list = raw_qualities_to_histogram(quals)
            cdf_bins.append(cdfs_list)
        else:
            logger.debug('Mean quality bin #%s of length < 1. Skipping' % i)
            cdf_bins.append([])
        i += 1
    return cdf_bins


def raw_qualities_to_histogram(qualities):
    """Approximate the distribution of base quality at each position in a read
    using a pseudo 2d kernel density estimation

    Generate cumulative distribution functions

    Args:
        qualities (list): raw count of all phred scores

    Returns:
        list: list of cumulative distribution functions. One cdf per base. The
            list has the size of the read length
    """
    logger = logging.getLogger(__name__)

    # moved this in quality_bins_to_histogram to try parallelization
    # quals = util.split_list([i for i in zip(*qualities)], n_parts=cpus)
    cdfs_list = []
    for q in qualities:
        numpy_log_handler = np.seterrcall(util.nplog)
        with np.errstate(under='ignore', divide='call'):
            try:
                kde = stats.gaussian_kde(q, bw_method=0.2 / np.std(q, ddof=1))
            except np.linalg.linalg.LinAlgError as e:
                # if np.std of array is 0, we modify the array slightly to be
                # able to divide by ~ np.std
                # this will print a FloatingPointError in DEBUG mode
                # logger.debug('np.std of quality array is zero: %s' % e)
                q = list(q)
                q[-1] += 1
                kde = stats.gaussian_kde(q, bw_method=0.2 / np.std(q, ddof=1))
            kde = kde.evaluate(range(41))
        cdf = np.cumsum(kde)
        cdf = cdf / cdf[-1]
        cdfs_list.append(cdf)
    return cdfs_list


def dispatch_subst(base, read, read_has_indels):
    """Return the x and y position of a substitution to be inserted in the
    substitution matrix.

    The substitution matrix is a 2D array of size 301 * 16
    The x axis (301) corresponds to the position in the read, while
    the y axis (16) represents the match or substitution (see the dispatch
    dict in the function). Positions 0, 4, 8 and 12 are matches, other
    positions are substitutions

    The size of x axis is 301 because we haven't calculated the read length yet

    Args:
        base (tuple): one base from an aligmnent object. According to the
            pysam documentation: an alignment is a list of tuples: aligned read
            (query) and reference positions. the parameter with_seq adds the
            ref sequence as the 3rd element of the tuples.
            substitutions are lower-case.
        read (read): a read object, from which the alignment comes from
        read_has_indels (bool): a boolean flag to keep track if the read has
            an indel or not

    Returns:
        tuple: x and y position for incrementing the substitution matrix and a
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
    """Transform a substitution matrix into probabilties of substitutions for
    each base and at every position

    From the raw mismatches at one position, returns a dictionary with
    probabilties of substitutions

    Args:
        substitution_matrix (np.array): the substitution matrix is a 2D array
            of size read_length * 16. fhe x axis (read_length) corresponds to
            the position in the read, while the y axis (16) represents the
            match or substitution. Positions 0, 4, 8 and 12 are matches, other
            positions are substitutions
        read_length (int): read length

    Returns:
        list: list of dictionaries representing
            the substitution probabilities for a collection of reads
    """
    logger = logging.getLogger(__name__)

    nucl_choices_list = []
    for pos in range(read_length):
        sums = {
            'A': np.sum(substitution_matrix[pos][1:4]),
            'T': np.sum(substitution_matrix[pos][5:8]),
            'C': np.sum(substitution_matrix[pos][9:12]),
            'G': np.sum(substitution_matrix[pos][13:])
        }
        # we want to avoid 'na' in the data so we raise FloatingPointError
        # if we try to divide by 0 (no count data for that nucl at that pos)
        # we assume equal rate of substitution
        with np.errstate(all='raise'):
            nucl_choices = {}
            try:
                A = (
                    ['T', 'C', 'G'],
                    [count / sums['A'] for
                        count in substitution_matrix[pos][1:4]])
            except FloatingPointError as e:
                logger.debug(e, exc_info=True)
                A = (['T', 'C', 'G'], [1/3, 1/3, 1/3])
            try:
                T = (
                    ['A', 'C', 'G'],
                    [count / sums['T'] for
                        count in substitution_matrix[pos][5:8]])
            except FloatingPointError as e:
                logger.debug(e, exc_info=True)
                T = (['A', 'C', 'G'], [1/3, 1/3, 1/3])
            try:
                C = (
                    ['A', 'T', 'G'],
                    [count / sums['C'] for
                        count in substitution_matrix[pos][9:12]])
            except FloatingPointError as e:
                logger.debug(e, exc_info=True)
                C = (['A', 'T', 'G'], [1/3, 1/3, 1/3])
            try:
                G = (
                    ['A', 'T', 'C'],
                    [count / sums['G'] for
                        count in substitution_matrix[pos][13:]])
            except FloatingPointError as e:
                logger.debug(e, exc_info=True)
                G = (['A', 'T', 'C'], [1/3, 1/3, 1/3])

            nucl_choices['A'] = A
            nucl_choices['T'] = T
            nucl_choices['C'] = C
            nucl_choices['G'] = G
        nucl_choices_list.append(nucl_choices)
    return nucl_choices_list


def dispatch_indels(read):
    """Return the x and y position of a insertion or deletion to be inserted in
    the indel matrix.

    The substitution matrix is a 2D array of size 301 * 9
    The x axis (301) corresponds to the position in the read, while
    the y axis (9) represents the match or indel (see the dispatch
    dict in the function). Positions 0 is match or substitution, other
    positions in 'N1' are insertions, 'N2 are deletions'

    The size of x axis is 301 because we haven't calculated the read length yet

    Args:
        read (read): an aligned read object

    Yields:
        tuple: a tuple with the x, y position for dispatching the indel in the
        indel matrix
    """
    logger = logging.getLogger(__name__)

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
            try:
                indel = dispatch_indels[insertion]
                dispatch_tuple = (position, indel)
                position += cigar_length
            except KeyError as e:  # we avoid ambiguous bases
                # logger.debug(
                #     '%s not in dispatch: %s' % (insertion, e), exc_info=True)
                position += cigar_length
                continue
        elif cigar_type == 2:  # deletion
            ref_base = read.query_alignment_sequence[position]
            deletion = ref_base.upper() + '2'
            try:
                indel = dispatch_indels[deletion]
                dispatch_tuple = (position, indel)
                position -= cigar_length
            except KeyError as e:  # we avoid ambiguous bases
                # logger.debug(
                #     '%s not in dispatch: %s' % (deletion, e), exc_info=True)
                position -= cigar_length
                continue
        else:
            logger.debug("CIGAR %s. Skipping read." % cigar_type)
            continue
        yield dispatch_tuple


def indel_matrix_to_choices(indel_matrix, read_length):
    """Transform an indel matrix into probabilties of indels for
    at every position

    From the raw indel count at one position, returns a dictionary with
    probabilties of indel

    Args:
        indel_matrix (np.array): the substitution matrix is a 2D array of
            size read_length * 16. fhe x axis (read_length) corresponds to the
            position in the read, while the y axis (9) represents the match or
            indel. Positions 0 is match or substitution, other positions in
            'N1' are insertions, 'N2 are deletions'
        read_length (int): read length

    Returns:
        tuple: tuple containing two lists of dictionaries representing the
        insertion or deletion probabilities for a collection of reads
    """
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
