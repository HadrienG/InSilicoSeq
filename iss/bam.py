#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy import stats

import pysam
import numpy as np


def get_mismatches(bam_file, read_length):
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

    nucl_choices_f = []
    nucl_choices_r = []
    ins_f = []
    ins_r = []
    del_f = []
    del_r = []

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

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not read.is_unmapped:
                # alignment is a list of tuples: aligned read (query) and
                # reference positions. the parameter with_seq adds the ref
                # sequence as the 3rd element of the tuples.
                # substitutions are lower-case.
                alignment = read.get_aligned_pairs(
                    matches_only=True,
                    with_seq=True
                    )
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
                # once we've counted the substitutions, we count the indels
                # looking at the cigar
                if has_indels == 1:
                    position = 0
                    for (cigar_type, cigar_length) in read.cigartuples:
                        if cigar_type == 0:  # match
                            position += cigar_length
                        elif cigar_type == 1:  # insertion
                            insertion = query_base.upper() + '1'
                            if read.is_read1:
                                indel_array_f[
                                    position,
                                    dispatch_indels[insertion]] += 1
                            elif read.is_read2:
                                indel_array_r[
                                    position,
                                    dispatch_indels[insertion]] += 1
                            position += cigar_length
                        elif cigar_type == 2:  # deletion
                            deletion = ref_base.upper() + '2'
                            if read.is_read1:
                                indel_array_f[
                                    position,
                                    dispatch_indels[deletion]] += 1
                            elif read.is_read2:
                                indel_array_r[
                                    position,
                                    dispatch_indels[deletion]] += 1
                            position -= cigar_length
                        else:
                            print('error')

    for position in range(read_length):
        # update base count in indel arrays
        indel_array_f[position][0] = sum(subst_array_f[position][::4])
        indel_array_r[position][0] = sum(subst_array_r[position][::4])

        nucl_choices_f.append(subst_matrix_to_choices(subst_array_f[position]))
        nucl_choices_r.append(subst_matrix_to_choices(subst_array_r[position]))
        ins_f.append(indel_rate(indel_array_f[position])[0])
        ins_r.append(indel_rate(indel_array_r[position])[0])
        del_f.append(indel_rate(indel_array_f[position])[1])
        del_r.append(indel_rate(indel_array_r[position])[1])

    return nucl_choices_f, nucl_choices_r, ins_f, ins_r, del_f, del_r


def subst_matrix_to_choices(substitutions_array):
    """from the raw mismatches at one position, returns nucleotides
    and probabilties of state change (substitutions)"""
    sums = {
        'A': np.sum(substitutions_array[1:4]),
        'T': np.sum(substitutions_array[5:8]),
        'C': np.sum(substitutions_array[9:12]),
        'G': np.sum(substitutions_array[13:])
    }
    nucl_choices = {
        'A': (
            ['T', 'C', 'G'],
            [count / sums['A'] for count in substitutions_array[1:4]]
            ),
        'T': (
            ['A', 'C', 'G'],
            [count / sums['T'] for count in substitutions_array[5:8]]
            ),
        'C': (
            ['A', 'T', 'G'],
            [count / sums['C'] for count in substitutions_array[9:12]]
            ),
        'G': (
            ['A', 'T', 'C'],
            [count / sums['G'] for count in substitutions_array[13:]]
            )
    }
    return nucl_choices


def indel_rate(indels_array):
    """from the raw deletion rates at one position, returns nucleotides
    and probabilties of deletion"""
    insertions = {
        'A': indels_array[1] / indels_array[0],
        'T': indels_array[2] / indels_array[0],
        'C': indels_array[3] / indels_array[0],
        'G': indels_array[4] / indels_array[0]
    }
    deletions = {
        'A': indels_array[5] / indels_array[0],
        'T': indels_array[6] / indels_array[0],
        'C': indels_array[7] / indels_array[0],
        'G': indels_array[8] / indels_array[0]
    }
    return insertions, deletions


def quality_distribution(model, bam_file):
    """Generate numpy histograms for each position of the input mapped reads.
    A histogram contains the distribution of the phred scores for one position
    in all the reads. Returns a numpy array of the histograms for each position

    Arguments:
        bam_file (:obj:`str`): input bam file

    Returns:
        tuple: (histograms_forward, histograms_reverse)
    """
    # deal with the forward reads
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        array_gen_f = (np.array(
            read.query_qualities) for read in bam.fetch()
                if not read.is_unmapped and read.is_read1)
        if model == 'cdf':
            histograms_forward = [np.histogram(
                i, bins=range(0, 41)) for i in zip(*array_gen_f)]
        elif model == 'kde':
            quals_forward = [i for i in zip(*array_gen_f)]

    # deal with the reverse reads
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        array_gen_r = (np.array(
            read.query_qualities) for read in bam.fetch()
                if not read.is_unmapped and read.is_read2)
        if model == 'cdf':
            histograms_reverse = [np.histogram(
                i, bins=range(0, 41)) for i in zip(*array_gen_r)]
        elif model == 'kde':
            quals_reverse = [i for i in zip(*array_gen_r)]

    # calculate weights and indices
    if model == 'cdf':
        weights_forward = []
        weights_reverse = []
        for hist in histograms_forward:
            values, indices = hist
            weights = values / np.sum(values)
            weights_forward.append((indices, weights))
        for hist in histograms_reverse:
            values, indices = hist
            weights = values / np.sum(values)
            weights_reverse.append((indices, weights))
        return weights_forward, weights_reverse

    if model == 'kde':
        cdfs_forward = []
        cdfs_reverse = []
        for x in quals_forward:
            # print(x)
            kde = stats.gaussian_kde(x, bw_method=0.2 / np.std(x, ddof=1))
            kde = kde.evaluate(range(41))
            cdf = np.cumsum(kde)
            cdf = cdf / cdf[-1]
            cdfs_forward.append(cdf)
        for x in quals_reverse:
            kde = stats.gaussian_kde(x, bw_method=0.2 / np.std(x, ddof=1))
            kde = kde.evaluate(range(41))
            cdf = np.cumsum(kde)
            cdf = cdf / cdf[-1]
            cdfs_reverse.append(cdf)
        return cdfs_forward, cdfs_reverse


def get_insert_size(bam_file):
    """Get the mean insert size give mapped reads in a bam file

    Arguments:
        bam_file (:obj:`str`): input bam file

    Returns:
        int: mean insert size"""
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        i_size_dist = [
            abs(read.template_length) for read in bam.fetch()
            if not read.is_unmapped and read.is_proper_pair]
        i_size = np.mean(i_size_dist)
    return int(i_size)


def write_to_file(read_length, hist_f, hist_r, sub_f, sub_r, ins_f,
                  ins_r, del_f, del_r, i_size, output):
    """write variables to a .npz file"""
    np.savez_compressed(
        output,
        read_length=read_length,
        insert_size=i_size,
        quality_hist_forward=hist_f,
        quality_hist_reverse=hist_r,
        subst_choices_forward=sub_f,
        subst_choices_reverse=sub_r,
        ins_forward=ins_f,
        ins_reverse=ins_r,
        del_forward=del_f,
        del_reverse=del_r
    )
