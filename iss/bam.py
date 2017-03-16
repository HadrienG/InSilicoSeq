#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pysam
import numpy as np


def substitutions(bam_file, read_length):
    array_f = np.empty([read_length, 16])
    array_r = np.empty([read_length, 16])

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

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if not read.is_unmapped:
                # a list of tuples: aligned read (query) and reference
                # positions. the parameter with_seq adds the ref sequence as
                # the 3rd element of the tuples. substitutions are lower-case
                alignment = read.get_aligned_pairs(
                    matches_only=True,
                    with_seq=True
                    )
                for base in alignment:
                    if read.seq[base[0]] != 'N':  # let's not deal with Ns
                        query_pos = base[0]
                        query_base = read.seq[query_pos]
                        ref_base = base[2]
                        dispatch_key = ref_base + query_base
                        if '^' in read.get_tag('MD'):
                            continue
                            # how to handle deletions?
                            # cannot base ourselves of positions
                            # since you get a deletion, the bases are off :(
                        if read.is_read1:
                            array_f[
                                query_pos,
                                dispatch_dict[dispatch_key]] += 1
                        elif read.is_read2:
                            array_r[
                                query_pos,
                                dispatch_dict[dispatch_key]] += 1
    return array_f, array_r


def quality_distribution(bam_file):
    """Generate numpy histograms for each position of the input mapped reads.
    A histogram contains the distribution of the phred scores for one position
    in all the reads. Returns a numpy array of the histograms for each position

    Arguments:
    bam_file -- the input bam file containing mapped reads
    """
    # deal with the forward reads
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        array_gen_f = (np.array(
            read.query_qualities) for read in bam.fetch()
                if not read.is_unmapped and read.is_read1)
        histograms_forward = [np.histogram(
            i, bins=range(0, 41)) for i in zip(*array_gen_f)]

    # deal with the reverse reads
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        array_gen_r = (np.array(
            read.query_qualities) for read in bam.fetch()
                if not read.is_unmapped and read.is_read2)
        histograms_reverse = [np.histogram(
            i, bins=range(0, 41)) for i in zip(*array_gen_r)]

    return histograms_forward, histograms_reverse


def get_insert_size(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        i_size_dist = [
            abs(read.template_length) for read in bam.fetch()
            if not read.is_unmapped and read.is_proper_pair]
        i_size = np.mean(i_size_dist)
    return int(i_size)


def write_to_file(read_length, hist_f, hist_r, sub_f, sub_r, i_size, output):
    np.savez_compressed(
        output,
        read_length=read_length,
        insert_size=i_size,
        quality_hist_forward=hist_f,
        quality_hist_reverse=hist_r,
        subst_matrix_forward=sub_f,
        subst_matrix_reverse=sub_r
    )
