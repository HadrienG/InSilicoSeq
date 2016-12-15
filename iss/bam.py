#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pysam
import numpy as np


def substitutions(bam_file):
    array = np.empty([125, 16])

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
                alignment = read.get_aligned_pairs(
                    matches_only=True,
                    with_seq=True
                    )
                for base in alignment:
                    if read.seq[base[0]] != 'N':  # let's not deal with Ns
                        query_pos = base[0]
                        ref_base = base[2]
                        query_base = read.seq[query_pos]
                        dispatch_key = ref_base + query_base
                        if '^' in read.get_tag('MD'):
                            continue
                            # how to handle deletions?
                            # cannot base ourselves of positions
                            # since you get a deletion, the bases are off :(
                        array[query_pos, dispatch_dict[dispatch_key]] += 1
    print(array)


def quality_distribution(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        n_reads = bam.count()
        array_gen = (np.array(
            read.query_qualities) for read in bam.fetch()
                if not read.is_unmapped)

        histograms = [np.histogram(i, bins=37) for i in zip(*array_gen)]
    return np.array(histograms)


def write_to_file(histograms, output):
    np.save(output, histograms)

    # TODO: THAT WILL GO TO THE ERROR_MODEL
    # phred_list = []
    # for hist in histograms:
    #     values, indices = hist
    #     weights = values / np.sum(values)
    #     random_quality = np.random.choice(
    #         indices[1:], p=weights
    #     )
    #     phred_list.append(round(random_quality))
    # print(phred_list)
