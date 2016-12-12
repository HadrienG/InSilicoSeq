#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pysam
import numpy as np


def parse_bam(bam_file):
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
        # print(dir(bam))
        for read in bam.fetch():
            # print(dir(read))
            # print(read.is_paired)
            if not read.is_unmapped:
                md_tag = read.get_aligned_pairs(
                    matches_only=True,
                    with_seq=True
                    )
                for base in md_tag:
                    if read.seq[base[0]] != 'N':
                        query_pos = base[0]
                        ref_base = base[2]
                        query_base = read.seq[query_pos]
                        dispatch_key = ref_base + query_base
                        print(read.is_reverse, read.mate_is_reverse)
                        array[query_pos, dispatch_dict[dispatch_key]] += 1
    print(array)
