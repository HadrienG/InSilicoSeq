#!/usr/bin/env python
# -*- coding: utf-8 -*-

from scipy import stats

import sys
import pysam
import numpy as np


def read_bam(bam_file):
    """read a bam file and yield read"""
    try:
        bam = pysam.AlignmentFile(bam_file, "rb")
    except OSError as e:
        print('Error:', e)
        print('Couldn\'t read bam file')
        sys.exit(1)
    except ValueError as e:
        print('Error:', e)
        print('Bad bam file')
        sys.exit(1)
    else:
        with bam:
            for read in bam.fetch():
                if not read.is_unmapped():
                    yield read


def to_model(bam_path, model):
    """from a bam file, return all variables needed for modelling reads"""
    insert_size_dist = []
    qualities_forward = []
    qualities_reverse = []
    subst_array_f = np.zeros([read_length, 16])
    subst_array_r = np.zeros([read_length, 16])
    indel_array_f = np.zeros([read_length, 9])
    indel_array_r = np.zeros([read_length, 9])

    for read in read_bam(bam_path):
        # get insert size distribution
        if read.is_proper_pair:
            template_length = abs(read.template_length)
            i_size = template_length - (2 * len(read.seq))
            insert_size_dist.append(i_size)

        # get qualities
        if read.is_read1:
            qualities_forward.append(read.query_qualities)
        elif read.is_read2:
            qualities_reverse.append(read.query_qualities)

        # get mismatches
        alignment = read.get_aligned_pairs(
            matches_only=True,
            with_seq=True
            )
        has_indels = False
        for base in alignment:
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

    insert_size = int(np.mean(i_size_dist))
    hist_f = raw_qualities_to_histogram(qualities_forward, model)
    hist_r = raw_qualities_to_histogram(qualities_reverse, model)
    read_length = len(hist_f)


def write_to_file(model, read_length, hist_f, hist_r, sub_f, sub_r, ins_f,
                  ins_r, del_f, del_r, i_size, output):
    """write variables to a .npz file"""
    np.savez_compressed(
        output,
        model=model,
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
