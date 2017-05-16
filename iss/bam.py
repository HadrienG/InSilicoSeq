#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import modeller
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
                if not read.is_unmapped:
                    yield read


def write_to_file(model, read_length, hist_f, hist_r, sub_f, sub_r, ins_f,
                  ins_r, del_f, del_r, i_size, output):
    """write variables to a .npz file"""
    try:
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
    except PermissionError as e:
        print('Error:', e)
        print('Permission denied')
        sys.exit(1)


def to_model(bam_path, model, output):
    """from a bam file, write all variables needed for modelling reads in
    a .npz model file"""
    insert_size_dist = []
    qualities_forward = []
    qualities_reverse = []
    subst_matrix_f = np.zeros([301, 16])  # we do no know the len of the reads
    subst_matrix_r = np.zeros([301, 16])  # yet
    indel_matrix_f = np.zeros([301, 9])
    indel_matrix_r = np.zeros([301, 9])

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
        read_has_indels = False
        for base in alignment:
            pos, subst, read_has_indels = modeller.dispatch_subst(
                base, read, read_has_indels)
            if read.is_read1 and subst is not None:
                subst_matrix_f[pos, subst] += 1
            elif read.is_read2 and subst is not None:
                subst_matrix_r[pos, subst] += 1
        if read_has_indels:
            for pos, indel in modeller.dispatch_indels(read):
                if read.is_read1:
                    indel_matrix_f[pos, indel] += 1
                elif read.is_read2:
                    indel_matrix_r[pos, indel] += 1

    insert_size = int(np.mean(insert_size_dist))
    hist_f = modeller.raw_qualities_to_histogram(qualities_forward, model)
    hist_r = modeller.raw_qualities_to_histogram(qualities_reverse, model)
    read_length = len(hist_f)
    # now we can resize the substitution and indel matrices before
    # doing operations on them
    subst_matrix_f.resize([read_length, 16])
    subst_matrix_r.resize([read_length, 16])
    indel_matrix_f.resize([read_length, 9])
    indel_matrix_r.resize([read_length, 9])

    subst_f = modeller.subst_matrix_to_choices(subst_matrix_f, read_length)
    subst_r = modeller.subst_matrix_to_choices(subst_matrix_r, read_length)

    # update the base count in indel matrices
    for position in range(read_length):
        indel_matrix_f[position][0] = sum(subst_matrix_f[position][::4])
        indel_matrix_r[position][0] = sum(subst_matrix_r[position][::4])

    ins_f, del_f = modeller.indel_matrix_to_choices(
        indel_matrix_f, read_length)
    ins_r, del_r = modeller.indel_matrix_to_choices(
        indel_matrix_r, read_length)

    write_to_file(
        model,
        read_length,
        hist_f,
        hist_r,
        subst_f,
        subst_r,
        ins_f,
        ins_r,
        del_f,
        del_r,
        insert_size,
        output + '.npz')
