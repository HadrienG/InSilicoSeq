#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util
from iss import modeller

from scipy import stats
from random import random

import sys
import pysam
import logging
import numpy as np


def read_bam(bam_file, n_reads=1000000):
    """Bam file reader. Select random mapped reads from a bam file

    Args:
        bam_file (string): path to a bam file

    Yields:
        read: a pysam read object
    """
    logger = logging.getLogger(__name__)

    try:
        lines = pysam.idxstats(bam_file).splitlines()
        total_records = sum([int(l.split("\t")[2])
                            for l in lines if not l.startswith("#")])
        # total_records = sum(1 for _ in bam.fetch() if not _.is_unmapped)
        random_fraction = n_reads / total_records
        bam = pysam.AlignmentFile(bam_file, 'rb')  # reopen the file

    except (IOError, ValueError,
            ZeroDivisionError, pysam.utils.SamtoolsError) as e:
        logger.error('Failed to read bam file: %s' % e)
        sys.exit(1)
    else:
        logger.info('Reading bam file: %s' % bam_file)
        c = 0
        with bam:
            for read in bam.fetch():
                if not read.is_unmapped and random() < random_fraction:
                    c += 1
                    if logger.getEffectiveLevel() == 10:
                        print(
                            'DEBUG:iss.bam:Subsampling %s / %s reads' % (
                                c, n_reads),
                            end='\r')
                    yield read
                elif c >= n_reads:
                    break


def write_to_file(model, read_length, mean_f, mean_r, hist_f, hist_r,
                  sub_f, sub_r, ins_f, ins_r, del_f, del_r, i_size, output):
    """Write variables to a .npz file

    Args:
        model (string): the type of error model
        read_length (int): read length of the dataset
        mean_f (list): list of mean bin sizes
        mean_r (list): list of mean bin sizes
        hist_f (list): list of cumulative distribution functions for the
            forward read quality
        hist_r (list): list of cumulative distribution functions for the
            reverse read quality
        sub_f (list): list of dictionaries representing the substitution
            probabilities for the forward reads
        sub_r (list): list of dictionaries representing the substitution
            probabilities for the reverse reads
        ins_f (list): list of dictionaries representing the insertion
            probabilities for the forward reads
        ins_r (list): list of dictionaries representing the insertion
            probabilities for the reverse reads
        del_f (list): list of dictionaries representing the deletion
            probabilities for the forward reads
        del_r (list): list of dictionaries representing the deletion
            probabilities for the reverse reads
        i_size (int): distribution of insert size for the aligned reads
        output (string): prefix of the output file
    """
    logger = logging.getLogger(__name__)

    try:
        logger.info('Writing model to file: %s' % output)
        np.savez_compressed(
            output,
            model=model,
            read_length=read_length,
            insert_size=i_size,
            mean_count_forward=mean_f,
            mean_count_reverse=mean_r,
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
        logger.error('Failed to open output file: %s' % e)
        sys.exit(1)


def to_model(bam_path, output):
    """from a bam file, write all variables needed for modelling reads in
    a .npz model file

    For a brief description of the variables that will be written to the
        output file, see the bam.write_to_file function

    Args:
        bam_path (string): path to a bam file
        output (string): prefix of the output file
    """
    logger = logging.getLogger(__name__)

    insert_size_dist = []
    qualities_forward = []
    qualities_reverse = []
    subst_matrix_f = np.zeros([301, 16])  # we dont know the len of the reads
    subst_matrix_r = np.zeros([301, 16])  # yet. we will find out from the
    indel_matrix_f = np.zeros([301, 9])   # len of the quality lists
    indel_matrix_r = np.zeros([301, 9])

    # read the bam file and extract info needed for modelling
    for read in read_bam(bam_path):
        # get insert size distribution
        if read.is_proper_pair:
            template_length = abs(read.template_length)
            i_size = template_length - (2 * len(read.seq))
            insert_size_dist.append(i_size)

        # get qualities
        if read.is_read1:
            # get mean quality too
            quality_means = []
            read_quality = read.query_qualities
            mean_quality = np.mean(read_quality)
            if read.is_reverse:
                read_quality = read_quality[::-1]  # reverse the list

            quality_plus_mean = [
                (quality, mean_quality) for quality in read_quality]
            qualities_forward.append(np.asarray(quality_plus_mean))
            # qualities_forward.append(read.query_qualities)
        elif read.is_read2:
            # get mean quality too
            quality_means = []
            read_quality = read.query_qualities
            mean_quality = np.mean(read_quality)
            if read.is_reverse:
                read_quality = read_quality[::-1]  # reverse the list

            quality_plus_mean = [
                (quality, mean_quality) for quality in read_quality]
            qualities_reverse.append(np.asarray(quality_plus_mean))
            # qualities_reverse.append(read.query_qualities)

        # get mismatches
        alignment = read.get_aligned_pairs(
            matches_only=True,
            with_seq=True
            )
        read_has_indels = False
        for base in alignment:  # dispatch mismatches in matrix
            pos, subst, read_has_indels = modeller.dispatch_subst(
                base, read, read_has_indels)
            if read.is_read1 and subst is not None:
                subst_matrix_f[pos, subst] += 1
            elif read.is_read2 and subst is not None:
                subst_matrix_r[pos, subst] += 1
        if read_has_indels:  # dispatch indels in matrix
            for pos, indel in modeller.dispatch_indels(read):
                if read.is_read1:
                    indel_matrix_f[pos, indel] += 1
                elif read.is_read2:
                    indel_matrix_r[pos, indel] += 1

    logger.info('Calculating insert size distribution')
    # insert_size = int(np.mean(insert_size_dist))
    hist_insert_size = modeller.insert_size(insert_size_dist)

    logger.info('Calculating mean and base quality distribution')
    quality_bins_f = modeller.divide_qualities_into_bins(qualities_forward)
    quality_bins_r = modeller.divide_qualities_into_bins(qualities_reverse)

    # getting distribution of mean sequence quality
    mean_f = [len(quality_bin) for quality_bin in quality_bins_f]
    mean_r = [len(quality_bin) for quality_bin in quality_bins_r]

    hists_f = modeller.quality_bins_to_histogram(quality_bins_f)
    hists_r = modeller.quality_bins_to_histogram(quality_bins_r)

    # modern illumina instruments return reads of the same length
    # in case our bam file contains aligned reads of different length,
    # we coerce the model's read length to the smallest read of the bam file
    length_forward = min((len(x) for x in hists_f if len(x) > 1))
    length_reverse = min((len(x) for x in hists_r if len(x) > 1))
    read_length = min(length_forward, length_reverse)

    # now we can resize the substitution and indel matrices before
    # doing operations on them
    subst_matrix_f.resize([read_length, 16], refcheck=False)
    subst_matrix_r.resize([read_length, 16], refcheck=False)
    indel_matrix_f.resize([read_length, 9], refcheck=False)
    indel_matrix_r.resize([read_length, 9], refcheck=False)

    logger.info('Calculating substitution rate')
    subst_f = modeller.subst_matrix_to_choices(subst_matrix_f, read_length)
    subst_r = modeller.subst_matrix_to_choices(subst_matrix_r, read_length)

    logger.info('Calculating indel rate')
    # update the base count in indel matrices
    for position in range(read_length):
        indel_matrix_f[position][0] = sum(subst_matrix_f[position][::4])
        indel_matrix_r[position][0] = sum(subst_matrix_r[position][::4])

    ins_f, del_f = modeller.indel_matrix_to_choices(
        indel_matrix_f, read_length)
    ins_r, del_r = modeller.indel_matrix_to_choices(
        indel_matrix_r, read_length)

    write_to_file(
        'kde',
        read_length,
        mean_f,
        mean_r,
        hists_f,
        hists_r,
        subst_f,
        subst_r,
        ins_f,
        ins_r,
        del_f,
        del_r,
        hist_insert_size,
        output + '.npz')
