#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import util

import sys
import random
import logging
import numpy as np


class ErrorModel(object):
    """Main ErrorModel Class

    This class is used to create inheriting classes and contains all
    the functions that are shared by all ErrorModel classes
    """
    @property
    def logger(self):
        component = "{}.{}".format(type(self).__module__, type(self).__name__)
        return logging.getLogger(component)

    def load_npz(self, npz_path, model):
        """load the error profile .npz file

        Args:
            npz_path (string): path to the npz file
            model (string): type of model. Could be 'cdf' or 'kde'. 'cdf' has
                been deprecated and is no longer available

        Returns:
            ndarray: numpy object containg variables necessary
                for error model construction
        """
        try:
            error_profile = np.load(npz_path, allow_pickle=True)
            assert error_profile['model'] == model
        except (OSError, IOError) as e:
            self.logger.error('Failed to read ErrorModel file: %s' % e)
            sys.exit(1)
        except AssertionError as e:
            self.logger.error(
                'Trying to load a %s ErrorModel in %s mode' % (
                    error_profile['model'], model))
            sys.exit(1)
        else:
            self.logger.debug('Loaded ErrorProfile: %s' % npz_path)
        return error_profile

    def introduce_error_scores(self, record, orientation):
        """Add phred scores to a SeqRecord according to the error_model

        Args:
            record (SeqRecord): a read record
            orientation (string): orientation of the read. Can be 'forward' or
                'reverse'

        Returns:
            SeqRecord: a read record with error scores
        """
        if orientation == 'forward':
            record.letter_annotations["phred_quality"] = self.gen_phred_scores(
                self.quality_forward, 'forward')
        elif orientation == 'reverse':
            record.letter_annotations["phred_quality"] = self.gen_phred_scores(
                self.quality_reverse, 'reverse')
        return record

    def mut_sequence(self, record, orientation):
        """Introduce substitution errors to a sequence

        If a random probability is higher than the probability of the basecall
        being correct, introduce a substitution error

        Args:
            record (SeqRecord): a read record with error scores
            orientation (string): orientation of the read. Can be 'forward' or
                'reverse'

        Returns:
            Seq: a sequence
        """

        # get the right subst_matrix
        if orientation == 'forward':
            nucl_choices = self.subst_choices_for
        elif orientation == 'reverse':
            nucl_choices = self.subst_choices_rev

        mutable_seq = record.seq.tomutable()
        quality_list = record.letter_annotations["phred_quality"]
        position = 0
        for nucl, qual in zip(mutable_seq, quality_list):
            if random.random() > util.phred_to_prob(qual) \
                    and nucl.upper() not in 'RYWSMKHBVDN':
                mutable_seq[position] = str(np.random.choice(
                    nucl_choices[position][nucl.upper()][0],
                    p=nucl_choices[position][nucl.upper()][1]))
            position += 1
        return mutable_seq.toseq()

    def adjust_seq_length(self, mut_seq, orientation, full_sequence, bounds):
        """Truncate or Extend reads to make them fit the read length

        When insertions or deletions are introduced to the reads, their length
        will change. This function takes a (mutable) read and a reference
        sequence, and extend or truncate the read if it has had an insertion
        or a deletion

        Args:
            mut_seq (MutableSeq): a mutable sequence
            orientation (string): orientation of the read. Can be 'forward' or
                'reverse'
            full_sequence (Seq): the reference sequence from which mut_seq
                comes from
            bounds (tuple): the position of the read in the full_sequence

        Returns:
            Seq: a sequence fitting the ErrorModel
        """
        read_start, read_end = bounds
        if len(mut_seq) == self.read_length:
            return mut_seq.toseq()
        elif len(mut_seq) > self.read_length:
            while len(mut_seq) > self.read_length:
                mut_seq.pop()
            return mut_seq.toseq()
        else:  # len smaller
            to_add = self.read_length - len(mut_seq)
            if orientation == 'forward':
                for i in range(to_add):
                    if read_end + i >= len(full_sequence):
                        nucl_to_add = 'A'
                    else:
                        nucl_to_add = str(full_sequence[read_end + i])
                    mut_seq.append(nucl_to_add)
            elif orientation == 'reverse':
                for i in range(to_add):
                    if read_end + i >= len(full_sequence):
                        nucl_to_add = 'A'
                    else:
                        nucl_to_add = util.rev_comp(
                            full_sequence[read_end + i])
                    mut_seq.append(nucl_to_add)
            return mut_seq.toseq()

    def introduce_indels(self, record, orientation, full_seq, bounds):
        """Introduce insertions or deletions in a sequence

        Introduce insertion and deletion errors according to the probabilities
        present in the indel choices list


        Args:
            record (SeqRecord): a sequence record
            orientation (string): orientation of the read. Can be 'forward' or
                'reverse'
            full_seq (Seq): the reference sequence from which mut_seq
                comes from
            bounds (tuple): the position of the read in the full_sequence

        Returns:
            Seq: a sequence with (eventually) indels
        """

        # get the right indel arrays
        if orientation == 'forward':
            insertions = self.ins_for
            deletions = self.del_for
        elif orientation == 'reverse':
            insertions = self.ins_rev
            deletions = self.del_rev

        mutable_seq = record.seq.tomutable()
        position = 0
        for nucl in range(self.read_length - 1):
            try:
                # skip ambiguous nucleotides
                if mutable_seq[nucl].upper() in 'RYWSMKHBVDN':
                    position += 1
                    continue
                for nucl_to_insert, prob in insertions[position].items():
                    if random.random() < prob:
                        # we want to insert after the base read
                        mutable_seq.insert(position + 1, str(nucl_to_insert))
                if random.random() < deletions[position][mutable_seq[nucl].upper()]:
                    mutable_seq.pop(position)
                position += 1
            except IndexError as e:
                continue

        seq = self.adjust_seq_length(
            mutable_seq, orientation, full_seq, bounds)
        return seq
