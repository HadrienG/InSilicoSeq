#!/usr/bin/env python
# -*- coding: utf-8 -*-

from iss import generator

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import random


def test_reads():
    random.seed(42)
    ref_genome = SeqRecord(
        Seq(str('CGTTTCAACC' * 40),
            IUPAC.unambiguous_dna
            ),
        id='my_genome',
        description='test genome'
        )
    read_gen = generator.reads(ref_genome, 100, 1)
    big_read = ''.join(str(read.seq) for read in read_gen)
    assert big_read == 'ACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCC\
GTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCATTTCAACCCGTTTCAACCCGTTTCAACC\
CGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGCG\
TTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTC\
AACCCGTTTCAACCCGTTTCAACCCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTT\
CAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTTCAACCCGTTT'
