#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from iss import download


def test_to_fasta():
    ref_genome = SeqRecord(
        Seq(str('ATATA' * 100),
            IUPAC.unambiguous_dna
            ),
        id='ref_genome',
        description='test reference'
        )
    genome_list = [ref_genome, ref_genome]
    download.to_fasta(genome_list, 'test_genomes.fasta')


def test_ncbi():
    genome_list = download.ncbi('bacteria', 2)
    assert len(genome_list) == 2
