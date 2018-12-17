#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord

from iss import download
from iss.util import cleanup

from nose.tools import with_setup


def setup_function():
    output_file_prefix = 'data/.test'


def teardown_function():
    cleanup(['data/test_download.fasta'])


@with_setup(setup_function, teardown_function)
def download_to_fasta():
    ftp_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/737/615/GCF_000737615.1_ASM73761v1/GCF_000737615.1_ASM73761v1_genomic.fna.gz'
    download.download_to_fasta(ftp_url, 'data/test_download.fasta')


@with_setup(setup_function, teardown_function)
def test_ncbi():
    genome_list = download.ncbi('bacteria', 2, 'data/test_download.fasta')
