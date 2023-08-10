#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pytest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from iss import download
from iss.util import cleanup

def setup_function():
    output_file_prefix = 'data/.test'


def teardown_cleanup():
    cleanup(['data/test_download.fasta'])


@pytest.fixture
def setup_and_teardown():
    setup_function()
    yield
    teardown_cleanup()


def download_to_fasta(setup_and_teardown):
    ftp_url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/737/615/GCF_000737615.1_ASM73761v1/GCF_000737615.1_ASM73761v1_genomic.fna.gz'
    download.assembly_to_fasta(ftp_url, 'data/test_download.fasta')


def test_ncbi(setup_and_teardown):
    genome_list = download.ncbi('bacteria', 2, 'data/test_download.fasta')
