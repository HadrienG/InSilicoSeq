#!/usr/bin/env python2
# -*- coding: utf-8 -*-

import pytest

from iss import download
from iss.util import cleanup


@pytest.fixture
def setup_and_teardown():
    yield
    cleanup(["data/test_download.fasta"])


def download_to_fasta(setup_and_teardown):
    ftp_url = (
        "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/737/615/GCF_000737615.1_ASM73761v1/"
        "GCF_000737615.1_ASM73761v1_genomic.fna.gz"
    )
    download.assembly_to_fasta(ftp_url, "data/test_download.fasta")


def test_ncbi(setup_and_teardown):
    download.ncbi("bacteria", 2, "data/test_download.fasta")
