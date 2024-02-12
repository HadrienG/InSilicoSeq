#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import random

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from iss import generator
from iss.error_models import basic, kde
from iss.util import cleanup


@pytest.fixture
def setup_and_teardown():
    yield
    cleanup(
        [
            "data/.test.iss.tmp.my_genome.0_R1.fastq",
            "data/.test.iss.tmp.my_genome.0_R2.fastq",
            "data/.test.iss.tmp.my_genome.tsv",
        ]
    )


def test_cleanup_fail():
    with pytest.raises(SystemExit):
        cleanup("data/does_not_exist")


@pytest.fixture
def with_seed():
    random.seed(42)
    np.random.seed(42)


def test_simulate_and_save(setup_and_teardown):
    err_mod = basic.BasicErrorModel(451, 0)
    ref_genome = SeqRecord(Seq(str("AAAAACCCCC" * 100)), id="my_genome", description="test genome")
    forward_handle = open("data/.test.iss.tmp.my_genome.0_R1.fastq", "w")
    reverse_handle = open("data/.test.iss.tmp.my_genome.0_R2.fastq", "w")
    mutations_handle = open("data/.test.iss.tmp.my_genome.tsv", "w")
    generator.simulate_reads(
        ref_genome, err_mod, 1000, 0, forward_handle, reverse_handle, mutations_handle, "metagenomics", gc_bias=True
    )


def test_simulate_and_save_short(setup_and_teardown):
    err_mod = basic.BasicErrorModel()
    ref_genome = SeqRecord(Seq(str("AACCC" * 100)), id="my_genome", description="test genome")
    forward_handle = open("data/.test.iss.tmp.my_genome.0_R1.fastq", "w")
    reverse_handle = open("data/.test.iss.tmp.my_genome.0_R2.fastq", "w")
    mutations_handle = open("data/.test.iss.tmp.my_genome.tsv", "w")
    generator.simulate_reads(
        ref_genome, err_mod, 1000, 0, forward_handle, reverse_handle, mutations_handle, "metagenomics", gc_bias=True
    )


def test_small_input():
    with pytest.raises(AssertionError):
        err_mod = kde.KDErrorModel("data/ecoli.npz")
        ref_genome = SeqRecord(Seq(str("AAAAACCCCC")), id="my_genome", description="test genome")

        generator.simulate_read(ref_genome, err_mod, 1, 0, "metagenomics")


def test_basic(with_seed):
    err_mod = basic.BasicErrorModel(450, 0)
    ref_genome = SeqRecord(Seq(str("AAAAACCCCC" * 100)), id="my_genome", description="test genome")
    read_tuple = generator.simulate_read(ref_genome, err_mod, 1, 0, "metagenomics")
    big_read = "".join(str(read_tuple[0].seq) + str(read_tuple[1].seq))
    assert big_read[-15:] == "TTTTGGGGGTTTTTG"


def test_kde(with_seed):
    err_mod = kde.KDErrorModel("data/ecoli.npz")
    ref_genome = SeqRecord(Seq(str("CGTTTCAACC" * 400)), id="my_genome", description="test genome")
    read_tuple = generator.simulate_read(ref_genome, err_mod, 1, 0, "metagenomics")
    big_read = "".join(str(read_tuple[0].seq) + str(read_tuple[1].seq))
    assert big_read[:15] == "CCGTTTCAACCCGTT"


def test_kde_short(with_seed):
    err_mod = kde.KDErrorModel("data/ecoli.npz", 1000, 10)
    ref_genome = SeqRecord(Seq(str("AAACC" * 100)), id="my_genome", description="test genome")
    read_tuple = generator.simulate_read(ref_genome, err_mod, 1, 0, "metagenomics")
    big_read = "".join(str(read_tuple[0].seq) + str(read_tuple[1].seq))
    assert big_read == "ACCAAACCAAACCAAACCAAGGTTTGGTTTGGTTTGGTAT"


def test_amp(with_seed):
    err_mod = kde.KDErrorModel(os.path.join(os.path.dirname(__file__), "../profiles/MiSeq"), 1000, 10)
    ref_genome = SeqRecord(
        Seq(
            (
                "TTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGG"
                "CCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
                "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTT"
            )
        ),
        id="my_amplicon",
        description="test amplicon",
    )
    read_tuple = generator.simulate_read(ref_genome, err_mod, 1, 0, "amplicon")
    assert len(read_tuple[0].seq) == 301
    assert read_tuple[0].seq.startswith("TTTAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA")
    assert len(read_tuple[1].seq) == 301
    assert read_tuple[1].seq.startswith("AAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCT")
