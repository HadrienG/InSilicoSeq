#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio import Entrez

import io
import sys
import time
import zlib
import random
import logging

import requests


class BadRequestError(Exception):
    """Exception to raise when a http request does not return 200
    """

    def __init__(self, url, status_code):
        super().__init__('%s returned %d' % (url, status_code))


def ncbi(kingdom, n_genomes, output):
    """download random genomes sequences from ncbi genomes with entrez eutils
        and requests.

    Args:
        kingdom (string): the kingdom from which the sequences are from
        n_genomes (int): the number of genomes to download

    Returns:
        str: the output file
    """
    logger = logging.getLogger(__name__)
    Entrez.email = 'hadrien.gourle@slu.se'
    Entrez.tool = 'InSilicoSeq'
    Entrez.api_key = 'd784b36672ca73601f4a19c3865775a17207'
    full_id_list = Entrez.read(Entrez.esearch(
        'assembly',
        term='%s[Organism] AND "latest refseq"[filter] AND "complete genome"[filter]'
        % kingdom, retmax=100000))['IdList']
    n = 0
    logger.info('Searching for %s to download' % kingdom)
    while n < n_genomes:
        ident = random.choice(full_id_list)
        genome_info = Entrez.read(
            Entrez.esummary(
                db='assembly',
                id=ident))["DocumentSummarySet"]["DocumentSummary"][0]
        if genome_info['FtpPath_RefSeq']:
            url = genome_info['FtpPath_RefSeq']
            url = "%s/%s_%s_genomic.fna.gz" \
                % (genome_info['FtpPath_RefSeq'],
                   genome_info['AssemblyAccession'],
                   genome_info['AssemblyName'])
            logger.info('Downloading %s' % genome_info['AssemblyAccession'])
            try:
                assembly_to_fasta(url, output)
            except BadRequestError as e:
                logger.debug('Could not download %s' %
                             genome_info['AssemblyAccession'])
                logger.debug('Skipping and waiting two seconds')
                time.sleep(2)
            else:
                n += 1
                full_id_list.remove(ident)
    return output


def assembly_to_fasta(url, output, chunk_size=1024):
    """download an assembly from the ncbi ftp and append
    the chromosome sequence to a fasta file.

    This function discards the plsamid sequences!

    Args:
        url (string): an url to a fasta file
        output (string): the output file name

    Returns:
        str: the file name
    """
    logger = logging.getLogger(__name__)
    if url.startswith("ftp://"):  # requests doesnt support ftp
        url = url.replace("ftp://", "https://")
    if url:
        request = requests.get(url)
        if request.status_code == 200:
            request = zlib.decompress(
                request.content, zlib.MAX_WBITS | 32).decode()
        else:
            raise BadRequestError(url, request.status_code)

    with io.StringIO(request) as fasta_io:
        seq_handle = SeqIO.parse(fasta_io, 'fasta')
        chromosome = filter_plasmids(seq_handle)

        try:
            f = open(output, 'a')
        except (IOError, OSError) as e:
            logger.error('Failed to open output file: %s' % e)
            sys.exit(1)
        else:
            logger.debug('Writing genome to %s' % output)
            with f:
                SeqIO.write(chromosome, f, 'fasta')
    return output


def filter_plasmids(handle):
    """returns the largest sequence from a sequence handle
    """
    n = 0
    for record in handle:
        if len(record) > n:
            n = len(record)
            largest = record
    return largest
