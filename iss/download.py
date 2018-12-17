#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
from Bio import Entrez

import zlib
import random
import logging

import requests


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
    genomes = []
    n = 0
    logger.info('Searching for %s to download' % kingdom)
    while n <= n_genomes:
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
            download_to_fasta(url, output)
            n += 1
    return output


def download_to_fasta(url, output, chunk_size=1024):
    """download an url and append to a fasta file

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
        request = zlib.decompress(request.content, zlib.MAX_WBITS | 32)

    try:
        f = open(output, 'ab')
    except (IOError, OSError) as e:
        logger.error('Failed to open output file: %s' % e)
        sys.exit(1)
    else:
        logger.debug('Writing genome to %s' % output)
        with f:
            f.write(request)
    return output
