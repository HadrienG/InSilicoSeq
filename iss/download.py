#!/usr/bin/env python

import random

from Bio import SeqIO
from Bio import Entrez


def ncbi(kingdom, n_genomes):
    """download random genomes sequences from ncbi genomes with entrez eutils.

    Args:
        kingdom (string): the kingdom from which the sequences are from
        n_genomes (int): the number of genomes to download

    Returns:
        list: list of handles
    """
    Entrez.email = ''
    Entrez.tool = 'InSilicoSeq'
    full_id_list = Entrez.read(Entrez.esearch(
        'genome', term='%s[Organism]' % kingdom, retmax=100000))['IdList']
    genomes = []
    n = 0
    while n < n_genomes:
        ident = random.choice(full_id_list)
        genome_info = Entrez.read(
            Entrez.esummary(db='genome', id=ident))[0]
        if genome_info['Assembly_Accession']:
            nucleotide_id = Entrez.read(Entrez.elink(
                dbfrom='genome',
                db='nucleotide',
                id=ident))[0]['LinkSetDb'][0]['Link'][0]['Id']
            nucleotide_info = Entrez.read(
                Entrez.esummary(db='nucleotide', id=nucleotide_id))[0]
            if nucleotide_info['AccessionVersion'].startswith('NC'):
                print('found %s' % genome_info['Assembly_Accession'])
                genome_record = Entrez.efetch(
                    'nucleotide',
                    id=nucleotide_id,
                    rettype='fasta',
                    retmode='txt')
                genomes.append(genome_record)
                n += 1
            else:
                continue
        else:
            continue

    return genomes


def to_fasta(genomes, output):
    """Write genomes to fasta

    Take the genomes from the ncbi function and write them to a fasta file:
    output_genomes.fasta

    Args:
        genomes (list): list of genome handles
        output (string): the output file prefix

    Returns:
        str: the file name
    """
    # define name of output files
    output_genomes = output + '_genomes.fasta'
    try:
        f = open(output_genomes, 'a')
    except PermissionError as e:
        logger.error('Failed to open output file: %s' % e)
        sys.exit(1)
    else:
        with f:
            for genome in genomes:
                record = SeqIO.read(genome, 'fasta')
                SeqIO.write(record, f, 'fasta')
    return output_genomes
