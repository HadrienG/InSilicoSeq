.. _model:

Creating an Error Model
=======================

If you do not wish to use the pre-computed error models provided with
InSilicoSeq, it is possible to create your own.

InSilicoSeq creates error models from .bam files. The input bam file should be
a set of reads aligned against a reference genome.

Given you have two read files, `reads_R1.fastq.gz` and`reads_R2.fastq.gz`,
and a referene genome `genome.fasta`:

Align you reads against the reference
-------------------------------------

.. code-block:: bash

    bowtie2-build genomes.fasta genomes
    bowtie2 -x genomes -1 reads_R1.fastq.gz \
        -2 reads_R2.fastq.gz | samtools view -bS | samtools sort -o genomes.bam
    samtools index genomes.bam

Build the model
---------------

.. code-block:: bash

    iss model -b genomes.bam -o genomes

which will create a `genomes.npz` file containing your newly built model

Additional options
------------------

--model
^^^^^^^

Error model to build. If not specified, using kernel density estimation
(default: kde). Can be 'kde' or 'cdf'

--quiet
^^^^^^^

Disable info logging

--debug
^^^^^^^

Enable debug logging
