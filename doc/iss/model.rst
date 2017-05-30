.. _model:

Creating an Error Model
=======================

If you do not wish to use the pre-computed error models provided with
InSilicoSeq, it is possible to create your own.

InSilicoSeq creates error models from .bam files. The input bam file should be
a set of reads aligned against a reference genome.

Given you have two read files, `reads_R1.fastq.gz` and`reads_R2.fastq.gz`,
and a referene genome `genome.fasta`:

Align you reads against the reference:

.. code-block:: bash

    bowtie2-build genome.fasta genome
    bowtie2 -x genome -1 reads_R1.fastq.gz \
        -2 reads_R2.fastq.gz | samtools view -bS | samtools sort -o genome.bam
    samtools index genome.bam

then build the model:

.. code-block:: bash

    iss model -b genome.bam -o genome

which will create a `genome.npz` file containing your newly built model
