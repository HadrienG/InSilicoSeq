.. _generate:

Generating reads
================

InSilicoSeq comes with a set of pre-computed error models to allow the user
to easily generate reads from:

- HiSeq 2500
- MiSeq (300bp)

Per example generate 1 million MiSeq 300bp reads from a set of input genomes:

.. code-block:: bash

    iss generate --genomes genomes.fasta --model_file MiSeq \
    --output MiSeq_reads

This will create 2 files, `MiSeq_reads_R1.fastq` and `MiSeq_reads_R2.fastq` in
your current directory

If you have created your custom model, change `--model_file MiSeq` to your
custom model file:

.. code-block:: bash

    iss generate --genomes genomes.fasta --model_file my_model.npz \
    --output my_model_reads


Required input files
--------------------

By default, InSilicoSeq only requires 1 file in order to start generating
reads: 1 (multi-)fasta files containing your input genome(s).

If you don't want to use a multi-fasta file or don't have one at hand but are
equipped with an Internet connection, you can download random genomes from the
ncbi:

.. code-block:: bash

    iss generate --ncbi bacteria --n_genomes 10 --model_file MiSeq \
    --output MiSeq_ncbi

In addition the the 2 fastq files, the downloaded genomes will be in the file
`MiSeq_ncbi_genomes.fasta` in your current directory.

*Note: If possible, I recommend using InSilicoSeq with a fasta file as input.*
*The eutils utilities from the ncbi can be slow and quirky.*


Abundance distribution
----------------------

The abundance of the input genomes is determined (by default) by a log-normal
distribution.

Alternatively, you can use other distribution with the `--abundance` parameter:
`uniform`, `halfnormal`, `exponential` or `zero-inflated-lognormal`

If you wish to fine-tune the distribution of your genomes, InSilicoSeq also
accepts an abundance file:

.. code-block:: bash

    iss generate --genomes genomes.fasta --abundance_file abundance.txt \
    --model_file HiSeq2500 --output HiSeq_reads

Example abundance file for a multi-fasta containing 2 genomes: genome_A and
genome_B.

.. code-block:: bash

    genome_A    0.2
    genome_B    0.8


For the abundance to make sense, the total abundance in your abundance file
must equal 1.


Full list of options
--------------------

--genomes
^^^^^^^^^

Input genome(s) from where the reads will originate

--ncbi
^^^^^^

Download input genomes from RefSeq instead of using --genomes. Requires
--n_genomes option. Can be bacteria, viruses or archaea.

--n_genomes
^^^^^^^^^^^

How many genomes will be downloaded from the ncbi.
Required if --ncbi is set.

--abundance
^^^^^^^^^^^

abundance distribution (default: lognormal). Can be uniform, halfnormal,
exponential, lognormal or zero-inflated-lognormal.

--abundance_file
^^^^^^^^^^^^^^^^

abundance file for coverage calculations (default: None).

--n_reads
^^^^^^^^^

Number of reads to generate (default: 1000000)

--model
^^^^^^^

Error model. If not specified, using kernel density estimation (default: kde).
Can be 'kde', 'cdf' or 'basic'

--model_file
^^^^^^^^^^^^

Error model file. If not specified, using a basic error model instead
(default: None). Use 'HiSeq2500' or 'MiSeq' for a pre-computed error model
provided with the software.

--quiet
^^^^^^^

Disable info logging

--debug
^^^^^^^

Enable debug logging

--output
^^^^^^^^

Output file prefix (Required)
