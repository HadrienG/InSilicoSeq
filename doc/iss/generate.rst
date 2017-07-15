.. _generate:

Generating reads
================

InSilicoSeq comes with a set of pre-computed error models, to allow the user
to easily generate reads from:

- HiSeq 2500
- MiSeq (250bp)
- MiSeq (300bp)

Per example generate 1 million MiSeq 300bp reads from a set of input genomes
with an abundance file:

.. code-block:: bash

    iss generate --genomes genomes.fasta --abundance abundance_file.txt \
        --model_file MiSeq300bp --output MiSeq_reads


If you have created your custom model, just give the .npz file to --model_file
instead of one of the pre-built ones.


Required input files
--------------------

InSilicoSeq requires 2 files in order to start generating reads:
1 (multi-)fasta files containing your input genome(s) and 1 abundance file,
containing the abudance of the genomes in the sample.

Example genome file:

.. code-block:: bash

    >genome_A
    ATGC...
    >genome_B
    CCGT...
    ...

Example abundance file:

(The total abundance must be 1!)

.. code-block:: bash

    genome_A    0.2
    genome_B    0.4
    ...


Additional options
------------------

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
