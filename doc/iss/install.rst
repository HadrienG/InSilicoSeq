.. _install:

Installing InSilicoSeq
======================

.. _using_pip:

Using pip
---------

To install InSilicoSeq, type the following in your terminal:

.. code-block:: bash

    pip install InSilicoSeq

It will install InSilicoSeq as well as the following dependencies:

* biopython
* numpy
* pysam
* scipy

.. _using_docker:

Using docker
-----------

If you wish to use InSilicoSeq using docker

.. code-block:: bash

    docker pull hadrieng/insilicoseq:0.5.0

To use InSilicoSeq with docker, you need to provide a `volume` to the
``docker run`` command. Given with the ``-v`` option, the volume is your way
to exchanging data (in this case, your input and output files) with the docker
container.

.. code-block:: bash

    docker run -v /Users/hadrien/data:/mnt/data -it --rm \
        hadrieng/insilicoseq:0.5.0 iss generate \
        --genomes /mnt/data/ecoli.fasta -f MiSeq \
        -o /mnt/data/reads_ecoli_miseq

The above command will mount the local folder ``/Users/hadrien/data`` onto
``/mnt/data`` on the docker side. The output reads will be located in
``/Users/hadrien/data`` when InSilicoSeq has finished running.
