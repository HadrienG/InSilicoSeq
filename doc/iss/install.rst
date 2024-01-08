.. _install:

Installing InSilicoSeq
======================

.. _using_pip:

Using pip
---------

InSilicoSeq requires python>=3.5.
To install InSilicoSeq, type the following in your terminal:

.. code-block:: bash

    pip install InSilicoSeq

It will install InSilicoSeq as well as the following dependencies:

* biopython
* numpy
* pysam
* scipy

Other installation options
^^^^^^^^^^^^^^^^^^^^^^^^^^

* With conda:

.. code-block:: bash

    conda install -c bioconda insilicoseq

* Upgrading InSilicoSeq to the latest version:

.. code-block:: bash

    pip install --upgrade InSilicoSeq

* If you don't have administration rights on your machine:

.. code-block:: bash

    pip3 install --user InSilicoSeq

* if you wish to install InSilicoSeq at a custom location (i.e with a module system):

.. code-block:: bash

    prefix="/path/to/install/prefix"
    pip install --install-option="--prefix=/$prefix" InSilicoSeq

then add ``$prefix/bin`` to your ``PATH``, and
``$prefix/lib/python3.X/site-packages`` to your ``PYTHONPATH`` (replacing
python3.X with your python version)

.. _using_docker:

Using docker
-----------

If you wish to use InSilicoSeq using docker

.. code-block:: bash

    docker pull hadrieng/insilicoseq:latest

To use InSilicoSeq with docker, you need to provide a `volume` to the ``docker run`` command.
Given with the ``-v`` option, the volume is your way of exchanging data (in this case, your input and output files) with the docker container.

.. code-block:: bash

    docker run -v /Users/hadrien/data:/mnt/data -it --rm \
        hadrieng/insilicoseq iss generate \
        --genomes /mnt/data/genomes.fasta -m miseq \
        -o /mnt/data/reads

The above command will mount the local folder ``/Users/hadrien/data`` onto ``/mnt/data`` on the docker side.
The output reads will be located in ``/Users/hadrien/data`` when InSilicoSeq has finished running.
