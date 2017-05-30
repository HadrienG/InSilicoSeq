.. InSilicoSeq documentation master file, created by
   sphinx-quickstart on Tue May 30 11:45:01 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

InSilicoSeq
============

InSilicoSeq (or iss) is a sequencing simulator producing (relatively) realistic
Illumina reads primarily intended for simulating metagenomic samples, although
it can be used to produce sequencing data from a single genome.

InSilicoSeq is written in python, and use a kernel density estimation model to
model the read quality of real sequencing data.

InSilicoSeq support substitution, insertion and deletion errors. If you don't
have the use for insertion and deletion error a basic error model is provided.

Contents
--------

.. toctree::
   :maxdepth: 2
   :glob:

   iss/install
   iss/generate
   iss/model
   iss/iss

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |date| date::
.. |time| date:: %H:%M

*This documentation was generated on* |date| *at* |time|.
