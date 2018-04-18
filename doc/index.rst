.. InSilicoSeq documentation master file, created by
   sphinx-quickstart on Tue May 30 11:45:01 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

InSilicoSeq
============

InSilicoSeq is a sequencing simulator producing realistic Illumina reads.
Primarily intended for simulating metagenomic samples, it can also be used to produce sequencing data from a single genome.

InSilicoSeq is written in python, and use kernel density estimators to model the read quality of real sequencing data.

InSilicoSeq support substitution, insertion, deletion errors, and models gc bias and insert size distribution.

Details
-------

**Authors:** Hadrien Gourlé, Juliette Hayer, Oskar E. Karlsson and Erik Bomgcam-Rudloff
**Contact:** `hadrien.gourle@slu.se<hadrien.gourle@slu.se>`_
**GitHub:** `https://github.com/HadrienG/InSilicoSeq<https://github.com/HadrienG/InSilicoSeq>`_
**License:** MIT


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
