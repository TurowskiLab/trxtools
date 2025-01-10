.. trxtools documentation master file, created by
   sphinx-quickstart on Tue Oct 18 13:14:23 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

trxtools
========

version 0.3.0

trxtools is set of bioinformatic tools and methods to facilitate analysis of transcriptional data. 

Major focus is put on postprocessing of high-troughput sequencing data originating from the followin methods:

* CRAC/CLIP
* RNA-seq
* tRNA-seq
* NET-seq
* ChIP-seq

Some functions in the package use a specific file naming system:

``'expID', 'expDate', 'protein', 'condition1', 'condition2', 'condition3'``

i.e. ``C123_JK050420_POLR2A-FLAG_wt_arsenite30min`` will be used for an experiment ``C123`` conducted by a person ``JK`` on ``020420``. This is puldown with ``POLR2A-FLAG`` in ``wt`` cells after ``arsenite30min``

Reporting bugs
--------------
If you encounter a bug, please report it on the `GitHub issue tracker <https://github.com/TurowskiLab/trxtools/issues>`_.

Contents
--------

.. toctree::
   :maxdepth: 4

   scripts


.. toctree::
   :maxdepth: 4
   
   trxtools


.. toctree::
   :maxdepth: 4
   
   release

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
