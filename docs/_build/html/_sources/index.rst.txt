.. trxtools documentation master file, created by
   sphinx-quickstart on Thu Jun 30 11:29:59 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

trxtools
========

trxtools was initiated as a storage for functions and scripts using in Turowski Lab to conduct various bioinformatics.
There is no intention to create a comprehensive package that can be operated by non-experienced user, but some scripts are easy to run for everyone. 
Importantly the package provides useful functions for skilled individuals.

Major focus is put on postprocessing of high-troughput sequencing data originating from the followin methods:

* CRAC/CLIP
* RNA-seq
* tRNA-seq
* NET-seq
* ChIP-seq

Some functions in the package use a specific file naming system:

``'expID', 'expDate', 'protein', 'condition1', 'condition2', 'condition3'``

i.e. ``C123_JK050420_POLR2A-FLAG_wt_arsenite30min`` will be used for an experiment ``C123`` conducted by a person ``JK`` on ``020420``. This is puldown with ``POLR2A-FLAG`` in ``wt`` cells after ``arsenite30min``

.. toctree::
   :maxdepth: 4

   scripts


.. toctree::
   :maxdepth: 4

   trxtools



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
