***********
mim-tRNAseq
***********

`Behrens <https://doi.org/10.1016/j.molcel.2021.01.028>`_ et al., 2021, Molecular Cell 81, 1–14


.. image:: ../img/globular_multi.png
   :align: center
   :scale: 35%

:Author: Drew Behrens

:Version: 0.3.4

Modification-induced misincorporation tRNA sequencing.

This package is a semi-automated analysis pipeline for the quantification and analysis of tRNA expression. Given trimmed sequencing reads in fastq format, this pipeline will:

* Cluster tRNAs, index modifications, and perform SNP-tolerant read alignment with GSNAP_.
* Calculate coverage information and plots (useful for QC).
* Quantify expression.
* Calculate tRNA differential expression with DESeq2_.
* Analyze functional tRNA pools and tRNA completeness via 3'-CCA analysis.
* Comprehensive modification quantification and misincorporation signature analysis.

.. _GSNAP: http://research-pub.gene.com/gmap/
.. _DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html


Index
=====

.. toctree::
   :maxdepth: 2

   intro.rst
   start.rst
   output.rst
   contact.rst

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
