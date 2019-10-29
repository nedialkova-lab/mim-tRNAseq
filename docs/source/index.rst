***********
mim-tRNAseq
***********

.. image:: /img/logo.png
   :align: center

:Author: Drew Behrens

:Version: 0.1

Modification-induced misincorporation based sequencing of tRNAs using high-throughput RNA sequencing datasets.

This package is a semi-automated analysis pipeline for the quantification and analysis of tRNA expression. Given trimmed sequencing reads in fastq format, this pipeline will:

* Cluster tRNAs, index modifications, and perform SNP-tolerant read alignment with [GSNAP](http://research-pub.gene.com/gmap/)
* Calculate coverage information and plots (useful for QC)
* Quantify expression
* Calculate tRNA differential expression with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
* Analyse functional tRNA pools and tRNA completeness via 3'-CCA analysis
* Comprehensive modifcation quantification and misincorporation signature analysis


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
