***********
mim-tRNAseq
***********

**Original paper**: `Behrens <https://doi.org/10.1016/j.molcel.2021.01.028>`_ et al., 2021, Molecular Cell 81, 1–14

**Protocols paper** : `Behrens <https://doi.org/10.1016/j.xpro.2022.101579>`_ and Nedialkova, 2022, Experimental and computational workflow for the analysis of tRNA pools from eukaryotic cells by mim-tRNAseq. STAR Protocols. 3, 101579

.. image:: ../img/globular_multi.png
   :align: center
   :scale: 35%

:Author: Drew Behrens

:Version: 1.2

Modification-induced misincorporation tRNA sequencing.

This package is a semi-automated analysis pipeline for the quantification and analysis of tRNA expression. Given trimmed sequencing reads in fastq format, this pipeline will:

* Cluster tRNAs, index modifications, and perform SNP-tolerant read alignment with GSNAP_.
* Calculate coverage information and plots (useful for QC).
* Quantify expression.
* Calculate tRNA differential expression with DESeq2_.
* Analyze functional tRNA pools and tRNA completeness via 3'-CCA analysis.
* Comprehensive modification quantification and misincorporation signature analysis.
* Detect coordination between pairs of modifications and modification-aminoacylation with SLAC_ (SingLe-read Analysis of Crosstalks).

.. _GSNAP: http://research-pub.gene.com/gmap/
.. _DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
.. _SLAC: https://doi.org/10.1093/nar/gkac1185


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
