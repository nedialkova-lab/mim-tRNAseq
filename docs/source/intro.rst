mim-tRNAseq introduction
========================

Method strategy
^^^^^^^^^^^^^^^

The driving force behind this pipeline is the unique alignment strategy it uses to accurately place tRNA-seq reads during alignment. The method is based off the ability of a reverse transcriptase (RT), TGIRT, to misincorportate incorrect nucleotides at modified tRNA positions. 

Due to the abundance of modifications to tRNA residues as well as high sequence conservation between tRNA genes and isodecoders, the generation and alignment of tRNA sequencing datasets is faced with two problems: 1) the modifications cause normal RTs to fail during the cDNA synthesis step of RNA-seq library preparation, and 2) even when long enough reads are sequenced, they cannot be uniquely placed during alignment because of the high degree of sequeunce similarity between tRNAs. This results in high rates of multi-mapping reads which are difficult to use for downstream analysis.

Using mim-tRNAseq and TGIRT overcomes these problems. First, tRNA genes sharing anticodons are clustered according to a user-defined ID threshold to limit read placement ambiguity during alignment. Next, the misincorporation by TGIRT at modified nucleotides allows read-through of modifications producing longer reads and libraries with less biased tRNA coverage, further improving alignment statistics and quantification. Additionally, indexed modification data is used in SNP-tolerant alignments by GSNAP (Wu and Nacu, 2010) to account for misincorporations, and thereby guide more accurate read placement. Collectively, mim-tRNAseq improves tRNA gene coverage from sequencing data and thereby reduces bias, enables more data from libraries to be used by reducing multi-mapping, and overall improves estimation of tRNA expression.

Detailed methodology is shown in the image below, and described in Behrens et al. (2020)

.. image:: ../img/method.png


Dependencies
^^^^^^^^^^^^

Please install all dependencies below before running mim-tRNAseq. In most cases, newer versions of the packages should be fine, but if you encounter any errors when running, first try to install the exact verisons of dependencies listed below.

**Unix command line dependencies:**

+-----------------+-------------------+-----------+
|Tool             | Version >=        | Link      |
+=================+===================+===========+
| Subread package | 1.6.2             | subread_  |
+-----------------+-------------------+-----------+
| GMAP-GSNAP      | 2019-02-26        | GSNAP_    |
+-----------------+-------------------+-----------+
| samtools        | 1.7               | samtools_ |
+-----------------+-------------------+-----------+
| usearch         | 10.0.240          | usearch_  |
+-----------------+-------------------+-----------+
| bedtools        | 2.26.0            | bedtools_ |
+-----------------+-------------------+-----------+
| R               | 3.5.2             | R_        |
+-----------------+-------------------+-----------+
| INFERNAL        | 1.1.2 (July 2016) | INFERNAL_ |
+-----------------+-------------------+-----------+
| BLAST           | 2.7.1             | BLAST_    |
+-----------------+-------------------+-----------+

.. _subread: http://subread.sourceforge.net/
.. _GSNAP: http://research-pub.gene.com/gmap/
.. _samtools: http://www.htslib.org/
.. _usearch: https://www.drive5.com/usearch/
.. _bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html
.. _R: https://www.r-project.org/
.. _INFERNAL: http://eddylab.org/infernal/
.. _BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

**Required R packages:**

+----------------+------------+----------------------+
| Package        | Version >= | Link                 |
+================+============+======================+
| DESeq2         | 1.22.2     | DESeq2_              |
+----------------+------------+----------------------+
| RColorBrewer   | 1.1.2      | RColorBrewer_        |
+----------------+------------+----------------------+
| pheatmap       | 1.0.12     | pheatmap_            |
+----------------+------------+----------------------+
| calibrate      | 1.7.2      | calibrate_           |
+----------------+------------+----------------------+
| gridExtra      | 2.3        | gridExtra_           |
+----------------+------------+----------------------+
| plyr           | 1.8.4      | plyr_                |
+----------------+------------+----------------------+
| reshape2       | 1.4.3      | reshape2_            |
+----------------+------------+----------------------+
| circlize       | 0.4.7      | circlize_            |
+----------------+------------+----------------------+
| tidyverse      | 1.2.1      | tidyverse_           |
+----------------+------------+----------------------+
| ggpol          | 0.0.5      | \* ggpol_            |
+----------------+------------+----------------------+
| ComplexHeatmap | 1.99.5     | \*\* ComplexHeatmap_ |
+----------------+------------+----------------------+

.. _DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
.. _RColorBrewer: https://www.rdocumentation.org/packages/RColorBrewer/versions/1.1-2
.. _pheatmap: https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12
.. _calibrate: https://cran.r-project.org/web/packages/calibrate/index.html
.. _gridExtra: https://cran.r-project.org/web/packages/gridExtra/index.html
.. _plyr: https://www.rdocumentation.org/packages/plyr/versions/1.8.4
.. _reshape2: https://cran.r-project.org/web/packages/reshape2/index.html
.. _circlize: https://cran.r-project.org/web/packages/circlize/index.html
.. _tidyverse: https://www.tidyverse.org/packages/
.. _ggpol: https://github.com/erocoar/ggpol
.. _ComplexHeatmap: https://github.com/jokergoo/ComplexHeatmap

\* To install ggpol, please get development version from github:

::
    if (!require(devtools)) {
    install.packages('devtools')
    }
    devtools::install_github('erocoar/ggpol')

\*\* To install ComplexHeatmap, please get development version from github:

::
    if (!require(devtools)) {
    install.packages('devtools')
    }
    devtools::install_github('jokergoo/ComplexHeatmap')	

**Required Python packages:**

+------------+------------+-------------+
| Package    | Version >= | Link        |
+============+============+=============+
| Biopython  | 1.70       | Biopython_  |
+------------+------------+-------------+
| pyfiglet   | 0.7.5      | pyfiglet_   |
+------------+------------+-------------+
| pybedtools | 0.7.10     | pybedtools_ |
+------------+------------+-------------+
| pysam      | 0.14.1     | pysam_      |
+------------+------------+-------------+
| pandas     | 0.22.0     | pandas_     |
+------------+------------+-------------+
| numpy      | 1.14.2     | NumPy_      |
+------------+------------+-------------+

.. _Biopython: https://biopython.org/
.. _pyfiglet: https://pypi.org/project/pyfiglet/0.7/
.. _pybedtools: https://daler.github.io/pybedtools/
.. _pysam: https://pysam.readthedocs.io/en/latest/api.html
.. _pandas: https://pandas.pydata.org/
.. _NumPy: https://numpy.org/







