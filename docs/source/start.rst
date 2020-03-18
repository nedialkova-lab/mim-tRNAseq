Quick-start guide
=================

Installation
^^^^^^^^^^^^

To use mim-tRNAseq, please clone this git repository (:code:`git clone https://github.com/nedialkova-lab/mim-tRNAseq`, or download zip and extract) and run the mim-seq.py script in the `scripts/` folder.
::

	./scripts/mim-seq.py

This will display the package usage and help page. Note the REQUIRED arguments and inputs. 
The package also comes with a `data/` folder which has the required tRNAscan-SE input files (and mitochondrial tRNA inputs where available) for a few species. Note that data folders containing "eColitK" in the name contain the E. coli Lys-TTT reference used as a spike-in in the paper. Using this reference in an experiment without this spike-in should not effect the results.

**Note:** plans for future versions include interfacing R code from within Python with rpy2 and packaging the Python package on PyPI and conda.
This will significantly improve the installation and usage of mim-tRNAseq.

Dependencies
^^^^^^^^^^^^

Please install all dependencies below before running mim-tRNAseq. In most cases, newer versions of the packages should be fine, but if you encounter any errors when running, first try to install the exact verisons of dependencies listed below.

**Unix command line dependencies:**

+-----------------+-------------------+-----------+
|Tool             | Version >=        | Link      |
+=================+===================+===========+
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
| pysam      | 0.14.1     | pysam_      |
+------------+------------+-------------+
| pandas     | 0.22.0     | pandas_     |
+------------+------------+-------------+
| numpy      | 1.14.2     | NumPy_      |
+------------+------------+-------------+

.. _Biopython: https://biopython.org/
.. _pyfiglet: https://pypi.org/project/pyfiglet/0.7/
.. _pysam: https://pysam.readthedocs.io/en/latest/api.html
.. _pandas: https://pandas.pydata.org/
.. _NumPy: https://numpy.org/


Usage
^^^^^

An example command to run mim-tRNAseq may look as follows. This will run an analysis between iPSC and HEK293T cells on an example dataset included in the package:
::

	./scripts/mim-seq.py -t data/hg19_eColitK/hg19_eColitK.fa -o data/hg19_eColitK/hg19_eschColi-tRNAs.out 
	-m data/hg19_eColitK/hg19-mitotRNAs.fa --snp-tolerance --cluster --cluster-id 0.97 --threads 15 
	--min-cov 2000 --max-mismatches 0.1 --control-condition HEK293T --cca-analysis -n hg19_mix 
	--out-dir hg19_HEKvsK562 --max-multi 4 --remap --remap-mismatches 0.05 sampleData_HEKvsK562.txt

The run should take around 15 minutes on a server using 15 processors (--threads 15).


Input files
^^^^^^^^^^^

Note: mim-tRNAseq does not require an input from Modomics_ for modification indexing, but automatically connexts to the Modomics servers and retrieves this information. Therefore an **internet connection is required** to run mim-tRNAseq.

mim-tRNAseq requires a few input files depending on the species of interest. Data for some of these species is already present in the `data/` folder. If not here, you may be able to obtain the required files from the gtRNAdb_. Failing this, the input files can be generated using tRNAscanSE_ on a genome reference file. Input files include:

* Genomic tRNA sequences: DNA sequences of tRNA loci in genome of interest in fasta format, including introns but excluding trailer and leader sequences.
* tRNA ".out" file: contains important info about tRNA introns.
* Experiment sample file: User-generated tab-delimited file with 2 columns. The first is the absolute path to trimmed tRNAseq reads. The second is the condition name, used to group replicates (e.g. WT or knock-out etc)
* OPTIONAL mitochondrial tRNA sequences: Can be obtained from the mitotRNAdb_ if available. First, find the organism of interest in the "Search Database" tab, select all sequences for organism, choose "Send FASTA" in the drop-down at the bottom of the results, and click "Submit".

`additionalMods.txt` is automatically read in by mim-tRNAseq to add additional modifications to the modification index that may not be in Modomics yet. Some important modifications have already been added for certain species, mainly based on Clark et al. tRNA base methylation identification and quantification via high-throughput sequencing (2016), and Rafels-Ybern et al. Codon adaptation to tRNAs with Inosine modification at position 34 is widespread among Eukaryotes and present in two Bacterial phyla (2018).

.. _Modomics: http://modomics.genesilico.pl/
.. _gtRNAdb: http://gtrnadb.ucsc.edu/
.. _tRNAscanSE: http://trna.ucsc.edu/tRNAscan-SE/
.. _mitotRNAdb: http://mttrna.bioinf.uni-leipzig.de/mtDataOutput/