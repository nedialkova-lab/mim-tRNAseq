Quick-start guide
=================

Installation
^^^^^^^^^^^^

To use mim-tRNAseq, it is recommended to install the package using `conda`, preferably in its own environment. Significant time improvements can be made to installing mimseq using mamba which we will use within the mimseq environment:
::
	conda create -n mimseq python=3.7
	conda activate mimseq
	conda install -c conda-forge mamba
	mamba install -c bioconda mimseq

usearch needs to be acquired and installed. Please do the following:
::
	wget https://drive5.com/downloads/usearch10.0.240_i86linux32.gz
	gunzip usearch10.0.240_i86linux32.gz
	chmod +x usearch10.0.240_i86linux32
	mv usearch10.0.240_i86linux32 usearch
	cp usearch /usr/local/bin

For this last `cp` command, root access is required. However, if this is not possible, please add the path of the usearch binary to your PATH:
::
	export PATH=$PATH:full/path/to/usearch

Alternatively, mim-tRNAseq can be installed with `pip`, in which case all additional non-python package dependencies (see below) will also need to be installed.
::
	pip install mimseq

The source code is also available on GitHub_

.. _GitHub: https://github.com/nedialkova-lab/mim-tRNAseq

Once installed, mim-tRNAseq should be executable and help displayed, by running
::
	mimseq --help

The package also comes with a `data/` folder which has the required tRNAscan-SE input files (and mitochondrial tRNA inputs where available) for a few species. Note that data folders containing "eColitK" in the name contain the E. coli Lys-TTT reference used as a spike-in in the paper. Using this reference in an experiment without this spike-in should not affect the results. Therefore, default inputs when using the `--species` are the references including E. coli Lys-TTT sequence.


Dependencies
^^^^^^^^^^^^

When using `conda` to install mim-tRNAseq, all dependencies below are managed and installed automatically. We therefore strongly recommend using conda to install mim-tRNAseq in a separate environment (see Installation above).
If you install from source or PyPi, please install all dependencies below before running mim-tRNAseq. In most cases, newer versions of the packages should be fine, but if you encounter any errors when running, first try to install the exact versions of dependencies listed below.

**Unix command line dependencies:**

+-----------------+-------------------+-----------+
|Tool             | Version           | Link      |
+=================+===================+===========+
| GMAP-GSNAP      | 2019-02-26        | GSNAP_    |
+-----------------+-------------------+-----------+
| samtools        | 1.7               | samtools_ |
+-----------------+-------------------+-----------+
| usearch         | 10.0.240          | usearch_  |
+-----------------+-------------------+-----------+
| bedtools        | 2.26.0            | bedtools_ |
+-----------------+-------------------+-----------+
| INFERNAL        | 1.1.2 (July 2016) | INFERNAL_ |
+-----------------+-------------------+-----------+
| BLAST           | 2.7.1             | BLAST_    |
+-----------------+-------------------+-----------+
| gcc             | 4.8.5             | gcc_      |
+-----------------+-------------------+-----------+

.. _GSNAP: http://research-pub.gene.com/gmap/
.. _samtools: http://www.htslib.org/
.. _usearch: https://www.drive5.com/usearch/
.. _bedtools: https://bedtools.readthedocs.io/en/latest/content/installation.html
.. _INFERNAL: http://eddylab.org/infernal/
.. _BLAST: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download
.. _gcc: https://gcc.gnu.org/

**Required R packages:**

+----------------+------------+----------------------+
| Package        | Version    | Link                 |
+================+============+======================+
| R base         | <3.7.0a0   | R_                   |
+----------------+------------+----------------------+
| DESeq2         | 1.26.0     | DESeq2_              |
+----------------+------------+----------------------+
| RColorBrewer   | 1.1.2      | RColorBrewer_        |
+----------------+------------+----------------------+
| pheatmap       | 1.0.12     | pheatmap_            |
+----------------+------------+----------------------+
| calibrate      | 1.7.5      | calibrate_           |
+----------------+------------+----------------------+
| gridExtra      | 2.3        | gridExtra_           |
+----------------+------------+----------------------+
| plyr           | 1.8.6      | plyr_                |
+----------------+------------+----------------------+
| dplyr          | 1.0.0      | dplyr_               |
+----------------+------------+----------------------+
| reshape2       | 1.4.3      | reshape2_            |
+----------------+------------+----------------------+
| circlize       | 0.4.8      | circlize_            |
+----------------+------------+----------------------+
| tidyverse      | 1.3.0      | tidyverse_           |
+----------------+------------+----------------------+
| ComplexHeatmap | 2.2.0      | ComplexHeatmap_      |
+----------------+------------+----------------------+
| devtools       | 2.2.2      | devtools_            |
+----------------+------------+----------------------+
| ggplot2        | =3.2.1     | ggplot2_             |
+----------------+------------+----------------------+

.. _R: https://cran.r-project.org/
.. _DESeq2: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
.. _RColorBrewer: https://www.rdocumentation.org/packages/RColorBrewer/versions/1.1-2
.. _pheatmap: https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12
.. _calibrate: https://cran.r-project.org/web/packages/calibrate/index.html
.. _gridExtra: https://cran.r-project.org/web/packages/gridExtra/index.html
.. _plyr: https://www.rdocumentation.org/packages/plyr/versions/1.8.4
.. _dplyr: https://cran.r-project.org/web/packages/dplyr/index.html
.. _reshape2: https://cran.r-project.org/web/packages/reshape2/index.html
.. _circlize: https://cran.r-project.org/web/packages/circlize/index.html
.. _tidyverse: https://www.tidyverse.org/packages/
.. _ComplexHeatmap: https://bioconductor.org/packages/release/bioc/html/ComplexHeatmap.html
.. _devtools: https://cran.r-project.org/web/packages/devtools/index.html

**Required Python packages:**

+------------+------------+-------------+
| Package    | Version    | Link        |
+============+============+=============+
| Python     | >=3        | Python_     |
+------------+------------+-------------+
| Biopython  | 1.70       | Biopython_  |
+------------+------------+-------------+
| pyfiglet   | 0.7.5      | pyfiglet_   |
+------------+------------+-------------+
| pysam      | 0.15.3     | pysam_      |
+------------+------------+-------------+
| pandas     | 0.22.0     | pandas_     |
+------------+------------+-------------+
| numpy      | 1.14.2     | NumPy_      |
+------------+------------+-------------+
| seaborn    | 0.10.1     | seaborn_    |
+------------+------------+-------------+
| pybedtools | 0.8.1      | pybedtools_ |
+------------+------------+-------------+
| requests   | 2.23.0     | requests_   |
+------------+------------+-------------+

.. _Python: https://www.python.org/
.. _Biopython: https://biopython.org/
.. _pyfiglet: https://pypi.org/project/pyfiglet/0.7/
.. _pysam: https://pysam.readthedocs.io/en/latest/api.html
.. _pandas: https://pandas.pydata.org/
.. _NumPy: https://numpy.org/
.. _seaborn: https://seaborn.pydata.org/
.. _pybedtools: https://daler.github.io/pybedtools/
.. _requests: https://requests.readthedocs.io/en/master/


Usage
^^^^^

An example command to run mim-tRNAseq may look as follows. This will run an analysis between HEK293T and K562 cells on an example dataset included in the package:
::

	mimseq --species Hsap38 --cluster-id 0.95 --threads 15 --min-cov 2000 --max-mismatches 0.1 --control-condition HEK293T -n hg38_test --out-dir hg38_HEK239vsK562 --max-multi 4 --remap --remap-mismatches 0.075 sampleData_HEKvsK562.txt

The run should take around 15 minutes on a server using 15 processors (`--threads 15`: please adjust according to your server capabilities).


Input files
^^^^^^^^^^^

Note: mim-tRNAseq does not require an input from Modomics_ for modification indexing, but automatically connects to the Modomics server and retrieves this information. Therefore an **internet connection is required** to run mim-tRNAseq. However, there is an offline copy of Modomics so that mim-tRNAseq can still run without connection, or if the Modomics database is offline.

mim-tRNAseq requires a few input files depending on the species of interest. Data for some of these species is already present in the `data/` folder and can be specified easily with the `--species` parameter. If not here, you may be able to obtain the required files from the gtRNAdb_. Failing this, the input files can be generated using tRNAscanSE_ on a genome reference file. Input files include:

* Genomic tRNA sequences: DNA sequences of tRNA loci in genome of interest in fasta format, including introns but excluding trailer and leader sequences.
* tRNA ".out" file: contains important info about tRNA introns.
* Experiment sample file: User-generated tab-delimited file with 2 columns. The first is the absolute path to trimmed tRNAseq reads. The second is the condition name, used to group replicates (e.g. WT or knock-out etc)
* OPTIONAL mitochondrial tRNA sequences: Can be obtained from the mitotRNAdb_ if available. First, find the organism of interest in the "Search Database" tab, select all sequences for organism, choose "Send FASTA" in the drop-down at the bottom of the results, and click "Submit".

`additionalMods.txt` is automatically read in by mim-tRNAseq to add additional modifications to the modification index that may not be in Modomics yet. Some important modifications have already been added for certain species, mainly based on Clark et al. tRNA base methylation identification and quantification via high-throughput sequencing (2016), and Rafels-Ybern et al. Codon adaptation to tRNAs with Inosine modification at position 34 is widespread among Eukaryotes and present in two Bacterial phyla (2018).

.. _Modomics: http://modomics.genesilico.pl/
.. _gtRNAdb: http://gtrnadb.ucsc.edu/
.. _tRNAscanSE: http://trna.ucsc.edu/tRNAscan-SE/
.. _mitotRNAdb: http://mttrna.bioinf.uni-leipzig.de/mtDataOutput/