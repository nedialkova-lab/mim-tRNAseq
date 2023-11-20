Quick-start guide
=================

Installation
^^^^^^^^^^^^

To use mim-tRNAseq, it is recommended to install the package using `conda`, preferably in its own environment. Significant time and dependency-related improvements can be made to using conda for managing environment and installing mimseq using the Miniforge_ version of conda which oncludes optional use for Mamba. We recommend installing Miniforge and then following the steps below:
::
	conda create -n mimseq python=3.7
	conda activate mimseq
	mamba install -c bioconda mimseq

usearch needs to be acquired and installed. Please do the following:
::
	wget https://drive5.com/downloads/usearch10.0.240_i86linux32.gz
	gunzip usearch10.0.240_i86linux32.gz
	chmod +x usearch10.0.240_i86linux32
	mv usearch10.0.240_i86linux32 usearch
	cp usearch /usr/local/bin

For this last `cp` command, root access is required. However, if this is not possible, please add the path of the usearch binary to your PATH (replace `full/path/to/usearch` with location of your usearch binary from above):
::
	export PATH=$PATH:full/path/to/usearch

Alternatively, mim-tRNAseq can be installed with `pip`, in which case all additional non-python package dependencies (including `usearch` as above, `BLAST`, `infernal`, `GMAP/GSNAP`, and all required R packages) will also need to be installed.
::
	pip install mimseq

The source code is also available on GitHub_

.. _GitHub: https://github.com/nedialkova-lab/mim-tRNAseq
.. _Miniforge: https://github.com/conda-forge/miniforge

Once installed, mim-tRNAseq should be executable and help displayed, by running
::
	mimseq --help

The package also comes with a `data/` folder which has the required tRNAscan-SE input files (and mitochondrial/plastid tRNA inputs where available) for a few species. Note that data folders containing "eColitK" in the name contain the E. coli Lys-TTT reference used as a spike-in in the paper. Using this reference in an experiment without this spike-in should not affect the results. Therefore, default inputs when using the `--species` are the references including E. coli Lys-TTT sequence.


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
| samtools        | >=1.11            | samtools_ |
+-----------------+-------------------+-----------+
| usearch         | 10.0.240          | usearch_  |
+-----------------+-------------------+-----------+
| bedtools        | >=2.30.0          | bedtools_ |
+-----------------+-------------------+-----------+
| INFERNAL        | >=1.1.4           | INFERNAL_ |
+-----------------+-------------------+-----------+
| BLAST           | 2.10.1            | BLAST_    |
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
| R base         | >=4        | R_                   |
+----------------+------------+----------------------+
| DESeq2         | >=1.26.0   | DESeq2_              |
+----------------+------------+----------------------+
| RColorBrewer   | 1.1.2      | RColorBrewer_        |
+----------------+------------+----------------------+
| pheatmap       | >=1.0.12   | pheatmap_            |
+----------------+------------+----------------------+
| calibrate      | >=1.7.7    | calibrate_           |
+----------------+------------+----------------------+
| gridExtra      | >=2.3      | gridExtra_           |
+----------------+------------+----------------------+
| plyr           | >=1.8.6    | plyr_                |
+----------------+------------+----------------------+
| dplyr          | >=1.0.6    | dplyr_               |
+----------------+------------+----------------------+
| reshape2       | >=1.4.3    | reshape2_            |
+----------------+------------+----------------------+
| circlize       | >=0.4.8    | circlize_            |
+----------------+------------+----------------------+
| tidyverse      | >=1.3.0    | tidyverse_           |
+----------------+------------+----------------------+
| ComplexHeatmap | >=2.2.0    | ComplexHeatmap_      |
+----------------+------------+----------------------+
| devtools       | >=2.4.1    | devtools_            |
+----------------+------------+----------------------+
| ggplot2        | >=3.3.5    | ggplot2_             |
+----------------+------------+----------------------+
| ggpol          | >= 0.0.7   | ggpol_               |
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
.. _ggplot2: https://ggplot2.tidyverse.org/
.. _ggpol: https://github.com/erocoar/ggpol

**Required Python packages:**

+------------+------------+-------------+
| Package    | Version    | Link        |
+============+============+=============+
| Python     | =3.7       | Python_     |
+------------+------------+-------------+
| Biopython  | >=1.79     | Biopython_  |
+------------+------------+-------------+
| pyfiglet   | >=0.8.post1| pyfiglet_   |
+------------+------------+-------------+
| pysam      | >=0.16.0.1 | pysam_      |
+------------+------------+-------------+
| pandas     | >=1.3.1    | pandas_     |
+------------+------------+-------------+
| numpy      | >=1.21.1   | NumPy_      |
+------------+------------+-------------+
| seaborn    | >=0.11.1   | seaborn_    |
+------------+------------+-------------+
| pybedtools | >=0.8.2    | pybedtools_ |
+------------+------------+-------------+
| requests   | >=2.26.0   | requests_   |
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

	mimseq --species Hsap --cluster-id 0.97 --threads 15 --min-cov 0.0005 --max-mismatches 0.075 --control-condition HEK293T -n hg38_test --out-dir hg38_HEK239vsK562 --max-multi 4 --remap --remap-mismatches 0.05 sampleData_HEKvsK562.txt

The run should take around 15 minutes on a server using 15 processors (`--threads 15`: please adjust according to your server capabilities).


Input files
^^^^^^^^^^^

Note: mim-tRNAseq does not require an input from Modomics_ for modification indexing, but automatically connects to the Modomics server and retrieves this information. Therefore an **internet connection is required** to run mim-tRNAseq. However, there is an offline copy of Modomics so that mim-tRNAseq can still run without connection, or if the Modomics database is offline.

mim-tRNAseq requires a few input files depending on the species of interest. Data for some of these species is already present in the `data/` folder and can be specified easily with the `--species` parameter (see :ref:`Pre-built references` below for available references). If not here, you may be able to obtain the required files from the GtRNAdb_, or request new predictions from the maintainers if your species of interest is not there. Failing this, the input files can be generated using tRNAscanSE_ on a genome reference file, but the annotation and naming of tRNAs becomes crucial for mim-tRNAseq functioning. Information on the tRNAscan-SE ID given in parantheses in the fasta file must match entries in the ".out" file for proper processing. This kind of manual prediction, annotation, and input into mim-tRNAseq can conceivably create many issues, as mim-tRNAseq expects files and annotations as thos formatted in GtRNADB files. This functionality has also not been extensively tested. 

Input files include:

* Genomic tRNA sequences: DNA sequences of tRNA loci in genome of interest in fasta format, including introns but excluding trailer and leader sequences.
* tRNA ".out" file: contains important info about tRNA introns.
* Experiment sample file: User-generated tab-delimited file with 2 columns. The first is the absolute path to trimmed tRNAseq reads. The second is the condition name, used to group replicates (e.g. WT or knock-out etc)
* OPTIONAL mitochondrial and/or plastid (in case of plant species) tRNA sequences: Can be obtained from the mitotRNAdb_ if available. First, find the organism of interest in the "Search Database" tab, select all sequences for organism, choose "Send FASTA" in the drop-down at the bottom of the results, and click "Submit". Or, for plant species, obtain sequences from PtRNAdb_ by going to "Search", choosing "Mitochondrial" and/or Plastid" in "Search by Genome", enabling "Search by Plant Name:" and searching for your species of interest. Download the results, and then reformat them to the correct format using the example `convertPtRNAdbSearch.py` script in the *Arabidopsis thaliana* data_ folder, making sure to change the file names in the script before running. Mitochondrial sequences can be specified to mim-tRNAseq with the `-m` or `--mito-trnas` parameter. Plastid sequences can be specified to mim-tRNAseq with the `-p` or `--plastid-trnas` parameter.

`additionalMods.txt` is automatically read in by mim-tRNAseq to add additional modifications to the modification index that may not be in Modomics yet. Some important modifications have already been added for certain species, mainly based on Clark et al. tRNA base methylation identification and quantification via high-throughput sequencing (2016), and Rafels-Ybern et al. Codon adaptation to tRNAs with Inosine modification at position 34 is widespread among Eukaryotes and present in two Bacterial phyla (2018).


Pre-built references
^^^^^^^^^^^^^^^^^^^^

`mimseq` contains a few pre-built references which available to specify at runtime with `--species`. All of these references include the *E. coli* `tRNA-Lys-TTT <http://gtrnadb.ucsc.edu/genomes/bacteria/Esch_coli_K_12_MG1655/genes/tRNA-Lys-TTT-1-1.html>`_ spike-in sequence as detailed in the original method (`Behrens <https://doi.org/10.1016/j.molcel.2021.01.028>`_ et al., 2021). Details on these references are given below:

* Hsap: *H. sapiens* hg38 (GRCh38) `tRNA predictions <http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/>`_.  
* Hsap19: *H. sapiens* hg19 (GRCh37) `tRNA predictions <http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi19/>`_.
* Mmus: *M. musculus* mm39 (GRCm39) `tRNA predictions <http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Mmusc39/>`_.
* Rnor: *R. norvegicus* rn7 (mRatBN7.2) `tRNA predictions <http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Rnorv7/>`_.
* Scer: *S. cerevisiae* S228C sacCer3 `tRNA predictions <http://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/>`_.
* Spom: *S. pombe* 972h- ASM294v2 `tRNA predictions <http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Schi_pomb_972h/>`_.
* Dmel: *D. melanogaster* dm6 (Aug. 2014 BDGP Release 6) `tRNA predictions <http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Dmela6/>`_.
* Drer: *D. rerio* danRer11 (GRCz11) `tRNA predictions <http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Dreri11/>`_.
* Ecol: *E. coli* str. K-12 substr. MG1655 ASM584v2 `tRNA predictions <http://gtrnadb.ucsc.edu/genomes/bacteria/Esch_coli_K_12_MG1655/>`_.
* Atha: *A. thaliana* TAIR10 (Feb 2011) `tRNA predictions <http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Athal10/>`_.

**Note:**

The Hsap, Hsap19, and Mmus references were built using the bed file supplied in the GtRNAdb downloads, which can be obtained from "Download tRNAscan-SE Results" on a species page. This bed file represents the "High Confidence Set and Top 30 Hits in Each Isotype of Filtered Sets" according to GtRNAdb (for hg38 example see `here <http://gtrnadb.ucsc.edu/GtRNAdb2/genomes/eukaryota/Hsapi19/Hsapi19-displayed-gene-list.html>`_). These predictions are reached by simply clicking "tRNA Predictions" on the left panel on a species page. We opted for this set of sequences to represent a less stringent set of tRNAs that might show expression despite filtering by tRNAScan-SE, thus allowing mimseq to filter unexpressed genes instead (using *--min-cov*). 

To create these references (since the fasta file is not directly supplied by GtRNAdb for this set of tRNAs), we extracted the sequence from the genome using bedtools, and subsequently renamed and reformatted the sequence headers with a custom script, *FastaHeadersforMimseq.py*. This analysis can be recreated for another species or genome by following the `README <https://github.com/nedialkova-lab/mim-tRNAseq/blob/master/mimseq/data/mm39-eColitK/README.txt>`_ for mouse mm39 as an example. Be sure to edit *FastaHeadersforMimseq.py* to suite your needs.

.. _Modomics: http://modomics.genesilico.pl/
.. _gtRNAdb: http://gtrnadb.ucsc.edu/
.. _tRNAscanSE: http://trna.ucsc.edu/tRNAscan-SE/
.. _mitotRNAdb: http://mttrna.bioinf.uni-leipzig.de/mtDataOutput/
.. _PtRNAdb: http://14.139.61.8/PtRNAdb/index.php
.. _data: https://github.com/nedialkova-lab/mim-tRNAseq/blob/master/mimseq/data/araTha1-eColitK/convertPtRNAdbSearch.py
