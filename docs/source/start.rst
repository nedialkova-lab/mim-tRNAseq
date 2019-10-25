Quick-start guide
=================

Installation
^^^^^^^^^^^^

To use mim-tRNAseq, please clone this git repository (`git clone https://github.com/drewjbeh/mim-tRNAseq.git`, or download zip and extract) and run the mim-seq.py script in the scripts/ folder.

::
    ./scripts/mim-seq.py

This will display the package usage and help page. Note the REQUIRED arguments and inputs. 
The package also comes with a data/ folder which has the required tRNAscan-SE input files for a few species. Note that data folders containing "eColitK" in the name contain the E. coli Lys-TTT reference used as a spike-in in the paper. Using this reference in an experiment without this spike-in should not effect the results.

**Note:** plans for future versions include interfacing R code from within Python with rpy2 and packaging the Python package on PyPI and conda
This will significantly improve the installation and usage of mim-tRNAseq

Usage
^^^^^

An example command to run mim-tRNAseq may look as follows:
::
    ./scripts/mim-seq.py -t data/hg19_eColitK/hg19_eColitK.fa -o data/hg19_eColitK/hg19_eschColi-tRNAs.out 
    -m data/hg19_eColitK/hg19-mitotRNAs.fa --snp-tolerance --cluster --cluster-id 0.97 --threads 15 
    --min-cov 1000 --max-mismatches 0.1 --control-condition kiPS --cca-analysis -n hg19_mix 
    --out-dir hg19_all_0.1_remap0.05_ID0.97 --max-multi 6 --remap --remap-mismatches 0.05 sampleData_hg19_all.txt

Please see our CodeOcean container for an example run of mim-tRNAseq on some sample data


Input files
^^^^^^^^^^^

Note: mim-tRNAseq does not require an input from Modomics_ for modification indexing, but automatically connexts to the Modomics servers and retrieves this information. Therefore an **internet connection is required** to run mim-tRNAseq.

mim-tRNAseq requires a few input files depending on the species of interest. Data for some of these species is already present in the data/ folder. If not here, you may be able to obtain the required files from the gtRNAdb_. Failing this, the input files can be generated using tRNAscanSE_ on a genome reference file. Input files include:

* Genomic tRNA sequences: DNA sequences of tRNA loci in genome of interest in fasta format, including introns but excluding trailer and leader sequences.
* tRNA ".out" file: contains important info about tRNA introns.
* Experiment sample file: User-generated tab-delimited file with 2 columns. The first is the absolute path to trimmed tRNAseq reads. The second is the condition name, used to group replicates (e.g. WT or knock-out etc)
* OPTIONAL mitochondrial tRNA sequences: Can be obtained from the mitotRNAdb_ if available. First, find the organism of interest in the "Search Database" tab, select all sequences for organism, choose "Send FASTA" in the drop-down at the bottom of the results, and click "Submit".

.. _Modomics: http://modomics.genesilico.pl/
.. _gtRNAdb: http://gtrnadb.ucsc.edu/
.. _tRNAscanSE: http://trna.ucsc.edu/tRNAscan-SE/
.. _mitotRNAdb: http://mttrna.bioinf.uni-leipzig.de/mtDataOutput/