# mim-tRNAseq
### Modification-induced misincorporation based sequencing of tRNAs using high-throughput RNA sequencing datasets.

This package is a semi-automated analysis pipeline for the quantification and analysis of tRNA expression. Given trimmed sequencing reads in fastq format, this pipeline will:
* Cluster tRNAs, index modifications, and perform SNP-tolerant read alignment with [GSNAP](http://research-pub.gene.com/gmap/)
* Calculate coverage information and plots (useful for QC)
* Quantify expression
* Calculate tRNA differential expression with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
* Analyse functional tRNA pools and tRNA completeness via 3'-CCA analysis
* Comprehensive modifcation quantification and misincorporation signature analysis

## Method strategy

The driving force behind this pipeline is the unique alignment strategy it uses to accurately place tRNA-seq reads during alignment. The method is based off the ability of a reverse transcriptase (RT), TGIRT, to misincorportate incorrect nucleotides at modified tRNA positions. 

Due to the abundance of modifications to tRNA residues as well as high sequence conservation between tRNA genes and isodecoders, the generation and alignment of tRNA sequencing datasets is faced with two problems: 1) the modifications cause normal RTs to fail during the cDNA synthesis step of RNA-seq library preparation, and 2) even when long enough reads are sequenced, they cannot be uniquely placed during alignment because of the high degree of sequeunce similarity between tRNAs. This results in high rates of multi-mapping reads which are difficult to use for downstream analysis.

Using mim-tRNAseq and TGIRT overcomes these problems. First, tRNA genes sharing anticodons are clustered according to a user-defined ID threshold to limit read placement ambiguity during alignment. Next, the misincorporation by TGIRT at modified nucleotides allows read-through of modifications producing longer reads and libraries with less biased tRNA coverage, further improving alignment statistics and quantification. Additionally, indexed modification data is used in SNP-tolerant alignments by GSNAP (Wu and Nacu, 2010) to account for misincorporations, and thereby guide more accurate read placement. Collectively, mim-tRNAseq improves tRNA gene coverage from sequencing data and thereby reduces bias, enables more data from libraries to be used by reducing multi-mapping, and overall improves estimation of tRNA expression.

Detailed methodology is shown in the image below, and described in Behrens et al. (2020)

![methods](/img/method.png)

## Dependencies

Please install all dependencies below before running mim-tRNAseq. In most cases, newer versions of the packages should be fine, but if you encounter any errors when running, first try to install the exact verisons of dependencies listed below.

Unix command line dependencies:

Tool | Version >= | Link
-----|------------|-----
Subread package | 1.6.2 | [subread package](http://subread.sourceforge.net/)
GMAP-GSNAP | 2019-02-26 | [GMAP/GSNAP](http://research-pub.gene.com/gmap/)
samtools | 1.7 | [samtools](http://www.htslib.org/)
usearch | 10.0.240 | [usearch](https://www.drive5.com/usearch/)
bedtools | 2.26.0 | [bedtools](https://bedtools.readthedocs.io/en/latest/content/installation.html)
R | 3.5.2 | [R](https://www.r-project.org/)
INFERNAL | 1.1.2 (July 2016) | [INFERNAL](http://eddylab.org/infernal/)
BLAST | 2.7.1 | [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

Required R packages:

Package | Version >= | Link
--------|------------|-----
DESeq2 | 1.22.2 | [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
RColorBrewer | 1.1.2 | [RColorBrewer](https://www.rdocumentation.org/packages/RColorBrewer/versions/1.1-2)
pheatmap | 1.0.12 | [pheatmap](https://www.rdocumentation.org/packages/pheatmap/versions/1.0.12)
calibrate | 1.7.2 | [calibrate](https://cran.r-project.org/web/packages/calibrate/index.html)
gridExtra | 2.3 | [girExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)
plyr | 1.8.4 | [plyr](https://www.rdocumentation.org/packages/plyr/versions/1.8.4)
reshape2 | 1.4.3 | [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)
circlize | 0.4.7 | [circlize](https://cran.r-project.org/web/packages/circlize/index.html)
tidyverse | 1.2.1 | [Install tidyverse](https://www.tidyverse.org/packages/)
ggpol | 0.0.5 | * [ggpol](https://github.com/erocoar/ggpol)
ComplexHeatmap | 1.99.5 | ** [ComplexHeatmap](https://github.com/jokergoo/ComplexHeatmap)

\* To install ggpol, please get development version from github:
```R
if (!require(devtools)) {
	install.packages('devtools')
    }
devtools::install_github('erocoar/ggpol')
```
\*\* To install ComplexHeatmap, please get development version from github:

```R
if (!require(devtools)) {
	install.packages('devtools')
    }
devtools::install_github('jokergoo/ComplexHeatmap')	
```
Required Python packages:

Package | Version >= | Link
--------|------------|-----
biopython | 1.70 | [Biopython](https://biopython.org/)
pyfiglet | 0.7.5 | [pyfiglet](https://pypi.org/project/pyfiglet/0.7/)
pybedtools | 0.7.10 | [pybedtools](https://daler.github.io/pybedtools/)
pysam | 0.14.1 | [pysam](https://pysam.readthedocs.io/en/latest/api.html)
pandas | 0.22.0 | [pandas](https://pandas.pydata.org/)
numpy | 1.14.2 | [NumPy](https://numpy.org/)
 
## Installation and usage

To use mim-tRNAseq, please clone this git repository (`git clone https://github.com/drewjbeh/mim-tRNAseq.git`, or download zip and extract) and run the mim-seq.py script in the scripts/ folder.
```bash
./scripts/mim-seq.py
```
This will display the package usage and help page. Note the REQUIRED arguments and inputs. 
The package also comes with a data/ folder which has the required tRNAscan-SE input files for a few species. Note that data folders containing "eColitK" in the name contain the E. coli Lys-TTT reference used as a spike-in in the paper. Using this reference in an experiment without this spike-in should not effect the results.

An example command to run mim-tRNAseq may look as follows:
```bash
./scripts/mim-seq.py -t data/hg19_eColitK/hg19_eColitK.fa -o data/hg19_eColitK/hg19_eschColi-tRNAs.out 
-m data/hg19_eColitK/hg19-mitotRNAs.fa --snp-tolerance --cluster --cluster-id 0.97 --threads 15 
--min-cov 1000 --max-mismatches 0.1 --control-condition kiPS --cca-analysis -n hg19_mix 
--out-dir hg19_all_0.1_remap0.05_ID0.97 --max-multi 6 --remap --remap-mismatches 0.05 sampleData_hg19_all.txt
```

## Input formatting

Note: mim-tRNAseq does not require an input from [Modomics](http://modomics.genesilico.pl/) for modification indexing, but automatically connexts to the Modomics servers and retrieves this information. Therefore an **internet connection is required** to run mim-tRNAseq.

mim-tRNAseq requires a few input files depending on the species of interest. Data for some of these species is already present in the data/ folder. If not here, you may be able to obtain the required files from the [gtRNAdb](http://gtrnadb.ucsc.edu/). Failing this, the input files can be generated using [tRNAscan-SE](http://trna.ucsc.edu/tRNAscan-SE/) on a genome reference file. Input files include:
* Genomic tRNA sequences: DNA sequences of tRNA loci in genome of interest in fasta format, including introns but excluding trailer and leader sequences.
* tRNA ".out" file: contains important info about tRNA introns.
* Experiment sample file: User-generated tab-delimited file with 2 columns. The first is the absolute path to trimmed tRNAseq reads. The second is the condition name, used to group replicates (e.g. WT or knock-out etc)
* OPTIONAL mitochondrial tRNA sequences: Can be obtained from the [mitotRNAdb](http://mttrna.bioinf.uni-leipzig.de/mtDataOutput/) if available. First, find the organism of interest in the "Search Database" tab, select all sequences for organism, choose "Send FASTA" in the drop-down at the bottom of the results, and click "Submit".

## Outputs

mim-tRNAseq automatically generates many data and figure output files, some are dependent on paramter options specified when running mim-tRNAseq.
The output are split into various subdirectories described below:

### indices

This folder contains the indices required for GSNAP alignment.

* `tRNAgenome\` contains the index of mature, processed tRNA transcripts
* `snp_index\` contains the SNP index generated from modified positions, needed by GSNAP for SNP-tolerant alignment. Only generated if --snp-tolerance is enabled.
* \*.log files contain log info from index generation

### annotation

This folder contains various files describing the tRNA trasncriptome of the genome of interest.

* 


## Contact

Drew Behrens: abehrens@biochem.mpg.de

Danny Nedialkova: nedialkova@biochem.mpg.de

Nedialkova laboratory: https://www.biochem.mpg.de/nedialkova


## Cite

Behrens et al., High-resolution quantitative profiling of tRNA pools by mim-tRNAseq (2020)

