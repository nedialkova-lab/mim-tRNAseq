# mim-tRNAseq
### Modification-induced misincorporation based sequencing of tRNAs using high-throughput RNA sequencing datasets.

This package is a semi-automated analysis pipeline for the quantification and analysis of tRNA expression. Given trimmed sequencing reads in fastq format, this pipeline will:
* Cluster tRNAs, index modifications, and perform SNP-tolerant read alignment with [GSNAP](http://research-pub.gene.com/gmap/)
* Calculate coverage information and plots (useful for QC)
* Quantify expression
* Calculate tRNA differential expression with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
* Analyse functional tRNA pools and tRNA completeness via 3'-CCA analysis
* Comprehensive modifcation quantification and misincorporation signature analysis

## Alignment strategy

The driving force behind this pipeline is the unique alignment strategy it uses to accurately place tRNA-seq reads during alignment. The method is based off the ability of a reverse transcriptase (RT), TGIRT, to misincorportate incorrect nucleotides at modified tRNA positions. 

Due to the abundance of modifications to tRNA residues as well as high sequence conservation between tRNA genes and isodecoders, the generation and alignment of tRNA sequencing datasets is faced with two problems: 1) the modifications cause normal RTs to fail during the cDNA synthesis step of RNA-seq library preparation, and 2) even when long enough reads are sequenced, they cannot be uniquely placed during alignment because of the high degree of sequeunce similarity between tRNAs. This results in high rates of multi-mapping reads which are difficult to use for downstream analysis.

Using mim-tRNAseq and TGIRT overcomes these problems. First, tRNA genes sharing anticodons are clustered according to a user-defined ID threshold to limit read placement ambiguity during alignment. Next, the misincorporation by TGIRT at modified nucleotides allows read-through of modifications producing longer reads and libraries with less biased tRNA coverage, further improving alignment statistics and quantification. Additionally, indexed modification data is used in SNP-tolerant alignments by GSNAP (Wu and Nacu, 2010) to account for misincorporations, and thereby guide more accurate read placement. Collectively, mim-tRNAseq improves tRNA gene coverage from sequencing data and thereby reduces bias, enables more data from libraries to be used by reducing multi-mapping, and overall improves estimation of tRNA expression.

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
ggpol | 0.0.5 | *
ComplexHeatmap | 1.99.5 | **

'* To install ggpol, please get development version from github:


Required Python packages:

Package | Version >=
--------|---------
biopython | 1.70
pyfiglet | 0.7.5
pybedtools | 0.7.10
pandas | 0.22.0
numpy | 1.14.0
 
## Usage

## Input formatting

## Outputs

## Contact



## Cite

## 
