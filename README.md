<p align="center">
	<img src="docs/img/globular_multi.png" width="450" height="450">
</p>

# mim-tRNAseq
[![Documentation Status](https://readthedocs.org/projects/mim-trnaseq/badge/?version=latest)](https://mim-trnaseq.readthedocs.io/en/latest/?badge=latest)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
### Modification-induced misincorporation based sequencing of tRNAs using high-throughput RNA sequencing datasets.

This package is a semi-automated analysis pipeline for the quantification and analysis of tRNA expression. Given trimmed sequencing reads in fastq format, this pipeline will:
* Cluster tRNAs, index modifications, and perform SNP-tolerant read alignment with [GSNAP](http://research-pub.gene.com/gmap/)
* Calculate coverage information and plots (useful for QC)
* Quantify expression
* Calculate tRNA differential expression with [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).
* Analyse functional tRNA pools and tRNA completeness via 3'-CCA analysis
* Comprehensive modifcation quantification and misincorporation signature analysis

## Method strategy

Detailed methodology is shown in the image below, and described in Behrens et al. (2020)

![methods](/docs/img/method.png)

 
## Installation and usage

[![Documentation Status](https://readthedocs.org/projects/mim-trnaseq/badge/?version=latest)](https://mim-trnaseq.readthedocs.io/en/latest/?badge=latest)

Please see the full documentation for explanantions of dependencies, inputs formatting, and outputs.

To use mim-tRNAseq, please clone this git repository (`https://github.com/nedialkova-lab/mim-tRNAseq`, or download zip and extract) and run the mim-seq.py script in the scripts/ folder.
```bash
./scripts/mim-seq.py
```

**Note:** plans for future versions include interfacing R code from within Python with rpy2 and packaging the Python package on PyPI and conda.
This will significantly improve the installation and usage of mim-tRNAseq.

An example command to run mim-tRNAseq may look as follows. This will run an analysis between iPSC and HEK293T cells on an example dataset included in the package:
```bash
./scripts/mim-seq.py -t data/hg19_eColitK/hg19_eColitK.fa -o data/hg19_eColitK/hg19_eschColi-tRNAs.out 
-m data/hg19_eColitK/hg19-mitotRNAs.fa --snp-tolerance --cluster --cluster-id 0.97 --threads 15 
--min-cov 2000 --max-mismatches 0.1 --control-condition HEK293T --cca-analysis -n hg19_mix 
--out-dir hg19_HEKvsK562 --max-multi 4 --remap --remap-mismatches 0.05 sampleData_HEKvsK562.txt
```

The run should take around 15 minutes on a server using 15 processors (--threads 15).

## Contact

Drew Behrens: abehrens@biochem.mpg.de

Danny Nedialkova: nedialkova@biochem.mpg.de

Nedialkova laboratory: https://www.biochem.mpg.de/nedialkova


## Cite

Behrens et al., High-resolution quantitative profiling of tRNA pools by mim-tRNAseq (2020)

