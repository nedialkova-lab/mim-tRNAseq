mim-tRNAseq introduction
========================

Method strategy
^^^^^^^^^^^^^^^

The driving force behind this pipeline is the unique alignment strategy it uses to accurately place tRNA-seq reads during alignment. The method is based on the ability of a thermostable group II intron reverse transcriptase (TGIRT) to misincorportate nucleotides at some modified tRNA positions. 

Due to the abundance of modifications to tRNA residues, as well as high sequence similarity between tRNA genes and isodecoders, the generation and alignment of tRNA sequencing datasets is faced with two problems: 1) popular RTs often cannot readthrough Watson-Crick face modifications during the cDNA synthesis step of RNA-seq library preparation, resulting in hard stops to RT, and 2) even when full-length reads are sequenced, they cannot be uniquely placed during alignment because of the high degree of sequence similarity between tRNAs. This results in high rates of multi-mapping reads which are difficult to use for downstream analysis.

Using mim-tRNAseq and TGIRT overcomes these problems. First, tRNA genes sharing anticodons are clustered according to a user-defined ID threshold to limit read placement ambiguity during alignment. Next, misincorporation by TGIRT at modified nucleotides allows read-through of modifications producing longer reads and libraries with less biased tRNA coverage, further improving alignment statistics and quantification. Additionally, indexed modification data is used in SNP-tolerant alignments by GSNAP (Wu and Nacu, 2010) to account for misincorporations, and thereby guide more accurate read placement. Finally, deconvolution of cluster-aligned reads to transcripts restores resolution. Collectively, mim-tRNAseq improves tRNA gene coverage from sequencing data and thereby reduces bias, enables more data from libraries to be used by reducing multi-mapping, and overall improves estimation of tRNA expression and modification abundance.

Further automated downstream analyses included in the pipeline allow a "one command" analysis of differential expression, detailed modification analysis and pairwise comparisons of their abundance between conditions, 3'-CCA completeness and coverage.

Detailed methodology is shown in the image below, and described in `Behrens <https://doi.org/10.1016/j.molcel.2021.01.028>`_ et al., 2021, Molecular Cell 81, 1–14

.. image:: ../img/method.png

Deconvolution
^^^^^^^^^^^^^

Since v0.3, the deconvolution algorithm has been updated to consider all mismatches as a set for each unique isodecoder (see image above). This theoretically allows the distinction and deconvolution of all unique tRNA sequences, since by definition each unique sequence has a distinct set of mismatches. 

Despite this, some clusters still have 3' coverage bias due to modifications that induce stops to RT. In these cases, complete deconvolution of reads for these sequences might not be possible if coverage at a mismatch defining that tRNA is significanly lower than at the 3' end of the tRNA (in this case, the number of reads assigned to this sequence will be underestimated). In order to overcome this, the *--deconv-cov-ratio* paramater can be used to set a threshold for this difference in coverage. Sequences not passing this threshold will be marked as not deconvoluted. 

Furthermore, some tRNAs may only be distinguishable from the parent cluster by positions that might also be modified sites. In these cases it becomes impossible to tell if a mismatch in a read is because it belongs to a member of the cluster, or if it is due to minsincorporations at the modification. Such tRNA transcripts (and the parent of the cluster) are also labelled as not devonvoluted.
See the image below for an overview of the handling of clusters that cannot be deconvoluted.

As of v0.4, transcripts that are not deconvoluted are renamed to provide details on which transcripts remain clustered. For example, if Ala-AGC-1 and Ala-AGC-2 are clustered and cannot be deconvoluted due to either insuffucient coverage, or mismatches at modified sites, these two transcripts will remain clustered, be renamed to Ala-AGC1/2 (the parent isodecoder number for the cluster is always listed first), and appear as a single entry for counts, modification analysis and coverage data/plots. 

.. image:: ../img/unsplitDeconv.png
