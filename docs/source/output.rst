Outputs
=======


mim-tRNAseq automatically generates many data and figure output files, some are dependent on parameter options specified when running mim-tRNAseq.
The outputs are split into various subdirectories described below:

`mim-tRNAseq\*.log`: Log file for the run.

**cov**

Outputs for coverage analysis.

* `\*coverage_byaa_norm.pdf`: Coverage per sample, stacked by isotype and normalised as a fraction of library size.
* `\*coverage_byaa_norm_scaled.pdf`: Normalized coverage as above, scaled relative to second last bin for comparability between samples.
* `sampleData\*.txt`: Sample data information with library size for coverage normalization.
* `coverage_bygene.txt`: Coverage and normalized coverage per tRNA/cluster per 4% bin by gene length.
* `coverage_byaa.txt`: Coverage and normalized coverage summed for isoacceptors per 4% bin by gene length.
* `coverage_bycluster_norm.pdf`: Plot of normalized coverage along gene length for each tRNA/cluster.
* `coverage_byaa_line.pdf`: Coverage per amino acid group in line graph.


**align**

Read alignment outputs.

* `mapping_stats.txt`: Alignment statistics from GSNAP alignment. Contains "NEW ALIGNMENT" when --remap is enabled.
* `align.log`: Log file from GSNAP alignment.
* `\*mult.bam`: bam file of multi-mapping read alignments per sample.
* `\*uniq.bam`: bam file of uniquely mapping read alignments per sample.
* `\*nomapping.bam`: bam file of unmapped reads per sample.
* `Primary_alignstats.pdf`: Plot of alignment stats from first alignment round.
* `Remap_alignstats.pdf`: Plot of alignment stats from read realignment if enabled.

**counts**

Read count outputs.

* `Anticodon_counts_raw.txt`: Raw read counts summed by tRNAs sharing anticodons.
* `Isodecoder_counts_raw.txt`: Raw deconvoluted isodecoder counts, parent info and size (i.e. number of identical sequences/gene copies). See methods in Behrens et al., 2021. Only produced if cluster-id < 1. The `Single_isodecoder` column indicates if this isodecoder was able to be deconvoluted from the parent or not based on coverage at the required mismatch. Counts for "False" sequences should be interpreted/used with caution. Such isodecoders and their parents are excluded from DESeq2 analysis.
* `*_counts_DESEqNormalized.txt`: Same as above, but read counts are DESeq2 counts normalized by library size factors. Same as counts appended to DESeq2 results files. Use these for as normalized expression counts for isodecoders and anticodon families.

**CCAanalysis**

Only generated if --cca-analysis flag is present. Contains data and plots for 3'-CCA analysis

* `\*ccaPlot.pdf`: Diverging bar plots indicating average proportions of 3'-CCA, 3'-CC, 3'-C and absent 3' ends for each condition. These proportions are calculated from uniquely aligned reads aligning to the 3' end of the reference transcript. If more than one condition is present in the sample file, there is one of these plots for each pairwise comparison. Otherwise, there is only one plot for the single condition. Percentages and vertical white line indicate average 3'-CCA proportions for the condition.
* `dinuc_plot`: Proportions of dinucleotide ends for *all* aligned reads for each alignment file.
* `CCAcounts.csv`: Data file used for plotting diverging bar plots. Counts of different 3' ends for each tRNA/cluster for each bam file.
* `CCAprops.csv`: Data file generated from `CCAcount.csv` in plotting diverging bar plots. Percentages of each 3' end type within each tRNA/cluster for each bam file.
* `AlignedDinucProportions.csv`: Data file for plotting dinuc_plot. Counts of dinucleotide ends for each bam file. 

**mods**

Only generated if --snp-tolerance is specified, or --max-mismatches is not 0.

* `*comb_heatmap.pdf`: Combined RT stop and misincorporation heatmaps for all tRNAs/clusters passing --min-cov threshold and successfully deconvoluted (see *mim-tRNAseq introduction*). Shows proportion of stops and misincorporations at canonical tRNA positions for all conditions. Available for mitochondrial/plastid tRNAs if -m is specified.
* `*misincProps.pdf`: Misincorporation proportions for each tRNA/cluster at selected conserved modified sites by identity of modified nucleotide. Available for mitochondrial/plastid tRNAs if -m is specified.
* `*misincSignatures_upstreamContext.pdf`: Signatures of misincorporated nucleotides at selected conserved modified sites, separated by identity of modified nucloetide and upstream nucleotide relative to RT direction. Available for mitochondrial/plastid tRNAs if -m is specified.
* `*misincSignatures_downstreamContext.pdf`: Signatures of misincorporated nucleotides at selected conserved modified sites, separated by identity of modified nucleotide and downstream nucleotide relative to RT direction. Available for mitochondrial/plastid tRNAs if -m is specified.
* `mismatchTable.csv`: Data table for all misincorporation analyses. Includes tRNA/cluster, misincorporation type, proportion of misincorporations, coverage at each position and canonical tRNA position information, among other useful information.
* `RTstopTable.csv`: Data table for RT stop analyses. Includes tRNA/cluster, canonical tRNA position and proportion of reads that stop at each position (normalized to total coverage of the reference sequence). This gives the relative frequency of reads stopping at all positions for a reference, the sum of which should equal 1.
* `readthroughTable.csv`: Data table similar to `RTstopTable.csv` except the proportion here represents 1 - proportion of reads that stop at a site, normalized to the total read coverage at that site only. This value will therefore estimates the proportion of reads per position that extend beyond that site
* `modContext.txt`: Nucleotide context information for selected modified positions. Note, `pos` here is not canonical position information but ungapped alignment positions for each tRNA/cluster.
* `allModsTable.csv`: All known modifications from Modomics and `additionalMods.txt` and newly detected mim-tRNAseq modifications for each cluster. 1-based numbering of position of modified base in mature transctipt sequence and canonical numbering. Used in new mods discovery.
* `predictedMods.csv`: Newly predicted modified sites per isodeocoder based on misinc-thresh.
* `modPos_totalMisincProp.csv`: Total misincrporation with coverage and nucleotide identity per conserved modified position for all transcripts.

**mods_logOR**

Differential modification analysis. Only generated if mods analysis is performed (see above) and if there are more than 1 condition in the experiment
Separate analyses for organellar and cytosolic tRNAs.

* `ConditionAvsConditionB_logOR.pdf`: Heatmaps of filtered, significant log10 odds ratios (logOR) for each tRNA at each position between condition A and B. Performed for all pairwise condition comparisons. Values are FDR adjusted chi-squared p-values <= 0.01, and filtered for known and newly detected modified sites in mimseq.
* 'ConditionAvsConditionB_logOR.csv': Data table for heatmap plotting.

**DESeq2**

Differential tRNA expression analyses using DESeq2. Cytosolic and organellar tRNA counts analyzed separately. Each analysis is split into two separate analyses performed on counts at the anticodon and isodecoder level, respectively. Both folders contain the same outputs.

* `qc-dispersions.png`: Diagnostic dispersion plot of isoacceptor/isodecoder expression dispersion before and after shrinkage of estimates towards the fitted estimates. See the DESeq2 analysis vignette here_ for details
* `qc-sampledists.png`: Sample distance heatmap based on variance stabilizing transformed counts. Hierarchical clustering used for clustering samples. Scale indicates sample distances.
* `qc-pca.png`: Principal component analysis plot for all samples according to normalized counts for isoacceptors/isodecoders.
* `vst-transformedCounts.csv`: Variace stabilizing transformed count data used for sample clustering. Also useful for comparing tRNA expression, although normalized counts are easier to understand for this purpose.
* `\*diffextr-countplot.pdf`: Count data plotted for each pairwise condition comparison. Significantly differentially expressed isoacceptors/isodecoders detected by DESeq2 (adjusted p-value < 0.05) indicated by coloured triangles. Note that isodecoders unable to be split from parents (as well as the corresponding parents) are excluded from DESeq2 analysis.
* `\*diffexpr-results.csv`: DESeq2 differential expression results for each pairwise condition comparison. Note, every pairwise comparison output also has normalized counts for *all* samples appended as the last set of columns.

.. _here: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#dispersion-plot-and-fitting-alternatives

**indices**

Indices required for GSNAP alignment.

* `tRNAgenome\` contains the index of mature, processed tRNA transcripts.
* `snp_index\` contains the SNP index generated from modified positions, needed by GSNAP for SNP-tolerant alignment. Only generated if --snp-tolerance is enabled.
* `\*.log`: files contain log info from index generation.

**annotation**

Various files describing the tRNA trascriptome of the genome of interest.

* `\*unsplitClusterInfo.txt`: details about cluster parents and members that were unable to be deconvoluted, including parent name, number of unsplit transcripts, the canonical tRNA position that prevented deconvolution, and the reason deconvolution was not possible.
* `\*tRNATranscripts.fa`: processed, intron spliced, 3'-CCA appended, and His 5'-G appended tRNA transcript sequences in fasta format.
* `\*modificationSNPs.txt`: SNP index information for each tRNA after matching to Modomics entries for species of interest.
* `\*isoacceptorInfo.txt`: Information on isoacceptor groups and their size in genome of interest.
* `\*maturetRNA.bed': bed6 file for mature tRNA transcripts - used for coverage calculations.
* `\*clusterTranscripts_align.stk`: Stockholm align file generated by INFERNAL cmalign for tRNA sequence and structural alignments. Used for metagene coverage plots.
* `cm.log`: log file for INFERNAL cmalign algorithm

	Parameter-dependent outputs:

	* `\*clusterTranscripts.fa`: Cluster parent transcript sequence if clustering is enabled.
	* `\*isodecoderInfo.txt`: Isodecoder representative gene with size of isodecoder group (i.e. number of identitical tRNA sequences). Only for cluster-id < 1
	* `\*clusters.bed`: bed6 file for cluster parents. Only if clustering is enabled.
	* `\*clusterInfo.txt`: Cluster parent-child relationship for every tRNA gene, with unique cluster number and size. Only if clustering is enabled.

**single_read_data**

Only generated if --crosstalks is specified. The analysis includes all modified sites based on misinc-thresh.

* `*crosstalks.tsv`: Data table for all tRNA crosstalk analyses by `SLAC <https://doi.org/10.1093/nar/gkac1185>`_. Includes tRNA/cluster, pair of crosstalking positions, Fisher exact test p-value, odds ratio, contingency table with read counts, FDR-corrected p-value, and canonical tRNA position information (NAs indicate low-coverage positions). The odds ratio informs whether two modifications/charging tend to appear together in the same read (OR > 1) or tend to be exclusive of one another (OR < 1).

