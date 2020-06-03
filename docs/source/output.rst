Outputs
=======


mim-tRNAseq automatically generates many data and figure output files, some are dependent on paramter options specified when running mim-tRNAseq.
The output are split into various subdirectories described below:

`mim-tRNAseq\*.log`: Log file for the run.

**cov**

Outputs for coverage analysis.

* `sampleData\*.txt`: Sample data information with library size for coverage normalization.
* `coverage_bygene.txt`: Coverage and normalized coverage per tRNA/cluster per 4% bin by gene length.
* `coverage_byaa.txt`: Coverage and normalized coverage summed for isoacceptors per 4% bin by gene length.
* `coverage_bycluster_norm.pdf`: Plot of normalized coverage along gene length for each tRNA/clustrer.
* `coverage_byaa_line.pdf`: Coverage per amino acid group in line graph.
* `\*coverage_byaa_norm.pdf`: Normalized coverage per amino acid and aligment file. Available for mitochondrial tRNAs if -m is given.
* `\*coverage_byaa_norm_scaled.pdf`: Normalized coverage scaled realtive to second last bin (94%) for comparability between samples. Available for mitochondrial tRNAs if -m is given.

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

* `counts.txt`: `featureCounts` generated counts per tRNA/cluster.
* `featureCounts.log`: Log file from `featureCounts` run.
* `counts.txt.summary`: Reap assignment statistics from featureCounts run.
* `Anticodon_counts.txt`: Read counts summed by tRNAs sharing anticodons.
* `Isodecoder_counts.txt`: Deconvoluted isodecoder counts calculated from mismatch proportions. See methods in Behrens et al., 2020. Only produced if cluster-id < 1.

**CCAanalysis**

Only generated if --cca-analysis flag is present. Contains data and plots for 3'-CCA analysis

* `\*ccaPlot.pdf`: Diverging bar plots indicating average proportions of 3'-CCA, 3'-CC, 3'-C and absent 3' ends for each condition. These proportions are calculated from uniquely aligned reads aligning to the 3' end of the reference transcript. If more than one condition is present in the sample file, there is one of these plots for each pairwise comparison. Otherwise, there is only one plot for the single condition. Percentages and vertical white line indicate average 3'-CCA proportions for the condition.
* `dinuc_plot`: Proportions of dinuleotide ends for *all* aligned reads for each alignment file.
* `CCAcounts.csv`: Data file used for plotting diverging bar plots. Counts of different 3' ends for each tRNA/cluster for each bam file.
* `AlignedDinucProportions.csv`: Data file for plotting dinuc_plot. Counts of dinucleotide ends for each bam file. 

**mods**

Only generated if --snp-tolerance is specified, or --max-mismatches is not 0.

* `\*comb_heatmap.pdf`: Combined RT stop and misincorporation heatmaps for all tRNAs/clusters passing --min-cov threshold, showing proportion of stops and misincorporations at canonical tRNA positions for all conditions. Available for mitochondrial tRNAs if -m is specified.
* `\*misincProps/pdf`: Misincorporation proportions for each tRNA/cluster at selected conserved modified sites by identity of modified nucleotide. Available for mitochondrial tRNAs if -m is specified.
* `\*misincSignatures_upstreamContext.pdf`: Signatures of misincorporated nucleotides at selected conserved modified sites, separated by identity of modified nucloetide and upstream nucleotide relative to RT direction. Available for mitochondrial tRNAs if -m is specified.
* `\*misincSignatures_downstreamContext.pdf`: Signatures of misincorporated nucleotides at selected conserved modified sites, separated by identity of modified nucloetide and downstream nucleotide relative to RT direction. Available for mitochondrial tRNAs if -m is specified.
* `mismatchTable.csv`: Data table for all misincorporation analyses. Includes tRNA/cluster, misincorporation type, proportion of misincorporations, coverage at each position and canonical tRNA position information, among other useful information.
* `RTstopTable.csv`: Data table for RT stop analyses. Includes tRNA/cluster, canonical tRNA position and proportion of reads that stop at each position (normalised to total coverage of the reference sequence). This gives the relative frequence of reads stopping at all positions for a reference, the sum of which should equal 1.
* `readthroughTable.csv`: Data table similar to `RTstopTable.csv` except the proportion here represents the proportion of reads that stop at a site normalised to the total read coverage at that site only. 1 - this value will therefore give the reads that contain readthrough at each position.
* `modContext.txt`: Nucleotide context information for selected modified positions. Note, `pos` here is not canonical position information but ungapped alignment positions for each tRNA/cluster.
* `allModsTable.csv`: All known modifications from Modomics and `additionalMods.txt` and newly detected mim-tRNAseq modifications for each cluster. Used in new mods discovery.
* `predictedMods.csv`: Newly predicted modified sites per tRNA/cluster based on misinc-thresh.

**DESeq2**

Differential tRNA expression analyses using DESeq2. Split into two separate analyses performed on counts at the isoacceptor and isodecoder level, respectively. Both folders contain the same outputs.

* `qc-dispersions.png`: Diagnostic dispersion plot of isoacceptor/isodecoder expression dispersion before and after shrinkage of estimates towards the fitted estimates. See the DESeq2 analysis vignette here_ for details
* `qc-sampledists.png`: Sample distance heatmap based on variance stabilizing transformed counts. Hierarchical clustering used for clusering samples. Scale indicates sample distances.
* `qc-pca.png`: Principal component analysis plot for all samples according to normalized counts for isoacceptors/isodecoders.
* `vst-transformedCounts.csv`: Variace stabilizing transfored count data used for sample clustering. Also useful for comparing tRNA expression, although normalized counts are easier to understand for this purpose.
* `\*diffextr-countplot.pdf`: Count data plotted for each pairwise condition comparison. Significantly differentially expressed isoacceptors/isodecoders detected by DESeq2 (adjusted p-value < 0.05) indicated by coloured triangles.
* `\*diffexpr-results.csv`: DESeq2 differential expression results for each pairwise condition comparison. Note, every pairwise comparison output also has normalized counts for *all* samples appended as the last set of columns.

.. _here: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#dispersion-plot-and-fitting-alternatives

**indices**

Indices required for GSNAP alignment.

* `tRNAgenome\` contains the index of mature, processed tRNA transcripts.
* `snp_index\` contains the SNP index generated from modified positions, needed by GSNAP for SNP-tolerant alignment. Only generated if --snp-tolerance is enabled.
* `\*.log`: files contain log info from index generation.

**annotation**

Various files describing the tRNA trasncriptome of the genome of interest.

* `\*tRNATranscripts.fa`: processed, intron spliced, 3'-CCA appeneded, and His 5'-G appended tRNA trancript sequences in fasta format.
* `\*modificationSNPs.txt`: SNP index information for each tRNA after matching to Modomics entries for species of interest.
* `\*isoacceptorInfo.txt`: Information on isoacceptor groups and their size in genome of interest.
* `\*maturetRNA.bed': bed6 file for mature tRNA transcripts - used for coverage calculations.
* `\*clusterTranscripts_align.stk`: Stockholm align file generated by INFERNAL cmalign for tRNA sequence and structural alignments. Used for metagene coverage plots.
* `cm.log`: log file for INFERNAL cmalign algorithm

	Parameter-dependent outputs:

	* `\*clusterTranscripts.fa`: Cluster parent transcript sequence if clustering is enabled.
	* `\*isodecoderInfo.txt`: Isodecoder representative gene with size of isodecoder group (i.e. number of identitical tRNA sequences). Onlu for cluster-id < 1
	* `\*clusters.bed`: bed6 file for cluster parents. Only if clustering is enabled.
	* `\*clusterInfo.txt`: Cluster parent-child relationship for every tRNA gene, with unique cluster number and size. Only if clustering is enabled.

