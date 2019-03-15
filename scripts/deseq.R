#!/usr/bin/Rscript

## RNA-seq analysis with DESeq2
## based on Stephen Turner's script, @genetics_blog
## https://gist.github.com/stephenturner/f60c1934405c127f09a6

suppressMessages(library(DESeq2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(ggplot2))
suppressMessages(library(calibrate))
suppressMessages(library(plyr))
suppressMessages(library(grid))

## Volcano plot with "significant" genes labeled
volcanoplot = function (res, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=FALSE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=19, main=main, xlab="log2 Fold Change", ylab="-log10(p-value)", ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=19, col="#d95f02", ...))
  if (labelsig) {
    require(calibrate)
    #with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
    with(subset(res, padj<sigthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR < ",sigthresh,sep="")), pch=19, col="#d95f02", bg="white")
}

# Function to order control/WT condition last in combinations so that contrasts for DE are always mutant/condition vs WT/control
lastlevel = function(f, control) {
    if (!is.factor(f)) stop("input for contrast combinations not a factor")
    orig_levels = levels(f)
    if (! control %in% orig_levels) stop("control must be a level of f")
    new_levels = c(setdiff(orig_levels, control), control)
    factor(f, levels = new_levels)
}

# for getting r-squared values back for plotting
lm_eqn <- function(linear_mod){
  l <- list(r2 = format(summary(linear_mod)$r.squared, digits = 2));
  eq <- substitute(~~italic(r)^2~"="~r2,l)
  
  as.character(as.expression(eq));                 
}

# Arguments
args = commandArgs(trailingOnly = TRUE)
outdir = args[1]
sampleData = args[2]
control_cond = args[3]
subdir = "DESeq2/"
subdir_cluster = "DESeq2/cluster"
subdir_anticodon = "DESeq2/anticodon"

# Output directory
dir.create(file.path(outdir, subdir), showWarnings = FALSE)
setwd(file.path(outdir))
dir.create(file.path(subdir, 'cluster'), showWarnings = FALSE)
dir.create(file.path(subdir, 'anticodon'), showWarnings = FALSE)

# Import data from featureCounts and sampleData
cluster_countdata = read.table("counts.txt", header=TRUE, row.names=1, check.names = FALSE)
anticodon_countdata = read.table("Anticodon_counts.txt", header=TRUE, row.names=1, check.names = FALSE)
coldata = read.table(paste(sampleData, sep=""), header=FALSE, sep = "\t", row.names=1)

# Remove first five columns (chr, start, end, strand, length)
cluster_countdata = cluster_countdata[ ,6:ncol(cluster_countdata), drop = FALSE]
coldata = data.frame(row.names=rownames(coldata), condition = coldata[,1])

# Remove .bam or .sam and outdir from column names
colnames(cluster_countdata) = gsub("\\.[sb]am$", "", colnames(cluster_countdata))
colnames(anticodon_countdata) = gsub("\\.[sb]am$", "", colnames(anticodon_countdata))
rownames(coldata) = gsub("\\.[sb]am$", "", rownames(coldata))
colnames(cluster_countdata) = substr(colnames(cluster_countdata),nchar(outdir)+1,nchar(colnames(cluster_countdata)))
colnames(anticodon_countdata) = substr(colnames(anticodon_countdata),nchar(outdir)+1,nchar(colnames(anticodon_countdata)))
rownames(coldata) = substr(rownames(coldata),nchar(outdir)+1,nchar(rownames(coldata)))
colnames(coldata) = "condition"

# Duplicate conditions with single samples to avoid errors - NOTE: DE results from these samples should not be interepretated either way!!
for (con in unique(coldata$condition)) {

  num = nrow(subset(coldata, condition == con))
  if (num == 1) {
    sample = rownames(coldata[which(coldata$condition == con), , drop = FALSE])
    new_sample = paste(sample, '2', sep = '')
    con = as.data.frame(con)
    rownames(con) = new_sample
    names(con) = 'condition'
    coldata = rbind(con, coldata)
    cluster_countdata[[new_sample]] = cluster_countdata[[sample]]
    anticodon_countdata[[new_sample]] = anticodon_countdata[[sample]]
  }
}

cluster_countdata = cluster_countdata[, order(colnames(cluster_countdata))]
anticodon_countdata = anticodon_countdata[, order(colnames(anticodon_countdata))]
coldata = coldata[order(rownames(coldata)), , drop = FALSE]

# Convert to matrix
cluster_countdata = as.matrix(cluster_countdata)
anticodon_countdata = as.matrix(anticodon_countdata)

# Get cluster and isoacceptor data
clusterFile = list.files(path="./", pattern="clusterInfo.txt", full.names=T)
if (length(clusterFile) == 1) {
  clusterInfo = read.table(clusterFile[1], header=T, row.names=1)
  clusterInfo = clusterInfo[ , 'cluster_size', drop=F]
  clusterInfo$rn = rownames(clusterInfo)
} else if (length(clusterFile == 0)) {
  clusterInfo = data.frame(cluster_size = 1, rn = rownames(dds))
  rownames(clusterInfo) = clusterInfo$rn
}

isoacceptorFile = list.files(path="./", pattern="isoacceptorInfo.txt", full.names=T)
isoacceptorInfo = read.table(isoacceptorFile[1], header=T, row.names=1)
isoacceptorInfo = isoacceptorInfo[ , 'size', drop=F]
isoacceptorInfo$rn = rownames(isoacceptorInfo)

# Analysis with DESeq2 ----------------------------------------------------

# Instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
# If only one condition, do not perform DE, but only generate normalised counts table
if (length(unique(coldata$condition)) == 1) {
  dds_cluster = DESeqDataSetFromMatrix(countData = cluster_countdata, colData = coldata, design = ~1)
  dds_anticodon = DESeqDataSetFromMatrix(countData = anticodon_countdata, colData = coldata, design = ~1)
  dds_cluster = estimateSizeFactors(dds_cluster)
  dds_anticodon = estimateSizeFactors(dds_anticodon)

  cluster_counts = as.data.frame(counts(dds_cluster, normalized = TRUE))
  anticodon_counts = as.data.frame(counts(dds_anticodon, normalized = TRUE))
  cluster_counts$rn = rownames(cluster_counts)
  anticodon_counts$rn = rownames(anticodon_counts)

  cluster_counts = merge(cluster_counts, clusterInfo, by = 'rn', type = 'left')
  anticodon_counts = merge(anticodon_counts, isoacceptorInfo, by = 'rn', type = 'left')

  names(cluster_counts)[1] = "Gene"
  names(anticodon_counts)[1] = "Gene"
  col_idx = grep("Gene", names(cluster_counts))
  resdata_cluster = cluster_counts[, c(col_idx, (1:ncol(cluster_counts))[-col_idx])]
  col_idx = grep("Gene", names(anticodon_counts))
  resdata_anticodon = anticodon_counts[, c(col_idx, (1:ncol(anticodon_counts))[-col_idx])]

  write.csv(resdata_cluster, file=paste(subdir_cluster, paste(toString(unique(coldata$condition)), "normCounts.csv", sep = '-'), sep="/"), row.names = FALSE)
  write.csv(resdata_anticodon, file=paste(subdir_anticodon, paste(toString(unique(coldata$condition)), 'normCounts.csv', sep="-"), sep="/"), row.names = FALSE)

} else {
  dds_cluster = DESeqDataSetFromMatrix(countData=cluster_countdata, colData=coldata, design=~condition)
  dds_anticodon = DESeqDataSetFromMatrix(countData=anticodon_countdata, colData=coldata, design=~condition)
  # Run the DESeq pipeline
  dds_cluster = DESeq(dds_cluster)
  dds_anticodon = DESeq(dds_anticodon)

  # count tables with mean per condition for dot plots below
  baseMeanPerLvl_cluster = as.data.frame(sapply(levels(dds_cluster$condition), function(lvl) rowMeans(counts(dds_cluster,normalized=TRUE)[,dds_cluster$condition == lvl] )))
  baseMeanPerLvl_anticodon = as.data.frame(sapply(levels(dds_anticodon$condition), function(lvl) rowMeans(counts(dds_anticodon,normalized=TRUE)[,dds_anticodon$condition == lvl] )))

  # Plot dispersions
  png(paste(subdir_cluster, "qc-dispersions.png", sep = "/"), 1000, 1000, pointsize=20)
  plotDispEsts(dds_cluster, main="Dispersion plot")
  dev.off()

  png(paste(subdir_anticodon, "qc-dispersions.png", sep = "/"), 1000, 1000, pointsize=20)
  plotDispEsts(dds_anticodon, main="Dispersion plot")
  dev.off()

  # Variance stabilizing transformation for clustering/heatmaps, etc
  vsd_cluster = varianceStabilizingTransformation(dds_cluster, blind=FALSE)
  write.csv(assay(vsd_cluster), file = paste(subdir_cluster, "vst-transformedCounts.csv", sep = "/"))

  vsd_anticodon = varianceStabilizingTransformation(dds_anticodon, blind=FALSE)
  write.csv(assay(vsd_anticodon), file = paste(subdir_anticodon, "vst-transformedCounts.csv", sep = "/"))

  # Sample distance heatmap
  sampleDists_cluster = dist(t(assay(vsd_cluster)))
  sampleDistMatrix_cluster <- as.matrix(sampleDists_cluster)
  rownames(sampleDistMatrix_cluster) <- vsd_cluster$condition
  colnames(sampleDistMatrix_cluster) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix_cluster,
           clustering_distance_rows=sampleDists_cluster,
           clustering_distance_cols=sampleDists_cluster,
           col=colors, filename=paste(subdir_cluster,"qc-sampledists.png",sep="/"))

  sampleDists_anticodon = dist(t(assay(vsd_anticodon)))
  sampleDistMatrix_anticodon <- as.matrix(sampleDists_anticodon)
  rownames(sampleDistMatrix_anticodon) <- vsd_anticodon$condition
  colnames(sampleDistMatrix_anticodon) <- NULL
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix_anticodon,
           anticodoning_distance_rows=sampleDists_anticodon,
           anticodoning_distance_cols=sampleDists_anticodon,
           col=colors, filename=paste(subdir_anticodon,"qc-sampledists.png",sep="/"))

  # Principal components analysis
  pcaData_cluster <- plotPCA(vsd_cluster, intgroup="condition", returnData=TRUE)
  percentVar <- round(100 * attr(pcaData_cluster, "percentVar"))
  ggplot(pcaData_cluster, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  ggsave(paste(subdir_cluster, "qc-pca.png", sep="/"), height = 7, width = 8)

  pcaData_anticodon <- plotPCA(vsd_anticodon, intgroup="condition", returnData=TRUE)
  percentVar <- round(100 * attr(pcaData_anticodon, "percentVar"))
  ggplot(pcaData_anticodon, aes(PC1, PC2, color=condition)) +
    geom_point(size=3) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
    coord_fixed()
  ggsave(paste(subdir_anticodon, "qc-pca.png", sep="/"), height = 7, width = 8)

  # Get combinations of coditions for various DE contrasts
  ordered_levels = levels(lastlevel(unique(dds_cluster$condition), control_cond))
  combinations = combn(ordered_levels, 2, simplify=FALSE)

  # For each contrast...
  for (i in 1:length(combinations)) {
    # Get differential expression results
    res_cluster = results(dds_cluster, contrast=c("condition",as.vector(combinations[[i]])))
    res_anticodon = results(dds_anticodon, contrast=c("condition",as.vector(combinations[[i]])))
    
    ## Order by adjusted p-value
    res_cluster = res_cluster[order(res_cluster$padj), ]
    res_anticodon = res_anticodon[order(res_anticodon$padj), ]
    res_cluster$rn = rownames(res_cluster)
    res_anticodon$rn = rownames(res_anticodon)
    count_df_cluster = as.data.frame(counts(dds_cluster, normalized=TRUE))
    count_df_anticodon = as.data.frame(counts(dds_anticodon, normalized=TRUE))
    count_df_cluster$rn = rownames(count_df_cluster)
    count_df_anticodon$rn = rownames(count_df_anticodon)
    
    ## Merge with normalized count data
    resdata_cluster = join_all(list(as.data.frame(res_cluster), count_df_cluster, clusterInfo), by="rn", type = 'left')
    resdata_anticodon = join_all(list(as.data.frame(res_anticodon), count_df_anticodon, isoacceptorInfo), by="rn", type = 'left')
    names(resdata_cluster)[7] = "Gene"
    names(resdata_anticodon)[7] = "Gene"
    col_idx = grep("Gene", names(resdata_cluster))
    resdata_cluster = resdata_cluster[, c(col_idx, (1:ncol(resdata_cluster))[-col_idx])]
    col_idx = grep("Gene", names(resdata_anticodon))
    resdata_anticodon = resdata_anticodon[, c(col_idx, (1:ncol(resdata_anticodon))[-col_idx])]

    # Count plots
    # add significance to baseMean matrices for current contrast
    baseMeanPerLvl_cluster$sig = res_cluster[rownames(baseMeanPerLvl_cluster),'padj'] < 0.05
    baseMeanPerLvl_cluster$sig = !is.na(baseMeanPerLvl_cluster$sig) & baseMeanPerLvl_cluster$sig # deal with NAs turning them into FALSE
    baseMeanPerLvl_anticodon$sig = res_anticodon[rownames(baseMeanPerLvl_anticodon),'padj'] < 0.05
    baseMeanPerLvl_anticodon$sig = !is.na(baseMeanPerLvl_anticodon$sig) & baseMeanPerLvl_anticodon$sig # deal with NAs turning them into FALSE

    # add direction of DE to baseMean
    baseMeanPerLvl_cluster$direction = sign(res_cluster[rownames(baseMeanPerLvl_cluster), 'log2FoldChange'])
    baseMeanPerLvl_cluster[is.na(baseMeanPerLvl_cluster$direction),'direction'] = 0
    baseMeanPerLvl_anticodon$direction = sign(res_anticodon[rownames(baseMeanPerLvl_anticodon), 'log2FoldChange'])
    baseMeanPerLvl_anticodon[is.na(baseMeanPerLvl_anticodon$direction),'direction'] = 0

    baseMeanPerLvl_cluster$comb = paste(baseMeanPerLvl_cluster$sig, baseMeanPerLvl_cluster$direction, sep='')
    baseMeanPerLvl_cluster$comb[which(baseMeanPerLvl_cluster$comb %in% c('FALSE1','FALSE0','FALSE-1'))] = "FALSE"
    baseMeanPerLvl_anticodon$comb = paste(baseMeanPerLvl_anticodon$sig, baseMeanPerLvl_anticodon$direction, sep='')
    baseMeanPerLvl_anticodon$comb[which(baseMeanPerLvl_anticodon$comb %in% c('FALSE1','FALSE0','FALSE-1'))] = "FALSE"

    cluster_lin_mod = lm(baseMeanPerLvl_cluster[,combinations[[i]][1]] ~ baseMeanPerLvl_cluster[,combinations[[i]][2]])
    anticodon_lin_mod = lm(baseMeanPerLvl_anticodon[,combinations[[i]][1]] ~ baseMeanPerLvl_anticodon[,combinations[[i]][2]])

    cluster_dot = ggplot(subset(baseMeanPerLvl_cluster, select = c(combinations[[i]], "sig", "comb")), aes_string(x = combinations[[i]][2], y = combinations[[i]][1])) +
      geom_point(aes(color = comb, shape = comb, size = sig)) +
      scale_x_log10() + scale_y_log10() +
      #geom_smooth(method = 'lm', se = TRUE, alpha = 0.5, color = '#3182bd', fill = 'grey') +
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = '#3182bd', alpha = 0.8) + 
      scale_color_manual('Differential expression', labels = c("None", "Down", "Up"), values = c("darkgrey", "#f28f3b", "#4daf4a")) + 
      scale_shape_manual('Differential expression', labels = c("None", "Down", "Up"), values = c(19, 17, 17)) + 
      scale_size_manual(values = c(1,2), guide = FALSE) + theme_bw() + 
      labs(x = paste('log10', combinations[[i]][2], 'counts', sep = ' '), y = paste('log10', combinations[[i]][1], 'counts', sep = ' ')) +
      annotate("label", 0, Inf, hjust = 0, vjust = 1, label = lm_eqn(cluster_lin_mod), parse = TRUE)

    anticodon_dot = ggplot(subset(baseMeanPerLvl_anticodon, select = c(combinations[[i]], "sig", "comb")), aes_string(x = combinations[[i]][2], y = combinations[[i]][1])) +
      geom_point(aes(color = comb, shape = comb, size = sig)) +
      scale_x_log10() + scale_y_log10() + 
      #geom_smooth(method = 'lm', se = TRUE, alpha = 0.5, color = '#3182bd', fill = 'grey') +
      geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = '#3182bd', alpha = 0.8) + 
      scale_color_manual('Differential expression', labels = c("None", "Down", "Up"), values = c("darkgrey", "#f28f3b", "#4daf4a")) + 
      scale_shape_manual('Differential expression', labels = c("None", "Down", "Up"), values = c(19, 17, 17)) + 
      scale_size_manual(values = c(1,2), guide = FALSE) + theme_bw() +
      labs(x = paste('log10', combinations[[i]][2], 'counts', sep = ' '), y = paste('log10', combinations[[i]][1], 'counts', sep = ' ')) +
      annotate("label", 0, Inf, hjust = 0, vjust = 1, label = lm_eqn(anticodon_lin_mod), parse = TRUE)

    ggsave(paste(subdir_cluster, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-countplot.pdf", sep="_"), sep="/"), cluster_dot, height = 5, width = 8)
    ggsave(paste(subdir_anticodon, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-countplot.pdf", sep="_"), sep="/"), anticodon_dot, height = 5, width = 8)

    # Volcano plots
    # pdf(paste(subdir_cluster, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-volcanoplot.pdf", sep="_"), sep="/"), width=15, height=12, pointsize=20)
    # volcanoplot(resdata_cluster, main=paste(combinations[[i]], collapse=" vs "), sigthresh=0.05, textcx=.8, xlim=c(-max(abs(resdata_anticodon$log2FoldChange), na.rm=TRUE)-0.2, max(abs(resdata_anticodon$log2FoldChange), na.rm=TRUE)+0.2))
    # dev.off()

    # pdf(paste(subdir_anticodon, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-volcanoplot.pdf", sep="_"), sep="/"), width=15, height=12, pointsize=20)
    # volcanoplot(resdata_anticodon, main=paste(combinations[[i]], collapse=" vs "), sigthresh=0.05, textcx=.8, xlim=c(-max(abs(resdata_anticodon$log2FoldChange), na.rm=TRUE)-0.2, max(abs(resdata_anticodon$log2FoldChange), na.rm=TRUE)+0.2))
    # dev.off()
    
    ## Write results
    write.csv(resdata_cluster, file=paste(subdir_cluster, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-results.csv", sep="_"), sep="/"))
    write.csv(resdata_anticodon, file=paste(subdir_anticodon, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-results.csv", sep="_"), sep="/"))
  }

  ## MA plot
  png(paste(subdir_cluster, "diffexpr-maplot.png", sep="/"), 1500, 1000, pointsize=20)
  plotMA(dds_cluster)
  dev.off()

  png(paste(subdir_anticodon, "diffexpr-maplot.png", sep="/"), 1500, 1000, pointsize=20)
  plotMA(dds_anticodon)
  dev.off()
}

