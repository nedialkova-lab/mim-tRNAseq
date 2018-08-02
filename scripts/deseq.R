#!/usr/bin/Rscript

## RNA-seq analysis with DESeq2
## based on Stephen Turner's script, @genetics_blog
## https://gist.github.com/stephenturner/f60c1934405c127f09a6

library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(calibrate)

# Arguments
args = commandArgs(trailingOnly = TRUE)
outdir = args[1]
sampleData = args[2]
subdir = "DESeq2"

# Output directory
dir.create(file.path(outdir, subdir), showWarnings = FALSE)
setwd(file.path(outdir))

message("\n+----------------------------------------------+
| Differential expression analysis with DESeq2 |
+----------------------------------------------+\n")

# Import data from featureCounts and sampleData
countdata = read.table(paste("counts.txt", sep=""), header=TRUE, row.names=1)
coldata = read.table(paste(sampleData, sep=""), header=FALSE, sep = "\t", row.names=1)

# Remove first five columns (chr, start, end, strand, length)
countdata = countdata[ ,6:ncol(countdata)]
coldata = data.frame(row.names=rownames(coldata), condition = coldata[,1])

# Remove .bam or .sam and outdir from filenames
colnames(countdata) = gsub("\\.[sb]am$", "", colnames(countdata))
rownames(coldata) = gsub("\\.[sb]am$", "", rownames(coldata))
colnames(countdata) = substr(colnames(countdata),nchar(outdir)+1,nchar(colnames(countdata)))
rownames(coldata) = substr(rownames(coldata),nchar(outdir)+1,nchar(rownames(coldata)))
colnames(coldata) = "condition"

# Convert to matrix
countdata = as.matrix(countdata)

# Analysis with DESeq2 ----------------------------------------------------

# Instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
dds = DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

# Run the DESeq pipeline
dds = DESeq(dds)

# Plot dispersions
png(paste(subdir, "qc-dispersions.png", sep = "/"), 1000, 1000, pointsize=20)
plotDispEsts(dds, main="Dispersion plot")
dev.off()

# Variance stabilizing transformation for clustering/heatmaps, etc
vsd = varianceStabilizingTransformation(dds, blind=FALSE)
write.csv(assay(vsd), file = paste(subdir, "vst-transformedCounts.csv", sep = "/"))

# Sample distance heatmap
sampleDists = dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, filename=paste(subdir,"qc-sampledists.png",sep="/"))

# Principal components analysis
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()
ggsave(paste(subdir, "qc-pca.png", sep="/"), height = 7, width = 8)

# Get combinations of coditions for various DE contrasts
combinations = combn(unique(dds$condition), 2, simplify=FALSE)

# For each contrast...
for (i in 1:length(combinations)) {
  # Get differential expression results
  res = results(dds, contrast=c("condition",as.vector(combinations[[i]])))
  ## Order by adjusted p-value
  res = res[order(res$padj), ]
  ## Merge with normalized count data
  resdata = merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
  names(resdata)[1] = "Gene"
  ## Write results
  write.csv(resdata, file=paste(subdir, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-results.csv", sep="_"), sep="/"))
}

## MA plot
png(paste(subdir, "diffexpr-maplot.png", sep="/"), 1500, 1000, pointsize=20)
plotMA(dds)
dev.off()

## Volcano plot with "significant" genes labeled
volcanoplot = function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}
png(paste(subdir, "diffexpr-volcanoplot.png", sep="/"), 1200, 1000, pointsize=20)
volcanoplot(resdata, lfcthresh=1, sigthresh=0.05, textcx=.8, xlim=c(-2.3, 2))
dev.off()


message("DESeq2 outputs located in: ", getwd())