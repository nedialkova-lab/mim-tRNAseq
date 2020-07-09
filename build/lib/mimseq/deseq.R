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
suppressMessages(library(dplyr))

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
    if (! control %in% orig_levels) stop("Control must be a level of data - potentially removed due to single replicate. Control condition must have replicates!")
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
cluster_id = as.numeric(args[4])
subdir = "DESeq2/"
subdir_isodecoder = "DESeq2/isodecoder"
subdir_anticodon = "DESeq2/anticodon"

# Set output directory
setwd(file.path(outdir))

# Import data from counts and sampleData
# Filter out mito counts
# Import data from counts and sampleData
# Filter out mito counts
anticodon_countdata = read.table("Anticodon_counts.txt", header=TRUE, row.names=1, check.names = FALSE)
anticodon_countdata = anticodon_countdata[!grepl("mito", rownames(anticodon_countdata)), ,drop = FALSE]
isodecoder_countdata = read.table("Isodecoder_counts.txt", header=TRUE, row.names=1, check.names = FALSE)
isodecoder_countdata = isodecoder_countdata[grepl("True", isodecoder_countdata$Single_isodecoder),]
isodecoder_countdata$Single_isodecoder = NULL
isodecoder_countdata = isodecoder_countdata[!grepl("mito", rownames(isodecoder_countdata)), ,drop = FALSE]
coldata = read.table(paste(sampleData, sep=""), header=FALSE, sep = "\t", row.names=1)

coldata = data.frame(row.names=rownames(coldata), condition = coldata[,1])

anticodon_countdata = anticodon_countdata[, order(colnames(anticodon_countdata))]
isodecoder_countdata = isodecoder_countdata[, order(colnames(isodecoder_countdata))]
coldata = coldata[order(rownames(coldata)), , drop = FALSE]

# Convert to matrix
anticodon_countdata = as.matrix(anticodon_countdata)
isodecoder_countdata = as.matrix(isodecoder_countdata)

isoacceptorFile = list.files(path="./", pattern="isoacceptorInfo.txt", full.names=T)
isoacceptorInfo = read.table(isoacceptorFile[1], header=T, row.names=1)
isoacceptorInfo = isoacceptorInfo[ , 'size', drop=F]
isoacceptorInfo$rn = rownames(isoacceptorInfo)

# cluster_id == 1 means that mim-seq clusters are already representative of isodecoders, therefore get isodecoder info directly from custerInfo
# cluster_id == '' means clustering is disabled
if (cluster_id == 1){
  isodecoderFile = list.files(path="./", pattern="clusterInfo.txt", full.names=T)
  isodecoderInfo = read.table(isodecoderFile[1], header=T)
  isodecoderInfo = isodecoderInfo[!duplicated(isodecoderInfo$parent),]
  isodecoderInfo = isodecoderInfo[, 'cluster_size', drop = F]
  colnames(isodecoderInfo) = 'size'
  isodecoderInfo$rn = rownames(isodecoderInfo)
} else {
  isodecoderFile = list.files(path="./", pattern="isodecoderInfo.txt", full.names=T)
  isodecoderInfo = read.table(isodecoderFile[1], header=T, row.names=1)
  isodecoderInfo = isodecoderInfo[ , 'size', drop=F]
  isodecoderInfo$rn = rownames(isodecoderInfo)
} 


if (nrow(coldata) == 1) {
  print("Too few samples to continue with DESeq analysis! Use raw count data in counts/ folder.")
} else{
  
  # Remove .bam or .sam and outdir from column names
  colnames(anticodon_countdata) = gsub("\\.[sb]am$", "", colnames(anticodon_countdata))
  colnames(isodecoder_countdata) = gsub("\\.[sb]am$", "", colnames(isodecoder_countdata))
  rownames(coldata) = gsub("\\.[sb]am$", "", rownames(coldata))
  colnames(anticodon_countdata) = substr(colnames(anticodon_countdata),nchar(outdir)+1,nchar(colnames(anticodon_countdata)))
  colnames(isodecoder_countdata) = substr(colnames(isodecoder_countdata),nchar(outdir)+1,nchar(colnames(isodecoder_countdata)))
  rownames(coldata) = substr(rownames(coldata),nchar(outdir)+1,nchar(rownames(coldata)))
  colnames(coldata) = "condition"
  
  # Make subdirectories
  dir.create(file.path(subdir), showWarnings = FALSE)
  dir.create(file.path(subdir, 'anticodon'), showWarnings = FALSE)
  dir.create(file.path(subdir, 'isodecoder'), showWarnings = FALSE)
  
  # Instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
  dds_anticodon = DESeqDataSetFromMatrix(countData = anticodon_countdata, colData = coldata, design = ~1)
  dds_isodecoder = DESeqDataSetFromMatrix(countData = isodecoder_countdata, colData = coldata, design = ~1)
  dds_anticodon = estimateSizeFactors(dds_anticodon)
  dds_isodecoder = estimateSizeFactors(dds_isodecoder)
  
  anticodon_counts = data.frame(row.names = rownames(anticodon_countdata))
  isodecoder_counts = data.frame(row.names = rownames(isodecoder_countdata))
  
  single_reps = FALSE
  # Remove conditions with single samples so DE analysis is not performed
  for (con in unique(coldata$condition)) {
    
    num = nrow(subset(coldata, condition == con))
    if (num == 1) {
      single_reps = TRUE
      a = rownames(coldata[coldata$condition == con, ,drop = FALSE])
      # build data frames of normalized counts for single replicte conditions
      anticodon_counts$a = counts(dds_anticodon, normalized = TRUE)[,a]
      isodecoder_counts$a = counts(dds_isodecoder, normalized = TRUE)[,a]
      # rename columns
      names(anticodon_counts)[names(anticodon_counts) == 'a'] = a
      names(isodecoder_counts)[names(isodecoder_counts) == 'a'] = a
      # remove these samples from the countdata
      isodecoder_countdata = isodecoder_countdata[,-grepl(a, colnames(isodecoder_countdata)), drop = FALSE]
      anticodon_countdata = anticodon_countdata[,-grepl(a, colnames(anticodon_countdata)), drop = FALSE]
      coldata = coldata[!rownames(coldata) == a, , drop = FALSE]
      
    }
  }
  
  if (single_reps == TRUE) {
    # Merge normalized count data with genen copy number info
    anticodon_counts$rn = rownames(anticodon_counts)
    isodecoder_counts$rn = rownames(isodecoder_counts)
    
    anticodon_counts = merge(anticodon_counts, isoacceptorInfo, by = 'rn', type = 'left')
    isodecoder_counts = merge(isodecoder_counts, isodecoderInfo, by = 'rn', type = 'left')
    
    names(anticodon_counts)[1] = "Gene"
    names(isodecoder_counts)[1] = "Gene"
    col_idx = grep("Gene", names(anticodon_counts))
    resdata_anticodon = anticodon_counts[, c(col_idx, (1:ncol(anticodon_counts))[-col_idx])]
    col_idx = grep("Gene", names(isodecoder_counts))
    resdata_isodecoder = isodecoder_counts[, c(col_idx, (1:ncol(isodecoder_counts))[-col_idx])]
    
    write.csv(resdata_anticodon, file=paste(subdir_anticodon, 'singleRep-normCounts.csv', sep="/"), row.names = FALSE)
    write.csv(resdata_isodecoder, file=paste(subdir_isodecoder, 'singleRep-normCounts.csv', sep="/"), row.names = FALSE)
  }
  
  if (dim(anticodon_countdata)[2] == 0) {
    print("DESeq2 analysis not run as there are no sample replicates!") 
  } else {
    
    # Analysis with DESeq2 ----------------------------------------------------
    
    # If only one condition, do not perform DE, but only generate normalised counts table
    if (length(unique(coldata$condition)) == 1) {
      print("Only one condition with replicates! Producing normalized counts only.")
      
      # reintialize dds in case samples were removed above
      dds_anticodon = DESeqDataSetFromMatrix(countData = anticodon_countdata, colData = coldata, design = ~1)
      dds_isodecoder = DESeqDataSetFromMatrix(countData = isodecoder_countdata, colData = coldata, design = ~1)
      dds_anticodon = estimateSizeFactors(dds_anticodon)
      dds_isodecoder = estimateSizeFactors(dds_isodecoder)
      
      anticodon_counts = as.data.frame(counts(dds_anticodon, normalized = TRUE))
      isodecoder_counts = as.data.frame(counts(dds_isodecoder, normalized = TRUE))
      anticodon_counts$rn = rownames(anticodon_counts)
      isodecoder_counts$rn = rownames(isodecoder_counts)
      
      anticodon_counts = merge(anticodon_counts, isoacceptorInfo, by = 'rn', type = 'left')
      isodecoder_counts = merge(isodecoder_counts, isodecoderInfo, by = 'rn', type = 'left')
      
      names(anticodon_counts)[1] = "Gene"
      names(isodecoder_counts)[1] = "Gene"
      col_idx = grep("Gene", names(anticodon_counts))
      resdata_anticodon = anticodon_counts[, c(col_idx, (1:ncol(anticodon_counts))[-col_idx])]
      col_idx = grep("Gene", names(isodecoder_counts))
      resdata_isodecoder = isodecoder_counts[, c(col_idx, (1:ncol(isodecoder_counts))[-col_idx])]
      
      write.csv(resdata_anticodon, file=paste(subdir_anticodon, paste(toString(unique(coldata$condition)), 'normCounts.csv', sep="-"), sep="/"), row.names = FALSE)
      write.csv(resdata_isodecoder, file=paste(subdir_isodecoder, paste(toString(unique(coldata$condition)), 'normCounts.csv', sep="-"), sep="/"), row.names = FALSE)
      
      
    } else {
      dds_anticodon = DESeqDataSetFromMatrix(countData=anticodon_countdata, colData=coldata, design=~condition)
      dds_isodecoder = DESeqDataSetFromMatrix(countData=isodecoder_countdata, colData=coldata, design=~condition)
      
      # Run the DESeq pipeline
      dds_anticodon = DESeq(dds_anticodon)
      dds_isodecoder = DESeq(dds_isodecoder)
      
      # count tables with mean per condition for dot plots below
      baseMeanPerLvl_anticodon = as.data.frame(sapply(levels(dds_anticodon$condition), function(lvl) rowMeans(counts(dds_anticodon,normalized=TRUE)[,dds_anticodon$condition == lvl] )))
      baseMeanPerLvl_isodecoder = as.data.frame(sapply(levels(dds_isodecoder$condition), function(lvl) rowMeans(counts(dds_isodecoder,normalized=TRUE)[,dds_isodecoder$condition == lvl] )))
      
      png(paste(subdir_anticodon, "qc-dispersions.png", sep = "/"), 1000, 1000, pointsize=20)
      plotDispEsts(dds_anticodon, main="Dispersion plot")
      dev.off()
      
      png(paste(subdir_isodecoder, "qc-dispersions.png", sep = "/"), 1000, 1000, pointsize=20)
      plotDispEsts(dds_isodecoder, main="Dispersion plot")
      dev.off()
      
      # Variance stabilizing transformation for clustering/heatmaps, etc
      
      vsd_anticodon = varianceStabilizingTransformation(dds_anticodon, blind=FALSE)
      write.csv(assay(vsd_anticodon), file = paste(subdir_anticodon, "vst-transformedCounts.csv", sep = "/"))
      
      vsd_isodecoder = varianceStabilizingTransformation(dds_isodecoder, blind=FALSE)
      write.csv(assay(vsd_isodecoder), file = paste(subdir_isodecoder, "vst-transformedCounts.csv", sep = "/"))
      
      sampleDists_anticodon = dist(t(assay(vsd_anticodon)))
      sampleDistMatrix_anticodon <- as.matrix(sampleDists_anticodon)
      rownames(sampleDistMatrix_anticodon) <- vsd_anticodon$condition
      colnames(sampleDistMatrix_anticodon) <- NULL
      colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
      pheatmap(sampleDistMatrix_anticodon,
               anticodoning_distance_rows=sampleDists_anticodon,
               anticodoning_distance_cols=sampleDists_anticodon,
               col=colors, filename=paste(subdir_anticodon,"qc-sampledists.png",sep="/"))
      
      sampleDists_isodecoder = dist(t(assay(vsd_isodecoder)))
      sampleDistMatrix_isodecoder <- as.matrix(sampleDists_isodecoder)
      rownames(sampleDistMatrix_isodecoder) <- vsd_isodecoder$condition
      colnames(sampleDistMatrix_isodecoder) <- NULL
      colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
      pheatmap(sampleDistMatrix_isodecoder,
               isodecoder_distance_rows=sampleDists_isodecoder,
               isodecoder_distance_cols=sampleDists_isodecoder,
               col=colors, filename=paste(subdir_isodecoder,"qc-sampledists.png",sep="/"))
      
      # Principal components analysis
      
      pcaData_anticodon <- plotPCA(vsd_anticodon, intgroup="condition", returnData=TRUE)
      percentVar <- round(100 * attr(pcaData_anticodon, "percentVar"))
      ggplot(pcaData_anticodon, aes(PC1, PC2, color=condition)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()
      ggsave(paste(subdir_anticodon, "qc-pca.png", sep="/"), height = 7, width = 8)
      
      pcaData_isodecoder <- plotPCA(vsd_isodecoder, intgroup="condition", returnData=TRUE)
      percentVar <- round(100 * attr(pcaData_isodecoder, "percentVar"))
      ggplot(pcaData_isodecoder, aes(PC1, PC2, color=condition)) +
        geom_point(size=3) +
        xlab(paste0("PC1: ",percentVar[1],"% variance")) +
        ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
        coord_fixed()
      ggsave(paste(subdir_isodecoder, "qc-pca.png", sep="/"), height = 7, width = 8)
      
      # Gene expression heatmaps
      
      cols = as.data.frame(colData(dds_isodecoder)[,'condition'])
      rownames(cols) <- colnames(dds_isodecoder)
      names(cols) <- "Condition"
      isodecoder_mat <- assay(vsd_isodecoder)
      isodecoder_mat = (isodecoder_mat - rowMeans(isodecoder_mat))
      isodecoder_mat[is.na(isodecoder_mat)] = 0
      isodecoder_mat = isodecoder_mat[,order(cols$Condition)]
      pheatmap(isodecoder_mat, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, 
               annotation_col = cols, filename = paste(subdir_isodecoder, "Isodecoder_vst_hm.png", sep="/"), 
               color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), 
               show_colnames = FALSE)
      
      cols = as.data.frame(colData(dds_anticodon)[,'condition'])
      rownames(cols) = colnames(dds_anticodon)
      names(cols) <- "Condition"
      anticodon_mat <- assay(vsd_anticodon)
      anticodon_mat = (anticodon_mat - rowMeans(anticodon_mat))
      rownames(anticodon_mat) = paste(unlist(strsplit(rownames(anticodon_mat), "-"))[3*(1:nrow(anticodon_mat))-1],unlist(strsplit(rownames(anticodon_mat), "-"))[3*(1:nrow(anticodon_mat))], sep = "-")
      anticodon_mat[is.na(anticodon_mat)] = 0
      anticodon_mat = anticodon_mat[,order(cols$Condition)]
      pheatmap(anticodon_mat, cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = cols, 
               filename = paste(subdir_anticodon, "Anticodon_vst_hm.png", sep="/"), 
               color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100), show_colnames = FALSE)
      
      # Get combinations of coditions for various DE contrasts
      ordered_levels = levels(lastlevel(unique(dds_anticodon$condition), control_cond))
      combinations = combn(ordered_levels, 2, simplify=FALSE)
      
      # For each contrast...
      for (i in 1:length(combinations)) {
        # Get differential expression results
        res_anticodon = results(dds_anticodon, contrast=c("condition",as.vector(combinations[[i]])))
        res_isodecoder = results(dds_isodecoder, contrast=c("condition",as.vector(combinations[[i]])))
        
        ## Order by adjusted p-value
        res_anticodon = res_anticodon[order(res_anticodon$padj), ]
        res_isodecoder = res_isodecoder[order(res_isodecoder$padj), ]
        res_anticodon$rn = rownames(res_anticodon)
        res_isodecoder$rn = rownames(res_isodecoder)
        count_df_anticodon = as.data.frame(counts(dds_anticodon, normalized=TRUE))
        count_df_isodecoder = as.data.frame(counts(dds_isodecoder, normalized=TRUE))
        count_df_anticodon$rn = rownames(count_df_anticodon)
        count_df_isodecoder$rn = rownames(count_df_isodecoder)
        
        ## Merge with normalized count data
        resdata_anticodon = join_all(list(as.data.frame(res_anticodon), count_df_anticodon, isoacceptorInfo), by="rn", type = 'left')
        resdata_isodecoder = join_all(list(as.data.frame(res_isodecoder), count_df_isodecoder, isodecoderInfo), by="rn", type = 'left')
        names(resdata_anticodon)[7] = "Gene"
        names(resdata_isodecoder)[7] = "Gene"
        col_idx = grep("Gene", names(resdata_anticodon))
        resdata_anticodon = resdata_anticodon[, c(col_idx, (1:ncol(resdata_anticodon))[-col_idx])]
        col_idx = grep("Gene", names(resdata_isodecoder))
        resdata_isodecoder = resdata_isodecoder[, c(col_idx, (1:ncol(resdata_isodecoder))[-col_idx])]
        
        # Count plots
        # add significance to baseMean matrices for current contrast
        baseMeanPerLvl_anticodon$sig = res_anticodon[rownames(baseMeanPerLvl_anticodon),'padj'] < 0.05
        baseMeanPerLvl_anticodon$sig = !is.na(baseMeanPerLvl_anticodon$sig) & baseMeanPerLvl_anticodon$sig # deal with NAs turning them into FALSE
        baseMeanPerLvl_isodecoder$sig = res_isodecoder[rownames(baseMeanPerLvl_isodecoder),'padj'] < 0.05
        baseMeanPerLvl_isodecoder$sig = !is.na(baseMeanPerLvl_isodecoder$sig) & baseMeanPerLvl_isodecoder$sig # deal with NAs turning them into FALSE
        
        # add direction of DE to baseMean
        baseMeanPerLvl_anticodon$direction = sign(res_anticodon[rownames(baseMeanPerLvl_anticodon), 'log2FoldChange'])
        baseMeanPerLvl_anticodon[is.na(baseMeanPerLvl_anticodon$direction),'direction'] = 0
        baseMeanPerLvl_isodecoder$direction = sign(res_isodecoder[rownames(baseMeanPerLvl_isodecoder), 'log2FoldChange'])
        baseMeanPerLvl_isodecoder[is.na(baseMeanPerLvl_isodecoder$direction),'direction'] = 0
        
        baseMeanPerLvl_anticodon$comb = paste(baseMeanPerLvl_anticodon$sig, baseMeanPerLvl_anticodon$direction, sep='')
        baseMeanPerLvl_anticodon$comb[which(baseMeanPerLvl_anticodon$comb %in% c('FALSE1','FALSE0','FALSE-1'))] = "FALSE"
        baseMeanPerLvl_isodecoder$comb = paste(baseMeanPerLvl_isodecoder$sig, baseMeanPerLvl_isodecoder$direction, sep='')
        baseMeanPerLvl_isodecoder$comb[which(baseMeanPerLvl_isodecoder$comb %in% c('FALSE1','FALSE0','FALSE-1'))] = "FALSE"
        
        anticodon_lin_mod = lm(baseMeanPerLvl_anticodon[,combinations[[i]][1]] ~ baseMeanPerLvl_anticodon[,combinations[[i]][2]])
        isodecoder_lin_mod = lm(baseMeanPerLvl_isodecoder[,combinations[[i]][1]] ~ baseMeanPerLvl_isodecoder[,combinations[[i]][2]])
        anticodon_cor = format(cor(baseMeanPerLvl_anticodon[,combinations[[i]][1]], baseMeanPerLvl_anticodon[,combinations[[i]][2]]), digits = 2)
        isodecoder_cor = format(cor(baseMeanPerLvl_isodecoder[,combinations[[i]][1]], baseMeanPerLvl_isodecoder[,combinations[[i]][2]]), digits = 2)
        
        # Convert x and y variables into symbols  to handle str and numeric variable names
        x = sym(combinations[[i]][2])
        y = sym(combinations[[i]][1])
        anticodon_data = subset(baseMeanPerLvl_anticodon, select = c(combinations[[i]], "sig", "comb"))
        isodecoder_data = subset(baseMeanPerLvl_isodecoder, select = c(combinations[[i]], "sig", "comb"))
        
        anticodon_dot = ggplot(anticodon_data, aes_string(x, y)) +
          geom_point(aes(color = comb, shape = comb, size = sig)) +
          scale_x_log10() + scale_y_log10() + 
          #geom_smooth(method = 'lm', se = TRUE, alpha = 0.5, color = '#3182bd', fill = 'grey') +
          geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = '#3182bd', alpha = 0.8) + 
          scale_color_manual('Differential expression', labels = c("None", "Down", "Up"), values = c("darkgrey", "#f28f3b", "#4daf4a")) + 
          scale_shape_manual('Differential expression', labels = c("None", "Down", "Up"), values = c(19, 17, 17)) + 
          scale_size_manual(values = c(1,2), guide = FALSE) + theme_bw() +
          labs(x = paste('log10', combinations[[i]][2], 'counts', sep = ' '), y = paste('log10', combinations[[i]][1], 'counts', sep = ' ')) +
          annotate("label", 0, Inf, hjust = 0, vjust = 1, label = paste("italic(r) == ", anticodon_cor), parse = TRUE)
        
        isodecoder_dot = ggplot(isodecoder_data, aes_string(x, y)) +
          geom_point(aes(color = comb, shape = comb, size = sig)) +
          scale_x_log10() + scale_y_log10() + 
          #geom_smooth(method = 'lm', se = TRUE, alpha = 0.5, color = '#3182bd', fill = 'grey') +
          geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = '#3182bd', alpha = 0.8) + 
          scale_color_manual('Differential expression', labels = c("None", "Down", "Up"), values = c("darkgrey", "#f28f3b", "#4daf4a")) + 
          scale_shape_manual('Differential expression', labels = c("None", "Down", "Up"), values = c(19, 17, 17)) + 
          scale_size_manual(values = c(1,2), guide = FALSE) + theme_bw() +
          labs(x = paste('log10', combinations[[i]][2], 'counts', sep = ' '), y = paste('log10', combinations[[i]][1], 'counts', sep = ' ')) +
          annotate("label", 0, Inf, hjust = 0, vjust = 1, label = paste("italic(r) == ", isodecoder_cor), parse = TRUE)
        
        ggsave(paste(subdir_anticodon, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-countplot.pdf", sep="_"), sep="/"), anticodon_dot, height = 5, width = 8)
        ggsave(paste(subdir_isodecoder, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-countplot.pdf", sep="_"), sep="/"), isodecoder_dot, height = 5, width = 8)
        
        ## Write results
        write.csv(resdata_anticodon, file=paste(subdir_anticodon, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-results.csv", sep="_"), sep="/"))
        write.csv(resdata_isodecoder, file=paste(subdir_isodecoder, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-results.csv", sep="_"), sep="/"))
      }
      
      # ## MA plot
      
      # png(paste(subdir_anticodon, "diffexpr-maplot.png", sep="/"), 1500, 1000, pointsize=20)
      # plotMA(dds_anticodon)
      # dev.off()
      
      # png(paste(subdir_isodecoder, "diffexpr-maplot.png", sep="/"), 1500, 1000, pointsize=20)
      # plotMA(dds_isodecoder)
      # dev.off()
    }
  }
}
