#!/usr/bin/Rscript

## RNA-seq analysis with DESeq2
## based on Stephen Turner's script, @genetics_blog
## https://gist.github.com/stephenturner/f60c1934405c127f09a6

suppressMessages(library(DESeq2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(ggplot2))
suppressMessages(library(calibrate))
suppressMessages(library(plyr))
suppressMessages(library(grid))
suppressMessages(library(dplyr))
suppressMessages(library(circlize))

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
p_adj = as.numeric(args[5])
subdir = "DESeq2/"
subdir_isodecoder = "DESeq2/isodecoder"
subdir_anticodon = "DESeq2/anticodon"

# Set output directory
setwd(file.path(outdir))

# Import data from counts and sampleData
# Filter out mito counts
# Import data from counts and sampleData
# Filter out mito counts
anticodon_countdata = read.table("Anticodon_counts_raw.txt", header=TRUE, row.names=1, check.names = FALSE)
anticodon_countdata = anticodon_countdata[!grepl("mito", rownames(anticodon_countdata)), ,drop = FALSE]
isodecoder_countdata = read.table("Isodecoder_counts_raw.txt", header=TRUE, row.names=1, check.names = FALSE)
isodecoder_countdata = isodecoder_countdata[grepl("True", isodecoder_countdata$Single_isodecoder),]
isodecoder_countdata$Single_isodecoder = NULL
isodecoder_countdata$size = NULL
isodecoder_countdata$parent = NULL
isodecoder_countdata = isodecoder_countdata[!grepl("mito", rownames(isodecoder_countdata)), ,drop = FALSE]
coldata = read.table(paste(sampleData, sep=""), header=FALSE, sep = "\t", row.names=1)

coldata = data.frame(row.names=rownames(coldata), condition = coldata[,1])

anticodon_countdata = anticodon_countdata[, order(colnames(anticodon_countdata))]
isodecoder_countdata = isodecoder_countdata[, order(colnames(isodecoder_countdata))]
coldata = coldata[order(rownames(coldata)), , drop = FALSE]

# Convert to matrix
anticodon_countdata = as.matrix(anticodon_countdata)
isodecoder_countdata = as.matrix(isodecoder_countdata)

isoacceptorFile = list.files(path="./", pattern = "isoacceptorInfo.txt", full.names = T)
isoacceptorInfo = read.table(isoacceptorFile[1], header = T, row.names = 1)
isoacceptorInfo = isoacceptorInfo[ , 'size', drop=F]
isoacceptorInfo$rn = rownames(isoacceptorInfo)

# cluster_id == 1 means that mim-seq clusters are already representative of isodecoders, therefore get isodecoder info directly from custerInfo
# cluster_id == '' means clustering is disabled
if (cluster_id == 1){
  isodecoderFile = list.files(path="./", pattern = "clusterInfo.txt", full.names = T)
  isodecoderInfo = read.table(isodecoderFile[1], header = T)
  isodecoderInfo = isodecoderInfo[!duplicated(isodecoderInfo$parent),]
  isodecoderInfo = isodecoderInfo[, 'cluster_size', drop = F]
  colnames(isodecoderInfo) = 'size'
  isodecoderInfo$rn = rownames(isodecoderInfo)
} else {
  isodecoderFile = list.files(path="./", pattern = "isodecoderInfo.txt", full.names = T)
  isodecoderInfo = read.table(isodecoderFile[1], header = T, row.names = 1)
  isodecoderInfo = isodecoderInfo[ , 'size', drop = F]
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
    # Merge normalized count data with gene copy number info
    anticodon_counts$rn = rownames(anticodon_counts)
    isodecoder_counts$rn = rownames(isodecoder_counts)
    
    anticodon_counts = merge(anticodon_counts, isoacceptorInfo, by = 'rn', type = 'left')
    isodecoder_counts = merge(isodecoder_counts, isodecoderInfo, by = 'rn', type = 'left')
    
    names(anticodon_counts)[1] = "Gene"
    names(isodecoder_counts)[1] = "Gene"
    col_idx = grep("Gene", names(anticodon_counts))
    resdata_anticodon = anticodon_counts[, c(col_idx, (seq_len(ncol(anticodon_counts)))[-col_idx])]
    col_idx = grep("Gene", names(isodecoder_counts))
    resdata_isodecoder = isodecoder_counts[, c(col_idx, (seq_len(ncol(isodecoder_counts)))[-col_idx])]
    
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
      resdata_anticodon = anticodon_counts[, c(col_idx, (seq_len(ncol(anticodon_counts)))[-col_idx])]
      col_idx = grep("Gene", names(isodecoder_counts))
      resdata_isodecoder = isodecoder_counts[, c(col_idx, (seq_len(ncol(isodecoder_counts)))[-col_idx])]
      
      write.csv(resdata_anticodon, file=paste(subdir_anticodon, paste(toString(unique(coldata$condition)), 'normCounts.csv', sep="-"), sep="/"), row.names = FALSE)
      write.csv(resdata_isodecoder, file=paste(subdir_isodecoder, paste(toString(unique(coldata$condition)), 'normCounts.csv', sep="-"), sep="/"), row.names = FALSE)
      
      
    } else {
      dds_anticodon = DESeqDataSetFromMatrix(countData=anticodon_countdata, colData=coldata, design=~condition)
      dds_isodecoder = DESeqDataSetFromMatrix(countData=isodecoder_countdata, colData=coldata, design=~condition)
      
      # Run the DESeq pipeline
      dds_anticodon = DESeq(dds_anticodon)
      dds_isodecoder = DESeq(dds_isodecoder)
      
      # Normalized counts
      count_df_anticodon = as.data.frame(counts(dds_anticodon, normalized=TRUE))
      count_df_isodecoder = as.data.frame(counts(dds_isodecoder, normalized=TRUE))
      count_df_anticodon$rn = rownames(count_df_anticodon)
      count_df_isodecoder$rn = rownames(count_df_isodecoder)
      
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
      
      # Get combinations of coditions for various DE contrasts
      ordered_levels = levels(lastlevel(unique(dds_anticodon$condition), control_cond))
      combinations = combn(ordered_levels, 2, simplify = FALSE)
      
      # initialise lists for DE info from all contrasts vs control, for DE heatmaps
      DEcounts_list_isodecoder = list()
      DEcounts_list_anticodon = list()
      # For each contrast...
      for (i in seq_len(length(combinations))) {
        # Get differential expression results
        res_anticodon = results(dds_anticodon, contrast=c("condition",as.vector(combinations[[i]])))
        res_isodecoder = results(dds_isodecoder, contrast=c("condition",as.vector(combinations[[i]])))
        
        ## Order by adjusted p-value
        res_anticodon = res_anticodon[order(res_anticodon$padj), ]
        res_isodecoder = res_isodecoder[order(res_isodecoder$padj), ]
        res_anticodon$rn = rownames(res_anticodon)
        res_isodecoder$rn = rownames(res_isodecoder)
        
        ## Merge with normalized count data
        resdata_anticodon = join_all(list(as.data.frame(res_anticodon), count_df_anticodon, isoacceptorInfo), by="rn", type = 'left')
        resdata_isodecoder = join_all(list(as.data.frame(res_isodecoder), count_df_isodecoder, isodecoderInfo), by="rn", type = 'left')
        names(resdata_anticodon)[7] = "Gene"
        names(resdata_isodecoder)[7] = "Gene"
        col_idx = grep("Gene", names(resdata_anticodon))
        resdata_anticodon = resdata_anticodon[, c(col_idx, (seq_len(ncol(resdata_anticodon)))[-col_idx])]
        col_idx = grep("Gene", names(resdata_isodecoder))
        resdata_isodecoder = resdata_isodecoder[, c(col_idx, (seq_len(ncol(resdata_isodecoder)))[-col_idx])]
        
        # Count plots
        # add significance to baseMean matrices for current contrast
        baseMeanPerLvl_anticodon$sig = res_anticodon[rownames(baseMeanPerLvl_anticodon),'padj'] <= p_adj
        baseMeanPerLvl_anticodon$sig = !is.na(baseMeanPerLvl_anticodon$sig) & baseMeanPerLvl_anticodon$sig # deal with NAs turning them into FALSE
        baseMeanPerLvl_isodecoder$sig = res_isodecoder[rownames(baseMeanPerLvl_isodecoder),'padj'] <= p_adj
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
          scale_color_manual(paste('Differential expression\n(adj-p <=',p_adj, ')'), labels = c("None", "Down", "Up"), values = c("darkgrey", "#f28f3b", "#4daf4a")) + 
          scale_shape_manual(paste('Differential expression\n(adj-p <=',p_adj, ')'), labels = c("None", "Down", "Up"), values = c(19, 17, 17)) + 
          scale_size_manual(values = c(1,2), guide = FALSE) + theme_bw() +
          labs(x = paste('log10', combinations[[i]][2], 'counts', sep = ' '), y = paste('log10', combinations[[i]][1], 'counts', sep = ' ')) +
          annotate("label", 0, Inf, hjust = 0, vjust = 1, label = paste("italic(r) == ", anticodon_cor), parse = TRUE)
        
        isodecoder_dot = ggplot(isodecoder_data, aes_string(x, y)) +
          geom_point(aes(color = comb, shape = comb, size = sig)) +
          scale_x_log10() + scale_y_log10() + 
          #geom_smooth(method = 'lm', se = TRUE, alpha = 0.5, color = '#3182bd', fill = 'grey') +
          geom_abline(intercept = 0, slope = 1, linetype = 'dashed', color = '#3182bd', alpha = 0.8) + 
          scale_color_manual(paste('Differential expression\n(adj-p <=',p_adj, ')'), labels = c("None", "Down", "Up"), values = c("darkgrey", "#f28f3b", "#4daf4a")) + 
          scale_shape_manual(paste('Differential expression\n(adj-p <=',p_adj, ')'), labels = c("None", "Down", "Up"), values = c(19, 17, 17)) + 
          scale_size_manual(values = c(1,2), guide = FALSE) + theme_bw() +
          labs(x = paste('log10', combinations[[i]][2], 'counts', sep = ' '), y = paste('log10', combinations[[i]][1], 'counts', sep = ' ')) +
          annotate("label", 0, Inf, hjust = 0, vjust = 1, label = paste("italic(r) == ", isodecoder_cor), parse = TRUE)
        
        ggsave(paste(subdir_anticodon, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-countplot.pdf", sep="_"), sep="/"), anticodon_dot, height = 5, width = 8)
        ggsave(paste(subdir_isodecoder, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-countplot.pdf", sep="_"), sep="/"), isodecoder_dot, height = 5, width = 8)
        
        ## Write results
        write.csv(resdata_anticodon, file=paste(subdir_anticodon, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-results.csv", sep="_"), sep="/"))
        write.csv(resdata_isodecoder, file=paste(subdir_isodecoder, paste(paste(combinations[[i]], collapse="vs"), "diffexpr-results.csv", sep="_"), sep="/"))
        
        # Build DE tables for heatmaps only if the current contrast is to the control condition
        if (control_cond %in% combinations[[i]]) {
          temp_lfc_isodecoder = as.data.frame(res_isodecoder[,c("log2FoldChange","padj")])
          temp_lfc_isodecoder = subset(temp_lfc_isodecoder, padj <= p_adj)
          colnames(temp_lfc_isodecoder) = c(paste(paste(combinations[[i]], collapse = "vs"), "_l2FC", sep = ""), paste(paste(combinations[[i]], collapse = "vs"), "_padj", sep = ""))
          temp_lfc_isodecoder = tibble::rownames_to_column(temp_lfc_isodecoder, var = "isodecoder")
          temp_lfc_isodecoder = temp_lfc_isodecoder[!grepl("Escherichia_coli", temp_lfc_isodecoder$isodecoder),]
          DEcounts_list_isodecoder[[paste(combinations[[i]], collapse = "vs")]] = temp_lfc_isodecoder
          
          temp_lfc_anticodon = as.data.frame(res_anticodon[,c("log2FoldChange","padj")])
          temp_lfc_anticodon = subset(temp_lfc_anticodon, padj <= p_adj)
          colnames(temp_lfc_anticodon) = c(paste(paste(combinations[[i]], collapse = "vs"), "_l2FC", sep = ""), paste(paste(combinations[[i]], collapse = "vs"), "_padj", sep = ""))
          temp_lfc_anticodon = tibble::rownames_to_column(temp_lfc_anticodon, var = "Anticodon")
          temp_lfc_anticodon = temp_lfc_anticodon[!grepl("Escherichia_coli", temp_lfc_anticodon$Anticodon),]
          DEcounts_list_anticodon[[paste(combinations[[i]], collapse = "vs")]] = temp_lfc_anticodon
        }
      }
      
      ## Write normalized counts to /counts dir
      count_df_isodecoder_out = left_join(count_df_isodecoder, isodecoderInfo, by="rn")
      count_df_anticodon_out = left_join(count_df_anticodon, isoacceptorInfo, by="rn")
      names(count_df_anticodon_out)[names(count_df_anticodon_out) == "rn"] = "Anticodon"
      names(count_df_isodecoder_out)[names(count_df_isodecoder_out) == "rn"] = "isodecoder"
      col_idx = grep("Anticodon", names(count_df_anticodon_out))
      count_df_anticodon_out = count_df_anticodon_out[, c(col_idx, (seq_len(ncol(count_df_anticodon_out)))[-col_idx])]
      col_idx = grep("isodecoder", names(count_df_isodecoder_out))
      count_df_isodecoder_out = count_df_isodecoder_out[, c(col_idx, (seq_len(ncol(count_df_isodecoder_out)))[-col_idx])]
      write.table(count_df_isodecoder_out, file = "Isodecoder_counts_DESEqNormalized.csv", sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(count_df_anticodon_out, file = "Anticodon_counts_DESEqNormalized.csv", sep = "\t", quote = FALSE, row.names = FALSE)
      print("DESeq2 normalized counts saved to counts/*_normalized.txt")
      
      # Build combined filtered tables for DE heatmaps and scale counts
      comb_isodecoder = DEcounts_list_isodecoder %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="isodecoder"), .)
      basemean_isodecoder = as.data.frame(res_isodecoder) %>% dplyr::select(baseMean) %>% tibble::rownames_to_column(var = "isodecoder")
      normcounts_isodecoder = subset(count_df_isodecoder_out, select = -c(size))
      comb_isodecoder = list(comb_isodecoder, basemean_isodecoder, normcounts_isodecoder) %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="isodecoder"), .)
      comb_isodecoder$isodecoder = sub(".*?_.*?_tRNA-","",comb_isodecoder$isodecoder)
      comb_isodecoder = comb_isodecoder[!grepl("tRX", comb_isodecoder$isodecoder),]
      
      scaled_counts_isodecoder = as.matrix(t(scale(t(comb_isodecoder[,!grepl("isodecoder|l2FC|padj|baseMean", colnames(comb_isodecoder))]))), rownames.force = TRUE)
      scaled_counts_isodecoder[is.na(scaled_counts_isodecoder)] = 0
      rownames(scaled_counts_isodecoder) = comb_isodecoder$isodecoder
      
      comb_anticodon = DEcounts_list_anticodon %>% Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="Anticodon"), .)
      basemean_anticodon = as.data.frame(res_anticodon) %>% dplyr::select(baseMean) %>% tibble::rownames_to_column(var = "Anticodon")
      normcounts_anticodon = subset(count_df_anticodon_out, select = -c(size))
      comb_anticodon = list(comb_anticodon, basemean_anticodon, normcounts_anticodon) %>% Reduce(function(dtf1,dtf2) left_join(dtf1,dtf2,by="Anticodon"), .)
      comb_anticodon$Anticodon = sub(".*?_.*?_tRNA-","",comb_anticodon$Anticodon)
      comb_anticodon = comb_anticodon[!grepl("tRX", comb_anticodon$Anticodon),]
      
      scaled_counts_anticodon = as.matrix(t(scale(t(comb_anticodon[,!grepl("Anticodon|l2FC|padj|baseMean", colnames(comb_anticodon))]))), rownames.force = TRUE)
      scaled_counts_anticodon[is.na(scaled_counts_anticodon)] = 0
      rownames(scaled_counts_anticodon) = comb_anticodon$Anticodon
      
      # Heatmaps
      # note that all tables and plots only built if DE tables are not empty
      col_fun = colorRamp2(c(-3, 0, 3), c("#CC5803", "#f7f7f7", "#36682B"))
      
      # Annotation for baseMean counts
      if (nrow(comb_isodecoder) != 0) {
        baseMeanAnno_iso = rowAnnotation(base_mean = anno_lines(comb_isodecoder$baseMean, 
                                                                gp = gpar(lwd = 1.5, col = "#084081")), 
                                         annotation_name_rot = 90,
                                         width = unit("0.8", "cm"))
        
        hm_iso = Heatmap(scaled_counts_isodecoder,
                         col = col_fun,
                         row_title = paste('n = ', nrow(scaled_counts_isodecoder), sep = ""),
                         show_row_names = FALSE,
                         border = "gray20",
                         cluster_columns = TRUE,
                         heatmap_legend_param = list(
                            title = paste("Scaled expression\nlog2 fold-change (adj-p <= ", p_adj,")", sep = "")
                          )
        )
      }
      
      if (nrow(comb_anticodon) != 0) {
        baseMeanAnno_anti = rowAnnotation(base_mean = anno_lines(comb_anticodon$baseMean, 
                                                                 gp = gpar(lwd = 1.5, col = "#084081")), 
                                          annotation_name_rot = 90,
                                          width = unit("0.8", "cm"))
        
        hm_anti = Heatmap(scaled_counts_anticodon,
                          col = col_fun,
                          row_title = paste('n = ', nrow(scaled_counts_anticodon), sep = ""),
                          row_names_side = "left",
                          row_names_gp = gpar(fontsize = 6),
                          border = "gray20",
                          cluster_columns = TRUE,
                          heatmap_legend_param = list(
                            title = paste("Scaled expression\nlog2 fold-change (adj-p <= ", p_adj,")", sep = "")
                          )
        )
      }
      
      for (i in seq_len(length(combinations))){
        if (control_cond %in% combinations[[i]]) {
          lfc = paste(paste(combinations[[i]], collapse="vs"), "l2FC", sep="_")
          if (nrow(comb_isodecoder) != 0) {
            lfc_iso = Heatmap(comb_isodecoder[lfc],
                              col = col_fun,
                              width = unit(0.5, "cm"), na_col = "white",
                              border = "gray20",
                              show_heatmap_legend = FALSE)
            hm_iso = hm_iso + lfc_iso
          }
          if (nrow(comb_anticodon) != 0) {
            lfc_anti = Heatmap(comb_anticodon[lfc],
                               col = col_fun,
                               width = unit(0.5, "cm"), na_col = "white",
                               border = "gray20",
                               show_heatmap_legend = FALSE)
          hm_anti = hm_anti + lfc_anti
          }
        }
      }
      
      if (nrow(comb_isodecoder) != 0) {
        hm_iso = hm_iso + baseMeanAnno_iso
        
        pdf(paste(subdir_isodecoder,"DE_isodecodersScaled_hm.pdf",sep="/"))
        draw(hm_iso)
        dev.off()
      }
      
      if (nrow(comb_anticodon) != 0) {
        hm_anti = hm_anti + baseMeanAnno_anti
        
        
        pdf(paste(subdir_anticodon,"DE_anticodonScaled_hm.pdf",sep="/"))
        draw(hm_anti)
        dev.off()
      }
    }
  }
}