#!/usr/bin/Rscript
#!/usr/bin/env Rscript

suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(gridExtra))
suppressMessages(library(plyr))
suppressMessages(library(tidyverse))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))

# Plotting of modification and stops heatmaps, and misincorporation scatter plots

args = commandArgs(trailingOnly = TRUE)
out = args[1]
mod_sites = args[2]
mod_sites = unlist(strsplit(mod_sites, "_"))
col_fun = colorRamp2(c(0, 0.5, 1), c("#f7fcf0", "#7bccc4", "#084081"))
cols = brewer.pal(9, "GnBu")[-(1:2)]
mito_trnas = args[4]
cons_pos = args[3]
cons_pos = unlist(strsplit(cons_pos, "_"))

# read in mods and aggregate for total misinc. (sum of all types) and by condition (mean)
mods = read.table(paste(out, "mods/mismatchTable.csv", sep = ''), header=T, sep = "\t", quote = '')
mods$proportion[is.na(mods$proportion)] = 0
mods$cluster = ifelse(grepl("mito", mods$cluster), sub(".*_mito_tRNA-","mito",mods$cluster), sub(".*_tRNA-","",mods$cluster))
mods$cluster = ifelse(mods$cluster == 'eColiLys-TTT-1-1', 'eColiLys', mods$cluster)
mods_agg = aggregate(mods$proportion, by = list(cluster=mods$cluster, pos=mods$pos, bam=mods$bam, struct=mods$struct, condition=mods$condition, canon_pos=mods$canon_pos), FUN = sum)
mods_agg = aggregate(mods_agg$x, by = list(cluster=mods_agg$cluster, pos=mods_agg$pos, struct=mods_agg$struct, condition=mods_agg$condition, canon_pos=mods_agg$canon_pos), FUN = mean)

# read in stops table and process as above for mods
stops = read.table(paste(out, "mods/RTstopTable.csv", sep = ''), header = T, sep = "\t", quote = '')
stops$proportion[is.na(stops$proportion)] = 0
stops$cluster = ifelse(grepl("mito", stops$cluster), sub(".*_mito_tRNA-","mito", stops$cluster), sub(".*_tRNA-","",stops$cluster))
stops$cluster = ifelse(stops$cluster == 'eColiLys-TTT-1-1', 'eColiLys', stops$cluster)

#stops = stops[-5396, ] ## NB!! This cluster had info for pos 122 which doesn't exist - removed manually. Must be changed for other libraries ##
stops_agg = aggregate(stops$proportion, by = list(cluster=stops$cluster, pos=stops$pos, condition=stops$condition, struct=stops$struct, canon_pos=stops$canon_pos), FUN = mean)

# read in context info created by ssAlign module
context_info = read.table(paste(out, "mods/modContext.txt", sep = ''), header = TRUE)
context_info$cluster = ifelse(grepl("mito", context_info$cluster), sub(".*_mito_tRNA-","mito", context_info$cluster), sub(".*_tRNA-","", context_info$cluster))
context_info$cluster = ifelse(context_info$cluster == 'eColiLys-TTT-1-1', 'eColiLys', context_info$cluster)

# for each condition make a misincorporation and stops heatmap as a combined figure using ComplexHeatmap
# ... make a scatter plot of misincorporation rates faceted by positions in cons_mods (selected known mod sites of interest) and by identity of nucleotide
for (i in unique(mods_agg$condition)) {
  
  # mods
  sub_mods_agg = subset(mods_agg, condition == i)
  sub_mods_wide = dcast(sub_mods_agg[,c('cluster','pos', 'x')], list(.(cluster), .(pos)), value.var = 'x', fun.aggregate = mean)
  sub_mods_wide[is.na(sub_mods_wide)] = 0
  rownames(sub_mods_wide) = sub_mods_wide$cluster
  sub_mods_wide = sub_mods_wide[, -1]
  sub_mods_mat = as.matrix(sub_mods_wide)
  col_anno = HeatmapAnnotation(Mean = anno_barplot(aggregate(sub_mods_agg$x, by = list(pos = sub_mods_agg$pos), FUN = mean)$x, height = unit(1.5, 'cm'),  gp = gpar(fill = '#C8553D')))
  count_mods = sub_mods_agg %>% group_by(cluster) %>% summarise(count = sum(x > 0.1))
  row_anno = rowAnnotation(Count = row_anno_barplot(count_mods$count, width = unit(1, 'cm'),  gp = gpar(fill = '#C8553D')))
  mods_hm = Heatmap(sub_mods_mat, column_labels = cons_pos, row_title = "Misincorporations", column_title = as.character(i), column_title_side = "top", cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, top_annotation = col_anno, right_annotation = row_anno, heatmap_legend_param = list(title = "Misincorporation proportion", direction = "horizontal"))
  
  # stops
  sub_stops_agg = subset(stops_agg, condition == i)
  sub_stops_wide = dcast(sub_stops_agg[,c('cluster','pos', 'x')], list(.(cluster), .(pos)), value.var = 'x', fun.aggregate = mean)
  sub_stops_wide[is.na(sub_stops_wide)] = 0
  rownames(sub_stops_wide) = sub_stops_wide$cluster
  sub_stops_wide = sub_stops_wide[, -1]
  sub_stops_mat = as.matrix(sub_stops_wide)
  col_anno = HeatmapAnnotation(Mean = anno_barplot(aggregate(sub_stops_agg$x, by = list(pos = sub_stops_agg$pos), FUN = mean)$x, height = unit(1.5, 'cm'),  gp = gpar(fill = '#C8553D')))
  count_stops = sub_stops_agg %>% group_by(cluster) %>% summarise(count = sum(x > 0.1))
  row_anno = rowAnnotation(Count = row_anno_barplot(count_stops$count, width = unit(1, 'cm'),  gp = gpar(fill = '#C8553D')))
  stops_hm = Heatmap(sub_stops_mat, column_labels = cons_pos, row_title = "RT stops", cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, top_annotation = col_anno, right_annotation = row_anno, heatmap_legend_param = list(title = "RT stop proportion", direction = "horizontal"))
  
  # combined heatmap
  heatmap_list = stops_hm %v% mods_hm
  pdf(paste(out, 'mods/', paste(i, "comb_heatmap.pdf", sep = "_"), sep = ''), width = 18, height = 16)
  draw(heatmap_list, ht_gap = unit(10, "mm"))
  dev.off()
  
  # scatter plots
  sub_mods_pos = subset(sub_mods_agg, canon_pos %in% mod_sites)
  sub_mods_pos[which(sub_mods_pos$x > 1), 'x'] = 1
  sub_mods_pos = merge(sub_mods_pos, context_info, by = c('cluster', 'pos'))
  names(sub_mods_pos)[names(sub_mods_pos) == 'x'] = 'Proportion'
  
  sub_mods_pos$canon_pos = factor(sub_mods_pos$canon_pos, levels = c('9','20', '20a', '26', '32','34','37','58'))
  
  mods_scatter = ggplot(subset(sub_mods_pos, !grepl("mito", sub_mods_pos$cluster)), aes(x=as.character(canon_pos), y = Proportion, color = Proportion)) + geom_jitter(width = 0.1, size = 3) +
    theme_bw() + facet_grid(identity~canon_pos, scales = "free_x", labeller = label_both) + scale_color_gradientn(colours = cols) +
    geom_hline(yintercept = 0.1, linetype = "dashed", alpha = 0.4) + 
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank()
    )
  
  ggsave(paste(out, "mods/", paste(i, 'misincProps.pdf', sep = '_'), sep = ''), mods_scatter, height=10, width=14)
  
  if (!is.na(mito_trnas)){
    mito_mods_scatter = ggplot(subset(sub_mods_pos, grepl("mito", sub_mods_pos$cluster)), aes(x=as.character(canon_pos), y = Proportion, color = Proportion)) + geom_jitter(width = 0.1, size = 3) +
      theme_bw() + facet_grid(identity~canon_pos, scales = "free_x", labeller = label_both) + scale_color_gradientn(colours = cols) +
      geom_hline(yintercept = 0.1, linetype = "dashed", alpha = 0.4) + 
      theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()
      )
    
    ggsave(paste(out, "mods/", paste('mito', i, 'misincProps.pdf', sep = '_'), sep = ''), mito_mods_scatter, height=10, width=14)
    
  }
  
  # Misinc signatures
  # create filter list of rows where total misinc. rate < 0.1 
  filter_0.1 = subset(sub_mods_agg, x < 0.1)
  # subset mods table for condition
  sub_mods_aggtype = subset(mods, condition == i)
  # use filter to filter rows from this table
  sub_mods_aggtype = anti_join(sub_mods_aggtype, filter_0.1, by=c("cluster","pos"))
  # add in context info
  sub_mods_aggtype = merge(sub_mods_aggtype, context_info, by = c("cluster","pos"))
  sub_mods_aggtype_cyt = subset(sub_mods_aggtype, !grepl("mito", sub_mods_aggtype$cluster))
  sub_mods_aggtype_cyt = aggregate(sub_mods_aggtype_cyt$proportion, by = list(identity = sub_mods_aggtype_cyt$identity, type = sub_mods_aggtype_cyt$type, upstream = sub_mods_aggtype_cyt$upstream, pos = sub_mods_aggtype_cyt$pos, canon_pos=sub_mods_aggtype_cyt$canon_pos), FUN = function(x) c(mean=mean(x), sd=sd(x)))
  sub_mods_aggtype_cyt = do.call("data.frame", sub_mods_aggtype_cyt)
  
  sub_mods_aggtype_cyt$canon_pos = factor(sub_mods_aggtype_cyt$canon_pos, levels = c('9', '20', '20a','26','32','34','37','58'))
  
  signature_plot = ggplot(sub_mods_aggtype_cyt, aes(x = type, y = x.mean, fill = type)) + 
    geom_bar(stat="identity", width = 0.8, position =position_dodge(width=0.9), alpha = 0.9) + 
    facet_grid(upstream~canon_pos+identity , scales = "free_x", labeller = label_both) + 
    geom_errorbar(aes(ymin = x.mean , ymax = x.mean + x.sd), width = 0.1, position = position_dodge(width=0.9)) + 
    theme_bw() +
    scale_fill_manual(values = c("#739FC2", "#7DB0A9", "#9F8FA9", "#C1B098"))
  
  ggsave(paste(out, "mods/", paste(i, 'misincSignatures.pdf', sep = '_'), sep = ''), signature_plot, height=10, width=14)
  
  if (!is.na(mito_trnas)){
    sub_mods_aggtype_mito = subset(sub_mods_aggtype, grepl("mito", sub_mods_aggtype$cluster))
    sub_mods_aggtype_mito = aggregate(sub_mods_aggtype_mito$proportion, by = list(identity = sub_mods_aggtype_mito$identity, type = sub_mods_aggtype_mito$type, upstream = sub_mods_aggtype_mito$upstream, pos = sub_mods_aggtype_mito$pos, canon_pos=sub_mods_aggtype_mito$canon_pos), FUN = function(x) c(mean=mean(x), sd=sd(x)))
    sub_mods_aggtype_mito = do.call("data.frame", sub_mods_aggtype_mito)
    sub_mods_aggtype_mito$canon_pos = factor(sub_mods_aggtype_mito$canon_pos, levels = c('9', '20', '20a', '26','32','34','37','58'))
    mito_signature_plot = ggplot(sub_mods_aggtype_mito, aes(x = type, y = x.mean, fill = type)) + 
      geom_bar(stat="identity", width = 0.8, position =position_dodge(width=0.9), alpha = 0.9) + 
      facet_grid(upstream~canon_pos+identity , scales = "free_x", labeller = label_both) + 
      geom_errorbar(aes(ymin = x.mean , ymax = x.mean + x.sd), width = 0.1, position = position_dodge(width=0.9)) + 
      theme_bw() +
      scale_fill_manual(values = c("#739FC2", "#7DB0A9", "#9F8FA9", "#C1B098"))
    
    ggsave(paste(out, "mods/", paste("mito", i, 'misincSignatures.pdf', sep = '_'), sep = ''), mito_signature_plot, height=10, width=14)
    
  }
  
}