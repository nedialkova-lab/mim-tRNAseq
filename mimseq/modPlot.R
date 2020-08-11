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
mito_trnas = args[5]
if (mito_trnas == ''){
  mito_trnas = NA
}
cons_pos = args[3]
cons_pos = unlist(strsplit(cons_pos, "_"))
cons_pos = cons_pos[!grepl("-",cons_pos)]
misinc_thresh = as.numeric(args[4])

# read in mods and aggregate for total misinc. (sum of all types) and by condition (mean)
mods = read.table(paste(out, "mods/mismatchTable.csv", sep = ''), header=T, sep = "\t", quote = '')
mods$proportion[is.na(mods$proportion)] = 0
mods$proportion[is.infinite(mods$proportion)] = 0
mods$isodecoder = sub(".*_mito_tRNA-","mito",mods$isodecoder)
mods$isodecoder = sub(".*_nmt_tRNA-","nmt",mods$isodecoder)
mods$isodecoder = sub(".*_tRNA-","",mods$isodecoder)
mods$isodecoder = ifelse(mods$isodecoder == 'eColiLys-TTT-1-1', 'eColiLys', mods$isodecoder)
mods = mods[!grepl("-", mods$canon_pos),]
mods_agg = aggregate(mods$proportion, by = list(isodecoder=mods$isodecoder, pos=mods$pos, bam=mods$bam, condition=mods$condition, canon_pos=mods$canon_pos), FUN = sum)
mods_agg = aggregate(mods_agg$x, by = list(isodecoder=mods_agg$isodecoder, pos=mods_agg$pos, condition=mods_agg$condition, canon_pos=mods_agg$canon_pos), FUN = mean)

# read in stops table and process as above for mods
stops = read.table(paste(out, "mods/RTstopTable.csv", sep = ''), header = T, sep = "\t", quote = '')
stops$proportion[is.na(stops$proportion)] = 0
stops$proportion[is.infinite(stops$proportion)] = 0
stops$isodecoder = sub(".*_mito_tRNA-","mito",stops$isodecoder)
stops$isodecoder = sub(".*_nmt_tRNA-","nmt",stops$isodecoder)
stops$isodecoder = sub(".*_tRNA-","",stops$isodecoder)
stops$isodecoder = ifelse(stops$isodecoder == 'eColiLys-TTT-1-1', 'eColiLys', stops$isodecoder)
stops = stops[!grepl("-", stops$canon_pos),]
stops_agg = aggregate(stops$proportion, by = list(isodecoder=stops$isodecoder, pos=stops$pos, condition=stops$condition, canon_pos=stops$canon_pos), FUN = mean)

# read in context info created by ssAlign module
context_info = read.table(paste(out, "mods/modContext.txt", sep = ''), header = TRUE)
colnames(context_info) = c("isodecoder", "pos", "identity", "upstream", "downstream")
context_info$isodecoder = sub(".*_mito_tRNA-","mito", context_info$isodecoder)
context_info$isodecoder = sub(".*_nmt_tRNA-","nmt", context_info$isodecoder)
context_info$isodecoder = sub(".*_tRNA-","", context_info$isodecoder)
context_info$isodecoder = ifelse(context_info$isodecoder == 'eColiLys-TTT-1-1', 'eColiLys', context_info$isodecoder)

# for each condition make a misincorporation and stops heatmap as a combined figure using ComplexHeatmap
# ... make a scatter plot of misincorporation rates faceted by positions in cons_mods (selected known mod sites of interest) and by identity of nucleotide
for (i in unique(mods_agg$condition)) {
  
  # cyto mods
  sub_mods_agg = mods_agg[mods_agg$condition == i & !grepl("mito", mods_agg$isodecoder) & !grepl("nmt", mods_agg$isodecoder),]
  sub_mods_wide = dcast(sub_mods_agg[,c('isodecoder','pos', 'x')], list(.(isodecoder), .(pos)), value.var = 'x', fun.aggregate = mean)
  sub_mods_wide[is.na(sub_mods_wide)] = 0
  rownames(sub_mods_wide) = sub_mods_wide$isodecoder
  sub_mods_wide = sub_mods_wide[, -1]
  sub_mods_mat = as.matrix(sub_mods_wide)
  col_anno = HeatmapAnnotation(Mean = anno_barplot(aggregate(sub_mods_agg$x, by = list(pos = sub_mods_agg$pos), FUN = mean)$x, height = unit(1.5, 'cm'),  gp = gpar(fill = '#C8553D')))
  count_mods = sub_mods_agg %>% group_by(isodecoder) %>% summarise(count = sum(x > misinc_thresh))
  row_anno = rowAnnotation(Count = row_anno_barplot(count_mods$count, width = unit(1, 'cm'),  gp = gpar(fill = '#C8553D')))
  cyto_mods_hm = Heatmap(sub_mods_mat, column_labels = cons_pos, row_title = "Misincorporations", column_title = as.character(i), column_title_side = "top", cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, top_annotation = col_anno, right_annotation = row_anno, heatmap_legend_param = list(title = "Misincorporation proportion", direction = "horizontal"))
  
  # cyto stops
  sub_stops_agg = stops_agg[stops_agg$condition == i & !grepl("mito", stops_agg$isodecoder) & !grepl("nmt", stops_agg$isodecoder),]
  sub_stops_wide = dcast(sub_stops_agg[,c('isodecoder','pos', 'x')], list(.(isodecoder), .(pos)), value.var = 'x', fun.aggregate = mean)
  sub_stops_wide[is.na(sub_stops_wide)] = 0
  rownames(sub_stops_wide) = sub_stops_wide$isodecoder
  sub_stops_wide = sub_stops_wide[, -1]
  sub_stops_mat = as.matrix(sub_stops_wide)
  col_anno = HeatmapAnnotation(Mean = anno_barplot(aggregate(sub_stops_agg$x, by = list(pos = sub_stops_agg$pos), FUN = mean)$x, height = unit(1.5, 'cm'),  gp = gpar(fill = '#C8553D')))
  count_stops = sub_stops_agg %>% group_by(isodecoder) %>% summarise(count = sum(x > misinc_thresh))
  row_anno = rowAnnotation(Count = row_anno_barplot(count_stops$count, width = unit(1, 'cm'),  gp = gpar(fill = '#C8553D')))
  cyto_stops_hm = Heatmap(sub_stops_mat, column_labels = cons_pos, row_title = "RT stops", cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, top_annotation = col_anno, right_annotation = row_anno, heatmap_legend_param = list(title = "RT stop proportion", direction = "horizontal"))
  
  # combined cyto heatmap
  heatmap_list = cyto_stops_hm %v% cyto_mods_hm
  pdf(paste(out, 'mods/', paste(i, "comb_heatmap.pdf", sep = "_"), sep = ''), width = 18, height = 16)
  draw(heatmap_list, ht_gap = unit(10, "mm"), column_title = "Cytoplasmic clusters")
  
  if (!is.na(mito_trnas)){
    # mito mods
    sub_mods_agg = mods_agg[mods_agg$condition == i & (grepl("mito", mods_agg$isodecoder) | grepl("nmt", mods_agg$isodecoder)),]
    
    if (nrow(sub_mods_agg) != 0) { 
      sub_mods_wide = dcast(sub_mods_agg[,c('isodecoder','pos', 'x')], list(.(isodecoder), .(pos)), value.var = 'x', fun.aggregate = mean)
      sub_mods_wide[is.na(sub_mods_wide)] = 0
      rownames(sub_mods_wide) = sub_mods_wide$isodecoder
      sub_mods_wide = sub_mods_wide[, -1]
      sub_mods_mat = as.matrix(sub_mods_wide)
      col_anno = HeatmapAnnotation(Mean = anno_barplot(aggregate(sub_mods_agg$x, by = list(pos = sub_mods_agg$pos), FUN = mean)$x, height = unit(1.5, 'cm'),  gp = gpar(fill = '#C8553D')))
      count_mods = sub_mods_agg %>% group_by(isodecoder) %>% summarise(count = sum(x > misinc_thresh))
      row_anno = rowAnnotation(Count = row_anno_barplot(count_mods$count, width = unit(1, 'cm'),  gp = gpar(fill = '#C8553D')))
      mito_mods_hm = Heatmap(sub_mods_mat, column_labels = cons_pos, row_title = "Misincorporations", column_title = as.character(i), column_title_side = "top", cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, top_annotation = col_anno, right_annotation = row_anno, heatmap_legend_param = list(title = "Misincorporation proportion", direction = "horizontal"))
    
      # mito stops
      sub_stops_agg = stops_agg[stops_agg$condition == i & (grepl("mito", stops_agg$isodecoder) | grepl("nmt", stops_agg$isodecoder)),]
      sub_stops_wide = dcast(sub_stops_agg[,c('isodecoder','pos', 'x')], list(.(isodecoder), .(pos)), value.var = 'x', fun.aggregate = mean)
      sub_stops_wide[is.na(sub_stops_wide)] = 0
      rownames(sub_stops_wide) = sub_stops_wide$isodecoder
      sub_stops_wide = sub_stops_wide[, -1]
      sub_stops_mat = as.matrix(sub_stops_wide)
      col_anno = HeatmapAnnotation(Mean = anno_barplot(aggregate(sub_stops_agg$x, by = list(pos = sub_stops_agg$pos), FUN = mean)$x, height = unit(1.5, 'cm'),  gp = gpar(fill = '#C8553D')))
      count_stops = sub_stops_agg %>% group_by(isodecoder) %>% summarise(count = sum(x > misinc_thresh))
      row_anno = rowAnnotation(Count = row_anno_barplot(count_stops$count, width = unit(1, 'cm'),  gp = gpar(fill = '#C8553D')))
      mito_stops_hm = Heatmap(sub_stops_mat, column_labels = cons_pos, row_title = "RT stops", cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, top_annotation = col_anno, right_annotation = row_anno, heatmap_legend_param = list(title = "RT stop proportion", direction = "horizontal"))
    
      # combined mito heatmap
      heatmap_list = mito_stops_hm %v% mito_mods_hm
      draw(heatmap_list, ht_gap = unit(10, "mm"), column_title = "Mitochondrial (and nuclear-encoded mito) clusters")
      dev.off()
    }
  }
  
  else {
    dev.off()
  }
  
  # scatter plots
  temp_mods = merge(mods, context_info, by = c('isodecoder', 'pos'))
  temp_mods = temp_mods %>% group_by(isodecoder, pos, bam, identity) %>% mutate(new_prop = proportion/sum(proportion))
  filter_proportions = temp_mods %>% group_by(isodecoder, pos, bam, identity) %>% filter((any(max(new_prop) > 0.95) & any(canon_pos != 34)) | (any(max(new_prop) > 0.95) & any(identity != 'A') & any(canon_pos == 34) & any(type != 'G')))
  sub_mods_agg = mods_agg[mods_agg$condition == i,]
  sub_mods_pos = sub_mods_agg[sub_mods_agg$canon_pos %in% mod_sites,]
  sub_mods_pos[which(sub_mods_pos$x > 1), 'x'] = 1
  sub_mods_pos = merge(sub_mods_pos, context_info, by = c('isodecoder', 'pos'))
  sub_mods_pos = anti_join(sub_mods_pos, filter_proportions, by = c('isodecoder', 'pos', 'identity'))
  names(sub_mods_pos)[names(sub_mods_pos) == 'x'] = 'Proportion'
  
  sub_mods_pos$canon_pos = factor(sub_mods_pos$canon_pos, levels = c('9','20', '20a', '26', '32','34','37','58'))
  
  mods_scatter = ggplot(sub_mods_pos[!grepl("mito", sub_mods_pos$isodecoder) & !grepl("nmt", sub_mods_pos$isodecoder), ], aes(x=as.character(canon_pos), y = Proportion, color = Proportion)) + geom_jitter(width = 0.1, size = 3) +
    theme_bw() + facet_grid(identity~canon_pos, scales = "free_x", labeller = label_both) + scale_color_gradientn(breaks = c(0.0, 0.25, 0.50, 0.75, 1.0), colours = cols) +
    geom_hline(yintercept = misinc_thresh, linetype = "dashed", alpha = 0.4) + 
    theme(
      axis.text.x=element_blank(),
      axis.title.x=element_blank(),
      axis.ticks.x=element_blank()
    )
  
  ggsave(paste(out, "mods/", paste(i, 'misincProps.pdf', sep = '_'), sep = ''), mods_scatter, height=10, width=14)
  
  if (!is.na(mito_trnas)){
    sub_mods_pos_mito = sub_mods_pos[grepl("mito", sub_mods_pos$isodecoder) | grepl("nmt", sub_mods_pos$isodecoder), ]

    if (nrow(sub_mods_pos_mito) != 0) {
      mito_mods_scatter = ggplot(sub_mods_pos_mito, aes(x=as.character(canon_pos), y = Proportion, color = Proportion)) + geom_jitter(width = 0.1, size = 3) +
        theme_bw() + facet_grid(identity~canon_pos, scales = "free_x", labeller = label_both) + scale_color_gradientn(breaks = c(0.0, 0.25, 0.50, 0.75, 1.0), colours = cols) +
        geom_hline(yintercept = misinc_thresh, linetype = "dashed", alpha = 0.4) + 
        theme(
          axis.text.x=element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x=element_blank()
        )
      ggsave(paste(out, "mods/", paste('mito', i, 'misincProps.pdf', sep = '_'), sep = ''), mito_mods_scatter, height=10, width=14)
    }
  }
  
  # Misinc signatures
  # create filter list of rows where total misinc. rate < misinc_thresh 
  filter_misincthresh = sub_mods_agg[sub_mods_agg$x < misinc_thresh, ]
  # subset mods table for condition
  sub_mods_aggtype = mods[mods$condition == i, ]
  # use filter to filter rows from this table
  sub_mods_aggtype = anti_join(sub_mods_aggtype, filter_misincthresh, by=c("isodecoder","pos"))
  # add in context info
  sub_mods_aggtype = merge(sub_mods_aggtype, context_info, by = c("isodecoder","pos"))
  sub_mods_aggtype$bam = sub(out, "", sub_mods_aggtype$bam)
  sub_mods_aggtype_cyt = sub_mods_aggtype[!grepl("mito", sub_mods_aggtype$isodecoder) & !grepl("nmt", sub_mods_aggtype$isodecoder), ]
  # renormalise by sum of misinc at each site for each isodecoder in each bam file - this makes sum all misinc types = 1, additionally filter all clusters at each pos where misinc of highest nucle > 0.95
  sub_mods_aggtype_cyt = sub_mods_aggtype_cyt %>% group_by(isodecoder, pos, bam, identity) %>% mutate(new_prop = proportion/sum(proportion)) %>% filter(any(max(new_prop) < 0.95) | (any(max(new_prop) >= 0.95 & any(identity == 'A') & any(canon_pos == 34) & any(type == 'G'))))
                                                                                                                                                                                  
  #sub_mods_aggtype_cyt_up = aggregate(sub_mods_aggtype_cyt$proportion, by = list(identity = sub_mods_aggtype_cyt$identity, type = sub_mods_aggtype_cyt$type, upstream = sub_mods_aggtype_cyt$upstream, pos = sub_mods_aggtype_cyt$pos, canon_pos=sub_mods_aggtype_cyt$canon_pos), FUN = function(x) c(mean=mean(x), sd=sd(x)))
  #sub_mods_aggtype_cyt_up = do.call("data.frame", sub_mods_aggtype_cyt_up)
  #sub_mods_aggtype_cyt_down = aggregate(sub_mods_aggtype_cyt$proportion, by = list(identity = sub_mods_aggtype_cyt$identity, type = sub_mods_aggtype_cyt$type, downstream = sub_mods_aggtype_cyt$downstream, pos = sub_mods_aggtype_cyt$pos, canon_pos=sub_mods_aggtype_cyt$canon_pos), FUN = function(x) c(mean=mean(x), sd=sd(x)))
  #sub_mods_aggtype_cyt_down = do.call("data.frame", sub_mods_aggtype_cyt_down)
  
  sub_mods_aggtype_cyt$canon_pos = factor(sub_mods_aggtype_cyt$canon_pos, levels = c('9', '20', '20a','26','32','34','37','58'))
  color_num = length(unique(sub_mods_aggtype_cyt$bam)) + 1
  dot_colors = brewer.pal(color_num, "Greys")[2:(color_num)]
  names(dot_colors) = unique(sub_mods_aggtype_cyt$bam)
  
  signature_plot_upstream = ggplot(sub_mods_aggtype_cyt, aes(x = type, y = new_prop, fill = type)) + 
    geom_jitter(aes(color = bam), alpha = 0.6, size = 0.7) +
    geom_boxplot(aes(color = type), lwd = 0.9, alpha = 0.4, outlier.shape = NA) +
    facet_grid(upstream~canon_pos+identity , scales = "free_x", labeller = label_both) + 
    theme_bw() +
    labs(y = "Proportion") + 
    theme(
      axis.title.x  = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values = c("A"="#739FC2", "C"="#7DB0A9", "G"="#9F8FA9", "T"="#C1B098", dot_colors)) +
    scale_fill_manual(values = c("#739FC2", "#7DB0A9", "#9F8FA9", "#C1B098")) +
    guides(color = "none", fill = guide_legend(override.aes = list(color = c("#739FC2", "#7DB0A9", "#9F8FA9", "#C1B098"))))
  
  ggsave(paste(out, "mods/", paste(i, 'misincSignatures_upstreamContext.pdf', sep = '_'), sep = ''), signature_plot_upstream, height=10, width=14, useDingbats=FALSE)
  
  signature_plot_downstream = ggplot(sub_mods_aggtype_cyt, aes(x = type, y = new_prop, fill = type)) + 
    geom_jitter(aes(color = bam), alpha = 0.6, size = 0.7) +
    geom_boxplot(aes(color = type), lwd = 0.9, alpha = 0.4, outlier.shape = NA) +
    facet_grid(downstream~canon_pos+identity , scales = "free_x", labeller = label_both) + 
    theme_bw() +
    labs(y = "Proportion") + 
    theme(
      axis.title.x  = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values = c("A"="#739FC2", "C"="#7DB0A9", "G"="#9F8FA9", "T"="#C1B098", dot_colors)) +
    scale_fill_manual(values = c("#739FC2", "#7DB0A9", "#9F8FA9", "#C1B098")) +
    guides(color = "none", fill = guide_legend(override.aes = list(color = c("#739FC2", "#7DB0A9", "#9F8FA9", "#C1B098"))))
  
  ggsave(paste(out, "mods/", paste(i, 'misincSignatures_downstreamContext.pdf', sep = '_'), sep = ''), signature_plot_downstream, height=10, width=14, useDingbats=FALSE)
  
  if (!is.na(mito_trnas)){
    sub_mods_aggtype_mito = sub_mods_aggtype[grepl("mito", sub_mods_aggtype$isodecoder) | grepl("nmt", sub_mods_aggtype$isodecoder), ]

    if(nrow(sub_mods_aggtype_mito) != 0) {
      # renormalise by sum of misinc at each site for each isodecoder in each bam file - this makes sum all misinc types = 1
      sub_mods_aggtype_mito = sub_mods_aggtype_mito %>% group_by(isodecoder, pos, bam) %>% mutate(new_prop = proportion/sum(proportion)) %>% filter(any(max(new_prop) < 0.95))
      #sub_mods_aggtype_mito = aggregate(sub_mods_aggtype_mito$proportion, by = list(identity = sub_mods_aggtype_mito$identity, type = sub_mods_aggtype_mito$type, upstream = sub_mods_aggtype_mito$upstream, downstream = sub_mods_aggtype_mito$downstream, pos = sub_mods_aggtype_mito$pos, canon_pos=sub_mods_aggtype_mito$canon_pos), FUN = function(x) c(mean=mean(x), sd=sd(x)))
      #sub_mods_aggtype_mito = do.call("data.frame", sub_mods_aggtype_mito)
      sub_mods_aggtype_mito$canon_pos = factor(sub_mods_aggtype_mito$canon_pos, levels = c('9', '20', '20a', '26','32','34','37','58'))

      # Reallocate bam dot colours to account for inconsistencies between mito and cyto bam numbers (usually due to very low count data)
      color_num = length(unique(sub_mods_aggtype_mito$bam)) + 1
      dot_colors = brewer.pal(color_num, "Greys")[2:(color_num)]
      names(dot_colors) = unique(sub_mods_aggtype_mito$bam)

      mito_signature_plot_upstream = ggplot(sub_mods_aggtype_mito, aes(x = type, y = new_prop, fill = type)) + 
        geom_jitter(aes(color = bam), alpha = 0.6, size = 0.7) +
        geom_boxplot(aes(color = type), lwd = 0.9, alpha = 0.4, outlier.shape = NA) +
        facet_grid(upstream~canon_pos+identity , scales = "free_x", labeller = label_both) + 
        theme_bw() +
        labs(y = "Proportion") + 
        theme(
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        scale_color_manual(values = c("A"="#739FC2", "C"="#7DB0A9", "G"="#9F8FA9", "T"="#C1B098", dot_colors)) +
        scale_fill_manual(values = c("#739FC2", "#7DB0A9", "#9F8FA9", "#C1B098")) +
        guides(color = "none", fill = guide_legend(override.aes = list(color = c("#739FC2", "#7DB0A9", "#9F8FA9", "#C1B098"))))
      
      ggsave(paste(out, "mods/", paste("mito", i, 'misincSignatures_upstreamContext.pdf', sep = '_'), sep = ''), mito_signature_plot_upstream, height=10, width=14)
      
      mito_signature_plot_downstream = ggplot(sub_mods_aggtype_mito, aes(x = type, y = new_prop, fill = type)) + 
        geom_jitter(aes(color = bam), alpha = 0.6, size = 0.7) +
        geom_boxplot(aes(color = type), lwd = 0.9, alpha = 0.4, outlier.shape = NA) +
        facet_grid(downstream~canon_pos+identity , scales = "free_x", labeller = label_both) + 
        theme_bw() +
        labs(y = "Proportion") + 
        theme(
          axis.title.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()
        ) +
        scale_color_manual(values = c("A"="#739FC2", "C"="#7DB0A9", "G"="#9F8FA9", "T"="#C1B098", dot_colors)) +
        scale_fill_manual(values = c("#739FC2", "#7DB0A9", "#9F8FA9", "#C1B098")) +
        guides(color = "none", fill = guide_legend(override.aes = list(color = c("#739FC2", "#7DB0A9", "#9F8FA9", "#C1B098"))))
      
      ggsave(paste(out, "mods/", paste("mito", i, 'misincSignatures_downstreamContext.pdf', sep = '_'), sep = ''), mito_signature_plot_downstream, height=10, width=14)
    }
  }
  
}