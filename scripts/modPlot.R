library(ggplot2)
library(reshape2)
library(gridExtra)
library(plyr)
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(RColorBrewer)

# Plotting of modification and stops heatmaps, and misincorporation scatter plots

args = commandArgs(trailingOnly = TRUE)
out = args[1]
cons_mods = args[2]

# read in mods and aggregate for total misinc. (sum of all types) and by condition (mean)
col_fun = colorRamp2(c(0, 0.5, 1), c("#f7fcf0", "#7bccc4", "#084081"))
mods = read.table(paste(out, "mods/mismatchTable.csv", sep = ''), header=T, sep = "\t", quote = '')
mods$proportion[is.na(mods$proportion)] = 0
mods$cluster = gsub(".*?_tRNA-", "", mods$cluster)
mods_agg = aggregate(mods$proportion, by = list(cluster=mods$cluster, pos=mods$pos, bam=mods$bam, struct=mods$struct, condition=mods$condition), FUN = sum)
mods_agg = aggregate(mods_agg$x, by = list(cluster=mods_agg$cluster, pos=mods_agg$pos, struct=mods_agg$struct, condition=mods_agg$condition), FUN = mean)

# read in stops table and process as above for mods
stops = read.table(paste(out, "mods/RTstopTable.csv", sep = ''), header = T, sep = "\t", quote = '')
stops$proportion[is.na(stops$proportion)] = 0
stops$cluster = gsub(".*?_tRNA-", "", stops$cluster)
#stops = stops[-5396, ] ## NB!! This cluster had info for pos 122 which doesn't exist - removed manually. Must be changed for other libraries ##
stops_agg = aggregate(stops$proportion, by = list(cluster=stops$cluster, pos=stops$pos, condition=stops$condition, struct=stops$struct), FUN = mean)

# read in context info created by ssAlign module
context_info = read.table(paste(out, "modContext.txt", sep = ''), header = TRUE)
context_info$cluster = sub(".*?_tRNA-","", context_info$cluster)

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
 	mods_hm = Heatmap(sub_mods_mat, row_title = "Misincorporations", column_title = as.character(i), column_title_side = "top", cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, top_annotation = col_anno, right_annotation = row_anno, heatmap_legend_param = list(title = "Misincorporation proportion", direction = "horizontal"))

 	# stops
	sub_stops_agg = subset(stops_agg, condition == i)
	sub_stops_wide = dcast(sub_stops_agg[,c('cluster','pos', 'x')], list(.(cluster), .(pos)), value.var = 'x', fun.aggregate = mean)
	sub_stops_wide[is.na(sub_stops_wide)] = 0
	rownames(sub_stops_wide) = sub_stops_wide$cluster
	sub_stops_wide = sub_stops_wide[, -1]
	sub_stops_mat = as.matrix(sub_stops_wide)
	col_anno = HeatmapAnnotation(Mean = anno_barplot(aggregate(sub_stops_agg$x, by = list(pos = sub_stops_agg$pos), FUN = mean)$x, height = unit(1.5, 'cm'),  gp = gpar(fill = '#C8553D')))
	count_stops = sub_stops_agg %>% group_by(cluster) %>% summarise(count = sum(x > 0.05))
	row_anno = rowAnnotation(Count = row_anno_barplot(count_stops$count, width = unit(1, 'cm'),  gp = gpar(fill = '#C8553D')))
	stops_hm = Heatmap(sub_stops_mat, row_title = "RT stops", cluster_columns = FALSE, cluster_rows = TRUE, col = col_fun, top_annotation = col_anno, right_annotation = row_anno, heatmap_legend_param = list(title = "RT stop proportion", direction = "horizontal"))

	# combined heatmap
	heatmap_list = mods_hm %v% stops_hm
	pdf(paste(out, 'mods/', paste(i, "comb_heatmap.pdf", sep = "_"), sep = ''), width = 16, height = 14)
	draw(heatmap_list, ht_gap = unit(10, "mm"))
	dev.off()

	# scatter plots
	sub_mods_pos = subset(sub_mods_agg, pos %in% cons_mods)
	sub_mods_pos = merge(sub_mods_pos, context_info, by = c('cluster', 'pos'))
	
	ggplot(sub_mods_agg, aes(x=as.character(pos), y = x, color = x)) + geom_jitter(width = 0.1, size = 3) +
	theme_bw() + facet_grid(identity~pos, scales = "free_x", labeller = label_both) + scale_color_gradientn(colours = cols) +
	geom_hline(yintercept = 0.1, linetype = "dashed", alpha = 0.4)


}





