#!/usr/bin/Rscript

# Plot 3' dinucleotide occurence of alingend reads
# Data table generated in mmQuant module

suppressMessages(library(ggplot2))

args = commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
	stop("At least one argument must be supplied (input file).n", call.=FALSE)

} else if (length(args)>0) {
	out = args[3]
	dinuc = read.table(args[1], header = TRUE, sep = "\t", na.strings="")
	dinuc = dinuc[!grepl("N", dinuc$dinuc),]
	dinuc = dinuc[-which(is.na(dinuc$dinuc)),]
	dinuc$color = "grey"
	dinuc$color[dinuc$dinuc == "CC"] = "red"
	dinuc$color[dinuc$dinuc == "CA"] = "green"

	dinuc_plot = ggplot(dinuc, aes(x=dinuc, y=proportion)) + geom_bar(stat="identity", aes(fill=color)) + 
	facet_wrap(~sample, ncol=3) + theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8), legend.position="none") + 
	scale_fill_manual(values = c("#1b9e77","#737373","#d95f02")) + xlab("3' dinucleotide") + ylab("Proportion")
	ggsave(paste(out, "dinuc_plot.pdf", sep = ''), dinuc_plot, height=5, width=10)

	cca = read.table(args[2], header = TRUE, sep = "\t")




}