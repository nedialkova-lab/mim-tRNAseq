#!/usr/bin/Rscript
#!/usr/bin/env Rscript

# Plot 3' dinucleotide occurence of aligned reads
# Data table generated in mmQuant module

suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(rlang))
suppressMessages(library(grid))
suppressMessages(library(gtable))

args = commandArgs(trailingOnly = TRUE)

# source facet_share.R
suppressMessages(source(args[5]))

if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
  
} else if (length(args)>0) {
  out = args[3]
  double_cca = args[4]
  if (double_cca == "True"){
    out_string = "double_ccaPlot.pdf"
  } else {
    out_string = "ccaPlot.pdf"
  }
  
  dinuc = read.table(args[1], header = TRUE, sep = "\t", na.strings="")
  dinuc = dinuc[!grepl("N", dinuc$dinuc),]
  dinuc$color = "grey"
  dinuc$color[dinuc$dinuc == "CC"] = "red"
  dinuc$color[dinuc$dinuc == "CA"] = "green"
  
  dinuc_plot = ggplot(dinuc, aes(x=dinuc, y=proportion)) +
               geom_bar(stat="identity", aes(fill=color)) +
               facet_wrap(~sample, ncol=3) +
               theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
                     legend.position="none") +
               scale_fill_manual(values = c("#66828D","#737373","#A05C45")) +
               xlab("3' dinucleotide") +
               ylab("Proportion")
  ggsave(paste(out, "dinuc_plot.pdf", sep = ''), dinuc_plot, height=5, width=10)
  
  cca_counts = read.table(args[2], header = TRUE, sep = "\t")
  cca_counts$gene = sub(".*_mito_tRNA-","mito",cca_counts$gene)
  cca_counts$gene = sub(".*_nmt_tRNA-","nmt",cca_counts$gene)
  cca_counts$gene = sub(".*_tRNA-","",cca_counts$gene)
  cca_counts$gene = ifelse(cca_counts$gene == 'eColiLys-TTT-1-1', 'eColiLys', cca_counts$gene)

  cca_prop = cca_counts %>% group_by(gene,sample) %>% 
    mutate(countT = sum(count)) %>% 
    group_by(end, .add = TRUE) %>% 
    mutate(per = round(100*count/countT,2))

  cca_summary = aggregate(cca_prop$per, by=list(gene = cca_prop$gene,
                                               end = cca_prop$end,
                                               condition = cca_prop$condition),
                          function(x) c(mean = mean(x), sd = sd(x)))
  cca_summary = do.call("data.frame", cca_summary)

  if (length(unique(cca_counts$condition)) > 1) {
    combinations = combn(unique(cca_counts$condition), 2, simplify = FALSE)

    # For each contrast...
    for (i in 1:length(combinations)) {
      cca_summary_sub = subset(cca_summary, condition %in% combinations[[i]])
      cca_summary_sub = cca_summary_sub %>%
                        mutate(x.mean = ifelse(condition == cca_summary_sub$condition[1],
                                                x.mean * -1,
                                                x.mean))
      
      cca_prop_sub = subset(cca_prop, condition %in% combinations[[i]])
      cca_prop_sub = cca_prop_sub %>%
                     mutate(per = ifelse(condition == cca_summary_sub$condition[1],
                                         per * -1,
                                         per))
      
      cca_prop_sub$bar_pos = NA
      cca_prop_sub$bar_pos[cca_prop_sub$end == "Absent"] = cca_prop_sub$per[cca_prop_sub$end == "Absent"]
      cca_prop_sub$bar_pos[cca_prop_sub$end == "C"] = cca_prop_sub$per[cca_prop_sub$end == "Absent"]+ 
        cca_prop_sub$per[cca_prop_sub$end == "C"]
      cca_prop_sub$bar_pos[cca_prop_sub$end == "CC"] = cca_prop_sub$per[cca_prop_sub$end == "Absent"]+ 
        cca_prop_sub$per[cca_prop_sub$end == "CC"] + cca_prop_sub$per[cca_prop_sub$end == "C"]
      cca_prop_sub$bar_pos[cca_prop_sub$end == "CA"] = cca_prop_sub$per[cca_prop_sub$end == "Absent"] + 
        cca_prop_sub$per[cca_prop_sub$end == "CA"] + cca_prop_sub$per[cca_prop_sub$end == "CC"] + 
        cca_prop_sub$per[cca_prop_sub$end == "C"]
      
      avg_cca = aggregate(cca_summary_sub$x.mean, by = list(condition = cca_summary_sub$condition, end = cca_summary_sub$end), mean)
      cca_summary_sub$end = factor(cca_summary_sub$end, levels = c('CA', 'CC', 'C', 'Absent'))
      
      cca_plot = ggplot(cca_summary_sub, aes(x = gene, y = x.mean, fill = end)) + 
        geom_bar(stat = 'identity', width = 0.8) +
        geom_hline(data = subset(avg_cca, end == 'CA'), aes(yintercept=x), color = "white", alpha = 0.9) + 
        geom_jitter(data = cca_prop_sub[cca_prop_sub$end == "CC",], aes(x = gene, y = bar_pos), size = 0.5, color = "#383D3B", alpha = 0.7) +
        geom_text(data = subset(avg_cca, end == 'CA'), aes(label = paste(abs(round(x,1)), '%'), x = Inf, y = x), size = 3.3, vjust = 1, color = '#3E606F', fontface='bold') +
        facet_share(~condition, dir = "h", scales = "free", reverse_num = TRUE) +
        coord_flip() + 
        scale_fill_manual(name = "", values = alpha(c(CA = "#F0F9ED", CC = "#427AA1", C = "#0D4282", Absent = "#133C55"), 0.8), labels = c("3'-CCA", "3'-CC", "3'-C", "Absent")) +
        scale_y_continuous(breaks = c(c(-100, -75, -50, -25, 0), c(0, 25, 50, 75, 100)))+
        scale_x_discrete(expand = c(0.03, 0)) +
        theme_minimal() + 
        theme(axis.title = element_blank(), 
              strip.text = element_text(face = "bold"), 
              axis.text.y = element_text(size = 9), 
              axis.text.x = element_text(face = 'bold'))
      
      ggsave(paste(out, paste(combinations[[i]][1], combinations[[i]][2], out_string, sep = '_'), sep = ''), cca_plot, height = 8, width = 9)
      
    }
    
  } else {
    cca_summary_sub = cca_summary
    cca_prop$bar_pos = NA
      cca_prop$bar_pos[cca_prop$end == "Absent"] = cca_prop$per[cca_prop$end == "Absent"]
      cca_prop$bar_pos[cca_prop$end == "C"] = cca_prop$per[cca_prop$end == "Absent"]+ 
        cca_prop$per[cca_prop$end == "C"]
      cca_prop$bar_pos[cca_prop$end == "CC"] = cca_prop$per[cca_prop$end == "Absent"]+ 
        cca_prop$per[cca_prop$end == "CC"] + cca_prop$per[cca_prop$end == "C"]
      cca_prop$bar_pos[cca_prop$end == "CA"] = cca_prop$per[cca_prop$end == "Absent"] + 
        cca_prop$per[cca_prop$end == "CA"] + cca_prop$per[cca_prop$end == "CC"] + 
        cca_prop$per[cca_prop$end == "C"]
         
    avg_cca = aggregate(cca_summary_sub$x.mean, by = list(condition = cca_summary_sub$condition, end = cca_summary_sub$end), mean)
    cca_summary_sub$end = factor(cca_summary_sub$end, levels = c('CA', 'CC', 'C', 'Absent'))
    
    cca_plot = ggplot(cca_summary_sub, aes(x = gene, y = x.mean, fill = end)) + 
        geom_bar(stat = 'identity', width = 0.8) +
        geom_hline(data = subset(avg_cca, end == 'CA'), aes(yintercept=x), color = "white", alpha = 0.9) + 
        geom_jitter(data = cca_prop[cca_prop$end == "CC",], aes(x = gene, y = bar_pos), size = 0.5, color = "#383D3B", alpha = 0.7) +
        geom_text(data = subset(avg_cca, end == 'CA'), aes(label = paste(abs(round(x,1)), '%'), x = Inf, y = x), size = 3.3, vjust = 1, color = '#3E606F', fontface='bold') +
        #facet_share(~condition, dir = "h", scales = "free", reverse_num = TRUE) +
        coord_flip() + 
        scale_fill_manual(name = "", values = alpha(c(CA = "#F0F9ED", CC = "#427AA1", C = "#0D4282", Absent = "#133C55"), 0.8), labels = c("3'-CCA", "3'-CC", "3'-C", "Absent")) +
        scale_y_continuous(breaks = c(c(-100, -75, -50, -25, 0), c(0, 25, 50, 75, 100)))+
        scale_x_discrete(expand = c(0.03, 0)) +
        theme_minimal() + 
        theme(axis.title = element_blank(), 
              strip.text = element_text(face = "bold"), 
              axis.text.y = element_text(size = 9), 
              axis.text.x = element_text(face = 'bold'))
   
    ggsave(paste(out, out_string, sep = ''), cca_plot, height = 8, width = 9)
    
  }
}
