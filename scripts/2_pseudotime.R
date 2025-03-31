#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#ARGS
#dir = "/Users/marialuisaratto/scripts/aggr/"
dir = as.character(args[1])
print(paste("Folder:", dir))
#file = cds
file = args[2]
print(paste("File:", file))

pseudotime = args[3]
print(paste("Pseudotime file:", pseudotime))

all = as.logical(args[4])
print(paste("All:", all))

print("Starting...")


suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(dplyr))
suppressMessages(library(monocle3))

if (file.exists(paste0(dir,"/", file)) == F) {
  stop(paste("Error: cds file does not exist at path:", paste0(dir,"/", file)))
}

if (file.exists(paste0(dir, "/", pseudotime)) == F) {
  stop(paste("Error: Pseudotime ile does not exist at path:", paste0(dir,"/", pseudotime)))
}


load(paste0(dir,"/", file))
pt = read.csv(paste0(dir, "/", pseudotime), row.names = 1)
analysis = names(pt)
names(pt) = c("pseudotime")

#CUMOLATIVE FREQ
pt$nomi = rownames(pt)
#AAACCAACAGGATTAA_30_AACCGGAG_SMAD2_TCTATT_1_SMAD2.2_78_KD_0_KD_KD_singleSample
pt = separate(pt, nomi, into = c("cellID", "UMI", "barcode", "gene", "UCI", "sample", "shRNA", "cellsxclone", "control", "replicate", "sampleC", "geneC", "SampleAnn"), sep = "_", remove = F, convert = T)
pt = mutate(pt, clone = paste(shRNA, UCI, sep = "_"))
pt$UMI = as.integer(pt$UMI)

if (length(unique(pt$gene) < 21)) {
  ggplot(pt, aes(pseudotime, colour = gene)) + stat_ecdf(geom = "step")+ 
    #scale_color_oi(order = c(9, 7))+
    theme_classic(base_size = 14)+
    ylab("Cumulative frequency")+
    xlab("monocle pseudotime")+
    ggtitle(paste("Cumulative frequency on", analysis))
  ggsave(filename = paste0(dir, "/cumulative_frequency_", analysis, "_gene.pdf"),
         width = 10, height = 14)
}else{
  print("Too many genes to plot cumulative_frequency_pseudotime_gene.pdf")
}


if (length(unique(pt$shRNA) < 21)) {
  ggplot(pt, aes(pseudotime, colour = shRNA)) + stat_ecdf(geom = "step")+ 
    #scale_color_oi(order = c(9, 7))+
    theme_classic(base_size = 14)+
    ylab("Cumulative frequency")+
    xlab("monocle pseudotime")+
    ggtitle(paste("Cumulative frequency on", analysis))+
    facet_grid(~factor(gene))
  ggsave(filename = paste0(dir, "/cumulative_frequency_", analysis, "_shRNA.pdf"),
         width = 14, height = 5)
}else{
  print("Too many shRNA to plot cumulative_frequency_pseudotime_shRNA.pdf")
}


if (length(unique(pt$clone) < 21)) {
  ggplot(pt, aes(pseudotime, colour = clone)) + stat_ecdf(geom = "step")+ 
    #scale_color_oi(order = c(9, 7))+
    theme_classic(base_size = 14)+
    ylab("Cumulative frequency")+
    xlab("monocle pseudotime")+
    ggtitle(paste("Cumulative frequency on", analysis))+
    facet_grid(~factor(shRNA))
  ggsave(filename = paste0(dir, "/cumulative_frequency_", analysis, "_clone.pdf"),
         width = 40, height = 10)
}else{
  print("Too many clones to plot cumulative_frequency_pseudotime_clone.pdf")
}
print("Finished plotting cumulative frequencies")

control_genes = unique(filter(pt, geneC == "control")$gene)
control_shRNA = unique(filter(pt, geneC == "control")$shRNA)
control_clones = unique(filter(pt, geneC == "control")$clone)

#KS TEST based on control samples 
if ("control" %in% pt$sampleC) {
  ks_stat = data.frame()
  
  for (c in setdiff(c(unique(pt$gene), unique(pt$shRNA), unique(pt$clone)), c(control_genes, control_shRNA, control_clones))) {
    tryCatch({
      print(c)
      data = filter(pt, gene == c | shRNA == c | clone == c)
      tmp_TET = filter(data, sampleC == "KD")
      tmp_CTR = filter(data, sampleC == "control")
      
      ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
      ks_pval_g = ks_g[["p.value"]]
      ks_stat_g = -ks_g$statistic[["D^+"]]
      print(ks_pval_g)
      
      ks_stat = rbind(ks_stat, data.frame(group = c, pval = ks_pval_g, statistic = ks_stat_g, analysis = "TETvsCTR", type = "greater"))
      
      ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
      ks_pval_l = ks_l[["p.value"]]
      ks_stat_l = ks_l$statistic[["D^-"]]
      print(ks_pval_l)
      ks_stat = rbind(ks_stat, data.frame(group = c, pval = ks_pval_l, statistic = ks_stat_l, analysis = "TETvsCTR", type = "less"))
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  
  #colnames(ks_stat) <- c("group", "pval", "analysis", "direction")
  
  ks_stat$padj <- p.adjust(ks_stat$pval, method = "BH")
  ks_stat$sig <- ks_stat$padj <= 0.05
  
  ks_stat_fil = ks_stat %>%
    group_by(group) %>%
    filter(padj == min(padj)) %>%
    ungroup() %>% 
    distinct()
  
  write.csv(ks_stat, file= paste0(dir, "/ks_statistics_vscontrol_samples_", analysis, ".csv"), row.names=FALSE)
  
  ks_stat_fil$gene <- sub("[._].*", "", ks_stat_fil$group)
  ggplot(data=ks_stat_fil, aes(x=statistic, y=log10(padj)*-1)) +
    geom_point(aes(colour = gene), size=2) +
    #scale_colour_manual(values=c("#009E73","#D55E00"))+
    ggrepel::geom_text_repel(
      data=ks_stat_fil,
      aes(label = group),
      size = 4,
      max.overlaps=10,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )+
    geom_hline(yintercept = log10(0.05)*-1,linetype="dotted", color="red")+
    geom_vline(xintercept = 0,linetype="dotted")+
    ggtitle("ks test statistics KD shRNA")+
    theme_minimal() +  # You can change the theme as needed
    labs(
      x = "KS statistic summary",
      y = "-Log10 adj p-value"
    ) 
  
  ggsave(filename = paste0(dir, "/volcano_ks_stats_vs_noTET_", analysis, ".pdf"),
         width = 14, height = 8)
}

#KS TEST based on control genes
#GENES
ks_stat = data.frame()
tmp2 = filter(pt, sampleC == "KD")

for (c in unique(pt$gene)) {
  tryCatch({
    print(c)
    tmp_TET = filter(tmp2, gene == c)
    
    if (all == T) {
      for (CG in control_genes) {
        tryCatch({
          print(CG)
          tmp_CTR = filter(tmp2, gene == CG)
          
          ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
          ks_pval_g = ks_g[["p.value"]]
          ks_stat_g = -ks_g$statistic[["D^+"]]
          print(ks_pval_g)
          ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_pval_g, statistic = ks_stat_g, analysis = paste0("vs", CG), type = "greater"))
          
          ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
          ks_pval_l = ks_l[["p.value"]]
          ks_stat_l = ks_l$statistic[["D^-"]]
          print(ks_pval_l)
          ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_pval_l, statistic = ks_stat_l, analysis = paste0("vs", CG), type = "less"))
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
    #control togheter
    print("ALL CONTROLS")
    tmp_CTR = filter(tmp2, gene %in% control_genes)
    
    ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
    ks_pval_g = ks_g[["p.value"]]
    ks_stat_g = -ks_g$statistic[["D^+"]]
    print(ks_pval_g)
    ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_pval_g, statistic = ks_stat_g, analysis = "vsControlGenes", type = "greater"))
    
    ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
    ks_pval_l = ks_l[["p.value"]]
    ks_stat_l = ks_l$statistic[["D^-"]]
    print(ks_pval_l)
    ks_stat = rbind(ks_stat, data.frame(gene = c, pval = ks_pval_l, statistic = ks_stat_l, analysis = "vsControlGenes", type = "less"))
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

ks_stat$padj <- p.adjust(ks_stat$pval, method = "BH")
ks_stat$sig <- ks_stat$padj <= 0.05

ks_stat_fil = ks_stat %>%
  group_by(gene, analysis) %>%
  filter(padj == min(padj)) %>%
  ungroup() %>% 
  distinct()

write.csv(ks_stat, file= paste0(dir, "/genes_ks_statistics_vscontrol_genes_all",all, "_", analysis, ".csv"), row.names=FALSE)

#ks_stat_fil$gene <- sub("[._].*", "", ks_stat_fil$group)
ggplot(data=ks_stat_fil, aes(x=statistic, y=log10(padj)*-1)) +
  geom_point(aes(colour = analysis), size=2) +
  #scale_colour_manual(values=c("#009E73","#D55E00"))+
  ggrepel::geom_text_repel(
    data=ks_stat_fil,
    aes(label = gene),
    size = 4,
    max.overlaps=10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  geom_hline(yintercept = log10(0.05)*-1,linetype="dotted", color="red")+
  geom_vline(xintercept = 0,linetype="dotted")+
  ggtitle("ks test statistics KD genes")+
  theme_minimal() +  # You can change the theme as needed
  labs(
    x = "KS statistic summary",
    y = "-Log10 adj p-value"
  ) 

ggsave(filename = paste0(dir, "/volcano_ks_stats_vs_control_genes_all",all, "_", analysis, ".pdf"),
       width = 14, height = 8)



# #SHRNAs
ks_stat = data.frame()
tmp2 = filter(pt, sampleC == "KD")

for (c in unique(pt$shRNA)) {
  tryCatch({
    #control genes
    print(c)
    tmp_TET = filter(tmp2, shRNA == c)
    
    if (all == T) {
      for (CG in c(control_genes, control_shRNA)) {
        tryCatch({
          tmp_CTR = filter(tmp2, gene == CG | shRNA == CG)
          ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
          ks_pval_g = ks_g[["p.value"]]
          ks_stat_g = -ks_g$statistic[["D^+"]]
          print(ks_pval_g)
          ks_stat = rbind(ks_stat, data.frame(shRNA = c, pval = ks_pval_g, statistic = ks_stat_g, analysis = paste0("vs", CG), type = "greater"))
          
          ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
          ks_pval_l = ks_l[["p.value"]]
          ks_stat_l = ks_l$statistic[["D^-"]]
          print(ks_pval_l)
          ks_stat = rbind(ks_stat, data.frame(shRNA = c, pval = ks_pval_l, statistic = ks_stat_l, analysis = paste0("vs", CG), type = "less"))
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
    }
    
    #controls together
    tmp_CTR = filter(tmp2, gene %in% control_genes)
    
    ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
    ks_pval_g = ks_g[["p.value"]]
    ks_stat_g = -ks_g$statistic[["D^+"]]
    print(ks_pval_g)
    ks_stat = rbind(ks_stat, data.frame(shRNA = c, pval = ks_pval_g, statistic = ks_stat_g, analysis = "vsControlGenes", type = "greater"))
    
    ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
    ks_pval_l = ks_l[["p.value"]]
    ks_stat_l = ks_l$statistic[["D^-"]]
    print(ks_pval_l)
    ks_stat = rbind(ks_stat, data.frame(shRNA = c, pval = ks_pval_l, statistic = ks_stat_l, analysis = "vsControlGenes", type = "less"))
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

ks_stat$padj <- p.adjust(ks_stat$pval, method = "BH")
ks_stat$sig <- ks_stat$padj <= 0.05

ks_stat_fil = ks_stat %>%
  group_by(shRNA, analysis) %>%
  filter(padj == min(padj)) %>%
  ungroup() %>% 
  distinct()

write.csv(ks_stat, file= paste0(dir, "/shRNA_ks_statistics_vscontrol_genes_all",all, "_", analysis, ".csv"), row.names=FALSE)

ks_stat_fil$gene <- sub("[._].*", "", ks_stat_fil$shRNA)
ggplot(data=ks_stat_fil, aes(x=statistic, y=log10(padj)*-1)) +
  geom_point(aes(colour = gene, shape = analysis), size=2) +
  #scale_colour_manual(values=c("#009E73","#D55E00"))+
  ggrepel::geom_text_repel(
    data=ks_stat_fil,
    aes(label = shRNA),
    size = 4,
    max.overlaps=10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  geom_hline(yintercept = log10(0.05)*-1,linetype="dotted", color="red")+
  geom_vline(xintercept = 0,linetype="dotted")+
  ggtitle("ks test statistics KD shRNA")+
  theme_minimal() +  # You can change the theme as needed
  labs(
    x = "KS statistic summary",
    y = "-Log10 adj p-value"
  ) 

ggsave(filename = paste0(dir, "/volcano_shRNA_ks_statistics_vscontrol_genes_all",all, "_", analysis, ".pdf"),
       width = 14, height = 8)

#CLONES
ks_stat = data.frame()
tmp2 = filter(pt, sampleC == "KD")

for (c in unique(pt$clone)) {
  tryCatch({
    #control genes
    print(c)
    tmp_TET = filter(pt, clone == c)
    
    # if (all == T) {
    #   for (CG in c(control_genes, control_shRNA, control_clones)) {
    #     tryCatch({
    #       tmp_CTR = filter(pt, gene == CG | shRNA == CG | clone == CG)
    #       ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
    #       ks_pval_g = ks_g[["p.value"]]
    #       ks_stat_g = -ks_g$statistic[["D^+"]]
    #       print(ks_pval_g)
    #       ks_stat = rbind(ks_stat, data.frame(clone = c, pval = ks_pval_g, statistic = ks_stat_g, analysis = paste0("vs", CG), type = "greater"))
    #       
    #       ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
    #       ks_pval_l = ks_l[["p.value"]]
    #       ks_stat_l = ks_l$statistic[["D^-"]]
    #       print(ks_pval_l)
    #       ks_stat = rbind(ks_stat, data.frame(clone = c, pval = ks_pval_l, statistic = ks_stat_l, analysis = paste0("vs", CG), type = "less"))
    #     }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    #   }
    # }
    # 
    #controls together
    tmp_CTR = filter(tmp2, gene %in% control_genes)
    
    ks_g = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="greater")
    ks_pval_g = ks_g[["p.value"]]
    ks_stat_g = -ks_g$statistic[["D^+"]]
    print(ks_pval_g)
    ks_stat = rbind(ks_stat, data.frame(clone = c, pval = ks_pval_g, statistic = ks_stat_g, analysis = "vsControlGenes", type = "greater"))
    
    ks_l = ks.test(tmp_TET$pseudotime, tmp_CTR$pseudotime, alternative="less")
    ks_pval_l = ks_l[["p.value"]]
    ks_stat_l = ks_l$statistic[["D^-"]]
    print(ks_pval_l)
    ks_stat = rbind(ks_stat, data.frame(clone = c, pval = ks_pval_l, statistic = ks_stat_l, analysis = "vsControlGenes", type = "less"))
    
    
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

ks_stat$padj <- p.adjust(ks_stat$pval, method = "BH")
ks_stat$sig <- ks_stat$padj <= 0.05

write.csv(ks_stat, file= paste0(dir,"/clones_ks_statistics_vscontrol_genes_", analysis, ".csv"), row.names=FALSE)

ks_stat_fil = ks_stat %>%
  group_by(clone, analysis) %>%
  filter(padj == min(padj)) %>%
  ungroup() %>% 
  distinct()

ks_stat_fil$gene <- sub("[._].*", "", ks_stat_fil$clone)
# get shRNA number
ks_stat_fil = separate(ks_stat_fil, col = clone, into = c("shRNA", "UCI"), sep = "_", remove = F)
ks_stat_fil = separate(ks_stat_fil, col = shRNA, into = c("G", "Number"), sep = "\\.", remove = F, fill = "right", extra = "merge")
ks_stat_fil$Number <- sub(".*\\.(.*)", "\\1", ks_stat_fil$Number)

p = ggplot(data=ks_stat_fil, aes(x=statistic, y=log10(padj)*-1)) +
  geom_point(aes(colour = Number), size=2) +
  #scale_colour_manual(values=c("#009E73","#D55E00"))+
  ggrepel::geom_text_repel(
    data=ks_stat_fil,
    aes(label = clone),
    size = 2,
    max.overlaps=10,
    box.padding = unit(0.35, "lines"),
    point.padding = unit(0.3, "lines")
  )+
  geom_hline(yintercept = log10(0.05)*-1,linetype="dotted", color="red")+
  geom_vline(xintercept = 0,linetype="dotted")+
  ggtitle("ks test statistics KD shRNA")+
  theme_minimal() +  # You can change the theme as needed
  labs(
    x = "KS statistic summary",
    y = "-Log10 adj p-value"
  ) +
  facet_wrap(vars(gene))

ggsave(p, filename = paste0(dir, "/volcano_clones_ks_statistics_vscontrol_genes_", analysis, ".pdf"),
       width = 20, height = 14)

#gene and pseudotime correlation
cds <- cds[
  ,
  colData(cds) %>%
    subset(
      sampleC == "KD"
    ) %>%
    row.names
]
pr_graph_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=1)
pr_graph_test_res = pr_graph_test_res[order(pr_graph_test_res$morans_I, decreasing = T), ]
write.csv(pr_graph_test_res, file= paste0(dir,"regression_pseudotime.csv"))
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

p = plot_cells(cds, genes=pr_deg_ids[1:30],
               show_trajectory_graph=FALSE,
               label_cell_groups=FALSE,
               label_leaves=FALSE)
ggsave(p, filename = paste0(dir, "/correlated_pseudotime_gene_exp", analysis, ".pdf"),
       width = 14, height = 14)

