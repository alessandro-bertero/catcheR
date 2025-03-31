#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
#ARGS
#dir = "/Users/marialuisaratto/scripts/AB03_test/"
dir = args[1]
#file = "processed_cds.RData"
file = args[2]
#meta = "cell_metadata.csv"
meta = args[3]
#timepoint = "PSC"
timepoint=args[4]
#control_gene = "SCR"
control_gene = args[5]


library(tidyverse)
library(ggplot2)
library(gtools)
library(dplyr)
library(RColorBrewer)

#install devtools if you don't have it already for easy installation
#install.packages("devtools")
suppressMessages(library(devtools))
suppressMessages(library(monocle3))
#install_github("chrisamiller/fishplot")
#library(fishplot)
#install.packages("see")
#library(see)

metadata = read.delim(paste0(dir, meta), header = T, sep = ",", row.names = 1)
metadata$cellsxclone = as.numeric(metadata$cellsxclone)
load(paste0(dir, file))

##in case of TET vs noTET
if ("control" %in% metadata$sampleC) {
  print("Control samples are present")
  
  #separate df based on condition
  #for (rep in unique(metadata$replicate)) {
  #print(rep)
  #data = filter(metadata, replicate == rep)
  data = distinct(metadata[c("clone", "cellsxclone", "sample", "replicate", "sampleC")])
  TET = filter(data, sampleC == "KD")
  tot = sum(TET$cellsxclone)
  TET = mutate(TET, percentage_TET = cellsxclone / tot * 100)
  TET = TET[c("clone", "percentage_TET")]
  CTR = filter(data, sampleC == "control")
  tot = sum(CTR$cellsxclone)
  CTR = mutate(CTR, percentage_CTR = cellsxclone / tot * 100)
  CTR = CTR[c("clone", "percentage_CTR")]
  #}
  
  data = merge(data, TET, by = "clone", all = T)
  data = merge(data, CTR, by = "clone", all = T)
  data[is.na(data)] <- 0
  
  frac.table = matrix(
    c(data$percentage_CTR, data$percentage_TET),
    ncol=2)
  data = separate(data, clone, c("gene", "barcode", "UCI"), sep = "_", remove = F)
  # parents =  c(rep(0, length(data$clone)), 1:length(data$clone))
  # fish = createFishObject(frac.table,parents = parents,timepoints=c(0,1), clone.annots =c(data$gene), clone.annots.cex = 0.4)
  # fish = setCol(fish, col = c(brewer.pal(length(unique(data$clone))-24, "Set3"),
  #                             brewer.pal(8, "Accent"),
  #                             brewer.pal(8, "Dark2"),
  #                             brewer.pal(8, "Set2")))
  # 
  # #calculate the layout of the drawing
  # fish = layoutClones(fish, separate.independent.clones = T)
  # 
  # #draw the plot, using the splining method (recommended)
  # #and providing both timepoints to label and a plot title
  # pdf('fish.pdf', width=8, height=5)
  # fishPlot(fish = fish, shape = "bezier", vlines = c(0,1), vlab = c("noTET", "TET"), 
  #          bg.type = "solid", bg.col = "white")
  # dev <- dev.off()
  
  #TOTAL NUMBER OF CELLS
  data_long = data %>% 
    select(-sampleC) %>%
    pivot_longer(cols = c(percentage_TET, percentage_CTR), names_to = "sampleC")
  
  p = ggplot(data = data_long, aes(x=sampleC, y=value, fill = gene)) +
    geom_col() +
    facet_wrap(vars(clone)) +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(filename = paste0(dir, "cellsxclone_TETvsnoTET.pdf"), plot = p, width = 20, height = 10)
  
  p = ggplot(data = data_long, aes(x=gene, y=value, fill = sampleC)) +
    #facet_grid(col =vars(gene)) +
    geom_col(position = "dodge") +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(filename = paste0(dir, "cellsxgene_TETvsnoTET.pdf"), plot = p, width = 20, height = 10)
  
  #FC 
  data = distinct(data)
  data = mutate(data, log2Percentage_CTR = log2(percentage_CTR))
  data = mutate(data, FC = percentage_TET / percentage_CTR)
  data = mutate(data, log2FC = log2(FC))
  #wide[is.na(wide)] <- 0
  p = ggplot(data = data, aes(x=clone, y= log2FC, fill = gene)) +
    geom_col() + 
    theme(axis.text.x = element_text(angle = 90))
  ggsave(filename = paste0(dir, "log2FC_TETvsnoTET.pdf"), plot = p, width = 15, height = 10)
  
  wide_big = filter(data, cellsxclone > 19)
  p = ggplot(data = wide_big, aes(x=clone, y= log2FC, fill = gene)) +
    geom_col() + 
    theme(axis.text.x = element_text(angle = 90))
  ggsave(filename = paste0(dir, "log2FC_bigclones_TETvsnoTET.pdf"), plot = p, width = 15, height = 10)
  
  p = ggplot(data = wide_big, aes(x=log2FC, y= log2Percentage_CTR, colour = gene)) + #colour = p.value
    geom_point() +
    geom_vline(xintercept = c(-1, 1), linetype = "dotted")
  ggsave(filename = paste0(dir, "volcano_plot_enrich_bigclones_TETvsnoTET.pdf"), plot = p)
  
  #p-value
  #WHICH IS THE COMPARISON FOR 2D???
  #wide = mutate(wide, !!paste0("p-value_", v) := fisher.test(!!sym(paste0("cellsxclone_", timepoint)), !!sym(paste0("cellsxclone_", v)), simulate.p.value = T)["p.value"][[1]])
  
  #}
  
  #STATS FOR CLUSTERS
  colData(cds)$clusters = clusters(cds)
  df = as.data.frame(colData(cds))
  
  ## GENE LEVEL
  condition = "gene"
  
  data_CTR = filter(df, sampleC == "control")
  # Calculate the number of cells per cluster and condition
  cluster_condition_counts_CTR <- data_CTR %>%
    group_by(clusters, gene) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate the percentage of KD and Control cells per cluster
  cluster_condition_percentage_CTR <- cluster_condition_counts_CTR %>%
    group_by(clusters) %>%
    mutate(total = sum(count)) %>%
    mutate(percentage = (count / total) * 100) %>%
    ungroup()
  
  data_TET = filter(df, sampleC == "KD")
  
  # Calculate the number of cells per cluster and condition
  cluster_condition_counts_TET <- data_TET %>%
    group_by(clusters, gene) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate the percentage of KD and Control cells per cluster
  cluster_condition_percentage_TET <- cluster_condition_counts_TET %>%
    group_by(clusters) %>%
    mutate(total = sum(count)) %>%
    mutate(percentage = (count / total) * 100) %>%
    ungroup()
  
  ## Combine the two data frames
  cluster_condition_percentage_TET$percentage_TET = cluster_condition_percentage_TET$percentage
  cluster_condition_percentage_TET$percentage = NULL
  cluster_condition_percentage_CTR$percentage_CTR = cluster_condition_percentage_CTR$percentage
  cluster_condition_percentage_CTR$percentage = NULL
  combined_df <- merge(
    cluster_condition_percentage_TET[,c("clusters", "gene", "percentage_TET")], 
    cluster_condition_percentage_CTR[,c("clusters", "gene", "percentage_CTR")])
  
  cluster_condition_percentage_TET$treat = "TET"
  cluster_condition_percentage_CTR$treat = "CTR"
  names(cluster_condition_percentage_CTR) = c("clusters","gene","count","total","percentage","treat")
  names(cluster_condition_percentage_TET) = c("clusters","gene","count","total","percentage","treat")
  combined_df = rbind(cluster_condition_percentage_CTR, cluster_condition_percentage_TET)
  # Plot
  p = ggplot(combined_df, aes(x = factor(clusters), y = percentage, fill=gene)) + geom_bar(stat = "identity", colour = "white", linewidth = 0.3)+
    #+theme(legend.position='none')
    #scale_fill_okabeito(reverse=T)+
    facet_grid(~treat)
  ggsave(p, filename = "genes_in_clusters_TETvsnoTET.pdf")
  
  
  ## shRNA level
  df = as.data.frame(colData(cds))
  condition = "shRNA"
  
  data_CTR = filter(df, sampleC == "control")
  # Calculate the number of cells per cluster and condition
  cluster_condition_counts_CTR <- data_CTR %>%
    group_by(clusters, shRNA) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate the percentage of KD and Control cells per cluster
  cluster_condition_percentage_CTR <- cluster_condition_counts_CTR %>%
    group_by(clusters) %>%
    mutate(total = sum(count)) %>%
    mutate(percentage = (count / total) * 100) %>%
    ungroup()
  
  data_TET = filter(df, sampleC == "KD")
  
  # Calculate the number of cells per cluster and condition
  cluster_condition_counts_TET <- data_TET %>%
    group_by(clusters, shRNA) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate the percentage of KD and Control cells per cluster
  cluster_condition_percentage_TET <- cluster_condition_counts_TET %>%
    group_by(clusters) %>%
    mutate(total = sum(count)) %>%
    mutate(percentage = (count / total) * 100) %>%
    ungroup()
  
  ## Combine the two data frames
  cluster_condition_percentage_TET$percentage_TET = cluster_condition_percentage_TET$percentage
  cluster_condition_percentage_TET$percentage = NULL
  cluster_condition_percentage_CTR$percentage_CTR = cluster_condition_percentage_CTR$percentage
  cluster_condition_percentage_CTR$percentage = NULL
  combined_df <- merge(
    cluster_condition_percentage_TET[,c("clusters", "shRNA", "percentage_TET")], 
    cluster_condition_percentage_CTR[,c("clusters", "shRNA", "percentage_CTR")])
  
  cluster_condition_percentage_TET$treat = "TET"
  cluster_condition_percentage_CTR$treat = "CTR"
  names(cluster_condition_percentage_CTR) = c("clusters","shRNA","count","total","percentage","treat")
  names(cluster_condition_percentage_TET) = c("clusters","shRNA","count","total","percentage","treat")
  combined_df = rbind(cluster_condition_percentage_CTR, cluster_condition_percentage_TET)
  # Plot
  p = ggplot(combined_df, aes(x = factor(clusters), y = percentage, fill=shRNA)) + geom_bar(stat = "identity", colour = "white", linewidth = 0.3)+
    #+theme(legend.position='none')
    #scale_fill_okabeito(reverse=T)+
    facet_grid(~treat)
  ggsave(p, filename = "shRNA_in_clusters_TETvsnoTET.pdf")
  
  
  ## CLUSTERS IN PERTURBATIONS
  ## GENE LEVEL
  condition = "gene"
  
  data_CTR = filter(df, sampleC == "control")
  # Calculate the number of cells per cluster and condition
  cluster_condition_counts_CTR <- data_CTR %>%
    group_by(clusters, gene) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate the percentage of KD and Control cells per cluster
  cluster_condition_percentage_CTR <- cluster_condition_counts_CTR %>%
    group_by(gene) %>%
    mutate(total = sum(count)) %>%
    mutate(percentage = (count / total) * 100) %>%
    ungroup()
  
  data_TET = filter(df, sampleC == "KD")
  
  # Calculate the number of cells per cluster and condition
  cluster_condition_counts_TET <- data_TET %>%
    group_by(clusters, gene) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate the percentage of KD and Control cells per cluster
  cluster_condition_percentage_TET <- cluster_condition_counts_TET %>%
    group_by(gene) %>%
    mutate(total = sum(count)) %>%
    mutate(percentage = (count / total) * 100) %>%
    ungroup()
  
  ## Combine the two data frames
  cluster_condition_percentage_TET$percentage_TET = cluster_condition_percentage_TET$percentage
  cluster_condition_percentage_TET$percentage = NULL
  cluster_condition_percentage_CTR$percentage_CTR = cluster_condition_percentage_CTR$percentage
  cluster_condition_percentage_CTR$percentage = NULL
  combined_df <- merge(
    cluster_condition_percentage_TET[,c("clusters", "gene", "percentage_TET")], 
    cluster_condition_percentage_CTR[,c("clusters", "gene", "percentage_CTR")])
  
  cluster_condition_percentage_TET$treat = "TET"
  cluster_condition_percentage_CTR$treat = "CTR"
  names(cluster_condition_percentage_CTR) = c("clusters","gene","count","total","percentage","treat")
  names(cluster_condition_percentage_TET) = c("clusters","gene","count","total","percentage","treat")
  combined_df = rbind(cluster_condition_percentage_CTR, cluster_condition_percentage_TET)
  # Plot
  p = ggplot(combined_df, aes(x = factor(treat), y = percentage, fill=clusters)) + 
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    #+theme(legend.position='none')
    #scale_fill_okabeito(reverse=T)+
    facet_grid(~gene) 
  ggsave(p, filename = paste0(dir, "clusters_in_genes_TETvsnoTET.pdf"), width = 20, height = 10)
  
  ##FISHER TEST
  complete_combinded <- combined_df %>%
    complete(clusters, gene, treat, fill = list(count = 0))
  complete_combinded_w <- complete_combinded %>%
    group_by(gene, treat, clusters) %>%
    summarise(count = sum(count), .groups = "drop") %>% 
    pivot_wider(names_from = clusters, values_from = count, values_fill = 0)
  res = data.frame(row.names = unique(complete_combinded$gene))
  
  for (g in unique(complete_combinded$gene)) {
    print(g)
    for (cluster in unique(complete_combinded$clusters)) {
      print(cluster)
      tmp = filter(complete_combinded, gene == g)
      ds = filter(tmp, clusters == cluster)[,c("treat", "count", "total")]
      ds = mutate(ds, out = total - count)
      ds$total = NULL
      rownames(ds) = ds$treat
      ds$treat = NULL
      # Check if both CTR and TET have at least 2 levels
      if (sum(is.na(ds)) > 0) {
        warning(paste("Skipping iteration", cluster, ": not enough rows"))
        res[g, cluster] = NA  
        next 
      }
      fisher = fisher.test(ds)
      p_value <- fisher[["p.value"]]
      print(p_value)
      res[g, cluster] = p_value
    }
  }
  
  res$gene = rownames(res)
  colnames(res) <- as.character(colnames(res))
  res <- res %>%
    pivot_longer(cols = c('1':as.character(ncol(res)-1)), names_to = "clusters", values_to = "p_value")
  
  # Add a q-value column using p.adjust
  res$padj <- p.adjust(res$p_value, method = "fdr")
  res$sig <- res$padj <= 0.05
  
  results = merge(res, complete_combinded)
  
  ## FC between CTR and TET
  # Aggregate the count for each gene and cluster, by treatment (TET vs CTR)
  fc_genes_clusters <- aggregate(count ~ gene + clusters + treat, data = results, sum)
  # Reshape the data so that you have counts for TET and CTR in separate columns
  fc_genes_clusters_wide <- reshape(fc_genes_clusters, idvar = c("gene", "clusters"), timevar = "treat", direction = "wide")
  # Calculate the Fold Change as TET/CTR
  fc_genes_clusters_wide$FC <- fc_genes_clusters_wide$count.TET / fc_genes_clusters_wide$count.CTR
  
  results = merge(results, fc_genes_clusters_wide)
  plotting = distinct(results[, c("gene", "clusters", "FC", "padj", "sig")])
  
  ## VOLCANO PLOT
  p = ggplot(data=plotting, aes(x=log2(as.numeric(FC)), y=log10(padj)*-1)) +
    geom_point(aes(colour = clusters), size=1.5) +
    ggrepel::geom_text_repel(
      data = subset(plotting, sig==T),
      aes(label = gene),
      size = 3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )+
    #geom_vline(xintercept = log2(1.5),linetype="dotted")+
    #geom_vline(xintercept = -log2(1.5),linetype="dotted")+
    geom_hline(yintercept = -log10(0.05),linetype="dotted")+
    ggtitle("fisher test statistics per gene in each cluster")+
    theme_minimal() +  
    labs(
      x = "Log2FC compared to CTR",
      y = "-Log10 adj p-value",
      color = "cluster"  # Legend title
    ) 
  
  ggsave(p, filename = paste0(dir, "gene_volcano_fisher_stats_TETvsnoTET.pdf"),
         width = 5, height = 4)
  
  
  
  ## shRNA LEVEL
  condition = "shRNA"
  
  data_CTR = filter(df, sampleC == "control")
  # Calculate the number of cells per cluster and condition
  cluster_condition_counts_CTR <- data_CTR %>%
    group_by(clusters, shRNA) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate the percentage of KD and Control cells per cluster
  cluster_condition_percentage_CTR <- cluster_condition_counts_CTR %>%
    group_by(shRNA) %>%
    mutate(total = sum(count)) %>%
    mutate(percentage = (count / total) * 100) %>%
    ungroup()
  
  data_TET = filter(df, sampleC == "KD")
  
  # Calculate the number of cells per cluster and condition
  cluster_condition_counts_TET <- data_TET %>%
    group_by(clusters, shRNA) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Calculate the percentage of KD and Control cells per cluster
  cluster_condition_percentage_TET <- cluster_condition_counts_TET %>%
    group_by(shRNA) %>%
    mutate(total = sum(count)) %>%
    mutate(percentage = (count / total) * 100) %>%
    ungroup()
  
  ## Combine the two data frames
  cluster_condition_percentage_TET$percentage_TET = cluster_condition_percentage_TET$percentage
  cluster_condition_percentage_TET$percentage = NULL
  cluster_condition_percentage_CTR$percentage_CTR = cluster_condition_percentage_CTR$percentage
  cluster_condition_percentage_CTR$percentage = NULL
  combined_df <- merge(
    cluster_condition_percentage_TET[,c("clusters", "shRNA", "percentage_TET")], 
    cluster_condition_percentage_CTR[,c("clusters", "shRNA", "percentage_CTR")])
  
  cluster_condition_percentage_TET$treat = "TET"
  cluster_condition_percentage_CTR$treat = "CTR"
  names(cluster_condition_percentage_CTR) = c("clusters","shRNA","count","total","percentage","treat")
  names(cluster_condition_percentage_TET) = c("clusters","shRNA","count","total","percentage","treat")
  combined_df = rbind(cluster_condition_percentage_CTR, cluster_condition_percentage_TET)
  # Plot
  p = ggplot(combined_df, aes(x = factor(treat), y = percentage, fill=clusters)) + 
    geom_bar(stat = "identity", colour = "white", linewidth = 0.3)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    #+theme(legend.position='none')
    #scale_fill_okabeito(reverse=T)+
    facet_grid(~shRNA) 
  ggsave(p, filename = paste0(dir, "clusters_in_shRNAs_TETvsnoTET.pdf"), width = 20, height = 7)
  
  
  ##FISHER TEST
  complete_combinded <- combined_df %>%
    complete(clusters, shRNA, treat, fill = list(count = 0))
  complete_combinded_w <- complete_combinded %>%
    group_by(shRNA, treat, clusters) %>%
    summarise(count = sum(count), .groups = "drop") %>% 
    pivot_wider(names_from = clusters, values_from = count, values_fill = 0)
  res = data.frame(row.names = unique(complete_combinded$shRNA))
  
  for (g in unique(complete_combinded$shRNA)) {
    print(g)
    for (cluster in unique(complete_combinded$clusters)) {
      print(cluster)
      tmp = filter(complete_combinded, shRNA == g)
      ds = filter(tmp, clusters == cluster)[,c("treat", "count", "total")]
      ds = mutate(ds, out = total - count)
      ds$total = NULL
      rownames(ds) = ds$treat
      ds$treat = NULL
      # Check if both CTR and TET have at least 2 levels
      if (sum(is.na(ds)) > 0) {
        warning(paste("Skipping iteration", cluster, ": not enough rows"))
        res[g, cluster] = NA  
        next 
      }
      fisher = fisher.test(ds)
      p_value <- fisher[["p.value"]]
      print(p_value)
      res[g, cluster] = p_value
    }
  }
  
  res$shRNA = rownames(res)
  colnames(res) <- as.character(colnames(res))
  res <- res %>%
    pivot_longer(cols = c('1':as.character(ncol(res)-1)), names_to = "clusters", values_to = "p_value")
  
  # Add a q-value column using p.adjust
  res$padj <- p.adjust(res$p_value, method = "fdr")
  res$sig <- res$padj <= 0.05
  
  results = merge(res, complete_combinded)
  
  ## FC between CTR and TET
  # Aggregate the count for each shRNA and cluster, by treatment (TET vs CTR)
  fc_shRNAs_clusters <- aggregate(count ~ shRNA + clusters + treat, data = results, sum)
  # Reshape the data so that you have counts for TET and CTR in separate columns
  fc_shRNAs_clusters_wide <- reshape(fc_shRNAs_clusters, idvar = c("shRNA", "clusters"), timevar = "treat", direction = "wide")
  # Calculate the Fold Change as TET/CTR
  fc_shRNAs_clusters_wide$FC <- fc_shRNAs_clusters_wide$count.TET / fc_shRNAs_clusters_wide$count.CTR
  
  results = merge(results, fc_shRNAs_clusters_wide)
  plotting = distinct(results[, c("shRNA", "clusters", "FC", "padj", "sig")])
  
  ## VOLCANO PLOT
  p = ggplot(data=plotting, aes(x=log2(as.numeric(FC)), y=log10(padj)*-1)) +
    geom_point(aes(colour = clusters), size=1.5) +
    ggrepel::geom_text_repel(
      data = subset(plotting, sig==T),
      aes(label = shRNA),
      size = 3,
      box.padding = unit(0.35, "lines"),
      point.padding = unit(0.3, "lines")
    )+
    #geom_vline(xintercept = log2(1.5),linetype="dotted")+
    #geom_vline(xintercept = -log2(1.5),linetype="dotted")+
    geom_hline(yintercept = -log10(0.05),linetype="dotted")+
    ggtitle("fisher test statistics per shRNA in each cluster")+
    theme_minimal() +  
    labs(
      x = "Log2FC compared to CTR",
      y = "-Log10 adj p-value",
      color = "cluster"  # Legend title
    ) 
  
  ggsave(p, filename = paste0(dir, "shRNA_volcano_fisher_stats_TETvsnoTET.pdf"),
         width = 5, height = 4)
  
} else {print("Control samples are not present")}


## VS TIMEPOINT
print(timepoint)
if (! is.na(timepoint)) {
  print(paste0("Starting timepoint is present. Comparing to ", timepoint, " and ", control_gene))
  data = distinct(metadata[c("clone", "cellsxclone", "sample", "replicate", "sample_ann")])
  data = separate(data, clone, c("gene", "barcode", "UCI"), sep = "_", remove = F)
  #data = split(data, data$sample)
  p = ggplot(data = data, aes(x=clone, y=cellsxclone, fill = sample_ann)) +
    #facet_grid(col =vars(gene)) +
    geom_col(position = "dodge") +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(filename = paste0(dir, "cellsxclone_samples.pdf"), plot = p, width = 20, height = 10)
  
  p = ggplot(data = data, aes(x=gene, y=cellsxclone, fill = sample_ann)) +
    #facet_grid(col =vars(gene)) +
    geom_col(position = "dodge") +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(filename = paste0(dir, "cellsxgene_samples.pdf"), plot = p, width = 20, height = 10)
  
  #CTR
  data = distinct(data[c("clone", "cellsxclone", "sample", "replicate", "sample_ann")])
  assign(timepoint, filter(data, sample_ann == timepoint))
  tot <- sum(get(timepoint)$cellsxclone)
  assign(timepoint, mutate(get(timepoint), !!paste0("cellsxclone_", timepoint) := cellsxclone))
  assign(timepoint, mutate(get(timepoint), !!paste0("percentage_", timepoint) := cellsxclone / tot * 100))
  assign(timepoint, get(timepoint)[c("clone", paste0("cellsxclone_", timepoint), paste0("percentage_", timepoint))])
  data <- merge(data, get(timepoint), by = "clone", all = TRUE)
  
  #FC 
  vec = unique(data$sample_ann)
  vec <- vec[vec != timepoint]
  
  for (v in vec) {
    print(v)
    assign(v, filter(data, sample_ann == v))
    assign(v, mutate(get(v), !!paste0("cellsxclone_", v) := cellsxclone))
    tot <- sum(get(v)$cellsxclone)
    assign(v, mutate(get(v), !!paste0("percentage_", v) := cellsxclone / tot * 100))
    assign(v, get(v)[c("clone", paste0("cellsxclone_", v), paste0("percentage_", v))])
    data <- merge(data, get(v), by = "clone", all = TRUE)
    data[is.na(data)] <- 0
    
    wide = mutate(data, !!paste0("log2Percentage_", timepoint) := log2(!!sym(paste0("percentage_", timepoint))))
    wide = mutate(wide, !!paste0("FC_", v) := !!sym(paste0("percentage_", v)) / !!sym(paste0("percentage_", timepoint)))
    wide = mutate(wide, !!paste0("log2FC_", v) := log2(!!sym(paste0("FC_", v))))
    wide[is.na(wide)] <- 0
    wide = separate(wide, clone, c("gene", "barcode", "UCI"), sep = "_", remove = F)
    
    wide_big <- filter(wide, wide[[paste0("cellsxclone_", v)]] > 19 | wide[[paste0("cellsxclone_", timepoint)]] > 19)
    
    # Plot both columns as y-values
    p <- ggplot(data = wide_big, aes(x = clone, y = .data[[paste0("log2FC_", v)]], fill = gene)) +
      geom_col(position = "dodge") +  # Use dodge to separate bars
      theme(axis.text.x = element_text(angle = 90))
    
    ggsave(filename = paste0(dir, "log2FC_bigclones_", v, ".pdf"), plot = p, width = 15, height = 10)
    
    
    #p-value
    #WHICH IS THE COMPARISON FOR 2D??? SCR
    CG = filter(wide, gene == control_gene)
    CG = CG[, c(paste0("cellsxclone_", timepoint), paste0("cellsxclone_", v))]
    CG = colSums(CG, na.rm = TRUE)
    pvalues = data.frame()
    for (c in unique(wide$clone)) {
      clone_vector = filter(wide, clone == c)
      clone_vector = clone_vector[, c(paste0("cellsxclone_", timepoint),paste0("cellsxclone_", v))]
      test = rbind(clone_vector, CG)
      pvalues[c,"pvalue"] = fisher.test(test)
    }
    
    pvalues$padj = p.adjust(pvalues$pvalue, method = "BH")
    pvalues$clone = rownames(pvalues)
    wide = merge(wide, pvalues, by = "clone")
    
    wide_big = filter(wide, if_any(matches("^cellsxclone"), ~ . > 19))
    p = ggplot(data = wide_big, aes(x=!!sym(paste0("log2FC_", v)), y= !!sym(paste0("log2Percentage_", timepoint)), colour = -log10(padj))) + #colour = p.value
      geom_point() +
      geom_vline(xintercept = c(-1, 1), linetype = "dotted")
    ggsave(filename = paste0(dir, "volcano_plot_enrich_", v, "_bigclones.pdf"), plot = p)
  }
  
  #STATS
  colData(cds)$clusters = clusters(cds)
  df = as.data.frame(colData(cds))
  
  ## GENE LEVEL VS SCR!!!
  # gene in clusters vs SCR in clusters
  assign(paste0("data_", timepoint), filter(df, sample_ann == timepoint))
  
  # Calculate the number of cells per cluster and condition
  assign(paste0("cluster_condition_counts_", timepoint), 
         get(paste0("data_", timepoint)) %>%
           group_by(clusters, gene) %>%
           summarise(count = n(), .groups = 'drop'))
  
  # Calculate the percentage of KD and Control cells per cluster
  assign(paste0("cluster_condition_percentage_", timepoint), 
         get(paste0("cluster_condition_counts_", timepoint)) %>%
           group_by(clusters) %>%
           mutate(total = sum(count)) %>%
           mutate(percentage = (count / total) * 100) %>%
           ungroup())
  
  for (v in vec) {
    assign(paste0("data_", v), filter(df, sample_ann == v))
    
    # Calculate the number of cells per cluster and condition
    assign(paste0("cluster_condition_counts_", v), 
           get(paste0("data_", v)) %>%
             group_by(clusters, gene) %>%
             summarise(count = n(), .groups = 'drop'))
    
    # Calculate the percentage of KD and Control cells per cluster
    assign(paste0("cluster_condition_percentage_", v), 
           get(paste0("cluster_condition_counts_", v)) %>%
             group_by(clusters) %>%
             mutate(total = sum(count)) %>%
             mutate(percentage = (count / total) * 100) %>%
             ungroup())
    
    # Store percentages with new column names
    assign(paste0("cluster_condition_percentage_", v), get(paste0("cluster_condition_percentage_", v)) %>%
             mutate(!!paste0("percentage_", v) := percentage) %>%
             select(-percentage))
    
    assign(paste0("cluster_condition_percentage_", timepoint), get(paste0("cluster_condition_percentage_", timepoint)) %>%
             mutate(!!paste0("percentage_", timepoint) := percentage) %>%
             select(-percentage))
    
    # Merge the two data frames
    combined_df <- merge(
      get(paste0("cluster_condition_percentage_", v))[, c("clusters", "gene", paste0("percentage_", v))], 
      get(paste0("cluster_condition_percentage_", timepoint))[, c("clusters", "gene", paste0("percentage_", timepoint))],
      all = T)
    combined_df[is.na(combined_df)] <- 0
    # Add treatment labels
    assign(paste0("cluster_condition_percentage_", v), get(paste0("cluster_condition_percentage_", v)) %>%
             mutate(treat = v))
    
    assign(paste0("cluster_condition_percentage_", timepoint), get(paste0("cluster_condition_percentage_", timepoint)) %>%
             mutate(treat = timepoint))
    
    
    
    # Standardize column names
    # Rename columns and reassign the object
    assign(paste0("cluster_condition_percentage_", v), 
           setNames(get(paste0("cluster_condition_percentage_", v)), 
                    c("clusters", "gene", "count", "total", "percentage", "treat")))
    
    #set names
    assign(paste0("cluster_condition_percentage_", timepoint), 
           setNames(get(paste0("cluster_condition_percentage_", timepoint)), 
                    c("clusters", "gene", "count", "total", "percentage", "treat")))
    # Combine datasets
    combined_df <- rbind(get(paste0("cluster_condition_percentage_", timepoint)), 
                         get(paste0("cluster_condition_percentage_", v)))
    
    # Plot
    p <- ggplot(combined_df, aes(x = factor(treat), y = percentage, fill = gene)) +
      geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
      facet_grid(~clusters)
    
    ggsave(p, filename = paste0(dir, "genes_in_clusters_", v, ".pdf"), width = 20, height = 10)
    
    
    ##FISHER TEST
    complete_combinded <- combined_df %>%
      complete(clusters, gene, treat, fill = list(count = 0))
    complete_combinded <- complete_combinded %>%
      group_by(clusters, treat) %>%
      mutate(total = ifelse(is.na(total), max(total, na.rm = F), total)) %>%
      ungroup()
    complete_combinded$gene[is.na(complete_combinded$gene)] <- "NA"
    complete_combinded[is.na(complete_combinded)] <- 0
    
    #calculate fraction of gene over control_gene
    complete_combinded = complete_combinded %>%
      group_by(clusters, treat) %>%
      mutate(X = unique(percentage[gene == control_gene]), 
             fraction_over_controlgene = percentage / X)
    
    
    complete_combinded_w <- complete_combinded %>%
      group_by(gene, treat, clusters) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      pivot_wider(names_from = clusters, values_from = count, values_fill = 0)
    
    # assign(paste0("for_pvalue_", v), 
    #        get(paste0("cluster_condition_percentage_", v)) %>% 
    #          complete(clusters, gene, treat, fill = list(count = 0))
    #        , envir = .GlobalEnv)
    # 
    # assign(paste0("for_pvalue_", v), 
    #        get(paste0("for_pvalue_", v)) %>%
    #   group_by(clusters) %>%
    #   mutate(total = ifelse(is.na(total), max(total, na.rm = TRUE), total)) %>%
    #   ungroup())
    # 
    # assign(paste0("for_pvalue_", timepoint), 
    #        get(paste0("cluster_condition_percentage_", timepoint)) %>% 
    #          complete(clusters, gene, treat, percentage, fill = list(count = 0))
    #        , envir = .GlobalEnv)
    # 
    # assign(paste0("for_pvalue_", timepoint), 
    #        get(paste0("for_pvalue_", timepoint)) %>%
    #          group_by(clusters) %>%
    #          mutate(total = ifelse(is.na(total), max(total, na.rm = TRUE), total)) %>%
    #          ungroup())
    complete_combinded$X = NULL
    list = c(unique(complete_combinded$gene))
    list = setdiff(list, control_gene)
    res = data.frame(row.names = list)
    for (g in list) {
      print(g)
      for (c in unique(complete_combinded$clusters)) {
        print(c)
        start = distinct(filter(complete_combinded, treat == v))
        start = filter(start, clusters == c)
        ds = filter(start, gene == g)[,c("treat", "count", "total")]
        ds = mutate(ds, out = total - count)
        ds$total = NULL
        #rownames(ds) = ds$treat
        ds$treat = NULL
        #compare to control_gene
        ctr = filter(start, gene == control_gene)[,c("treat", "count", "total")]
        ctr = mutate(ctr, out = total - count)
        ctr$total = NULL
        #rownames(ctr) = ctr$treat
        ctr$treat = NULL
        # Check if both CTR and TET have at least 2 levels
        ds = rbind(ds, ctr)
        #print(ds)
        #rownames(ds) = c(v, "ctr")
        if (sum(is.na(ds)) > 0) {
          warning(paste("Skipping iteration", c, ": not enough rows"))
          res[g, c] = NA  
          next 
        }
        fisher = fisher.test(ds)
        p_value <- fisher[["p.value"]]
        print(p_value)
        res[g, c] = p_value
      }
    }
    
    res$gene = rownames(res)
    colnames(res) <- as.character(colnames(res))
    res <- res %>%
      pivot_longer(cols = c('1':as.character(ncol(res)-1)), names_to = "clusters", values_to = "p_value")
    
    # Add a q-value column using p.adjust
    res$padj <- p.adjust(res$p_value, method = "BH")
    res$sig <- res$padj <= 0.05
    
    results = left_join(res, filter(complete_combinded, treat == v), by = c("gene", "clusters"))
    
    ## FC between CTR and TET
    # Aggregate the count for each gene and cluster, by treatment (TET vs CTR)
    fc_genes_clusters <- aggregate(count ~ gene + clusters + treat, data = results, sum)
    # Reshape the data so that you have counts for TET and CTR in separate columns
    fc_genes_clusters_wide <- reshape(fc_genes_clusters, idvar = c("gene", "clusters"), timevar = "treat", direction = "wide")
    write.table(fc_genes_clusters_wide, paste0(dir, "FC_results_", v, "vs", timepoint, ".txt"), quote = F)
    # # Calculate the Fold Change as TET/CTR
    # fc_genes_clusters_wide <- fc_genes_clusters_wide[fc_genes_clusters_wide[[paste0("count.", timepoint)]] != 0, ]
    # fc_genes_clusters_wide <- fc_genes_clusters_wide[fc_genes_clusters_wide[[paste0("count.", v)]] != 0, ]
    # fc_genes_clusters_wide$FC <- fc_genes_clusters_wide[[paste0("count.", v)]] / 
    #   fc_genes_clusters_wide[[paste0("count.", timepoint)]]
    # results = merge(results, fc_genes_clusters_wide)
    # plotting = distinct(results[, c("gene", "clusters", "FC", "padj", "sig")])
    
    results = mutate(results, log2_fraction = log2(fraction_over_controlgene))
    
    ## VOLCANO PLOT
    p = ggplot(data=results, aes(x=log2_fraction, y=log10(padj)*-1)) +
      geom_point(aes(colour = clusters), size=1.5) +
      ggrepel::geom_text_repel(
        data = subset(results, sig==T),
        aes(label = gene),
        size = 3,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      )+
      geom_vline(xintercept = 0,linetype="dotted")+
      #geom_vline(xintercept = -log2(1.5),linetype="dotted")+
      geom_hline(yintercept = -log10(0.05),linetype="dotted")+
      ggtitle("fisher test statistics per gene in each cluster")+
      theme_minimal() +  
      labs(
        x = paste0("Log2 fraction over ", control_gene),
        y = "-Log10 adj p-value",
        color = "cluster"  # Legend title
      ) 
    
    ggsave(p, filename = paste0(dir, "gene_volcano_fisher_stats_", v, ".pdf"),
           width = 5, height = 4)
    
  }
  
  
  
  
  #STATS
  colData(cds)$clusters = clusters(cds)
  df = as.data.frame(colData(cds))
  
  ## shRNA LEVEL VS SCR!!!
  # shRNA in clusters vs SCR in clusters
  assign(paste0("data_", timepoint), filter(df, sample_ann == timepoint))
  
  # Calculate the number of cells per cluster and condition
  assign(paste0("cluster_condition_counts_", timepoint), 
         get(paste0("data_", timepoint)) %>%
           group_by(clusters, shRNA) %>%
           summarise(count = n(), .groups = 'drop'))
  
  # Calculate the percentage of KD and Control cells per cluster
  assign(paste0("cluster_condition_percentage_", timepoint), 
         get(paste0("cluster_condition_counts_", timepoint)) %>%
           group_by(clusters) %>%
           mutate(total = sum(count)) %>%
           mutate(percentage = (count / total) * 100) %>%
           ungroup())
  
  for (v in vec) {
    assign(paste0("data_", v), filter(df, sample_ann == v))
    
    # Calculate the number of cells per cluster and condition
    assign(paste0("cluster_condition_counts_", v), 
           get(paste0("data_", v)) %>%
             group_by(clusters, shRNA) %>%
             summarise(count = n(), .groups = 'drop'))
    
    # Calculate the percentage of KD and Control cells per cluster
    assign(paste0("cluster_condition_percentage_", v), 
           get(paste0("cluster_condition_counts_", v)) %>%
             group_by(clusters) %>%
             mutate(total = sum(count)) %>%
             mutate(percentage = (count / total) * 100) %>%
             ungroup())
    
    # Store percentages with new column names
    assign(paste0("cluster_condition_percentage_", v), get(paste0("cluster_condition_percentage_", v)) %>%
             mutate(!!paste0("percentage_", v) := percentage) %>%
             select(-percentage))
    
    assign(paste0("cluster_condition_percentage_", timepoint), get(paste0("cluster_condition_percentage_", timepoint)) %>%
             mutate(!!paste0("percentage_", timepoint) := percentage) %>%
             select(-percentage))
    
    # Merge the two data frames
    combined_df <- merge(
      get(paste0("cluster_condition_percentage_", v))[, c("clusters", "shRNA", paste0("percentage_", v))], 
      get(paste0("cluster_condition_percentage_", timepoint))[, c("clusters", "shRNA", paste0("percentage_", timepoint))],
      all = T)
    combined_df[is.na(combined_df)] <- 0
    # Add treatment labels
    assign(paste0("cluster_condition_percentage_", v), get(paste0("cluster_condition_percentage_", v)) %>%
             mutate(treat = v))
    
    assign(paste0("cluster_condition_percentage_", timepoint), get(paste0("cluster_condition_percentage_", timepoint)) %>%
             mutate(treat = timepoint))
    
    
    
    # Standardize column names
    # Rename columns and reassign the object
    assign(paste0("cluster_condition_percentage_", v), 
           setNames(get(paste0("cluster_condition_percentage_", v)), 
                    c("clusters", "shRNA", "count", "total", "percentage", "treat")))
    
    #set names
    assign(paste0("cluster_condition_percentage_", timepoint), 
           setNames(get(paste0("cluster_condition_percentage_", timepoint)), 
                    c("clusters", "shRNA", "count", "total", "percentage", "treat")))
    # Combine datasets
    combined_df <- rbind(get(paste0("cluster_condition_percentage_", timepoint)), 
                         get(paste0("cluster_condition_percentage_", v)))
    
    # Plot
    p <- ggplot(combined_df, aes(x = factor(treat), y = percentage, fill = shRNA)) +
      geom_bar(stat = "identity", colour = "white", linewidth = 0.3) +
      facet_grid(~clusters)
    
    ggsave(p, filename = paste0(dir, "shRNAs_in_clusters_", v, ".pdf"), width = 20, height = 10)
    
    
    ##FISHER TEST
    complete_combinded <- combined_df %>%
      complete(clusters, shRNA, treat, fill = list(count = 0))
    complete_combinded <- complete_combinded %>%
      group_by(clusters, treat) %>%
      mutate(total = ifelse(is.na(total), max(total, na.rm = F), total)) %>%
      ungroup()
    complete_combinded$shRNA[is.na(complete_combinded$shRNA)] <- "NA"
    complete_combinded[is.na(complete_combinded)] <- 0
    
    #calculate fraction of shRNA over control_gene
    complete_combinded = complete_combinded %>%
      group_by(clusters, treat) %>%
      mutate(X = unique(percentage[shRNA == control_gene]), 
             fraction_over_controlshRNA = percentage / X)
    
    
    complete_combinded_w <- complete_combinded %>%
      group_by(shRNA, treat, clusters) %>%
      summarise(count = sum(count), .groups = "drop") %>%
      pivot_wider(names_from = clusters, values_from = count, values_fill = 0)
    
    
    complete_combinded$X = NULL
    list = c(unique(complete_combinded$shRNA))
    list <- list[!grepl(control_gene, list)]
    res = data.frame(row.names = list)
    for (g in list) {
      print(g)
      for (c in unique(complete_combinded$clusters)) {
        print(c)
        start = distinct(filter(complete_combinded, treat == v))
        start = filter(start, clusters == c)
        ds = filter(start, shRNA == g)[,c("treat", "count", "total")]
        ds = mutate(ds, out = total - count)
        ds$total = NULL
        #rownames(ds) = ds$treat
        ds$treat = NULL
        #compare to control_gene
        ctr = filter(start, shRNA == control_gene)[,c("treat", "count", "total")]
        ctr = mutate(ctr, out = total - count)
        ctr$total = NULL
        #rownames(ctr) = ctr$treat
        ctr$treat = NULL
        # Check if both CTR and TET have at least 2 levels
        ds = rbind(ds, ctr)
        #print(ds)
        #rownames(ds) = c(v, "ctr")
        if (sum(is.na(ds)) > 0) {
          warning(paste("Skipping iteration", c, ": not enough rows"))
          res[g, c] = NA  
          next 
        }
        fisher = fisher.test(ds)
        p_value <- fisher[["p.value"]]
        print(p_value)
        res[g, c] = p_value
      }
    }
    
    res$shRNA = rownames(res)
    colnames(res) <- as.character(colnames(res))
    res <- res %>%
      pivot_longer(cols = c('1':as.character(ncol(res)-1)), names_to = "clusters", values_to = "p_value")
    
    # Add a q-value column using p.adjust
    res$padj <- p.adjust(res$p_value, method = "BH")
    res$sig <- res$padj <= 0.05
    
    results = left_join(res, filter(complete_combinded, treat == v), by = c("shRNA", "clusters"))
    
    ## FC between CTR and TET
    # Aggregate the count for each shRNA and cluster, by treatment (TET vs CTR)
    fc_shRNAs_clusters <- aggregate(count ~ shRNA + clusters + treat, data = results, sum)
    # Reshape the data so that you have counts for TET and CTR in separate columns
    fc_shRNAs_clusters_wide <- reshape(fc_shRNAs_clusters, idvar = c("shRNA", "clusters"), timevar = "treat", direction = "wide")
    write.table(fc_shRNAs_clusters_wide, paste0(dir, "FC_results_", v, "vs", timepoint, ".txt"), quote = F)
    # Calculate the Fold Change as TET/CTR
    # fc_shRNAs_clusters_wide <- fc_shRNAs_clusters_wide[fc_shRNAs_clusters_wide[[paste0("count.", timepoint)]] != 0, ]
    # fc_shRNAs_clusters_wide <- fc_shRNAs_clusters_wide[fc_shRNAs_clusters_wide[[paste0("count.", v)]] != 0, ]
    # fc_shRNAs_clusters_wide$FC <- fc_shRNAs_clusters_wide[[paste0("count.", v)]] / 
    #   fc_shRNAs_clusters_wide[[paste0("count.", timepoint)]]
    # results = merge(results, fc_shRNAs_clusters_wide)
    # plotting = distinct(results[, c("shRNA", "clusters", "FC", "padj", "sig")])
    
    results = mutate(results, log2_fraction = log2(fraction_over_controlshRNA))
    
    ## VOLCANO PLOT
    p = ggplot(data=results, aes(x=log2_fraction, y=log10(padj)*-1)) +
      geom_point(aes(colour = clusters), size=1.5) +
      ggrepel::geom_text_repel(
        data = subset(results, sig==T),
        aes(label = shRNA),
        size = 3,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      )+
      geom_vline(xintercept = 0,linetype="dotted")+
      #geom_vline(xintercept = -log2(1.5),linetype="dotted")+
      geom_hline(yintercept = -log10(0.05),linetype="dotted")+
      ggtitle("fisher test statistics per shRNA in each cluster")+
      theme_minimal() +  
      labs(
        x = paste0("Log2 fraction over ", control_gene),
        y = "-Log10 adj p-value",
        color = "cluster"  # Legend title
      ) 
    
    ggsave(p, filename = paste0(dir, "shRNA_volcano_fisher_stats_", v, ".pdf"),
           width = 5, height = 4)
    
  }
} else {print("Starting timepoint is not present")}

