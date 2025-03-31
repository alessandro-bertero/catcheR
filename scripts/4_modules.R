#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
#ARGS
#dir = "/Users/marialuisaratto/scripts/aggr/"
dir = args[1]
#file = cds
file = args[2]
resolution = as.numeric(args[3])

if (file.exists(paste0(dir,"/", file)) == F) {
  stop(paste("Error: cds file does not exist at path:", paste0(dir,"/", file)))
}

print("Starting...")

suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(dplyr))
suppressMessages(library(monocle3))

load(paste0(dir,"/", file))

#REGRESSION
cds_KD <- cds[
  ,
  colData(cds) %>%
    subset(
      sampleC == "KD"
    ) %>%
    row.names
]
colData(cds_KD)$clusters = clusters(cds_KD)
print("regression...")
# # Fit models
# gene_fits <- fit_models(cds_KD, model_formula_str = paste0("~cluster + replicate"))
# print("genes fitted for cluster...")
# # Extract coefficient tables and adjust with BH
# fit_coefs <- coefficient_table(gene_fits)
# 
# regression = fit_coefs %>% filter(term != "(Intercept)")
# regression = regression %>% filter (q_value < 0.05) %>%
#   select(gene_short_name, term, q_value, estimate)
# write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_cluster.csv"), row.names=FALSE)
# 
# # Fit models for perturbations GENE
# gene_fits <- fit_models(cds_KD, model_formula_str = paste0("~gene + replicate"))
# print("genes fitted for genes...")
# # Extract coefficient tables and adjust with BH
# fit_coefs <- coefficient_table(gene_fits)
# 
# regression = fit_coefs %>% filter(term != "(Intercept)")
# regression = regression %>% filter (q_value < 0.05) %>%
#   select(gene_short_name, term, q_value, estimate)
# write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_gene.csv"), row.names=FALSE)
# 
# # Fit models for perturbations shRNA
# gene_fits <- fit_models(cds_KD, model_formula_str = paste0("~shRNA + replicate"))
# print("genes fitted for shRNAs...")
# # Extract coefficient tables and adjust with BH
# fit_coefs <- coefficient_table(gene_fits)
# 
# regression = fit_coefs %>% filter(term != "(Intercept)")
# regression = regression %>% filter (q_value < 0.05) %>%
#   select(gene_short_name, term, q_value, estimate)
# write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_shRNA.csv"), row.names=FALSE)
# 
# # Fit models for perturbations clone
# gene_fits <- fit_models(cds_KD, model_formula_str = paste0("~clone + replicate"))
# print("genes fitted for clones...")
# # Extract coefficient tables and adjust with BH
# fit_coefs <- coefficient_table(gene_fits)
# 
# regression = fit_coefs %>% filter(term != "(Intercept)")
# regression = regression %>% filter (q_value < 0.05) %>%
#   select(gene_short_name, term, q_value, estimate)
# write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_clone.csv"), row.names=FALSE)

#GRAPH AUTOCORRELATION
print("Modules for genes...")
pr_graph_test_res <- graph_test(cds_KD, neighbor_graph="knn", cores=2)
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
#gene module
gene_module_df <- find_gene_modules(cds_KD[pr_deg_ids,], resolution=resolution)
write.csv(gene_module_df, file= paste0(dir,"/gene_modules.csv"), row.names=FALSE)
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_KD)),
                                cell_group=colData(cds_KD)$gene)
agg_mat <- aggregate_gene_expression(cds_KD, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Gene ", colnames(agg_mat))
#plot
pdf(paste0(dir, "/pheatmap_gene_module.pdf"), width=7, height=10)
pheatmap::pheatmap(
  agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
  scale="column", clustering_method="ward.D2",
  fontsize=6)
dev.off()

#TOP modules 
#row diff
row_range <- apply(agg_mat, 1, function(x) max(x) - min(x))
sorted_rows_by_range <- agg_mat[order(-row_range), ]
top_n_rows_by_range <- head(sorted_rows_by_range, 10)
#plot
pdf(paste0(dir, "/pheatmap_gene_module_top10.pdf"), width=7, height=10)
pheatmap::pheatmap(
  top_n_rows_by_range, cluster_rows=TRUE, cluster_cols=TRUE,
  scale="column", clustering_method="ward.D2",
  fontsize=6)
dev.off()

#shRNA and module
print("Modules for shRNA...")
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_KD)),
                                cell_group=colData(cds_KD)$shRNA)
agg_mat <- aggregate_gene_expression(cds_KD, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("shRNA ", colnames(agg_mat))
#plot
pdf(paste0(dir, "/pheatmap_shRNA_module.pdf"), width=7, height=10)
pheatmap::pheatmap(
  agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
  scale="column", clustering_method="ward.D2",
  fontsize=6)
dev.off()

#TOP modules 
#row diff
row_range <- apply(agg_mat, 1, function(x) max(x) - min(x))
sorted_rows_by_range <- agg_mat[order(-row_range), ]
top_n_rows_by_range <- head(sorted_rows_by_range, 10)
#plot
pdf(paste0(dir, "/pheatmap_shRNA_module_top10.pdf"), width=7, height=10)
pheatmap::pheatmap(
  top_n_rows_by_range, cluster_rows=TRUE, cluster_cols=TRUE,
  scale="column", clustering_method="ward.D2",
  fontsize=6)
dev.off()

# #clone and module
# cell_group_df <- tibble::tibble(cell=row.names(colData(cds_KD)),
#                                 cell_group=colData(cds_KD)$clone)
# agg_mat <- aggregate_gene_expression(cds_KD, gene_module_df, cell_group_df)
# row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
# colnames(agg_mat) <- stringr::str_c("Clone ", colnames(agg_mat))
# #plot
# pdf(paste0(dir, "/pheatmap_clone_module.pdf"), width=7, height=10)
# pheatmap::pheatmap(
#   agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
#   scale="column", clustering_method="ward.D2",
#   fontsize=6)
# dev.off()
# 
# #TOP modules 
# #row diff
# row_range <- apply(agg_mat, 1, function(x) max(x) - min(x))
# sorted_rows_by_range <- agg_mat[order(-row_range), ]
# top_n_rows_by_range <- head(sorted_rows_by_range, 10)
# #plot
# pdf(paste0(dir, "/pheatmap_clone_module_top10.pdf"), width=7, height=10)
# pheatmap::pheatmap(
#   top_n_rows_by_range, cluster_rows=TRUE, cluster_cols=TRUE,
#   scale="column", clustering_method="ward.D2",
#   fontsize=6)
# dev.off()

#module x cell
print("Modules for cells...")
cell_group_df <- tibble::tibble(cell=row.names(colData(cds_KD)),
                                cell_group=row.names(colData(cds_KD)))
agg_mat <- aggregate_gene_expression(cds_KD, gene_module_df, cell_group_df)
agg_mat = as.data.frame(t(agg_mat))
dir.create(paste0(dir, "/modules_cells"))
for (mod in colnames(agg_mat)) {
  write.csv(as.data.frame(agg_mat[,mod, drop = FALSE]), paste0(dir, "/modules_cells/module_", mod, ".csv"), row.names = T)
}

#VS NO TET
colData(cds)$clusters = clusters(cds)
meta = as.data.frame(colData(cds)@listData)
if ("control" %in% meta$sampleC) {
  print("control sample analysis are present: running modules analyisis with paired controls...")
  
  colData(cds)$gene_sample <- paste(colData(cds)$gene, colData(cds)$sampleC, sep = "_")
  
  #GRAPH AUTOCORRELATION
  pr_graph_test_res <- graph_test(cds, neighbor_graph="knn", cores=2)
  pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
  #gene module
  #gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)),
                                  cell_group=colData(cds)$gene_sample)
  agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
  colnames(agg_mat) <- stringr::str_c("Gene ", colnames(agg_mat))
  #plot
  pdf(paste0(dir, "/pheatmap_vsNOTET_module.pdf"), width=7, height=10)
  pheatmap::pheatmap(
    agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
    scale="column", clustering_method="ward.D2",
    fontsize=6)
  dev.off()
}

