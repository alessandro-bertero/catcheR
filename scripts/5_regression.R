#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
#ARGS
#dir = "/Users/marialuisaratto/scripts/aggr/"
dir = args[1]
#file = cds
file = args[2]

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
# Fit models
gene_fits <- fit_models(cds_KD, model_formula_str = paste0("~cluster + replicate"))
print("genes fitted for cluster...")
# Extract coefficient tables and adjust with BH
fit_coefs <- coefficient_table(gene_fits)

regression = fit_coefs %>% filter(term != "(Intercept)")
regression = regression %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_cluster.csv"), row.names=FALSE)

# Fit models for perturbations GENE
gene_fits <- fit_models(cds_KD, model_formula_str = paste0("~gene + replicate"))
print("genes fitted for genes...")
# Extract coefficient tables and adjust with BH
fit_coefs <- coefficient_table(gene_fits)

regression = fit_coefs %>% filter(term != "(Intercept)")
regression = regression %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_gene.csv"), row.names=FALSE)

# Fit models for perturbations shRNA
gene_fits <- fit_models(cds_KD, model_formula_str = paste0("~shRNA + replicate"))
print("genes fitted for shRNAs...")
# Extract coefficient tables and adjust with BH
fit_coefs <- coefficient_table(gene_fits)

regression = fit_coefs %>% filter(term != "(Intercept)")
regression = regression %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_shRNA.csv"), row.names=FALSE)

# Fit models for perturbations clone
gene_fits <- fit_models(cds_KD, model_formula_str = paste0("~clone + replicate"))
print("genes fitted for clones...")
# Extract coefficient tables and adjust with BH
fit_coefs <- coefficient_table(gene_fits)

regression = fit_coefs %>% filter(term != "(Intercept)")
regression = regression %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
write.csv(regression, file= paste0(dir,"/regression_filtered_genes_per_clone.csv"), row.names=FALSE)
