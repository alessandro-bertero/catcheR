#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
print(args)
#ARGS
#dir = "/Users/marialuisaratto/scripts/aggr/"
dir = args[1]
#file = "filtered_annotated_silencing_matrix_complete_all_samples.csv"
file = args[2]
#control_genes = c("B2M", "SCR", NA)
control_genes = as.data.frame(read.delim(paste0(dir, "/", args[3]), header = F))
#control_samples = c(1,3)

# Directory path from args
file_control_samples <- paste0(dir, "/", args[4])
print(paste("Control samples file:", file_control_samples))
file_replicates <- paste0(dir, "/", args[5])
print(paste("Replicates file:", file_replicates))
file_sample_ann <- paste0(dir, "/", args[6])
print(paste("File samples annotation:", file_sample_ann))

# Load control samples file if it exists, else give a warning
if (file.exists(file_control_samples)) {
  control_samples = as.data.frame(read.delim(file_control_samples, header = F))
  print(paste("Control samples:", control_samples))
} else {
  warning(paste("Control samples not specified"))
}

# Load replicates file if it exists, else give a warning
if (file.exists(file_replicates)) {
  replicates = as.data.frame(read.delim(file_replicates, header = F))
  print(paste0("Replicates:", replicates))
} else {
  warning(paste("Replicates not specified"))
}

## ADD sample annotation
sample_ann = as.data.frame(read.csv(file_sample_ann, header = F))
names(sample_ann) = c("sample", "sample_ann")
print(paste("Sample annotation:", sample_ann))


#res = 5e-4
res = as.numeric(args[7])
#genes = "genelist.txt"
genes = args[8]

# #INSTALL MONOCLE
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.19")
# 
# BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
#                        'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
#                        'SummarizedExperiment', 'batchelor', 'HDF5Array',
#                        'terra', 'ggrastr'))
# 
# install.packages("devtools")
# devtools::install_github('cole-trapnell-lab/monocle3')
suppressMessages(library(devtools))
suppressMessages(library(monocle3))
suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(gtools))
suppressMessages(library(dplyr))
set.seed(1111)

#LOAD ANNOTATED GENE EXP
#exp = read.csv("/30tb/3tb/data/ratto/testing/annotated_silencing_matrix_complete_all_samples.csv", header = T)
exp = read.csv(paste0(dir,"/", file), header = T)
exp = exp[, colSums(exp != 0) > 0]
exp <- separate(
  exp,
  col="X",
  into=c("a","b"),
  sep = ":",
  remove = TRUE,
  convert = FALSE,
  extra = "warn",
  fill = "warn"
)
exp$a <- NULL
names <- make.unique(exp$b, sep=".")
rownames(exp) <- names
exp$b <- NULL

#CREATE ANNOTATION DATAFRAME
data = data.frame()
data = as.data.frame(colnames(exp))
colnames(data) = c("nomi")
rownames(data)=data$nomi
data = separate(data, nomi, into = c("cellID", "UMI", "barcode", "gene", "UCI", "sample"), sep = "_", remove = F, convert = T)
data = mutate(data, clone = paste(gene, barcode, UCI, sep = "_"))
data$UMI = as.integer(data$UMI)
data$sample = as.integer(data$sample)

#cellsxclone
# cellsxclone = distinct(data[,c("cellID", "clone")])
# cellsxclone_table = data.frame(table(cellsxclone$clone))
# names(cellsxclone_table) = c("clone", "cellsxclone")
# data = merge(data,cellsxclone_table)
# rownames(data) = data$nomi

cellsxclone = distinct(data[, c("sample", "cellID", "clone")])
cellsxclone_table = cellsxclone %>%
  group_by(sample, clone) %>%
  summarise(cellsxclone = n(), .groups = "drop")
data = merge(data, cellsxclone_table, by = c("sample", "clone"), all.x = TRUE)
rownames(data) = data$nomi

#Load genes correspondence and add gene column
gene_ass = read.delim(paste(dir,"/rc_barcodes_genes.csv",sep=""), header = F, sep = ",", col.names = c("barcode","shRNA"))
data = left_join(data, gene_ass, by = "barcode")
data = left_join(data, sample_ann, by = "sample")

#ADD CONTROLS AND REPLICATES
if (exists("control_samples")) {
  data = mutate(data, control = case_when(gene %in% control_genes$V1 | sample %in% control_samples$V1 ~ "control",
                                        .default = "KD"))

  data = mutate(data, sampleC = case_when(sample %in% control_samples$V1 ~ "control",
                                          .default = "KD"))
} else {
  data = mutate(data, control = case_when(gene %in% control_genes$V1 ~ "control",
                                          .default = "KD"))
  
  data$sampleC = "KD"
}

data = mutate(data, geneC = case_when(gene %in% control_genes$V1 ~ "control",
                                      .default = "KD"))

if (exists("replicates")) {
replicate = data.frame(sample = sort(unique(data$sample)), replicate = replicates$V1)
data = merge(data, replicate)
} else {
  data$replicate = 0
}

#update names
data = mutate(data, nomi = paste(nomi, shRNA, cellsxclone, control, replicate, sampleC, geneC, sample_ann, sep = "_"))
rownames(data) = data$nomi

#ORDER IN SAME WAY MATRIX AND ANNOTATION 
data_complete = exp[, order(colnames(exp))]
data = data[order(row.names(data)), ]
colnames(data_complete) = row.names(data)

#save data
write.csv(data_complete, paste(dir, "/espression_data.csv",sep=""))
write.csv(data, paste(dir, "/cell_metadata.csv",sep=""))

#CREATE MONOCLE DATASET
data_matrix = as.matrix(data_complete)
starting_cds <- new_cell_data_set(expression_data = data_matrix,
                                  cell_metadata = data,
                                  gene_metadata = NULL)
# Save dataframe as an RData file
save(starting_cds, file = paste0(dir, "/starting_cds.RData"))

#NORMALIZE for size and log
norm = normalized_counts(
  starting_cds,
  norm_method = "log",
  pseudocount = 1
)

cds <- new_cell_data_set(expression_data = norm,
                         cell_metadata = data,
                         gene_metadata = NULL)

# Save dataframe as an RData file
#save(cds, file = paste0(dir, "/norm_cds.RData"))

cds <- preprocess_cds(cds, num_dim = 100)
p = plot_pc_variance_explained(cds)
ggsave(p, filename = paste0(dir,"/PCA.pdf"),
       width = 5, height = 5)
cds <- align_cds(cds, alignment_group = "replicate", useNames = TRUE)
cds <- reduce_dimension(cds, reduction_method="UMAP", umap.fast_sgd = FALSE,cores=1,n_sgd_threads=1)

UMAP = as.data.frame(reducedDims(cds)$UMAP)
names(UMAP) = c("x", "y")
UMAP$nomi = rownames(UMAP)
UMAP = separate(UMAP, nomi, into = c("cellID", "UMI", "barcode", "gene", "UCI", "sample", "shRNA", "cellsxclone", "control", "replicate", "sampleC", "geneC", "sample_ann"), sep = "_", remove = F, convert = T)
p = ggplot(UMAP, aes(x, y, color = sample_ann, alpha = 0.1)) +
  geom_point(size = 0.4)
ggsave(p, filename = paste0(dir,"/UMAP.pdf"),
       width = 6, height = 6)

#PLOT GENE EXPRESSION
glist = unique(t(read.delim(paste0(dir, "/", genes), sep = ",", header = F)))
#put gene names in $gene_short_name to solve bug
rowData(cds)$gene_name <- rownames(cds)
rowData(cds)$gene_short_name <- rowData(cds)$gene_name

p = plot_cells(cds,
               genes = glist,
               label_cell_groups = T,
               show_trajectory_graph = FALSE)+ 
  labs(x = "UMAP 1", y = "UMAP 2", title = "")#+
#scale_color_viridis(option="E", discrete=F)

ggsave(p, filename = paste0(dir,"/UMAP_gene_expression.pdf"),
       width = 10, height = 10)

#CLUSTERING
cds <- cluster_cells(cds, resolution=res, random_seed = 42)
p = plot_cells(cds)
ggsave(p, filename = paste0(dir,"/UMAP_clustering.pdf"),
       width = 5, height = 5)
colData(cds)$clusters = clusters(cds)


#trajectories
print("Learning trajectories...")
cds <- learn_graph(cds)
p = plot_cells(cds,
               color_cells_by = "cluster",
               label_groups_by_cluster=FALSE,
               label_leaves=FALSE,
               label_branch_points=FALSE)
ggsave(p, filename = paste0(dir,"/UMAP_clustering_trajectories.pdf"),
       width = 5, height = 5)

# Save dataas an RData file
save(cds, file = paste0(dir, "/processed_cds.RData"))

#MARKER GENES
print("Finding clusters marker genes...")
marker_test_res <- top_markers(cds, group_cells_by="cluster",
                               reference_cells=1000, cores=8)

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))
#save data
write.csv(marker_test_res, paste(dir, "/gene_markers.csv",sep=""))
write.csv(top_specific_markers, paste(dir, "/top1_gene_markers.csv",sep=""))

print("Done!")