#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
# library(ggplot2)
# library(tidyverse)
#library(see)
args = commandArgs(trailingOnly=TRUE)

#Set working direcotry
#wd = "/2tb/torino/ratto/bertero/H001AS5/test_set/output_torino/CM2/"
#1 dir
dir = args[1]
#2 matrix
m = args[2]
#3 swaps
s = args[3]

#Load input matrix
matrix = read.csv(paste(dir,m,sep=""), row.names = 1, header = T)
swaps = read.table(paste(dir,s,sep=""), sep = "\t", header = T, row.names = 1)

#from names to barcodes
barcodes_genes = read.csv(paste(dir,"rc_barcodes_genes.csv",sep=""), header = F, col.names = c("barcode", "name"))
swaps = left_join(swaps, barcodes_genes)
swaps$gene = sub("\\..*", "", swaps$name)
swaps = mutate(swaps, clone = paste(barcode, gene, UCI, sep = "_"))

#from actual_name to actual_barcode
names(barcodes_genes) = c("actual_barcode","actual_name")
barcodes_genes <- barcodes_genes[!duplicated(barcodes_genes$actual_name), ]
swaps = left_join(swaps, barcodes_genes)
swaps$actual_gene = sub("\\..*", "", swaps$actual_name)
swaps = mutate(swaps, actual_clone = paste(actual_barcode, actual_gene, UCI, sep = "_"))
association = swaps[, c("clone", "actual_clone")]
  
#FIND COMMON CLONES
names = as.data.frame(colnames(matrix))
names(names) = c("name")
names = separate(names, name, into = c("cellID", "UMIxUCI", "barcode", "gene", "UCI"), sep = "_")
names = mutate(names, clone = paste(barcode, gene, UCI, sep = "_"))
names = left_join(names, association, by = "clone")

print("Swaps found:")
print(length(intersect(names$clone, swaps$clone)))

#NEW NAMES complete
names$actual_clone <- ifelse(is.na(names$actual_clone), names$clone, names$actual_clone)
names = mutate(names, new_names = paste(cellID, UMIxUCI, actual_clone, sep = "_"))

#change in matrix
colnames(matrix) = names$new_names
write.csv(matrix, paste(dir,"/silencing_matrix_updated.csv",sep=""))

