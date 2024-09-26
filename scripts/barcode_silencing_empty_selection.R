#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)

#Set working direcotry
args = commandArgs(trailingOnly=TRUE)
dir = args[1]

#Determine arguments
#2 matrix
matrix = args[2]

#3 threhsold
threshold = args[3]
#4 sample 
sample = args[4]
#5 ref
ref = args[5]

dir = paste0(dir, "/Results_", sample)


#Load from files obtained with bash script
with_cells_cellID = read.table(paste(dir, "/with_cells_cellID_", sample, ".txt", sep=""), sep = "\n")
with_cells_UMI = read.table(paste(dir,"/with_cells_UMI_", sample, ".txt", sep=""), sep = "\n")
with_cells_reference = read.table(paste(dir,"/with_cells_reference_", sample, ".txt", sep=""), sep = "\n")
with_cells_UCI = read.table(paste(dir,"/with_cells_UCI_", sample, ".txt", sep=""), sep = "\n")
with_cells_barcode = read.table(paste(dir,"/with_cells_barcode_", sample, ".txt", sep=""), sep = "\n")

#Create dataframe
complete_table = data.frame(with_cells_cellID, with_cells_UMI, with_cells_barcode, with_cells_UCI, with_cells_reference)
complete_table = distinct(complete_table)
names(complete_table) = c("cellID","UMI","barcode", "UCI", "reference")

#Load genes correspondence and add gene column
gene_ass = read.delim(paste(dir,"/rc_barcodes_genes.csv",sep=""), header = F, sep = ";", col.names = c("barcode","gene"))
complete_table = left_join(complete_table, gene_ass, by = "barcode")


#Select empty reads
tr_complete_table = filter(complete_table, reference == ref)

#UMIxreference
tmp = distinct(tr_complete_table[, c("cellID", "UMI")])
UMIxref = data.frame(table(tmp[,"cellID"]))
UMIxref = filter(UMIxref, Freq > threshold)
tr_complete_table = filter(tr_complete_table, cellID %in% UMIxref$Var1)

#Save complete table
write.csv(tr_complete_table, paste(dir,"/complete_table_empty_", sample, ".csv",sep=""))

#Load from files obtained from previous analysis
table = distinct(read.table(paste(dir,"/cells_with_at_least_1_UCI.csv",sep=""), sep = ",", header = T, row.names = 1))
names(table) = c("cellID")

empty = subset(tr_complete_table, !(cellID %in% table$cellID))
print(length(unique(tr_complete_table$cellID)))
print(paste("Empty cells:", length(unique(empty$cellID)), sep = ""))
empty = empty %>% mutate(Name = paste(cellID, "?", "empty", gene, "empty", sep = "_"))
empty = distinct(empty[,c("cellID", "Name")])
#single_UCI_table = single_UCI_table %>% unite(Name, c(cellID, UMIxUCI, barcode, gene, UCI), sep = "_", remove = FALSE)

#load transcriptome matrix
trans = read.table(paste(dir,matrix,sep="/"), sep = ",", header = T, row.names = 1)
names(trans) = sub("\\..*", "", names(trans))
trans = as.data.frame(t(trans))
tmp = trans %>% transmute(cellID = rownames(trans))
trans = data.frame(trans, tmp)

#change name
trans = inner_join(as.data.frame(trans),as.data.frame(empty), by = "cellID")
trans$cellID = NULL
trans$Name.y = NULL
rownames(trans) = trans$Name
trans$Name = NULL
trans = as.data.frame(t(trans))
head(trans)
#write file with new names
write.csv(trans, paste(dir,"/silencing_matrix_empty_", sample, ".csv",sep=""))
saveRDS(trans, file = paste(dir,"/silencing_matrix_empty_", sample, ".rds",sep=""))
bar = read.table(paste(dir,"/silencing_matrix_", sample, ".csv",sep=""), sep = ",", header = T, row.names = 1)
all = cbind(bar, trans)
write.csv(all, paste(dir,"/silencing_matrix_complete_", sample, ".csv",sep=""))
saveRDS(all, file = paste(dir,"/silencing_matrix_complete_", sample, ".rds",sep=""))

