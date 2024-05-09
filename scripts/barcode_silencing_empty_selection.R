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

#Load from files obtained with bash script
with_cells_cellID = read.table(paste(dir, "/with_cells_cellID.txt", sep=""), sep = "\n")
with_cells_UMI = read.table(paste(dir,"/with_cells_UMI.txt", sep=""), sep = "\n")
with_cells_reference = read.table(paste(dir,"/with_cells_reference.txt", sep=""), sep = "\n")
with_cells_UCI = read.table(paste(dir,"/with_cells_UCI.txt", sep=""), sep = "\n")
with_cells_barcode = read.table(paste(dir,"/with_cells_barcode.txt", sep=""), sep = "\n")

#Create dataframe
complete_table = data.frame(with_cells_cellID, with_cells_UMI, with_cells_barcode, with_cells_UCI, with_cells_reference)
complete_table = distinct(complete_table)
names(complete_table) = c("cellID","UMI","barcode", "UCI", "reference")

#Load genes correspondence and add gene column
gene_ass = read.delim(paste(dir,"/rc_barcodes_genes.csv",sep=""), header = F, sep = ";", col.names = c("barcode","gene"))
complete_table = left_join(complete_table, gene_ass, by = "barcode")


#Select empty reads
tr_complete_table = filter(complete_table, reference == "TACGCGTTCATCTGGGGGAGCCG")


# #Distribution of genes
# gene_dist = data.frame(table(tr_complete_table$gene))
# #ggplot(barcodes_dist, aes(x= Var1)) + geom_histogram()
# jpeg("gene_distribution_empty.jpg", width = 350, height = 350)
# barplot(gene_dist$Freq, names.arg=gene_dist$Var1, las = 2, col="#69b3a2")
# dev.off()

#Save complete table
write.csv(tr_complete_table, paste(dir,"/complete_table_empty.csv",sep=""))

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

#write file with new names
write.csv(trans, paste(dir,"/silencing_matrix_empty.csv",sep=""))
bar = read.table(paste(dir,"/silencing_matrix.csv",sep=""), sep = ",", header = T, row.names = 1)
all = cbind(bar, trans)
write.csv(all, paste(dir,"/silencing_matrix_complete.csv",sep=""))

