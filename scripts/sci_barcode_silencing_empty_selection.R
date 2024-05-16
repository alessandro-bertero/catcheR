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

#Select empty reads
tr_complete_table = filter(complete_table, reference == "TACGCGTTCATCTGGGGGAGCCG")

#change names to match
RT = read.table(paste(dir, "/sci-RNA-seq-8.RT.oligos", sep=""), sep = "\t")
RT$number = as.numeric(rownames(RT))-1
RT = data.frame(RT$V2, RT$number)
names(RT) = c("RTb", "RT")

tr_complete_table <- tr_complete_table %>%
  mutate(demux = substr(cellID, 1, 3),
         RTb = substr(cellID, 4, 13))

tr_complete_table = merge(tr_complete_table, RT, by = "RTb")
tr_complete_table$cellID = paste(tr_complete_table$demux, "__RT_", tr_complete_table$RT, sep = "")
tr_complete_table = mutate(tr_complete_table, UMI = paste(cellID, UMI, sep = "-"))

#UMIxreference
tmp = distinct(tr_complete_table[, c("cellID", "UMI")])
UMIxref = data.frame(table(tmp[,"cellID"]))
UMIxref = filter(UMIxref, Freq > threshold)
tr_complete_table = filter(tr_complete_table, cellID %in% UMIxref$Var1)
#Save complete table
write.csv(tr_complete_table, paste(dir,"/complete_table_empty.csv",sep=""))

#Load from files obtained from previous analysis
table = distinct(read.table(paste(dir,"/cells_with_at_least_1_UCI.csv",sep=""), sep = ",", header = T, row.names = 1))
names(table) = c("cellID")

empty = subset(tr_complete_table, !(cellID %in% table$cellID))
print(length(unique(tr_complete_table$cellID)))
print(paste("Empty cells:", length(unique(empty$cellID)), sep = ""))
empty = empty %>% mutate(Name = paste(cellID, "?", "empty", "empty", "empty", sep = "_"))
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

