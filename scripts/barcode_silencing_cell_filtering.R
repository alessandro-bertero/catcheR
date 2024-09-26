#!/usr/bin/env Rscript
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidyverse)
#library(see)
args = commandArgs(trailingOnly=TRUE)

#Set working direcotry
#wd = "/2tb/torino/ratto/bertero/H001AS5/test_set/output_torino/CM2/"
#1 dir
dir = args[1]
#2 percentage
percentage = as.numeric(args[2])
#3 UMI count
UMI_count = as.numeric(args[3])
#4 gene exp matrix
matrix = args[4]
#5 sample
sample = args[5]
dir = paste0(dir, "/Results_", sample)

#Load from files obtained from explorative analysis
table = read.table(paste(dir,"/complete_table_fin_", sample, ".csv",sep=""), sep = ",", header = T)
table[,1] = NULL

#Define thresholds
table = mutate(table, TorF = case_when(UMIpercentagexUCI > percentage & UMIxUCI > UMI_count ~ TRUE, UMIpercentagexUCI <= percentage | UMIxUCI <= UMI_count ~ F))


#plot 2D with thresholds
tmp = distinct(table[c("cellID", "UCI","UMIpercentagexUCI", "UMIxUCI")])
p = ggplot(tmp, aes(x = UMIpercentagexUCI, y = UMIxUCI)) + 
  geom_count(alpha = 0.2) + 
  scale_y_continuous(trans='log10') + 
  geom_vline(xintercept=percentage, colour = "red") + 
  geom_hline(yintercept=UMI_count, colour = "red") + 
  theme_minimal()

ggsave(filename = paste0(dir, "/2D_percentage_of_UMIxUCI_UMI_count_thresholds.jpg"),
       plot = p,
       width = 10, height = 6)

pdf(paste(dir,"/2D_percentage_of_UMIxUCI_UMI_count_thresholds.pdf",sep=""), width = 10, height = 6)
p
dev.off()

#Select only TRUE UCI
true_table = filter(table, TorF == "TRUE")

#cells with at least 1 UCI
write.csv(true_table[, "cellID"], paste(dir,"/cells_with_at_least_1_UCI.csv",sep=""))

#True UCI x cell
true_UCIxcellID = distinct(true_table[,c("cellID", "UCI")])
true_UCIxcellID_table = data.frame(table(true_UCIxcellID$cellID))
names(true_UCIxcellID_table) = c("cellID", "True_UCIxcellID")
write.csv(true_UCIxcellID_table, paste0(dir, "/table_true_UCIxcell.csv"))

#plot
p = ggplot(true_UCIxcellID_table, aes(x=True_UCIxcellID)) +geom_bar()+ 
  theme_minimal()
ggsave(filename = paste0(dir, "/true_UCIxcell.jpg"),
       plot = p,
       width = 5, height = 3)

pdf(paste(dir,"/true_UCIxcell.pdf",sep=""))
#hist(true_UCIxcellID_table$True_UCIxcellID, breaks = 5)
p
dev.off()

#Add column with true UCI x cell
true_table = left_join(true_UCIxcellID_table, true_table, by = "cellID")
table2 = merge(table, true_UCIxcellID_table, all = T)
single_UCI_table = data.frame(filter(true_table, True_UCIxcellID == 1))

#plot TRUE UCIs
table2[is.na(table2)] <- 0
tmp = distinct(table2[c("cellID", "UCI","UMIpercentagexUCI", "UMIxUCI", "True_UCIxcellID")])
p = ggplot(tmp, aes(x = UMIpercentagexUCI, y = UMIxUCI, colour = as.factor(True_UCIxcellID))) + 
  geom_count(alpha = 0.3) + 
  #ylim(0,500) + 
  scale_y_continuous(trans='log10') + 
  geom_vline(xintercept=percentage, colour = "red") + 
  geom_hline(yintercept=UMI_count, colour = "red") + 
  scale_color_manual(values = c("#C0C0C0", "#E69F00", "#0072B2","#0072B2","#0072B2","#0072B2", "#0072B2", "#0072B2", "#0072B2"))+
  theme_minimal()

ggsave(filename = paste0(dir, "/2D_percentage_of_UMIxUCI_UMI_count_trueorfalse.jpg"),
       plot = p,
       width = 10, height = 6)

pdf(paste(dir,"/2D_percentage_of_UMIxUCI_UMI_count_trueorfalse.pdf",sep=""), width = 10, height = 6)
p
dev.off()

#select valid cells = only 1 true UCI + top UCI UMIs / second top UCI UMIs > 5
tmp2 = left_join(true_UCIxcellID_table, table, by = "cellID")

ValidCell_df <- tmp2 %>%
  group_by(cellID) %>%
  mutate(
    UMI = NULL,
    True_UCIxcellID = ifelse(True_UCIxcellID == 1 & max(UMIxUCI) / max(UMIxUCI[UMIxUCI != max(UMIxUCI)]) >= 5 | 
                               (True_UCIxcellID == 1 & UMIpercentagexUCI == 100), TRUE, FALSE),
    ValidCell = True_UCIxcellID
  ) %>%
  distinct(cellID, ValidCell) %>%
  ungroup() #%>%
  #select(cellID, ValidCell)

true_table = merge(true_table, ValidCell_df)
table2 = merge(table2, ValidCell_df, all = T)

#cells x clone
cellsxUCI = distinct(true_table[,c("cellID", "UCI")])
cellsxUCI_table = data.frame(table(cellsxUCI$UCI))
names(cellsxUCI_table) = c("UCI", "cellsxUCI")


#plot
p = ggplot(cellsxUCI_table, aes(x=cellsxUCI)) +geom_bar() + theme_minimal()
ggsave(filename = paste0(dir,"/cellsxclone.jpg"),
       plot = p,
       width = 5, height = 3)

pdf(paste(dir,"/cellsxclone.pdf",sep=""), width = 9)
#hist(true_UCIxcellID_table$True_UCIxcellID, breaks = 5)
p
dev.off()

p = ggplot(cellsxUCI_table, aes(x=log10(cellsxUCI))) +geom_bar()+theme_minimal()
ggsave(filename = paste0(dir,"/cellsxclone_log10.jpg"),
       plot = p,
       width = 5, height = 3)

pdf(paste(dir,"/cellsxclone_log10.pdf",sep=""), width = 9)
p
dev.off()

#plot Valid cells
tmp = distinct(table2[c("cellID", "UCI","UMIpercentagexUCI", "UMIxUCI", "ValidCell")])
p = ggplot(tmp, aes(x = UMIpercentagexUCI, y = UMIxUCI, colour = as.factor(ValidCell))) + 
  geom_count(alpha = 0.3) + ylim(0,500) + 
  scale_y_continuous(trans='log10') + 
  geom_vline(xintercept=percentage, colour = "red") + 
  geom_hline(yintercept=UMI_count, colour = "red") + 
  scale_color_manual(values = c("#0072B2","#E69F00"))+
  theme_minimal()

ggsave(filename = paste0(dir,"/2D_percentage_of_UMIxUCI_UMI_count_ValidCells.jpg"),
        plot = p,
        width = 10, height = 6)

pdf(paste(dir,"/2D_percentage_of_UMIxUCI_UMI_count_ValidCells.pdf",sep=""), width = 10, height = 6)
p
dev.off()

#Combine UMIxUCI and gene with cellID for following sc analysis
single_UCI_table_old = single_UCI_table
single_UCI_table = filter(true_table, ValidCell == TRUE)
single_UCI_table = single_UCI_table %>% unite(Name, c(cellID, UMIxUCI, barcode, gene, UCI), sep = "_", remove = FALSE)

new_names = distinct(single_UCI_table[,1:2])
write.csv(new_names, paste(dir,"/new_cell_names.csv",sep=""))

#load transcriptome matrix
#trans = read.table(paste(dir,"y12.csv",sep="/"), sep = ",", header = T, row.names = 1)
trans = read.table(paste(dir,matrix,sep="/"), sep = ",", header = T, row.names = 1)
names(trans) = sub("\\..*", "", names(trans))
trans = as.data.frame(t(trans))
tmp3 = trans %>% transmute(cellID = rownames(trans))
trans = data.frame(trans, tmp3)

#change name
trans = inner_join(as.data.frame(trans),as.data.frame(new_names), by = "cellID")
trans$cellID = NULL
trans$Name.y = NULL
rownames(trans) = trans$Name
trans$Name = NULL
trans = as.data.frame(t(trans))
head(trans)
#write file with new names
print(paste(dir,"/silencing_matrix_", sample, ".csv",sep=""))
write.csv(trans, paste(dir,"/silencing_matrix_", sample, ".csv",sep=""))

fileConn<-file(paste(dir,"/log_part3.txt",sep=""))
writeLines(c(paste("Cells with at least 1 true UCI:", length(unique(true_table$cellID)), sep = ""), 
             paste("Cells with only 1 true UCI:", length(unique(single_UCI_table_old$cellID)), sep = ""),
             paste("Valid cells after filtering (UMIs of top UCI / UMIs of 2nd top UCI > 5):", length(unique(single_UCI_table$cellID)), sep = "")
             ), fileConn)
