#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#1 working dir 
dir = args[1]
print(dir)
#2 reference GGCGCGTTCATCTGGGGGAGCCG
ref = args[2]
print(ref)
#3 threshold mode
mode = args[3]
#4 sample
sample = args[4]
dir = paste0(dir, "/Results_", sample)
library(dplyr)
library(tidyr)
library(ggplot2)
library(quantmod)


#Redirect output to log2.txt
sink(file = paste(dir,"/log2.txt",sep = ""), append = T, type = "output", split = T)


#Load from files obtained with bash script
with_cells_cellID = read.table(paste(dir, "/with_cells_cellID_", sample, ".txt", sep=""), sep = "\n")
with_cells_UMI = read.table(paste(dir,"/with_cells_UMI_", sample, ".txt", sep=""), sep = "\n")
with_cells_reference = read.table(paste(dir,"/with_cells_reference_", sample, ".txt", sep=""), sep = "\n")
with_cells_UCI = read.table(paste(dir,"/with_cells_UCI_", sample, ".txt", sep=""), sep = "\n")
with_cells_barcode = read.table(paste(dir,"/with_cells_barcode_", sample, ".txt", sep=""), sep = "\n")
with_cells_reads_counts = read.table(paste(dir,"/with_cells_reads_counts_", sample, ".txt", sep=""), sep = "\n")

#Create dataframe
complete_table = data.frame(with_cells_cellID, with_cells_UMI, with_cells_barcode, with_cells_UCI, with_cells_reference, with_cells_reads_counts)
names(complete_table) = c("cellID","UMI","barcode", "UCI", "reference", "reads")
complete_table = mutate(complete_table, UMI = paste(cellID, UMI, sep = "-"))
complete_table = distinct(complete_table, across(-reads), .keep_all = TRUE)

#Load genes correspondence and add gene column
gene_ass = read.delim(paste(dir,"/rc_barcodes_genes.csv",sep=""), header = F, sep = ",", col.names = c("barcode","shRNA"))
complete_table = left_join(complete_table, gene_ass, by = "barcode")
remove_after_dot <- function(string) {
  gsub("\\.[^.]*$", "", string)
}
complete_table = mutate(complete_table, gene = remove_after_dot(shRNA))


#Select reads with reference
tr_complete_table = filter(complete_table, reference == ref)

print(paste("Cells:", length(unique(complete_table$cellID))))
print(paste("Transfected cells:", length(unique(tr_complete_table$cellID))))
print(paste("Transfected reads:", length(tr_complete_table$cellID)))

#Keep only known barcodes
tr_complete_table <- tr_complete_table[!(is.na(tr_complete_table$gene)), ]
tr_complete_table = distinct(tr_complete_table, across(-reads), .keep_all = TRUE)

print(paste("Transfected cells with known barcode:", length(unique(tr_complete_table$cellID))))
print(paste("Transfected reads with known barcode:", length(tr_complete_table$cellID)))
print(paste("Total UMIs:", length(tr_complete_table$UMI)))
print(paste("Unique UMIs of cells with known barcode:", length(unique(tr_complete_table$UMI))))
print(paste("Duplicated UMIs:", sum(duplicated(tr_complete_table$UMI))))

#remove duplicated UMIs
tr_complete_table = tr_complete_table[order(tr_complete_table[,'UMI'],-tr_complete_table[,'reads']),]
tr_complete_table = tr_complete_table[!duplicated(tr_complete_table$UMI),]
print(paste("Duplicated UMIs after filtering:", sum(duplicated(tr_complete_table$UMI))))

#Distribution of barcodes
barcodes_dist = data.frame(table(tr_complete_table$shRNA))
names(barcodes_dist) = c("shRNA", "UMIs")
p = ggplot(barcodes_dist, aes(x = shRNA, y = UMIs, fill = shRNA)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = NULL, y = "UMIs") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

if (file.exists(paste(dir, "/colors.csv", sep = ""))) {
  colors <- read.csv(paste(dir, "/colors.csv", sep = ""), header = F, col.names = c("shRNA", "colors"))
  color_mapping <- setNames(colors$colors, colors$shRNA)
  p = p+
    scale_fill_manual(values = color_mapping)

} else {
  print("Colors file does not exist.")
}

ggsave(filename = paste0(dir, "/barcode_distribution.jpg"),
       plot = p,
       width = 7, height = 9)

pdf(paste(dir,"/barcode_distribution.pdf", sep=""))
p
dev.off()

#Distribution of genes
gene_dist = data.frame(table(tr_complete_table$gene))
names(gene_dist) = c("Gene", "UMIs")
p = ggplot(gene_dist, aes(x = Gene, y = UMIs, fill = Gene)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  labs(x = NULL, y = "UMIs") + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(filename = paste0(dir, "/gene_distribution.jpg"),
       plot = p,
       width = 7, height = 9)

pdf(paste(dir,"/gene_distribution.pdf",sep=""))
p
dev.off()

#UMIxUCI
UMIxUCI = tr_complete_table[,c("cellID", "UCI", "UMI")]
UMIxUCI$cellID_UCI_tally <- paste(UMIxUCI$cellID, UMIxUCI$UCI, sep = "_")
UMIxUCI = UMIxUCI[,c("cellID_UCI_tally", "UMI")]
UMIxUCI = distinct(UMIxUCI)
UMIxUCI = data.frame(table(UMIxUCI[,"cellID_UCI_tally"]))
print("UMIxUCI")
U = length(UMIxUCI$Var1)
one = length(filter(UMIxUCI, Freq ==1)$Var1)
perc_1 = one / U *100
print(paste("Total UCIs:",U))
print(paste("UCIs with 1 UMI:", one))
print(paste("Percentage of UCIs with 1 UMI:",perc_1))
p = perc_1 / 100 * 1.35
print(paste("Percentage of UCIs to filter out:",p))
UMIxUCI <- UMIxUCI[order(-UMIxUCI$Freq), ]
n_rows_percent <- ceiling(nrow(UMIxUCI) * (1-p))
UMIxUCI_f <- UMIxUCI[1:n_rows_percent, ]
threshold_noise <- min(UMIxUCI_f$Freq)
print(paste("Suggested threshold of UMIsxUCI based on noise:",threshold_noise))

hist = as.data.frame(table(UMIxUCI$Freq))
threshold_bimodal = findValleys(hist$Freq, thresh=5)[1]
print(paste("Suggested threshold of UMIsxUCI based on bimodal distribution:",threshold_bimodal))

if (mode == "noise") {
  threshold = threshold_noise
  print("Mode: noise")
} else {
  threshold = threshold_bimodal
  print("Mode: bimodal")
}
write(threshold, paste(dir, "/UMI_threshold.txt", sep = ""))

#plot
p = ggplot(UMIxUCI, aes(x = Freq)) +
  geom_histogram(bins = 2000, alpha = 0.7) +
  labs(x = "UMIsxUCI")+
  geom_vline(xintercept = threshold, col = "red")+
  theme_minimal()

ggsave(filename = paste0(dir, "/UMIxUCI.jpg"),
       plot = p,
       width = 7, height = 9)

pdf(paste(dir,"/UMIxUCI.pdf",sep=""))
p
dev.off()

#plot zoom 
p <- p + coord_cartesian(xlim = c(0,100) ,ylim = c(0, 400)) 

ggsave(filename = paste0(dir, "/UMIxUCI_400_100.jpg"),
       plot = p,
       width = 7, height = 9)
pdf(paste(dir,"/UMIxUCI_400_100.pdf",sep=""))
p
dev.off()

#UCIxcell
UCIxcellID = tr_complete_table[,c("cellID","UCI")]
UCIxcellID = distinct(UCIxcellID)
UCIxcellID_table = as.data.frame(table(UCIxcellID[,"cellID"]))
#plot
p = ggplot(UCIxcellID_table, aes(x = Freq)) +
  geom_histogram(bins = 70) + 
  labs(x = "UCIsxcell") +        
  theme_minimal()
ggsave(filename = paste0(dir, "/UCIxcell.jpg"),
       plot = p,
       width = 7, height = 9)
pdf(paste(dir,"/UCIxcell.pdf",sep=""))
p
dev.off()

#Add UMIxcellID
UMIxcellID = data.frame(table(distinct(tr_complete_table[,c("cellID","UMI")])[,"cellID"]))
names(UMIxcellID) = c("cellID", "UMIxcellID")
tr_complete_table = left_join(UMIxcellID, tr_complete_table, by = "cellID")

#Add UMIxUCI
UMIxUCI = separate(UMIxUCI, col = "Var1", into = c("cellID","UCI"), sep = "_")
names(UMIxUCI) = c("cellID","UCI", "UMIxUCI")
tr_complete_table = left_join(UMIxUCI, tr_complete_table, by = c("cellID", "UCI"))

#Add percentage of UMI on total cell UMI supporting each UCI
tr_complete_table_fin <- tr_complete_table %>% mutate(UMIpercentagexUCI = ((UMIxUCI*100)/UMIxcellID))

#plot percentage
percentage_table = distinct(tr_complete_table_fin[c("UCI", "UMIpercentagexUCI")])
p = ggplot(percentage_table, aes(x = UMIpercentagexUCI)) +
  geom_histogram(bins = 100) + 
  labs(x = "UMIpercentagexUCI") +
  geom_vline(xintercept = 15,col = "red")+
  theme_minimal()

ggsave(filename = paste0(dir, "/percentage_of_UMIxUCI_dist.jpg"),
       plot = p,
       width = 7, height = 9)
pdf(paste(dir,"/percentage_of_UMIxUCI_dist.pdf",sep=""))
p
dev.off()

#2D plot percentage and UMI count x UCI
tmp = distinct(tr_complete_table_fin[c("cellID", "UCI","UMIpercentagexUCI", "UMIxUCI")])
p = ggplot(tmp, aes(x = UMIpercentagexUCI, y = UMIxUCI)) + 
  geom_count(alpha = 0.3) + 
  scale_y_continuous(trans='log10') +
  theme_minimal()

ggsave(filename = paste0(dir, "/2D_percentage_of_UMIxUCI_UMI_count.jpg"),
       plot = p,
       width = 7, height = 9)
pdf(paste(dir,"/2D_percentage_of_UMIxUCI_UMI_count.pdf",sep=""),width = 10, height = 6)
p
dev.off()

#Save complete table
write.csv(tr_complete_table_fin, paste(dir,"/complete_table_fin_", sample, ".csv",sep=""))
system("chmod 777 /data/scratch")