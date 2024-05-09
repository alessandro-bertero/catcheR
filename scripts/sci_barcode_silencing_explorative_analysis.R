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
print(mode)
#4 matrix
matrix = args[4]
file.exists(paste(dir,matrix,sep="/"))

library(quantmod, warn.conflicts = F)
library(Biostrings, warn.conflicts = F)
library(stringr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(tidyr, warn.conflicts = F)
library(ggplot2, warn.conflicts = F)

#Redirect output to log2.txt
sink(file = paste(dir,"/log2.txt",sep = ""), append = T, type = "output", split = T)


#Load from files obtained with bash script
with_cells_cellID = read.table(paste(dir, "/with_cells_cellID.txt", sep=""), sep = "\n")
with_cells_UMI = read.table(paste(dir,"/with_cells_UMI.txt", sep=""), sep = "\n")
with_cells_reference = read.table(paste(dir,"/with_cells_reference.txt", sep=""), sep = "\n")
with_cells_UCI = read.table(paste(dir,"/with_cells_UCI.txt", sep=""), sep = "\n")
with_cells_barcode = read.table(paste(dir,"/with_cells_barcode.txt", sep=""), sep = "\n")
#with_cells_reads_counts = read.table(paste(dir,"/with_cells_reads_counts.txt", sep=""), sep = "\n")

#Create dataframe
complete_table = data.frame(with_cells_cellID, with_cells_UMI, with_cells_barcode, with_cells_UCI, with_cells_reference)
names(complete_table) = c("cellID","UMI","barcode", "UCI", "reference")
complete_table = mutate(complete_table, UMI = paste(cellID, UMI, sep = "-"))
complete_table = distinct(complete_table)

#Load genes correspondence and add gene column
gene_ass = read.delim(paste(dir,"/rc_barcodes_genes.csv",sep=""), header = F, sep = ",", col.names = c("barcode","shRNA"))
complete_table = left_join(complete_table, gene_ass, by = "barcode")
remove_after_dot <- function(string) {
  gsub("\\.[^.]*$", "", string)
}
complete_table = mutate(complete_table, gene = remove_after_dot(shRNA))
#Select reads with reference
tr_complete_table = filter(complete_table, reference == ref)
head(tr_complete_table)


# #UMI x all cell
# length(unique(complete_table$cellID))
# cells_table = distinct(complete_table[,c("cellID", "UMI")])
# cells_table = as.data.frame(table(cells_table$cellID))
# cells_table = filter(cells_table, Freq > 500)
# length(unique(cells_table$Var1))

#UMI x tr cell
# length(unique(tr_complete_table$cellID))
# tr_cells_table = filter(complete_table, cellID %in% tr_complete_table$cellID)
# tr_cells_table = distinct(tr_cells_table[,c("cellID", "UMI")])
# tr_cells_table = as.data.frame(table(tr_cells_table$cellID))
# names(tr_cells_table) = c("cellID", "UMIxcell")
# tr_complete_table = left_join(tr_complete_table, tr_cells_table)
# length(unique(tr_cells_table$cellID))
# #tr_cells_table = filter(tr_cells_table, Freq > 1)
# jpeg(paste(dir,"/UMIxcell.jpg",sep=""), width = 700, height = 900)
# hist(tr_cells_table$UMIxcell, breaks = 1000)
# dev.off()
# tr_cells_table = filter(tr_cells_table, UMIxcell > 500)
# length(unique(tr_cells_table$cellID))

#check same lenght
#char_lengths <- nchar(tr_complete_table$cellID[1])
#all(char_lengths == char_lengths[1])
#get actual cell names and filter for sequences cells
RT = read.table(paste(dir, "/sci-RNA-seq-8.RT.oligos", sep=""), sep = "\t")
RT$number = as.numeric(rownames(RT))-1
RT = data.frame(RT$V2, RT$number)
names(RT) = c("RTb", "RT")

tr_complete_table <- tr_complete_table %>%
  mutate(demux = substr(cellID, 1, 3),
         RTb = substr(cellID, 4, 13))

tr_complete_table = merge(tr_complete_table, RT, by = "RTb")
tr_complete_table$cellID = paste(tr_complete_table$demux, "__RT_", tr_complete_table$RT, sep = "")
#tr_complete_table = filter(tr_complete_table, UMIxcell > 500)
trans_names = unique(tr_complete_table$cellID)
length(trans_names)
write.csv(trans_names, file = paste(dir, "trans_names.csv", sep = ""))

#Distribution of cells in rows and cols
#demux reads
data = as.data.frame(tr_complete_table$demux)
names(data) = "demux"
data <- data %>%
  mutate(letter = str_match(demux, "([A-Za-z]+)")[, 2],
         number = str_match(demux, "(\\d+)")[, 2])
letter_dist = data.frame(table(data$letter))
p=ggplot(letter_dist,aes(x= Var1,y = Freq))+geom_col()+
  labs(x = "Demultiplexing letter", y = "Unique reads")
ggsave(paste(dir,"/demux_letter_distribution.png", sep=""),p, device = "png")
ggsave(paste(dir,"/demux_letter_distribution.pdf", sep=""),p, device = "pdf")

number_dist = data.frame(table(data$number))
p=ggplot(number_dist,aes(x= Var1,y = Freq))+geom_col()+
  labs(x = "Demultiplexing number", y = "unique reads")
ggsave(paste(dir,"/demux_number_distribution.png", sep=""),p, device = "png")
ggsave(paste(dir,"/demux_number_distribution.pdf", sep=""),p, device = "pdf")

#demux cells
data = distinct(data.frame(tr_complete_table$demux, tr_complete_table$cellID))
data = as.data.frame(data$tr_complete_table.demux)
names(data) = "demux"
data <- data %>%
  mutate(letter = str_match(demux, "([A-Za-z]+)")[, 2],
         number = str_match(demux, "(\\d+)")[, 2])
letter_dist = data.frame(table(data$letter))
p=ggplot(letter_dist,aes(x= Var1,y = Freq))+geom_col()+
  labs(x = "Demultiplexing letter", y = "Cells")
ggsave(paste(dir,"/demux_letter_cell_distribution.png", sep=""),p, device = "png")
ggsave(paste(dir,"/demux_letter_cell_distribution.pdf", sep=""),p, device = "pdf")

number_dist = data.frame(table(data$number))
p=ggplot(number_dist,aes(x= Var1,y = Freq))+geom_col()+
  labs(x = "Demultiplexing number", y = "Cells")
ggsave(paste(dir,"/demux_number_cells_distribution.png", sep=""),p, device = "png")
ggsave(paste(dir,"/demux_number_cells_distribution.pdf", sep=""),p, device = "pdf")

#RT reads
RT_dist = data.frame(table(tr_complete_table$RT))
RT_dist$Var1 = as.numeric(RT_dist$Var1)
p=ggplot(RT_dist,aes(x= Var1,y = Freq))+geom_col()+
  labs(x = "RT number", y = "Unique reads")+
  geom_text(aes(label = signif(Var1)), nudge_y = 35)
ggsave(paste(dir,"/RT_distribution.png", sep=""),p, device = "png", width = 20, height = 10)
ggsave(paste(dir,"/RT_distribution.pdf", sep=""),p, device = "pdf", width = 20, height = 10)

#RT cells
RT_dist = distinct(data.frame(tr_complete_table$cellID, tr_complete_table$RT))
RT_dist = data.frame(table(RT_dist$tr_complete_table.RT))
RT_dist$Var1 = as.numeric(RT_dist$Var1)
p=ggplot(RT_dist,aes(x= Var1,y = Freq))+geom_col()+
  labs(x = "RT number", y = "Cells")+
  geom_text(aes(label = signif(Var1)), nudge_y = 2)
ggsave(paste(dir,"/RT_cells_distribution.png", sep=""),p, device = "png", width = 20, height = 10)
ggsave(paste(dir,"/RT_cells_distribution.pdf", sep=""),p, device = "pdf", width = 20, height = 10)

tr_complete_table$RTb = NULL
tr_complete_table$demux = NULL
tr_complete_table$RT = NULL
#tr = tr_complete_table

tr_complete_table$UMI <- sapply(strsplit(tr_complete_table$UMI, "-"), function(x) x[2])
tr_complete_table$UMI = paste(tr_complete_table$cellID, tr_complete_table$UMI, sep = "-")
#filter cell names
trans = read.table(paste(dir,matrix,sep="/"), sep = ",", header = T, row.names = 1)
names = colnames(trans)
#write.csv(names, file = paste(dir, "names.csv", sep = ""))

#names distribution
data = as.data.frame(names)
data <- data %>%
  mutate(letter = str_match(names, "([A-Za-z]+)")[, 2],
         number = str_match(names, "(\\d+)")[, 2],
         RT = str_match(names, ".*_(.*)")[, 2])

letter_dist = data.frame(table(data$letter))
p=ggplot(letter_dist,aes(x= Var1,y = Freq))+geom_col()+
  labs(x = "Demultiplexing letter", y = "Cells")
ggsave(paste(dir,"/demux_letter_cellnames_distribution.png", sep=""),p, device = "png")
ggsave(paste(dir,"/demux_letter_cellnames_distribution.pdf", sep=""),p, device = "pdf")

number_dist = data.frame(table(data$number))
p=ggplot(number_dist,aes(x= Var1,y = Freq))+geom_col()+
  labs(x = "Demultiplexing number", y = "Cells")
ggsave(paste(dir,"/demux_number_cellsnames_distribution.png", sep=""),p, device = "png")

RT_dist = data.frame(table(data$RT))
RT_dist$Var1 = as.numeric(RT_dist$Var1)
p=ggplot(RT_dist,aes(x= Var1,y = Freq))+geom_col()+
  labs(x = "RT number", y = "Cells")+
  geom_text(aes(label = signif(Var1)), nudge_y = 2)
ggsave(paste(dir,"/RT_cellsnames_distribution.png", sep=""),p, device = "png", width = 20, height = 10)
ggsave(paste(dir,"/RT_cellsnames_distribution.pdf", sep=""),p, device = "pdf", width = 20, height = 10)

#filtering for sequenced cells 
# Encoding(as.character(tr_complete_table$cellID))
# Encoding(names)
# tr_complete_table$cellID = lapply(FUN = iconv, X = tr_complete_table$cellID, to = "ASCII")
# names = lapply(FUN = iconv, names, to = "ASCII")
# length(intersect(tr_complete_table$cellID, names))
tr_complete_table= filter(tr_complete_table, cellID %in% names)
tr_complete_table = distinct(tr_complete_table)
#stats
print(paste("Total cell IDs (also not transfected or sequenced):", length(unique(complete_table$cellID))))
print(paste("Transfected cells:", length(unique(tr_complete_table$cellID))))
print(paste("Transfected reads:", length(tr_complete_table$cellID)))

#Keep only known barcodes
tr_complete_table <- tr_complete_table[!(is.na(tr_complete_table$gene)), ]
tr_complete_table = distinct(tr_complete_table)

print(paste("Transfected cells with known barcode:", length(unique(tr_complete_table$cellID))))
print(paste("Transfected reads with known barcode:", length(tr_complete_table$cellID)))
print(paste("Total UMIs:", length(tr_complete_table$UMI)))
print(paste("Unique UMIs of cells with known barcode:", length(unique(tr_complete_table$UMI))))
print(paste("Duplicated UMIs:", sum(duplicated(tr_complete_table$UMI))))

#remove duplicated UMIs
tr_complete_table = tr_complete_table[!duplicated(tr_complete_table$UMI),]
print(paste("Duplicated UMIs after filtering:", sum(duplicated(tr_complete_table$UMI))))
#tr = tr_complete_table

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
UMIxUCI$cellID_UCI_tally <- paste(UMIxUCI$cellID, UMIxUCI$UCI, sep = "/")
UMIxUCI = UMIxUCI[,c("cellID_UCI_tally", "UMI")]
#UMIxUCI = distinct(UMIxUCI)
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
UMIxUCI = separate(UMIxUCI, col = "Var1", into = c("cellID","UCI"), sep = "/")
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
write.csv(tr_complete_table_fin, paste(dir,"/complete_table_fin.csv",sep=""))

# length(unique(tr_complete_table$cellID))
# length(intersect(names, tr_complete_table$cellID))
# not_common = as.data.frame(c(setdiff(names, tr_complete_table_fin$cellID), 
#                            setdiff(tr_complete_table_fin$cellID,names)))

