#!/usr/bin/env Rscript
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

#1 dir
dir = args[1]
#2 threshold
#threshold = args[2]
#3 clones
clones = args[2]
if (clones == "NULL" | (clones == "")) {
  clones = NULL
}
print(clones)

#fastqc(group = "docker", data.folder = paste(dir, "../fastqc", sep = ""))

if (file.exists(paste(dir, "/colors.csv", sep = ""))) {
  colors <- read.csv(paste(dir, "/colors.csv", sep = ""), header = F, col.names = c("shRNA", "colors"))
  color_mapping <- setNames(colors$colors, colors$shRNA)
} else {
  print("Colors file does not exist.")
}

#Load from files obtained with bash script
UMI = fread(paste(dir,"/inter_UMI.txt", sep=""), sep = "\n", data.table = F)
UCI = fread(paste(dir,"/inter_UCI.txt", sep=""), sep = "\n", data.table = F)
barcode = fread(paste(dir,"/inter_BC.txt", sep=""), sep = "\n", data.table = F)
shRNA = fread(paste(dir,"/inter_shRNA.txt", sep=""), sep="\n", data.table=F)

#LOAD BARCODE NAMES
barcodes_names = distinct(read.table(paste(dir,"/rc_barcodes_genes.csv", sep=""), sep = ","))
names(barcodes_names) = c("barcode", "name")
print(head(barcodes_names))

#LOAD shRNA
expected = read.table(paste(dir,"/expected_shRNA_names.txt", sep=""), sep = " ")
names(expected) = c("actual_name", "shRNA")
print(head(expected))

#Create dataframe
complete_table = data.frame(UMI, barcode, UCI, shRNA)
names(complete_table) = c("UMI","barcode", "UCI", "shRNA")
complete_table = distinct(complete_table)
complete_table = mutate(complete_table, clone = paste(barcode, UCI, sep = "_"))

#distribution of clone
# clone_t = as.data.frame(table(complete_table[,c("clone")]))
# print("Total clones:")
# length(clone_t$Var1)
# print("Clones with only 1 read:")
# length(filter(clone_t, Freq == 1)$Var1)
# clone_t = separate(clone_t, Var1, into = c("barcode", "UCI"), sep = "_", remove = F)
# #add barcode name
# clone_t = left_join(clone_t, barcodes_names)
# p = ggplot(filter(clone_t, Freq > 1)) +
#   geom_col(aes(x = fct_reorder(Var1, Freq, .desc = T), y = Freq, fill = name))  +
#   theme_minimal() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# #labs(title = data$name) +
# #scale_fill_manual(values = c("#00AFBB", "#E7B800"), drop = F)
# pdf(paste(dir,"/clone_distribution_all.pdf",sep=""), width = 10, height = 6)
# p
# dev.off()
# ggsave(paste(dir,"/clone_distribution_all.jpg",sep=""),
#        plot = p,
#        width = 10, height = 15,
#        limitsize = F)
# ggsave(paste(dir,"/clone_distribution_all.pdf",sep=""),
#        plot = p,
#        width = 10, height = 15,
#        limitsize = F)
# 
# clone_t = filter(clone_t, Freq > threshold)
# p = ggplot(clone_t) +
#   geom_col(aes(x = fct_reorder(Var1, Freq, .desc = T), y = Freq, fill = name))  +
#   theme_minimal() +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
# 
# if (exists("color_mapping")) {
#   p = p+
#     scale_fill_manual(values = color_mapping)
# }
# #labs(title = data$name) +
# #scale_fill_manual(values = c("#00AFBB", "#E7B800"), drop = F)
# ggsave(paste(dir,"/clone_distribution_filter_", threshold, ".jpg",sep=""),
#        plot = p,
#        width = 10, height = 15,
#        limitsize = F)
# ggsave(paste(dir,"/clone_distribution_filter_", threshold, ".pdf",sep=""),
#        plot = p,
#        width = 10, height = 15,
#        limitsize = F)

#ADD NAME AND shRNA
complete_table = left_join(complete_table, barcodes_names, by = "barcode")
complete_table = separate(complete_table, name, into = c("gene", "trash"), remove = F)
complete_table$trash = NULL
complete_table = complete_table %>%
  mutate(shRNA = str_sub(shRNA, 5, -1))
expected = expected %>%
  mutate(shRNA = str_sub(shRNA, 1, 21))
complete_table = filter(complete_table, shRNA %in% expected$shRNA)
complete_table = left_join(complete_table, expected, by = "shRNA", relationship = "many-to-many")
complete_table = na.omit(complete_table)
complete_table$same = with(complete_table, ifelse(as.character(name) == as.character(actual_name), 'yes',
                                                  ifelse(as.character(name) != as.character(actual_name), 'no', NA)))
table = distinct(as.data.frame(table(complete_table[,c("UCI", "actual_name", "same", "name")])))

#complete_table = separate(complete_table, col = "BC", into = c("gene"), sep = "_", remove = F)

# distribution of swaps 
dir.create(paste(dir, "/shRNA_distribution", sep = ""))
for (seq in expected$actual_name) {
  data = filter(table, actual_name == seq)
  p = ggplot(data) + 
    geom_col(aes(x = name, y = Freq, fill = same))  +
    theme_minimal() + 
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
  ggsave(paste(dir,"/shRNA_distribution/", seq, ".jpg",sep=""), plot = p)
  ggsave(paste(dir,"/shRNA_distribution/", seq, ".pdf",sep=""), plot = p)
}

#distribution 
gene_dis = as.data.frame(table(distinct(complete_table[, c("UMI", "gene")])[, "gene"]))
names(gene_dis) = c("gene", "UMI")
p = ggplot(gene_dis) +
  geom_col(aes(x = gene, y = UMI, fill = gene)) +
  theme_minimal() + 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
ggsave(paste(dir,"/gene_distribution.jpg",sep=""), plot = p)
ggsave(paste(dir,"/gene_distribution.pdf",sep=""), plot = p)

name_dis = as.data.frame(table(distinct(complete_table[, c("UMI", "name")])[, "name"]))
names(name_dis) = c("name", "UMI")
#name_dis  = left_join(name_dis, complete_table[, c("name", "gene")])
p = ggplot(name_dis) +
  geom_col(aes(x = name, y = UMI, fill = name)) +
  theme_minimal()+ 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
if (exists("color_mapping")) {
  p = p+
    scale_fill_manual(values = color_mapping)
}
ggsave(paste(dir,"/name_distribution.jpg",sep=""), plot = p)
ggsave(paste(dir,"/name_distribution.pdf",sep=""), plot = p)


#total swaps
table = filter(table, Freq > 3)
table = mutate(table, clone = paste(UCI, name, sep = "_"))
print("Clones with at least 3 reads:")
length(unique(table$clone))

good_c = c()
res_table = data.frame()
for (c in unique(table$clone)) {
  df = filter(table, clone == c)
  if (max(df$Freq) > 30000) {
    #print(max(df$Freq))
    freqs = sort(df$Freq, decreasing = TRUE)
    #print(freqs[2])
    if (length(df$Freq) == 1) {
      good_c = append(good_c, c)
      res_table = rbind(res_table, filter(df, Freq == freqs[1]))
    }
    else if (freqs[1] / freqs[2] > 2) {
      good_c = append(good_c, c)
      res_table = rbind(res_table, filter(df, Freq == freqs[1]))
    }
  }
}

good_table = filter(table, clone %in% good_c)

swaps = filter(res_table, same == "no")
print(paste0("Swaps: ", length(swaps$clone)))
correct = filter(res_table, same == "yes")
print(paste0("Correct: ", length(correct$clone)))

write.table(good_table, paste(dir, "/reliable_clones_30000_2.txt", sep = ""), sep = "\t")
write.table(swaps, paste(dir, "/reliable_clones_swaps_30000_2.txt", sep = ""), sep = "\t")
write.table(correct, paste(dir, "/reliable_clones_confirmations_30000_2.txt", sep = ""), sep = "\t")

#correct dis
c_dis = as.data.frame(table(correct[, c("name")]))
names(c_dis) = c("name", "correct")
#swap dis
swap_dis = as.data.frame(table(swaps[, c("name")]))
names(swap_dis) = c("name", "swaps")
swap_dis = merge(swap_dis, c_dis)
swap_per = mutate(swap_dis, percentage = swaps / correct * 100)
swap_per[sapply(swap_per, is.infinite)] <- 100

p = ggplot(swap_dis) +
  geom_col(aes(x = name, y = swaps, fill = name)) +
  theme_minimal()+ 
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) 
if (exists("color_mapping")) {
  p = p+
    scale_fill_manual(values = color_mapping)
}
ggsave(paste(dir,"/swaps_distribution.jpg",sep=""), plot = p)
ggsave(paste(dir,"/swaps_distribution.pdf",sep=""), plot = p)


#OF INTEREST
if (is.null(clones) == F) {
  cs <- readLines(paste0(dir, "/", clones))
  ofint = filter(complete_table, clone %in% cs)
  length(unique(cs))
  #ofint = left_join(ofint, barcodes_names, by = "barcode")
  #ofint = separate(ofint, name, into = ("gene"), sep = "\\.", remove = F)
  #dist 
  clone_t = as.data.frame(table(ofint[,c("clone")]))
  length(clone_t$Var1)
  length(filter(clone_t, Freq == 1)$Var1)
  clone_t = separate(clone_t, Var1, into = c("barcode", "UCI"), sep = "_", remove = F)
  clone_t = left_join(clone_t, barcodes_names)
  clone_t = separate(clone_t, name, into = ("gene"), sep = "\\.", remove = F)
  clone_t = filter(clone_t, Freq > 1)
  p = ggplot(clone_t) +
    geom_col(aes(x = fct_reorder(Var1, Freq, .desc = T), y = Freq, fill = gene))  +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  ggsave(paste(dir,"/clone_distribution_ofint.jpg",sep=""),
         plot = p,
         width = 10, height = 15,
         limitsize = F)
  ggsave(paste(dir,"/clone_distribution_ofint.pdf",sep=""),
         plot = p,
         width = 10, height = 15,
         limitsize = F)
  
  # #modify shRNA to reach 21Ns and no AAAA
  # ofint = ofint %>%
  #   mutate(shRNA = str_sub(shRNA, 5, -1))
  # expected = read.table(paste(dir,"/expected_shRNA_names.txt", sep=""), sep = " ")
  # names(expected) = c("actual_name", "shRNA")
  # #reverse_complement <- function(sequence) paste0(sapply(rev(strsplit(toupper(sequence), NULL)[[1]]), function(n) c("T" = "A", "A" = "T", "C" = "G", "G" = "C")[n]), collapse = "")
  # #expected = mutate(expected, shRNA = reverse_complement(shRNA))
  # expected = expected %>%
  #   mutate(shRNA = str_sub(shRNA, 1, 21))
  #add expected name based on shRNA
  ofint = filter(ofint, shRNA %in% expected$shRNA)
  #ofint = left_join(ofint, expected, by = "shRNA")
  ofint = distinct(ofint)
  length(ofint$UMI)
  #dist
  clone_t = as.data.frame(table(ofint[,c("clone")]))
  clone_t = separate(clone_t, Var1, into = c("barcode", "UCI"), sep = "_", remove = F)
  clone_t = left_join(clone_t, barcodes_names)
  clone_t = separate(clone_t, name, into = ("gene"), sep = "\\.", remove = F)
  clone_t = filter(clone_t, Freq > 1)
  p = ggplot(clone_t) +
    geom_col(aes(x = fct_reorder(Var1, Freq, .desc = T), y = Freq, fill = gene))  +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
  
  ggsave(paste(dir,"/clone_distribution_ofint_shRNA.jpg",sep=""),
         plot = p,
         width = 10, height = 15,
         limitsize = F)
  ggsave(paste(dir,"/clone_distribution_ofint_shRNA.pdf",sep=""),
         plot = p,
         width = 10, height = 15,
         limitsize = F)
  #count reads
  t = as.data.frame(table(ofint[,c("shRNA", "clone", "gene", "name", "actual_name")]))
  t = filter(t, Freq > 0)
  #add same col
  t$name = as.character(t$name)
  t$actual_name = as.character(t$actual_name)
  t$same = with(t, ifelse(name == actual_name, 'yes',
                          ifelse(name != actual_name, 'no', NA)))
  write.table(t, paste(dir, "table_of_int.txt", sep = ""), sep = "\t")
  #plot
  dir.create(paste(dir, "/shRNA_distribution_UCI/", sep = ""))
  for (seq in unique(t$clone)) {
    data = filter(t, clone == seq)
    write.table(data, paste(dir, "/shRNA_distribution_UCI/table_", data$name[[1]], "_",seq, ".txt", sep = ""), sep = "\t")
    p = ggplot(data) + 
      geom_col(aes(x = actual_name, y = Freq, fill = same))  +
      labs(title = data$name) +
      scale_fill_manual(values = c("#00AFBB", "#E7B800"), drop = F) +
      theme_minimal() + 
      theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
    ggsave(paste(dir,"/shRNA_distribution_UCI/", seq, "_of_int.jpg",sep=""), plot = p)
    ggsave(paste(dir,"/shRNA_distribution_UCI/", seq, "_of_int.pdf",sep=""), plot = p)
  }
  
  #auto swaps
  good_c = c()
  res_table = data.frame()
  for (c in unique(t$clone)) {
    df = filter(t, clone == c)
    if (max(df$Freq) > 50) {
      #print(max(df$Freq))
      freqs = sort(df$Freq, decreasing = TRUE)
      #print(freqs[2])
      if (length(df$Freq) == 1) {
        good_c = append(good_c, c)
        res_table = rbind(res_table, filter(df, Freq == freqs[1]))
      }
      else if (freqs[1] / freqs[2] > 10) {
        good_c = append(good_c, c)
        res_table = rbind(res_table, filter(df, Freq == freqs[1]))
      }
    }
  }
  
  good_table = filter(table, clone %in% good_c)
  
  swaps = filter(res_table, same == "no")
  print(paste0("Swaps: ", length(swaps$clone)))
  correct = filter(res_table, same == "yes")
  print(paste0("Correct: ", length(correct$clone)))
  
  write.table(good_table, paste(dir, "/reliable_clones_ofint.txt", sep = ""), sep = "\t")
  write.table(swaps, paste(dir, "/reliable_clones_swaps_ofint.txt", sep = ""), sep = "\t")
  write.table(correct, paste(dir, "/reliable_clones_confirmations_ofint.txt", sep = ""), sep = "\t")
  
}
