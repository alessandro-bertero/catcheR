library(ggplot2)
library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

#1 dir
dir = args[1]
#dir = "/20tb/ratto/Bertero/plasmid analysis/MiSeq/"
#2 threshold percentage
DIs = args[2]
# #3 threshold
# threshold = args[3]
#4 clones of int
clones = args[3]
if (clones == "NULL" | (clones == "")) {
  clones = NULL
}
print(clones)

#fastqc(group = "docker", data.folder = paste(dir, "./fastqc", sep = ""))

UMI = fread(paste(dir,"/final_UMI.txt", sep=""), sep = "\n", data.table = F)
UCI = fread(paste(dir,"/final_UCI.txt", sep=""), sep = "\n", data.table = F)
reference = fread(paste(dir,"/final_reference.txt", sep=""), sep = "\n", data.table = F)
barcode = fread(paste(dir,"/final_BC.txt", sep=""), sep = "\n", data.table = F)
complete_table = data.frame(UMI, barcode, UCI)
names(complete_table) = c("UMI","barcode", "UCI")
complete_table = mutate(complete_table, clone = paste(barcode, UCI, sep = "_"))
write.csv(complete_table, paste(dir, "/complete_table_final_plasmid.csv", sep = ""))

#generate clone t
clone_t = as.data.frame(table(complete_table[,c("clone")]))
length(clone_t$Var1)
length(filter(clone_t, Freq == 1)$Var1)
clone_t = separate(clone_t, Var1, into = c("barcode", "UCI"), sep = "_", remove = F)
barcodes_names = distinct(read.table(paste(dir,"/rc_barcodes_genes.csv", sep=""), sep = ","))
names(barcodes_names) = c("barcode", "name")
clone_t = left_join(clone_t, barcodes_names, by = "barcode")
clone_t = na.omit(clone_t)

#PERCENTAGE
f300 = filter(clone_t, Freq > DIs)
tmp = filter(complete_table, clone %in% f300$Var1)
tmp = distinct(tmp[,c("barcode", "UCI")])
tmp = as.data.frame(table(tmp[,c("barcode")]))
sum(tmp$Freq)
clone_t_per = tmp
clone_t_per = mutate(clone_t_per, percentage = Freq / sum(Freq) * 100)

p = ggplot(clone_t_per, aes(x = percentage)) +
  geom_histogram() +
  theme_minimal()
ggsave(paste(dir,"/clone_percentage_all.jpg",sep=""), 
       plot = p, 
       width = 20, height = 15, 
       limitsize = F)
ggsave(paste(dir,"/clone_percentage_all.pdf",sep=""), 
       plot = p, 
       width = 20, height = 15, 
       limitsize = F)
write.csv(clone_t_per, paste(dir, "/percentages.csv", sep = ""))

clone_t = filter(clone_t, Freq > 1)
clone_t = separate(data = clone_t, col = name, into = c("gene"),sep = "\\.", remove = F)
# p = ggplot(clone_t) + 
#   geom_col(aes(x = fct_reorder(Var1, Freq, .desc = T), y = Freq, fill = name))  +
#   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
#   theme_minimal()
# #labs(title = data$name) +
# #scale_fill_manual(values = c("#00AFBB", "#E7B800"), drop = F) 
# if (file.exists(paste(dir, "/colors.csv", sep = ""))) {
#   colors <- read.csv(paste(dir, "/colors.csv", sep = ""), header = F, col.names = c("name", "colors"))
#   color_mapping <- setNames(colors$colors, colors$name)
#   p = p+
#     scale_fill_manual(values = color_mapping)
#   
# } else {
#   print("Colors file does not exist, using palette.")
# }
# 
# ggsave(paste(dir,"/clone_distribution_all.jpg",sep=""), 
#        plot = p, 
#        width = 20, height = 15, 
#        limitsize = F)
# ggsave(paste(dir,"/clone_distribution_all.pdf",sep=""), 
#        plot = p, 
#        width = 20, height = 15, 
#        limitsize = F)

write.csv(clone_t, paste(dir, "/distribution_all_clones.csv", sep = ""))

#filter
clone_t = filter(clone_t, Freq > DIs)
p = ggplot(clone_t) + 
  geom_col(aes(x = fct_reorder(Var1, Freq, .desc = T), y = Freq, fill = name))  +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  theme_minimal()
#labs(title = data$name) +
#scale_fill_manual(values = c("#00AFBB", "#E7B800"), drop = F) 
if (file.exists(paste(dir, "/colors.csv", sep = ""))) {
  colors <- read.csv(paste(dir, "/colors.csv", sep = ""), header = F, col.names = c("name", "colors"))
  color_mapping <- setNames(colors$colors, colors$name)
  p = p+
    scale_fill_manual(values = color_mapping)
  
} else {
  print("Colors file does not exist, using palette.")
}
ggsave(paste(dir,"/clone_distribution_filter_", DIs, ".jpg",sep=""), 
       plot = p, 
       width = 20, height = 15, 
       limitsize = F)
ggsave(paste(dir,"/clone_distribution_filter_", DIs, ".pdf",sep=""), 
       plot = p, 
       width = 20, height = 15, 
       limitsize = F)

#histogram 
# p = ggplot(clone_t, aes(x = Freq)) +
#   geom_histogram(aes(fill = name), binwidth = 200) + 
#   labs(x = "Reads per clone", 
#        y = "Log10 number of clones", 
#        title = "Histogram of reads per clone, no filtering, log10") +
#   scale_y_log10() + 
#   theme_minimal()
# if (file.exists(paste(dir, "/colors.csv", sep = ""))) {
#   colors <- read.csv(paste(dir, "/colors.csv", sep = ""), header = F, col.names = c("name", "colors"))
#   color_mapping <- setNames(colors$colors, colors$name)
#   p = p+
#     scale_fill_manual(values = color_mapping)
#   
# } else {
#   print("Colors file does not exist, using palette.")
# }
# ggsave(paste(dir,"/hist_all.jpg",sep=""), 
#        plot = p, 
#        width = 20, height = 15, 
#        limitsize = F)
# ggsave(paste(dir,"/hist_all.pdf",sep=""), 
#        plot = p, 
#        width = 20, height = 15, 
#        limitsize = F)

clone_t = filter(clone_t, Freq > DIs)
p = ggplot(clone_t) +
  geom_histogram(aes(y = after_stat(density), x = Freq, fill = name), binwidth = 1000) + 
  geom_line(aes(y = after_stat(density), x = Freq),stat = "density")+
  #geom_density(aes(y = after_stat(density), x = Freq))+
  labs(x = "Reads per clone", 
       y = "Number of clones", 
       title = paste("Histogram of reads per clone, Freq > ", DIs, sep = "")) +
  theme_minimal() #+
#stat_function(fun = dnorm, args = list(mean = mean(clone_t$Freq), sd = sd(clone_t$Freq)))
#coord_cartesian(ylim = c(0,5))
if (file.exists(paste(dir, "/colors.csv", sep = ""))) {
  colors <- read.csv(paste(dir, "/colors.csv", sep = ""), header = F, col.names = c("name", "colors"))
  color_mapping <- setNames(colors$colors, colors$name)
  p = p+
    scale_fill_manual(values = color_mapping)
  
} else {
  print("Colors file does not exist, using palette.")
}

ggsave(paste(dir,"/density_all_filtering_", DIs, ".jpg",sep=""), 
       plot = p, 
       width = 20, height = 15, 
       limitsize = F)
ggsave(paste(dir,"/density_all_filtering_", DIs, ".pdf",sep=""), 
       plot = p, 
       width = 20, height = 15, 
       limitsize = F)


#PIE CHART ALL

p = ggplot(clone_t, aes(x="", y=Freq, fill=name))+
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) + 
  theme_minimal()

if (file.exists(paste(dir, "/colors.csv", sep = ""))) {
  colors <- read.csv(paste(dir, "/colors.csv", sep = ""), header = F, col.names = c("name", "colors"))
  color_mapping <- setNames(colors$colors, colors$name)
  p = p+
    scale_fill_manual(values = color_mapping)
  
} else {
  print("Colors file does not exist, using palette.")
}

ggsave(paste(dir,"/pie_chart_all.jpg",sep=""), 
       plot = p, 
       width = 15, height = 15, 
       limitsize = F)
ggsave(paste(dir,"/pie_chart_all.pdf",sep=""), 
       plot = p, 
       width = 15, height = 15, 
       limitsize = F)

clone_t = separate(clone_t, name, into = c("gene"), sep = "\\.", remove = F)
p = ggplot(clone_t, aes(x="", y=Freq, fill=gene))+
  geom_bar(width = 1, stat = "identity") + 
  coord_polar("y", start=0) + 
  theme_minimal()
ggsave(paste(dir,"/pie_chart_all_gene.jpg",sep=""), 
       plot = p, 
       width = 15, height = 15, 
       limitsize = F)
ggsave(paste(dir,"/pie_chart_all_gene.pdf",sep=""), 
       plot = p, 
       width = 15, height = 15, 
       limitsize = F)


#only of int
if (is.null(clones) == F) {
  clones <- readLines(paste0(dir, "/clones.txt"))
  ofint = filter(complete_table, clone %in% clones)
  length(unique(ofint$clone))
  length(unique(clones))
  #distribution of clone
  clone_t = as.data.frame(table(ofint[,c("clone")]))
  length(clone_t$Var1)
  length(filter(clone_t, Freq == 1)$Var1)
  clone_t = separate(clone_t, Var1, into = c("barcode", "UCI"), sep = "_", remove = F)
  barcodes_names = distinct(read.table(paste(dir,"/rc_barcodes_genes.csv", sep=""), sep = ","))
  names(barcodes_names) = c("barcode", "name")
  clone_t = left_join(clone_t, barcodes_names)
  clone_t = filter(clone_t, Freq > DIs)
  p = ggplot(clone_t) + 
    geom_col(aes(x = fct_reorder(Var1, Freq, .desc = T), y = Freq, fill = name))  +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
    theme_minimal()
  #labs(title = data$name) +
  #scale_fill_manual(values = c("#00AFBB", "#E7B800"), drop = F) 
  ggsave(paste(dir,"/clone_distribution_ofint_filter_", DIs, ".jpg",sep=""), 
         plot = p, 
         width = 20, height = 15, 
         limitsize = F)
  ggsave(paste(dir,"/clone_distribution_ofint_filter_", DIs, ".pdf",sep=""), 
         plot = p, 
         width = 20, height = 15, 
         limitsize = F)
  
  #distribution of clone TOT
  clone_t = as.data.frame(table(ofint[,c("clone")]))
  length(clone_t$Var1)
  length(filter(clone_t, Freq == 1)$Var1)
  clone_t = separate(clone_t, Var1, into = c("barcode", "UCI"), sep = "_", remove = F)
  barcodes_names = distinct(read.table(paste(dir,"/rc_barcodes_genes.csv", sep=""), sep = ","))
  names(barcodes_names) = c("barcode", "name")
  clone_t = left_join(clone_t, barcodes_names)
  write.csv(clone_t, paste(dir, "/clones_final_plasmid.csv", sep = ""))
  #clone_t = filter(clone_t, Freq > 2000)
  p = ggplot(clone_t) + 
    geom_col(aes(x = fct_reorder(Var1, Freq, .desc = T), y = Freq, fill = name))  +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) + 
    theme_minimal()
  #labs(title = data$name) +
  #scale_fill_manual(values = c("#00AFBB", "#E7B800"), drop = F) 
  ggsave(paste(dir,"/clone_distribution_ofint_nofilter.jpg",sep=""), 
         plot = p, 
         width = 20, height = 15, 
         limitsize = F)
  ggsave(paste(dir,"/clone_distribution_ofint_nofilter.pdf",sep=""), 
         plot = p, 
         width = 20, height = 15, 
         limitsize = F)
  
  #PIE CHART
  p = ggplot(clone_t, aes(x="", y=Freq, fill=name))+
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0)+ 
    theme_minimal()
  
  if (file.exists(paste(dir, "/colors.csv", sep = ""))) {
    colors <- read.csv(paste(dir, "/colors.csv", sep = ""), header = F, col.names = c("name", "colors"))
    color_mapping <- setNames(colors$colors, colors$name)
    p = p+
      scale_fill_manual(values = color_mapping)
    
  } else {
    print("Colors file does not exist, using palette.")
  }
  
  ggsave(paste(dir,"/pie_chart_ofint.jpg",sep=""), 
         plot = p, 
         width = 15, height = 15, 
         limitsize = F)
  ggsave(paste(dir,"/pie_chart_ofint.pdf",sep=""), 
         plot = p, 
         width = 15, height = 15, 
         limitsize = F)
  
  clone_t = separate(clone_t, name, into = c("gene"), sep = "\\.", remove = F)
  p = ggplot(clone_t, aes(x="", y=Freq, fill=gene))+
    geom_bar(width = 1, stat = "identity") + 
    coord_polar("y", start=0)+ 
    theme_minimal()
  ggsave(paste(dir,"/pie_chart_ofint_gene.jpg",sep=""), 
         plot = p, 
         width = 15, height = 15, 
         limitsize = F)
  ggsave(paste(dir,"/pie_chart_ofint_gene.pdf",sep=""), 
         plot = p, 
         width = 15, height = 15, 
         limitsize = F)
  
}

#EMPTY
complete_table = data.frame(UMI, barcode, UCI, reference)
complete_table = distinct(complete_table)
names(complete_table) = c("UMI","barcode", "UCI", "reference")

#Select empty reads #TACGCGTTCATCTGGGGGAGCCG
tr_complete_table = filter(complete_table, reference == "TACGCGTTCATCTGGGGGAGCCG")
                                                         
#Save complete table
write.csv(tr_complete_table, paste(dir,"/complete_table_empty.csv",sep=""))


print(paste("Empty reads:", length(unique(tr_complete_table$UMI)), sep = ""))

