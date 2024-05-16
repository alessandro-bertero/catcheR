#load libraries 
library(rCASC)
library(tidyverse)
library(ggplot2)
library("SAVER") 
library(gtools)
library(phateR)
library(rlang)

#set timepoints. Directories for each timepoint are called "D" followed by the number of days. If there's only one timepoint, write just one number.
timepoints = c(".")

#set dir
home = "/20tb/ratto/AB03/Analysis/sc/"
scratch = "/scratch"
setwd(home)

# QC, normalization, annotation 
for (d in timepoints) {
  setwd(paste(home, d, sep = ""))
  mitoRiboUmi(
    group = "docker",
    scratch.folder="/scratch",
    file= paste(home, d, "/silencing_matrix_complete.csv", sep =""),
    separator=",",
    gtf.name="Homo_sapiens.GRCh38.109.gtf",
    bio.type="protein_coding",
    umiXgene=1
  )
}

d = "00"
tmp = read.table(paste(home,d, "/silencing_matrix_complete.csv", sep=""), sep=",", header=TRUE, row.names = 1)
colnames(tmp) = paste(colnames(tmp), d, sep = "_")
tmp00 = tmp
d = "02"
tmp = read.table(paste(home,d, "/silencing_matrix_complete.csv", sep=""), sep=",", header=TRUE, row.names = 1)
colnames(tmp) = paste(colnames(tmp), d, sep = "_")
tmp02 = tmp
d = "06"
tmp = read.table(paste(home,d, "/silencing_matrix_complete.csv", sep=""), sep=",", header=TRUE, row.names = 1)
colnames(tmp) = paste(colnames(tmp), d, sep = "_")
tmp06 = tmp
d = "12"
tmp = read.table(paste(home,d, "/silencing_matrix_complete.csv", sep=""), sep=",", header=TRUE, row.names = 1)
colnames(tmp) = paste(colnames(tmp), d, sep = "_")
tmp12 = tmp

raw.data <- cbind(tmp00, tmp02)
raw.data <- cbind(raw.data, tmp06)
raw.data <- cbind(raw.data, tmp12)

cortex.saver <- SAVER::saver(as.matrix(raw.data), estimates.only = T, ncores = 20) 
write.table(cortex.saver,paste(home,"SAVER_matrix.csv", sep=""),sep=",", row.names = T, col.names = T)


data <- data.frame(read.table(paste(home,"SAVER_matrix.csv", sep=""), sep = ",", header = T))
data <- separate(
  data,
  col="X",
  into=c("a","b"),
  sep = ":",
  remove = TRUE,
  convert = FALSE,
  extra = "warn",
  fill = "warn"
)
#data$b[1] <- "genes"
data$a <- NULL
names <- make.unique(data$a, sep=".")
#data[duplicated(data$b),1]
rownames(data) <- names
data$b <- NULL
write.table(data,paste(home,"SAVER_matrix_symbol.csv", sep=""),col.names=NA,sep=",")

scannobyGtf(
  group = "docker",
  file=paste(home, "SAVER_matrix.csv", sep = ""),
  gtf.name="Homo_sapiens.GRCh38.109.gtf",
  biotype = "protein_coding",
  mt = F,
  ribo.proteins = F,
  umiXgene = 3,
  riboStart.percentage = 0,
  riboEnd.percentage = 100,
  mitoStart.percentage = 0,
  mitoEnd.percentage = 0,
  thresholdGenes = 100
)


#create complete matrix
com_mat = function(timepoints, home, file, name){
  clust_fin = data.frame()
  for (day in timepoints) {
  #expression matrix
    exp = read.csv(paste(home,"/", day, "/",file, sep=""), header = T)
    data <- separate(
      exp,
      col="X",
      into=c("a","b"),
      sep = ":",
      remove = TRUE,
      convert = FALSE,
      extra = "warn",
      fill = "warn"
    )
    data$a <- NULL
    names <- make.unique(data$b, sep=".")
    rownames(data) <- names
    data$b <- NULL
  
    data = data.frame(t(data))
    rownames(data)= paste(rownames(data),day,name,sep = "_")
    #join
    clust_fin = rbind(clust_fin, data)
    return(clust_fin)
    write.table(clust_fin, paste(home, "/all_timepoints_expression.csv", sep = ""), sep = ",")
  }
}

#PHATE 
analysis_phate = function(data, home, clusters){
  dir.create(paste(home, "/PHATE_analysis", sep=""))
  setwd(paste(home, "/PHATE_analysis", sep=""))
  
  data_PHATE <- phate(data, mds.solver = "smacof")
  plot = ggplot(data_PHATE) +
    geom_point(aes(PHATE1, PHATE2))
  ggsave("PHATE_auto.pdf",plot = plot,device = "pdf", width = 30, height = 30, units = "cm")
  
  cluster_PHATE = cluster_phate(data_PHATE, k = clusters, seed = 1111)
  PHATE_dataframe = as.data.frame(data_PHATE)
  PHATE_dataframe$cluster = cluster_PHATE
  PHATE_dataframe$nomi = rownames(PHATE_dataframe)
  PHATE_dataframe = separate(PHATE_dataframe, nomi, into = c("cellID", "UMI", "barcode", "gene", "UCI", "day", "experiment"), sep = "_", remove = T, convert = T)
  PHATE_dataframe[is.na(PHATE_dataframe)] <- "empty"
  PHATE_dataframe = PHATE_dataframe %>% mutate(
    Control = case_when(gene == "SCR" | gene == "B2M" | is.na(gene) | gene == "empty" ~ "Control",
    TRUE ~ "KD"))
  
  plot = ggplot(PHATE_dataframe) +
    geom_point(aes(PHATE1, PHATE2, color = as.factor(cluster))) +
    labs(color="Cluster")
  ggsave("PHATE_clustering.pdf",plot = plot,device = "pdf", width = 30, height = 30, units = "cm")
  
  
  plot= ggplot(PHATE_dataframe) +
    geom_point(aes(PHATE1, PHATE2, color = as.factor(gene))) +
    labs(color="Gene") +
    facet_wrap(vars(UCI))
  ggsave("PHATE_facet_UCI.pdf",plot = plot,device = "pdf", width = 30, height = 30, units = "cm")
  
  
  plot=ggplot(PHATE_dataframe) +
    geom_point(aes(PHATE1, PHATE2, color = as.factor(gene))) +
    labs(color="Gene")
  ggsave("PHATE_gene.pdf",plot = plot,device = "pdf", width = 30, height = 30, units = "cm")
  
  
  plot=ggplot(PHATE_dataframe) +
    geom_point(aes(PHATE1, PHATE2, color = as.factor(Control), alpha = 0.4)) +
    labs(color="Control")
  ggsave("PHATE_control.pdf",plot = plot,device = "pdf", width = 30, height = 30, units = "cm")
  
  
  plot=ggplot(PHATE_dataframe) +
    geom_point(aes(PHATE1, PHATE2, color = as.factor(gene))) +
    labs(color="Gene")+
    facet_wrap(vars(Control))
  ggsave("PHATE_facet_control.pdf",plot = plot,device = "pdf", width = 60, height = 30, units = "cm")
  
  return(PHATE_dataframe)
}


#SELECT CLONES with >20 cells in at least 1 timepoint
big_c = function(data, n){
  data$clone <- paste(data$barcode,data$UCI, sep = "_")
  #cellsxclone
  cellsxclone = distinct(data[,c("cellID", "clone")])
  cellsxclone_table = data.frame(table(cellsxclone$clone))
  names(cellsxclone_table) = c("clone", "cellsxclone")
  data = merge(data,cellsxclone_table)
  keep = c()
  for (c in unique(data$clone)) {
    d = filter(data, clone == c)
    if (max(d$cellsxclone) >= n) {
      keep = append(c, keep)
    }
  }
  big_clones = subset(data, clone %in% keep)
  big_clones$gene = replace_na(big_clones$gene, "empty")
  big_clones$UMI = replace_na(big_clones$UMI, median(big_clones$UMI))
  big_clones$UCI = replace_na(big_clones$UCI, "empty")
  return(big_clones)
}

clust_fin5 = com_mat(timepoints = c("CM", "PSC"), home = "~/scripts/5", name = "exp5", file = "filtered_annotated_SAVER_matrix.csv")
clust_fin7 = com_mat(timepoints = c("D0","D2","D6","D12"), home = "~/scripts/7", name = "exp7", file = "filtered_annotated_SAVER_matrix_symbol.csv")
clust_fin_AB03 = com_mat(timepoints = ".", home = home, name = "AB03", file = "filtered_annotated_SAVER_matrix.csv")

data = clust_fin_AB03
PHATE_dataframe = analysis_phate(clust_fin_AB03, home = home, clusters = 2)
big_clones = big_c(PHATE_dataframe, 20)
sel = data[,c("GATA4", "SMAD2", "KMT2D", "B2M", "NKX2.5", "CHD7", "MYL7", "TTN", "TNNI1", "TNNI2")] #EGFP?
sel$nomi = rownames(sel)
sel = separate(sel, nomi, into = c("cellID"), sep = "_", remove = T, convert = T)
clust_fin = merge(sel, big_clones, by = "cellID")

view#PHATE
bg0 = filter(big_clones, day == "0") %>% merge(mytable_perc_day0[c("clone","side")], by = "clone")
bg2 = merge(filter(big_clones, day == "2"), mytable_perc_day0[c("clone","side")])
bg6 = merge(filter(big_clones, day == "6"), mytable_perc_day0[c("clone","side")])
bg12 = merge(filter(big_clones, day =="12"), mytable_perc_day0[c("clone","side")])

big_clones = merge(big_clones, mytable_perc_day0[c("clone","side")], by = "clone")

#THEME FOR GGPLOT
theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family = 'Helvetica')
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#D3D3D3"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
}

#VIOLIN PLOTS PDFs 
#KD GENE EXPRESSION
genes= c("GATA4", "SMAD2", "KMT2D", "B2M", "NKX2.5", "CHD7")
for (g in genes) {
  print(g)
  for (d in timepoints) {
    print(d)
    #for (s in unique(big_clones$side)) {
      f = filter(clust_fin, as.character(day) == as.character(d)) %>% filter(gene == g | gene == "SCR" | gene == "empty" | gene == "B2M") #%>% filter(side == s)
      gene_exp = sym(g)
      plot= ggplot(f) +
        geom_violin(aes(x = clone, !! gene_exp, colour=gene, fill = gene)) +
        stat_summary(aes(x = clone, !! gene_exp),fun = "median", geom = "crossbar", width = 0.5)+
        labs(title = paste(g," expression day ",d, sep = ""))+
        #facet_grid(col=vars(as.factor(side)))+
        theme_Publication()
      ggsave(paste("violin_",as.character(g),"_expression_day",d, ".pdf", sep=""),plot = plot,device = "pdf", width = 50, height = 20, units = "cm")
    #}
  }
}

#MARKERS GENES EXPRESSION
markers12 = c("MYL7", "TTN", "TNNI1", "TNNI2")
for (g in markers12) {
    #for (s in unique(clust_fin$side)) {
      f = filter(clust_fin) #, as.character(day) == as.character(12)) %>% filter(side == s)
      gene_exp = sym(g)
      plot= ggplot(f) +
        geom_violin(aes(x = clone, !! gene_exp, colour=gene, fill = gene)) +
        stat_summary(aes(x = clone, !! gene_exp),fun = "median", geom = "crossbar", width = 0.5)+
        labs(title = paste(g," expression day 12", sep = ""))+
        theme_Publication()
      ggsave(paste("violin_",as.character(g),"_expression_day12","_side", s, ".pdf", sep=""),plot = plot,device = "pdf", width = 50, height = 20, units = "cm")
  #}
}

markers6 = c("MYL7", "TTN", "TNNI1", "TNNI2")
for (g in markers6) {
    for (s in unique(big_clones$side)) {
      f = filter(big_clones, as.character(day) == as.character(6)) %>% filter(side == s)
      gene_exp = sym(g)
      plot= ggplot(f) +
        geom_violin(aes(x = clone, !! gene_exp, colour=gene, fill = gene)) +
        stat_summary(aes(x = clone, !! gene_exp),fun = "median", geom = "crossbar", width = 0.5)+
        labs(title = paste(g," expression day ",d, sep = ""))+
        theme_Publication()
      ggsave(paste("violin_",as.character(g),"_expression_day6","_side", s, ".pdf", sep=""),plot = plot,device = "pdf", width = 50, height = 20, units = "cm")
  }
}

markers2 = c("POU5F1")
for (g in markers2) {
    print(g)
      for (s in unique(big_clones$side)) {
        f = filter(big_clones, as.character(day) == as.character(2)) %>% filter(side == s)
        gene_exp = sym(g)
        plot= ggplot(f) +
          geom_violin(aes(x = clone, !! gene_exp, colour=gene, fill = gene)) +
          stat_summary(aes(x = clone, !! gene_exp),fun = "median", geom = "crossbar", width = 0.5)+
          labs(title = paste(g," expression day ",d, sep = ""))+
          theme_Publication()
        ggsave(paste("violin_",as.character(g),"_expression_day2","_side", s, ".pdf", sep=""),plot = plot,device = "pdf", width = 50, height = 20, units = "cm")
    }
}

#MARKER CM 
markersCM = c("MYL7", "TTN", "TNNI1", "TNNI2")
for (g in markersCM) {
  f = filter(big_clones, as.character(day) == as.character("CM"))
  gene_exp = sym(g)
  plot= ggplot(f) +
    geom_violin(aes(x = clone, !! gene_exp, colour=gene, fill = gene)) +
    stat_summary(aes(x = clone, !! gene_exp),fun = "median", geom = "crossbar", width = 0.5)+
    labs(title = paste(g," expression CM", sep = ""))+
    theme_Publication()
  ggsave(paste("violin_",as.character(g),"_expression_CM","_side", s, ".pdf", sep=""),plot = plot,device = "pdf", width = 50, height = 20, units = "cm")
}

markersPSC = c("POU5F1")
for (g in markersPSC) {
  print(g)
  f = filter(big_clones, as.character(day) == as.character("PSC")) 
  gene_exp = sym(g)
  plot= ggplot(f) +
    geom_violin(aes(x = clone, !! gene_exp, colour=gene, fill = gene)) +
    stat_summary(aes(x = clone, !! gene_exp),fun = "median", geom = "crossbar", width = 0.5)+
    labs(title = paste(g," expression PSC", sep = ""))+
    theme_Publication()
  ggsave(paste("violin_",as.character(g),"_expression_PSC","_side", s, ".pdf", sep=""),plot = plot,device = "pdf", width = 50, height = 20, units = "cm")
}

#CLONE depletion

#cellsxclone percentage
clones = distinct(big_clones[, c("clone", "gene", "cellsxclone", "barcode", "day")])
# tot_cells_PSC = distinct(filter(clones, day == "PSC")[, c("clone", "cellsxclone")])
# tot_cells_PSC = sum(tot_cells_PSC$cellsxclone)
# tot_cells_CM = distinct(filter(clones, day == "CM")[, c("clone", "cellsxclone")])
# tot_cells_CM = sum(tot_cells_CM$cellsxclone)
# clones_PSC = mutate(filter(clones, day == "PSC"), percentage_cellsxclone = (cellsxclone/tot_cells_PSC)*100)
# clones_CM = mutate(filter(clones, day == "CM"), percentage_cellsxclone = (cellsxclone/tot_cells_CM)*100)
# clones = rbind(clones_PSC, clones_CM)
# cl = clones_PSC[,c("clone","percentage_cellsxclone")]
# names(cl) = c("clone", "starting_percentage")
# clones = merge(cl, clones)
# plot=ggplot(clones)+
#   geom_col(aes(y=log2fc, x =clone, fill = gene))+
#   labs(y="Log2 foldchange of percentage of clone cells on total cells, day 12 vs day 0", 
#        x = "Clones", 
#        title = "Clone depletion")+
#   theme_Publication()
# ggsave(paste("barplot_foldchange_side", s, ".pdf", sep=""),plot = plot,device = "pdf", width = 40, height = 35, units = "cm")
# 
# plot=ggplot(filter(clones, day == "CM"))+
#   geom_point(aes(x=starting_percentage, y = percentage_cellsxclone, colour = gene, size = 2, alpha = 0.5))+
#   labs(x="Percentage of clone cells on total cells at PSC", 
#        y = "Percentage of clone cells on total cells at CM", 
#        title = "Clone depletion")+
#   xlim(c(0,21))+
#   ylim(c(0,21))+
#   geom_abline(intercept = 0, slope = 1, linetype="dotted", alpha=0.5, linewidth = 1, colour = "red")+
#   theme_Publication()
# ggsave(paste("2Dplot_foldchange.pdf", sep=""),plot = plot,device = "pdf", width = 40, height = 35, units = "cm")
# 
# clones = mutate(clones, log2fc = log2(as.numeric(percentage_cellsxclone/starting_percentage)))
# plot=ggplot(filter(clones, day == "CM"))+
#   geom_point(aes(x=log2fc, y =starting_percentage, colour = gene, size = 2, alpha = 0.3))+
#   geom_vline(xintercept = 1,linetype="dotted")+
#   annotate("text", x=1.1, y=10, label="Clones with doubled number of cells", angle=90)+
#   geom_vline(xintercept = -1,linetype="dotted")+
#   annotate("text", x=-1.1, y=10, label="Clones with halved number of cells", angle=90)+
#   xlim(c(-5,5))+
#   labs(y="Starting percentage of clone cells on total cells at PSC", 
#        x = "Log2 foldchange of percentage of clone cells on total cells, CM vs PSC", 
#        title = "Clone depletion vulcano plot")+
#   theme_Publication()
# ggsave(paste("vuncano_plot_depletion.pdf", sep=""),plot = plot,device = "pdf", width = 40, height = 35, units = "cm")




tot_cells_0 = distinct(filter(clones, day == 0)[, c("clone", "cellsxclone")])
tot_cells_0 = sum(tot_cells_0$cellsxclone)
tot_cells_2 = distinct(filter(clones, day == 2)[, c("clone", "cellsxclone")])
tot_cells_2 = sum(tot_cells_2$cellsxclone)
tot_cells_6 = distinct(filter(clones, day == 6)[, c("clone", "cellsxclone")])
tot_cells_6 = sum(tot_cells_6$cellsxclone)
tot_cells_12 = distinct(filter(clones, day == 12)[, c("clone", "cellsxclone")])
tot_cells_12 = sum(tot_cells_12$cellsxclone)

clones_0 = mutate(filter(clones, day == "0"), percentage_cellsxclone = (cellsxclone/tot_cells_0)*100)
clones_2 = mutate(filter(clones, day == "2"), percentage_cellsxclone = (cellsxclone/tot_cells_2)*100)
clones_6 = mutate(filter(clones, day == "6"), percentage_cellsxclone = (cellsxclone/tot_cells_6)*100)
clones_12 = mutate(filter(clones, day == "12"), percentage_cellsxclone = (cellsxclone/tot_cells_12)*100)
clones = rbind(clones_0, clones_2, clones_6, clones_12)
cl = clones_0[,c("clone","percentage_cellsxclone")]
names(cl) = c("clone", "starting_percentage")
clones = merge(cl, clones)

tmp = distinct(mytable_perc_day0[,c("clone","total")])
tot = sum(as.numeric(tmp$total))
mytable_perc_day0 = mutate(mytable_perc_day0, starting_percentage = (as.numeric(mytable_perc_day0$total) / tot)*100)

clones_12 = merge(clones_12, mytable_perc_day0, by = "clone")
clones_12 = mutate(clones_12, log2fc = log2(as.numeric(percentage_cellsxclone/starting_percentage)))

for (s in unique(clones_12$side)) {
  print(s)
  plot=ggplot(filter(clones_12, side == s))+
    geom_col(aes(y=log2fc, x =clone, fill = gene))+
    labs(y="Log2 foldchange of percentage of clone cells on total cells, day 12 vs day 0", 
         x = "Clones", 
         title = "Clone depletion")+
    theme_Publication()
  ggsave(paste("barplot_foldchange_side", s, ".pdf", sep=""),plot = plot,device = "pdf", width = 40, height = 35, units = "cm")
}


PHATE_dataframe = merge(sel, PHATE_dataframe, by = "cellID")
plot=ggplot(PHATE_dataframe) +
  geom_point(aes(PHATE1, PHATE2, color = GATA4)) +
  labs(color="GATA4")+
  theme_Publication()
ggsave("PHATE_facet_GATA4.pdf",plot = plot,device = "pdf", width = 60, height = 30, units = "cm")


plot=ggplot(PHATE_dataframe) +
  geom_point(aes(PHATE1, PHATE2, color = TNNI1)) +
  labs(color="TNNI1")+
  theme_Publication()+
  scale_color_viridis_c()
ggsave("PHATE_facet_TNNI1.pdf",plot = plot,device = "pdf", width = 60, height = 30, units = "cm")

plot=ggplot(PHATE_dataframe) +
  geom_point(aes(PHATE1, PHATE2, color = TTN)) +
  labs(color="TTN")+
  theme_Publication()+
  scale_color_viridis_c()
ggsave("PHATE_facet_TTN.pdf",plot = plot,device = "pdf", width = 60, height = 30, units = "cm")


plot=ggplot(PHATE_dataframe) +
  geom_point(aes(PHATE1, PHATE2, color = MYL7)) +
  labs(color="MYL7")+
  theme_Publication()+
  scale_color_viridis_c()
ggsave("PHATE_facet_MYL7.pdf",plot = plot,device = "pdf", width = 60, height = 30, units = "cm")

