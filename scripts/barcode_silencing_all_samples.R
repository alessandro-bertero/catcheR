#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)
library(ggplot2)

#Set working direcotry
args = commandArgs(trailingOnly=TRUE)
dir = args[1]
#2 samples
samples = args[2]


# Initialize the matrix with the first file
matrix = read.table(paste(dir, "Results_1/silencing_matrix_complete_1.csv", sep=""), sep = ",", header = TRUE, row.names = 1)
names(matrix) = paste(names(matrix), 1, sep = "_")

# Loop through the remaining files and bind them
for (s in 2:samples) {
  m = read.table(paste(dir, "Results_", s, "/silencing_matrix_complete_", s, ".csv", sep=""), sep = ",", header = TRUE, row.names = 1)
  names(m) = paste(names(m), s, sep = "_")
  matrix = cbind(matrix, m)
}

write.csv(matrix, paste(dir,"/silencing_matrix_complete_all_samples.csv",sep=""))
