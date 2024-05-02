# HEDGe_iPS2-seq
Pipelines and codes relative to the hPS2-seq screening platform
Data analysis pipeline from HEDGe lab. 

CatcheR functions curated by Maria Luisa Ratto
single cell RNAseq analysis curated by Elisa Balmas

See pubblication.

Repository contains 2 sets of scripts: 

- barcode pipeline: first step of the analysis, to assign perturbation to single cells. There are two different versions: one for 10X data and one for double indexing sci.
- single cell analysis: perturbation effect evaluation after clustering

# Prerequisites
The following R functions require docker, since each of them opens a docker, computes the analysis inside of it to ensure reproducibility and then closes it. 

# Barcode pipeline for 10X
Starting from single-cell gene expression matrix, the wrapper function catcheR_barcodes() produces a gene expression matrix where cell names contains annotations about the perturbation present in each cell. 
Inputs:
- read 1 fastq (or fastq.gz) containing barcodes sequencing
- read 2 fastq (or fastq.gz) containing barcodes sequencing
- gene expression matrix in csv format
- a file called rc_barcodes_gene.csv, containing the association between barcodes and shRNAs, comma separated. Use reverse complement! E.g.
  
      CTTCTTTC,CHD7.1
      GTACTCAA,CHD7.2
      TTCGTCAT,CHD7.3
      ATCTCTCA,CHD7.4
      AGGCGAGA,GATA4.1
      TCTTCAGC,GATA4.2
      CACAGATA,GATA4.3
      ACAATCTC,KMT2D.1
      TCGGAGCA,KMT2D.2
      ATCCGTAT,KMT2D.3
      GAGACCAT,KMT2D.4
      CTGCAGTA,NKX2.5.1
      CGTGATGC,NKX2.5.2
      TGATTCAG,NKX2.5.3
      CAAGAGCC,SMAD2.1
      AACCGGAG,SMAD2.2
      GAAGTTCG,SMAD2.3
      GCGGAACT,SMAD2.4
      TGGAACTG,B2M
      AGTAGGCT,EGFP
      GCCTGTGT,SCR
- (optional) colors.csv, a comma separated file indicating the colors to use for plotting of each shRNA data. E.g.

      CHD7.1,#AA0DFE
      CHD7.2,#3283FE
      CHD7.3,#85660D
      CHD7.4,#782AB6
      GATA4.1a,#565656
      ...

## catcheR_barcodes()

    catcheR_barcodes(group=c("docker","sudo"),folder, fastq.read1, fastq.read2, expression.matrix, reference = "GGCGCGTTCATCTGGGGGAGCCG", UCI.length = 6, threads = 2, percentage = 15, mode = "bimodal")



Wrapper function arguments: 

  - group: a character string. Two options: sudo or docker, depending to which group the user belongs. For docker running
  - folder: a character string indicating the path of the working folder containing the input files
  - fastq.read1: a character string indicating the filename of read 1 fastq (or fastq.gz). This library must have enrichment of the sequence of interest for perturbation deconvolution. 
  - fastq.read2: a character string indicating the filename of read 2 fastq (or fastq.gz). This library must have enrichment of the sequence of interest for perturbation deconvolution. 
  - expression.matrix: a character string indicating the filename of the gene expression matrix file (csv format). 
  - reference: a character string indicating the sequence to identify reads containing barcodes. Should be found at at the beginning of read2. Use reverse complement! Default is the sequence used in the iPS2-seq construct. 
  - UCI.length: integer indicating the length of Unique Clonal Identifier. Should be found on read2 after the reference. Default is 6, as in the iPS2-seq construct.
  - threads: integer number of threads to be used for parallelization
  - percentage: integer threshold of percentage of UMIs supporting a UCI over total UMIs supporting UCIs in the same cell, to consider the UCI valid. Suggested default is 15.
  - mode: a character string. Two options: "bimodal" or "noise". To evaluate a threshold number of UMIs to consider a UCI valid there are 2 options: "bimodal" (default) which sets the threshold at the valley of the UMIxUCI distribution, or "noise", which sets the threshold at 1.35 * number of UCI supported by a single UMI.

Example

    folder = "/20tb/ratto/catcheR/test_02_7/"
    catcheR_barcodes(group = "docker", 
                 folder = folder, 
                 fastq.read1 = list.files(folder, pattern = "R1"), 
                 fastq.read2 = list.files(folder, pattern = "R2"), 
                 expression.matrix = "matrix.csv", 
                 reference = "GGCGCGTTCATCTGGGGGAGCCG",
                 UCI.length = 6
                 threads = 12, 
                 percentage = 15,
                 mode = "noise")

