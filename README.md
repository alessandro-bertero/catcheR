# HEDGe_iPS2-seq
Pipelines and codes relative to the hPS2-seq screening platform
Data analysis pipeline from HEDGe lab. 

See preprint DOI: http://dx.doi.org/10.2139/ssrn.4854180

Balmas, Elisa and Ratto, Maria Luisa and Snijders, Kirsten E. and Calogero, Raffaele and Mendjan, Sasha and Bertero, Alessandro, Single Cell Transcriptional Perturbome in Pluripotent Stem Cell Models.
Available at SSRN: https://ssrn.com/abstract=4854180 or http://dx.doi.org/10.2139/ssrn.4854180

See associated protocol and CatcheR ![protocol](/DOCUMENTATION/catcheR.pdf) in DOCUMENTATION.

Repository contains 2 sets of scripts: 

- CatcheR barcode pipeline curated by Maria Luisa Ratto: first step of the analysis, to assign perturbation to single cells.
  There are two different versions: one for 10X data and one for double indexing sci-RNAseq. An overview of the workflow can be found in the DOCUMENTATION folder.
  
- single cell analysis curated by Elisa Balmas: Evaluate the perturbation effect of a gene or shRNA at the clonal or population level.
  Clustering is done with Monocle 3 and customized statistical methods have been employed to assess:
  (1) Cluster enrichment variation due to a perturbation; (2) Changes in Pseudotime or Module gene expression associated with a perturbation.
  Zenodo repository https://doi.org/10.5281/zenodo.11085619 contains the scratch folder to reproduce the analysis in the paper with the ![scripts](
/single_cell_analysis)

# CatcheR installation
Use the "install_github" function in the "devtools" package.

    library(devtools) 
    install_github("alessandro-bertero/catcheR")
    library(catcheR)

Alternatively, the repository can be download manually and loaded as follows: 

    library(devtools);
    load_all("."); # Working directory should be in the package catcheR

## Prerequisites
The following R functions require docker, since each of them opens a docker, computes the analysis inside of it to ensure reproducibility and then closes it. 

# Barcode pipeline for 10X
Starting from a single-cell gene expression matrix, the wrapper function catcheR_10Xcatch() produces a gene expression matrix where cell names contain annotations about the perturbation present in each cell. 
Inputs:
- read 1 fastq (or fastq.gz) containing the barcodes sequencing library (from cellranger mkfastq)
- read 2 fastq (or fastq.gz) containing the barcodes sequencing library (from cellranger mkfastq)
- gene expression matrix in csv format (transform the cellranger count matrix to a csv file with XXXXXXX)
- a file called rc_barcodes_gene.csv, containing the association between barcodes and shRNAs, comma separated. Use the reverse complement of the shRNA barcode. E.g.
  
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

## catcheR_10Xcatch()

    catcheR_10Xcatch(group=c("docker","sudo"),folder, fastq.read1, fastq.read2, expression.matrix, reference = "GGCGCGTTCATCTGGGGGAGCCG", UCI.length = 6, threads = 2, percentage = 15, mode = "bimodal")


Wrapper function arguments: 

  - group: a character string. Two options: sudo or docker, depending to which group the user belongs. 
  - folder: a character string indicating the path of the working folder containing the input files.
  - fastq.read1: a character string indicating the filename of read 1 fastq (or fastq.gz). This is the read 1 of the barcode library: this is the library created by enrichment of the sequence of interest with a specific iPS2seq primer and it is needed to deconvolute the perturbation. 
  - fastq.read2: a character string indicating the filename of read 2 fastq (or fastq.gz). This is the read 2 of the barcode library: this is the library created by enrichment of the sequence of interest with a specific iPS2seq primer and it is needed to deconvolute the perturbation. 
  - expression.matrix: a character string indicating the filename of the gene expression matrix file (csv format). 
  - reference: a character string indicating the reference sequence to identify reads containing the barcodes. This should be found at the beginning of read2. Use reverse complement! The default is the Tet repressor sequence used for the enrichment of the iPS2-seq construct (see the sequence in the example). 
  - UCI.length: integer indicating the length of the Unique Clonal Identifier. Should be found on read2 after the reference. Default is 6, as in the iPS2-seq construct.
  - threads: integer number of threads to be used for parallelization. The default is 2, but more threads are recommended to increase performance. 
  - percentage: integer threshold of the percentage of UMIs supporting a UCI over the total UMIs supporting UCIs in the same cell to consider the UCI valid. The suggested default is 15.
  - mode: a character string. Two options: "bimodal" or "noise". To evaluate a threshold number of UMIs to consider a UCI valid there are 2 options: "bimodal" (default) which sets the threshold at the valley of the UMIxUCI distribution, or "noise", which sets the threshold at 1.35 * number of UCI supported by a single UMI.

Example

    folder = "/20tb/ratto/catcheR/test_02_7/"
    catcheR_10Xcatch(group = "docker", 
                 folder = folder, 
                 fastq.read1 = list.files(folder, pattern = "R1"), 
                 fastq.read2 = list.files(folder, pattern = "R2"), 
                 expression.matrix = "matrix.csv", 
                 reference = "GGCGCGTTCATCTGGGGGAGCCG",
                 UCI.length = 6
                 threads = 12, 
                 percentage = 15,
                 mode = "noise")

Note. Check the UMIxUCI plots and the percentage_of_UMIxUCI_dist plots showing the distribution of shRNA and assess the noise of the dataset. Then re-run the previous analysis when, for example, to change the thresholds to custom values or the mode for threshold identification. In this case, run the catcheR_explorative function.

  catcheR_10XcatchQC(
                group = "docker", 
                folder = "/20tb/ratto/catcheR/test_CM5/", 
                reference = "GGCGCGTTCATCTGGGGGAGCCG", 
                mode = "noise")
                
  catcheR_filtercatch(group = "docker", 
                       folder = "/20tb/ratto/catcheR/sci_8/", 
                       expression.matrix = "exp_mat.csv", 
                       UMI.count = 5, 
                       percentage = 15)

