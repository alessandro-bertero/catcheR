# catcheR
## Clonality And Treatment Controlled sHrna Effect findeR

is a tool for iPS2-seq screening platform data from HEDGe lab (curated by Maria Luisa Ratto). 

See preprint DOI: http://dx.doi.org/10.2139/ssrn.4854180

Balmas, Elisa and Ratto, Maria Luisa and Snijders, Kirsten E. and Calogero, Raffaele and Mendjan, Sasha and Bertero, Alessandro, Single Cell Transcriptional Perturbome in Pluripotent Stem Cell Models.
Available at SSRN: https://ssrn.com/abstract=4854180 or http://dx.doi.org/10.2139/ssrn.4854180

See complete documentation at http://marialuisaratto.github.io/catcheRdocs

See associated protocol and CatcheR ![protocol](/DOCUMENTATION/catcheR.pdf) in DOCUMENTATION.

CatcheR analysis can be divided in two steps: single cell perturbation deconvolution and assignment, and perturbation effects statistical analysis (dimensionality reduction and clustering with Monocle3, pseudotime evaluation, enrichment / depletion analysis, genes modules analysis).

The folder manuscript_analysis_sc contains scripts used in the secondary analysis for the manuscript, curated by Elisa Balmas. The [Zenodo repository](https://doi.org/10.5281/zenodo.11085619) (https://doi.org/10.5281/zenodo.11085619) contains the scratch folders with the processed data and additional files to download and reproduce the analysis as in the paper with the scripts.

Processed and raw data files (fastq) are deposidet in the BioStudies repository:

iPS2-sci-seq - hiPSC-CMs (.fastq) E-MTAB-14102 https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14102
iPS2-10X-seq - monolayer differentiation (.fastq) E-MTAB-14065
iPS2-10X-seq - cardioids (.fastq) E-MTAB-14066
iPS2-seq - hiPSCs clonal drift (.fastq)	E-MTAB-15303
iPS2-multi-seq - hiPSCs (.fastq) E-MTAB-15332
iPS2-10X-seq - neural organoids (.fastq) E-MTAB-15308
iPS2-10X-seq - monoclonal cardioids (.fastq) E-MTAB-15307
iPS2-CITE-seq - polyclonal cardioids (.fastq) E-MTAB-15309



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

Library "rrundocker" is also required and can be installed with:

      install_github("Reproducible -Bioinformatics/rrundocker")
      library(rrundocker)

# Barcode pipeline for 10X
Starting from a single-cell gene expression matrix, the wrapper function catcheR_10Xcatch() produces a gene expression matrix where cell names contain annotations about the perturbation present in each cell. 
Inputs:
- read 1 fastq (or fastq.gz) containing the barcodes sequencing library (from cellranger mkfastq)
- read 2 fastq (or fastq.gz) containing the barcodes sequencing library (from cellranger mkfastq)
- gene expression matrix in csv format (transform the H5 cellranger count matrix to a csv file within cellranger or with rCASC h5tocsv (https://kendomaniac.github.io/rCASC/reference/h5tocsv.html)
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

To allow reproducibility of the scripts outside CatcheR we recommend pulling the Docker image below from the terminal:

    dockered

    docker pull hedgelab/rstudio-hedgelab

    docker run -d -itv /the/folder/you/want/to/share:/scratch --privileged=true -p 8787:8787 \\
    -e PASSWORD=<your_password> -e USER=rstudio --name=NAME_CONTAINER hedgelab/rstudio-hedgelab:iPS2seq

    docker exec -idt NAME_CONTAINER rstudio-server start

Replace <your_password> with your desired password. This command maps port 8787 on your host machine to port 8787 in the container, allowing you to access RStudio via your browser at http://localhost:8787. The USER=rstudio part ensures you'll log in as the rstudio user, and the PASSWORD variable sets the password you'll use to log in.
Then use Rstudio through a browser

