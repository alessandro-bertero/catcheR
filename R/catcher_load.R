#' @title catcheR_load
#' @description Data loading step of the data analysis pipeline for the iPS2-seq methods from HEDGe lab.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files and results of previous analysis
#' @param expression.matrix, a character string indicating the filename of the annotated gene expression matrix csv resulting from catcher catch ("filtered_annotated_silencing_matrix_complete_all_samples.csv")
#' @param control_genes, string indicating the file name of a file listing the genes which are used as controls
#' @param control_samples, string indicating the file name of a file listing the samples which are used as controls (no TET samples)
#' @param replicates, string indicating the file name of a file listing the name of the replicates and their order associated to samples (used for batch correction)
#' @param sample_names, string indicating the file name of a csv file listing the samples and the associated name
#' @param resolution, numeric indicating the resolution to be used for Monocle clustering 
#' @param genes, string indicating the file name of a file listing the genes of interest whose expressin will be plotted on the UMAP
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return a gene expression matrix with empty cells
#'
#' @examples
#'
#' catcheR_load(group="docker",folder="/3tb/data/ratto/aggr/test/", expression.matrix = "annotated_silencing_matrix_complete_all_samples.csv",control_genes =  "controls.txt",control_samples = "noTET.txt",replicates = "replicates.txt",sample_names = "samples.csv",resolution = 8e-4,genes = "genelist.txt")
#'
#' @export

catcheR_load <- function(
    group=c("docker","sudo"),
    folder, 
    expression.matrix,
    control_genes,
    control_samples = "NotSpecified",
    replicates = "NotSpecified",
    sample_names,
    resolution = 8e-4,
    genes = "NotSpecified"){
  
  #running time 1
  ptm <- proc.time()
  #setting the folder as working folder
  if (!file.exists(folder)){
    cat(paste("\nIt seems that the ",folder, " folder does not exist\n"))
    return(2)
  }
  
  #storing the position of the home folder
  home <- getwd()
  setwd(folder)
  #initialize status
  #system("echo 0 > ExitStatusFile")
  
  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    #system("echo 10 > ExitStatusFile")
    setwd(home)
    return(10)
  }
  # if(!test){
  #   cat("\nERROR: Docker seems not to be installed in your system\n")
  #   system("echo 10 >& ExitStatusFile")
  #   setwd(home)
  #   return(10)
  # }
  #if(logged){logged="TRUE"}else{logged="FALSE"}
  #executing the docker job
  run_in_docker(
    image_name = "docker.io/repbioinfo/catcher_sc_loc",
    volumes = list(
      c(folder, "/data/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/1_load_data.R",
      "/data/scratch",
      expression.matrix,
      control_genes,
      control_samples, 
      replicates,
      sample_names,
      resolution,
      genes
    )
  )
  #executing the docker job
  #docker run --platform linux/amd64 -v /20tb/ratto/catcheR/test_CM5/:/data/scratch repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_slicing.sh /data/scratch 1st2nd_hiPSC_CM_S5_R1_001.fastq 1st2nd_hiPSC_CM_S5_R2_001.fastq y12.csv GGCGCGTTCATCTGGGGGAGCCG 6 12
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data/scratch -d docker.io/repbioinfo/catcher_sc /home/1_load_data.R /data/scratch ", expression.matrix, " ", control_genes, " ", control_samples, " ", replicates, " ", sample_names, " ", resolution, " ", genes, " ", sep="")
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data -d docker.io/repbioinfo/desc.2018.01 Rscript /bin/top.R ", matrixName," ",format," ",separator, " ", logged, " ", threshold," ",type, sep="")
  setwd(home)
} 
