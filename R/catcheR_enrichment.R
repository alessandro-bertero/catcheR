#' @title catcheR_enrichment
#' @description Data loading step of the data analysis pipeline for the iPS2-seq methods from HEDGe lab.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files and results of previous analysis
#' @param file, a character string indicating the filename of the processed cds file resulting from catcheR_load
#' @param meta, string indicating the file name of the metadata file resulting from catcheR_load
#' @param timepoint, string indicating the starting timepoint to compare other timepoints to (e.g. day 0) (optional)
#' @param control_gene, string indicating the control gene to calculate statistics against in case internal control (no TET) is not present (e.g. "SCR") (optional)
#' @param min_cells_cluster, integer indicating the threshold for minimum number of cells per cluster to keep the cluster (optional) default: 70
#' @param min_cells_shRNA, integer indicating the threshold for minimum number of cells per shRNA to keep the cluster (optional) default: 40
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return plots and pvalues
#'
#' @examples
#'
#' catcheR_enrichment(group = "docker", folder = folder, file = "processed_cds.RData", meta = "cell_metadata.csv", timepoint = "PSC", control_gene = "SCR")
#'
#' @export

catcheR_enrichment <- function(
    group=c("docker","sudo"),
    folder, 
    file,
    meta,
    timepoint = NULL,
    control_gene = NULL,
    min_cells_cluster = 70,
    min_cells_shRNA = 40){
  
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
    image_name = "docker.io/repbioinfo/catcher_sc",
    volumes = list(
      c(folder, "/data/scratch/")
    ),
    additional_arguments = c(
      "Rscript /home/3_enrichment_depletion.R",
      "/data/scratch/",
      file,
      meta,
      timepoint,
      control_gene,
      min_cells_cluster,
      min_cells_shRNA
    )
  )
  #executing the docker job
  #docker run --platform linux/amd64 -v /20tb/ratto/catcheR/test_CM5/:/data/scratch repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_slicing.sh /data/scratch 1st2nd_hiPSC_CM_S5_R1_001.fastq 1st2nd_hiPSC_CM_S5_R2_001.fastq y12.csv GGCGCGTTCATCTGGGGGAGCCG 6 12
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data/scratch -d docker.io/repbioinfo/catcher_sc /home/1_load_data.R /data/scratch ", expression.matrix, " ", control_genes, " ", control_samples, " ", replicates, " ", sample_names, " ", resolution, " ", genes, " ", sep="")
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data -d docker.io/repbioinfo/desc.2018.01 Rscript /bin/top.R ", matrixName," ",format," ",separator, " ", logged, " ", threshold," ",type, sep="")
  setwd(home)
} 
