#' @title catcheR_scicatchQC
#' @description Explorative step of the data analysis pipeline for the iPS2-seq methods from HEDGe lab with SCI.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files
#' @param expression.matrix, a character string indicating the filename of the gene expression matrix csv
#' @param reference, a character string indicating the sequence to identify reads containing barcodes. Should be found at at the beginning of read2. Use reverse complement!
#' @param mode, a character string. Two options: "bimodal" or "noise". To evaluate a threshold number of UMIs to consider a UCI valid there are 2 options: "bimodal" (default) which sets the threshold at the valley of the UMIxUCI distribution, or "noise", which sets the threshold at 1.35 * number of UCI supported by a single UMI.
#' @param x, an integer indicating the x limit to zoom the UMIxUCI plot. Default is 100
#' @param y, an integer indicating the y limit to zoom the UMIxUCI plot. Default is 400
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return QC plots and intermediate data
#'
#' @examples
#'\dontrun{
#'
#' catcheR_scicatchQC(group = "docker", 
#'     folder = "/20tb/ratto/catcheR/sci_8/", 
#'     expression.matrix = "exp_mat.csv", 
#'     mode = "noise")
#'
#' @export


catcheR_scicatchQC <- function(
  group=c("docker","sudo"),
  folder, expression.matrix, reference = "GGCGCGTTCATCTGGGGGAGCCG", mode = "bimodal", x = 100, y = 400){ #noise or bimodal
  
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
  system("echo 0 > ExitStatusFile")
  
  #testing if docker is running
  test <- dockerTest()
  if(!test){
    cat("\nERROR: Docker seems not to be installed in your system\n")
    system("echo 10 > ExitStatusFile")
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
    image_name = "docker.io/repbioinfo/catcher_barcode_pipeline_update",
    volumes = list(
      c(folder, "/data/scratch/")
    ),
    additional_arguments = c(
      "Rscript /home/sci_barcode_silencing_explorative_analysis.R",
      "/data/scratch/",
      reference,
      mode,
      expression.matrix, 
      x, 
      y
    )
  )
  #executing the docker job
  #docker run --platform linux/amd64 -v /20tb/ratto/catcheR/test_CM5/:/data/scratch repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_slicing.sh /data/scratch 1st2nd_hiPSC_CM_S5_R1_001.fastq 1st2nd_hiPSC_CM_S5_R2_001.fastq y12.csv GGCGCGTTCATCTGGGGGAGCCG 6 12
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data/scratch -d docker.io/repbioinfo/catcher_barcode_pipeline_update /home/sci_barcode_silencing_explorative_analysis.R /data/scratch ", reference, " ", mode, " ", expression.matrix, " ", x, " ", y, sep="")
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data -d docker.io/repbioinfo/desc.2018.01 Rscript /bin/top.R ", matrixName," ",format," ",separator, " ", logged, " ", threshold," ",type, sep="")
  setwd(home)
} 
