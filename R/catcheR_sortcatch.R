#' @title catcheR_sortcatch
#' @description Modify annotated gene expression matrix from catcheR_10Xcatch or catcheR_scicatch based on the barcode swaps identified with catcheR_step1QC.
#' @param folder, a character string indicating the path of the working folder containing the input files
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param expression.matrix, a character string indicating the filename of the csv file of the annotated gene expression produced by catcheR_10Xcatch or catcheR_scicatch i.e. silencing_matrix.csv
#' @param swaps, a character string indicating the filename of the txt file indicating swaps to be corrected, output of catcheR_step1QC, i.e. ”reliable_clones_swaps.txt"
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return Updated gene expression matrix called ”silencing_matrix_updated.csv"
#'
#' @examples
#'\dontrun{
#'
#' catcheR_sortcatch(group="docker",folder = "/home/user/Documents/reassign/", expression.matrix="silencing_matrix.csv", swaps="reliable_clones_swaps_50.txt")
#'
#' @export


catcheR_sortcatch <- function(
    group=c("docker","sudo"),
    folder, expression.matrix, swaps){ #noise or bimodal
  
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
      c(folder, "/data/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/sortcatch.R",
      "/data/scratch",
      expression.matrix,
      swaps
    )
  )
  #executing the docker job
  #docker run --platform linux/amd64 -v /20tb/ratto/catcheR/test_CM5/:/data/scratch repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_slicing.sh /data/scratch 1st2nd_hiPSC_CM_S5_R1_001.fastq 1st2nd_hiPSC_CM_S5_R2_001.fastq y12.csv GGCGCGTTCATCTGGGGGAGCCG 6 12
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data/scratch -d docker.io/repbioinfo/catcher_barcode_pipeline Rscript /home/sortcatch.R /data/scratch/ ", expression.matrix, " ", swaps, sep="")
  setwd(home)
} 
