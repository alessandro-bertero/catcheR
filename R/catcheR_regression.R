#' @title catcheR_regression
#' @description Finds genes whose expression is correlated with a specific perturbation group of cell. Part of the data analysis pipeline for the iPS2-seq methods from HEDGe lab.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files and results of previous analysis
#' @param cds, a character string indicating the filename of the annotated gene expression matrix csv resulting from catcher catch ("filtered_annotated_silencing_matrix_complete_all_samples.csv")
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return tables of genes correlated with a specific perturbation group
#'
#' @examples
#'\dontrun{
#'
#' catcheR_regression(group="docker",folder="/3tb/data/ratto/aggr/test/", cds = "processed_cds.Rdata")
#'
#' @export

catcheR_regression <- function(
    group=c("docker","sudo"),
    folder, 
    cds){
  
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
  # test <- dockerTest()
  # if(!test){
  #   cat("\nERROR: Docker seems not to be installed in your system\n")
  #   #system("echo 10 > ExitStatusFile")
  #   setwd(home)
  #   return(10)
  # }

  #executing the docker job
  run_in_docker(
    image_name = "docker.io/repbioinfo/catcher_sc",
    volumes = list(
      c(folder, "/data/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/5_regression.R",
      "/data/scratch",
      cds
    )
  )
  setwd(home)
} 
