#' @title catcheR_10Xnocatch
#' @description Empty selection step of the data analysis pipeline for the iPS2-seq methods from HEDGe lab.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files
#' @param expression.matrix, a character string indicating the filename of the gene expression matrix csv
#' @param threshold, an integer indicating the minimum number of UMIs associated to the empty reference to consider a cell empty
#' @param samples, integer indicating the number of different samples present in the sample (aggregated with cell ranger aggr)
#' @param reference, the sequence indicating the empty plasmid
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return a gene expression matrix with empty cells
#'
#' @examples
#'\dontrun{
#'
#' catcheR_10Xnocatch(group = "docker", folder = folder, expression.matrix = "matrix.csv", threshold = 10)
#'
#' @export

catcheR_10Xnocatch <- function(
    group=c("docker","sudo"),
    folder, 
    expression.matrix,
    threshold,
    samples = 1,
    reference = "TACGCGTTCATCTGGGGGAGCCG"){
  
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
    image_name = "docker.io/repbioinfo/catcher_barcode_pipeline_update",
    volumes = list(
      c(folder, "/data/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/barcode_silencing_empty_selection.R",
      "/data/scratch",
      expression.matrix,
      threshold,
      samples,
      reference
    )
  )
  setwd(home)
} 
