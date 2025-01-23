#' @title catcheR_pseudotime
#' @description Data loading step of the data analysis pipeline for the iPS2-seq methods from HEDGe lab.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files and results of previous analysis
#' @param cds, a character string indicating the filename of the annotated gene expression matrix csv resulting from catcher catch ("filtered_annotated_silencing_matrix_complete_all_samples.csv")
#' @param pseudotime, a character string indicating the pseudotime csv file
#' @param all, logical indicating if KS test needs to be run against each control gene and shRNA separately (instead on only against them together). Default = F
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return statistics on pseudotime
#'
#' @examples
#'\dontrun{
#'
#' catcheR_pseudotime(group="docker",folder="/3tb/data/ratto/aggr/test/", cds = "processed_cds.Rdata",pseudotime = "pseudotime.csv")
#'
#' @export

catcheR_pseudotime <- function(
  group=c("docker","sudo"),
  folder, 
  cds,
  pseudotime,
  all = F){
  
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
      c(folder, "/data/scratch")
    ),
    additional_arguments = c(
      "Rscript /home/2_pseudotime.R",
      "/data/scratch",
      cds,
      pseudotime,
      all
    )
  )
  #executing the docker job
  #docker run --platform linux/amd64 -v /20tb/ratto/catcheR/test_CM5/:/data/scratch repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_slicing.sh /data/scratch 1st2nd_hiPSC_CM_S5_R1_001.fastq 1st2nd_hiPSC_CM_S5_R2_001.fastq y12.csv GGCGCGTTCATCTGGGGGAGCCG 6 12
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data/scratch -d docker.io/repbioinfo/catcher_sc /home/4_modules.R /data/scratch ", cds, " ", sep="")
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data -d docker.io/repbioinfo/desc.2018.01 Rscript /bin/top.R ", matrixName," ",format," ",separator, " ", logged, " ", threshold," ",type, sep="")
  #resultRun <- runDocker(group=group, params=params)
  setwd(home)
} 
