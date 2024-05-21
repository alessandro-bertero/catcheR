#' @title catcheR_10Xnocatch
#' @description Empty selection step of the data analysis pipeline for the iPS2-seq methods from HEDGe lab.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files
#' @param expression.matrix, a character string indicating the filename of the gene expression matrix csv
#' @param threshold, an integer indicating the minimum number of UMIs associated to the empty reference to consider a cell empty
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
    threshold){
  
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
  #docker run --platform linux/amd64 -v /20tb/ratto/catcheR/test_CM5/:/data/scratch repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_slicing.sh /data/scratch 1st2nd_hiPSC_CM_S5_R1_001.fastq 1st2nd_hiPSC_CM_S5_R2_001.fastq y12.csv GGCGCGTTCATCTGGGGGAGCCG 6 12
  params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data/scratch -d docker.io/repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_empty_selection.R /data/scratch ", expression.matrix, " ", threshold, sep="")
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data -d docker.io/repbioinfo/desc.2018.01 Rscript /bin/top.R ", matrixName," ",format," ",separator, " ", logged, " ", threshold," ",type, sep="")
  resultRun <- runDocker(group=group, params=params)
  
  #waiting for the end of the container work
  if(resultRun==0){
    cat("\nData filtering is finished\n")
  }
  
  #saving log and removing docker container
  container.id <- readLines(paste(folder,"/dockerID", sep=""), warn = FALSE)
  #system(paste("docker logs ", substr(container.id,1,12), " >& ",folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker logs ", substr(container.id,1,12), " > ",folder,"/", substr(container.id,1,12),".log 2>&1", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  
  
  #removing temporary folder
  cat("\n\nRemoving the temporary file ....\n")
  # system(paste("rm -R ",scrat_tmp.folder))
  #file.remove(paste0(folder,"out.info"))
  file.remove(paste0(folder,"dockerID"))
  #file.remove(paste0(folder,"tempFolderID"))
  #system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
} 
