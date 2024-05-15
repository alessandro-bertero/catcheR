#' @title catcheR_design
#' @description Fucntion to design oligos for iPS2-seq methods from HEDGe lab.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files
#' @param sequences, a character string indicating the filename of the csv file containing sequences of interest and corresponding barcodes
#' @param gibson.five, a character string indicating the sequence of the gibson overhang to put before the shRNA, at the beginning of the construct
#' @param gibson.three, a character string indicating the sequence of the gibson overhang to put at the end of the construct
#' @param fixed, a character string indicating the fixed sequence to put after the shRNA in the contruct, between poly T and UCI
#' @param restriction.sites, a character string indicating the filename of the txt file containing sequences of the restriction sites to avoid in the sequence 
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return a file called "output.txt" containing the sequences of the constructs. Restriction sites are in lower case. 
#'
#' @examples
#'\dontrun{
#'
#' catcheR_design(group = "docker", 
#' folder = "/20tb/ratto/catcheR/oligo/", 
#' sequences = "shRNAs.csv", 
#' restriction.sites = "restr.txt")
#'
#' @export

catcheR_design <- function(
  group=c("docker","sudo"),
  folder, sequences, gibson.five = "AGTTCCCTATCAGTGATAGAGATCCC", gibson.three = "GTAGCTCGCTGATCAGC", fixed = "GTCGACATTTAAATGGCGCGCC", restriction.sites){
  
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
  
  #check  if scratch folder exist
  #  if (!file.exists(scratch.folder)){
  #    cat(paste("\nIt seems that the ",scratch.folder, " folder does not exist\n"))
  #    system("echo 3 > ExitStatusFile 2>&1")
  setwd(folder)
  #   return(3)
  # }
  #  tmp.folder <- gsub(":","-",gsub(" ","-",date()))
  scrat_tmp.folder=folder
  writeLines(scrat_tmp.folder,paste(folder,"/tempFolderID", sep=""))
  # cat("\ncreating a folder in scratch folder\n")
  #dir.create(file.path(scrat_tmp.folder))
  #preprocess matrix and copying files
  if (!file.exists(paste(folder,"/",sequences,sep=""))){
    cat(paste("\n It Seems that fastaq read1 file is not in ",folder,"\n"))
    system("echo 3 > ExitStatusFile 2>&1 &")
    setwd(folder)
    return(3)
  }
  #executing the docker job
  params <- paste("--cidfile ",folder,"/dockerID -v ",folder,":/data/scratch -d docker.io/repbioinfo/catcher_barcode_pipeline /home/oligo.sh /data/scratch ", sequences, " " ,gibson.five, " ", fixed, " ", gibson.three, " ", restriction.sites, " ", sep="")
  
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
