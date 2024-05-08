#' @title catcheR_sci_tomatrix
#' @description From fastq to gene expression matrix for SCI data.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files
#' @param sample.name, a character string indicating the name of the experiment
#' @param UMI.cutoff, integer indicating the minimum number of UMI per cell to consider the cell valid
#'
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return a gene expression matrix + QC plots and stats
#'
#' @examples
#'\dontrun{
#'
#'catcheR_sci_tomatrix(
#'  group="docker",
#'  folder="/20tb/ratto/catcheR/tomatrix/", 
#'  sample.name="H001AS8", UMI.cutoff=500)
#'
#' @export

catcheR_sci_tomatrix <- function(
  group=c("docker","sudo"),
  folder, sample.name, UMI.cutoff){
  
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
  if (!dir.exists(paste(folder,"/fastq",sep=""))){
    cat(paste("\n It Seems that fastaq read1 file is not in ",folder,"\n"))
    system("echo 3 > ExitStatusFile 2>&1 &")
    setwd(folder)
    return(3)
  }
  
  #executing the docker job
  params <- paste("--cidfile ",folder,"/dockerID -v ",folder,":/data/scratch -d repbioinfo/sci_tomatrix /home/tomatrix.sh ", sample.name, " " ,UMI.cutoff, " ", sep="")
  
  resultRun <- runDocker(group=group, params=params)
  
  #waiting for the end of the container work
  if(resultRun==0){
    cat("\nData filtering is finished\n")
  }
  
  #running time 2
  ptm <- proc.time() - ptm
  dir <- dir(folder)
  dir <- dir[grep("run.info",dir)]
  if(length(dir)>0){
    con <- file("run.info", "r")
    tmp.run <- readLines(con)
    close(con)
    tmp.run[length(tmp.run)+1] <- paste("catcheR_barcodes user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("catcheR_barcodes system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("catcheR_barcodes elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("catcheR_barcodes run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("catcheR_barcodes system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("catcheR_barcodes elapsed run time mins ",ptm[3]/60, sep="")
    
    writeLines(tmp.run,"run.info")
  }
  
  #saving log and removing docker container
  container.id <- readLines(paste(folder,"/dockerID", sep=""), warn = FALSE)
  #system(paste("docker logs ", substr(container.id,1,12), " &> ",folder,"/", substr(container.id,1,12),".log", sep=""))
  system(paste("docker logs ", substr(container.id,1,12), " > ", folder,"/", substr(container.id,1,12),".log", " 2>&1", sep=""))
  system(paste("docker rm ", container.id, sep=""))
  #removing temporary folder
  #  cat("\n\nRemoving the temporary file ....\n")
  #  system(paste("rm -fR ",scrat_tmp.folder))
  system("rm -fR out.info")
  system("rm -fR dockerID")
  system("rm  -fR tempFolderID")
  system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",folder, sep=""))
  setwd(home)
} 

