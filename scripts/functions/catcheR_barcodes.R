catcheR_barcodes <- function(
  group=c("docker","sudo"),
  folder, fastq.read1, fastq.read2, expression.matrix, reference = "GGCGCGTTCATCTGGGGGAGCCG", UCI.length = 6, threads = 2, percentage = 15, mode = "bimodal"){
  
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
  if (!file.exists(paste(folder,"/",fastq.read1,sep=""))){
    cat(paste("\n It Seems that fastaq read1 file is not in ",folder,"\n"))
    system("echo 3 > ExitStatusFile 2>&1 &")
    setwd(folder)
    return(3)
  }
  if (!file.exists(paste(folder,"/",fastq.read2,sep=""))){
    cat(paste("\n It Seems that fastaq read2 file is not in ",folder,"\n"))
    system("echo 3 > ExitStatusFile 2>&1 &")
    setwd(folder)
    return(3)
  }
  #executing the docker job
  params <- paste("--cidfile ",folder,"/dockerID -v ",folder,":/data/scratch -d docker.io/repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_slicing.sh /data/scratch ", fastq.read1, " " ,fastq.read2, " ", expression.matrix, " ", reference, " ", UCI.length, " ", threads, " ", percentage, " ", mode, " ", sep="")

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
