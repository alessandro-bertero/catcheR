catcheR_inter_plasmids <- function(
  group=c("docker","sudo"),
  folder, fastq.read1, fastq.read2, clones = NULL){ 
  
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
  #docker run --platform linux/amd64 -v /20tb/ratto/catcheR/test_CM5/:/data/scratch repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_slicing.sh /data/scratch 1st2nd_hiPSC_CM_S5_R1_001.fastq 1st2nd_hiPSC_CM_S5_R2_001.fastq y12.csv GGCGCGTTCATCTGGGGGAGCCG 6 12
  params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data/scratch -d docker.io/repbioinfo/catcher_barcode_pipeline /home/plasmid_inter.sh /data/scratch ", fastq.read1, " ", fastq.read2, " ",clones, sep="")
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data -d docker.io/repbioinfo/desc.2018.01 Rscript /bin/top.R ", matrixName," ",format," ",separator, " ", logged, " ", threshold," ",type, sep="")
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
    tmp.run[length(tmp.run)+1] <- paste("catcheR_plasmids user run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("catcheR_plasmids system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("catcheR_plasmids elapsed run time mins ",ptm[3]/60, sep="")
    writeLines(tmp.run,"run.info")
  }else{
    tmp.run <- NULL
    tmp.run[1] <- paste("catcheR_plasmids run time mins ",ptm[1]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("catcheR_plasmids system run time mins ",ptm[2]/60, sep="")
    tmp.run[length(tmp.run)+1] <- paste("catcheR_plasmids elapsed run time mins ",ptm[3]/60, sep="")
    
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
