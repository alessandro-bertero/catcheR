#' @title catcheR_scicatch
#' @description Data analysis pipeline for the iPS2-seq methods from HEDGe lab, using combinatorial single cell indexing sequencing.
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files. Must contain a folder "fastq" with demultiplexed fastq files
#' @param expression.matrix, a character string indicating the filename of the gene expression matrix csv
#' @param reference, a character string indicating the sequence to identify reads containing barcodes. Should be found at at the beginning of read2. Use reverse complement!
#' @param UCI.legth, integer indicating the length of unique clonal identifier. Should be found after reference on read2
#' @param threads, integer number of threads to be used for parallelization 
#' @param percentage, integer threshold of percentage of UMIs supporting a UCI over total UMIs supporting UCIs in the same cell, to consider the UCI valid. Suggested default is 15.
#' @param ratio, select valid cells = only 1 true UCI and top UCI UMIs / second top UCI UMIs > ratio
#' @param mode, a character string. Two options: "bimodal" or "noise". To evaluate a threshold number of UMIs to consider a UCI valid there are 2 options: "bimodal" (default) which sets the threshold at the valley of the UMIxUCI distribution, or "noise", which sets the threshold at 1.35 * number of UCI supported by a single UMI.
#' @param x, an integer indicating the x limit to zoom the UMIxUCI plot. Default is 100
#' @param y, an integer indicating the y limit to zoom the UMIxUCI plot. Default is 400
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return a gene expression matrix where cell names contains annotations about the perturbation present in each cell + QC plots
#'
#' @examples
#'\dontrun{
#'
#' folder = "/20tb/ratto/catcheR/sci_8/"
#' catcheR_scicatch(group = "docker", 
#'                      folder = folder, 
#'                      expression.matrix = "exp_mat.csv", 
#'                      threads = 12,
#'                      ratio = 5, 
#'                      mode = "noise")
#'
#' @export

catcheR_scicatch <- function(
  group=c("docker","sudo"),
  folder, expression.matrix, reference = "GGCGCGTTCATCTGGGGGAGCCG", UCI.length = 6, threads = 2, percentage = 15, ratio = 5, mode = "bimodal", x = 100, y = 400){
  
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
    print("echo 10 > ExitStatusFile")
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
    print("echo 3 > ExitStatusFile 2>&1 &")
    setwd(folder)
    return(3)
  }
  
  #executing the docker job
  run_in_docker(
    image_name = "docker.io/repbioinfo/catcher_barcode_pipeline",
    volumes = list(
      c(folder, "/data/scratch/")
    ),
    additional_arguments = c(
      "/home/sci_barcode_silencing_slicing.sh",
      "/data/scratch/",
      expression.matrix,
      reference,
      UCI.length,
      threads,
      percentage,
      mode
    )
  )
  
  #waiting for the end of the container work
  # if(resultRun==0){
  #   cat("\nData filtering is finished\n")
  # }
  
  #saving log and removing docker container
  #container.id <- readLines(paste(folder,"/dockerID", sep=""), warn = FALSE)
  #system(paste("docker logs ", substr(container.id,1,12), " >& ",folder,"/", substr(container.id,1,12),".log", sep=""))
  #system(paste("docker logs ", substr(container.id,1,12), " > ",folder,"/", substr(container.id,1,12),".log 2>&1", sep=""))
  #system(paste("docker rm ", container.id, sep=""))
  
  
  #removing temporary folder
  #cat("\n\nRemoving the temporary file ....\n")
  # system(paste("rm -R ",scrat_tmp.folder))
  #file.remove(paste0(folder,"out.info"))
  #file.remove(paste0(folder,"dockerID"))
  #file.remove(paste0(folder,"tempFolderID"))
  #system(paste("cp ",paste(path.package(package="rCASC"),"containers/containers.txt",sep="/")," ",data.folder, sep=""))
  setwd(home)
} 

