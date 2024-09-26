#' @title catcheR_step2QC
#' @description For the analysis of the plasmids in their final form, i.e. where only the barcode is present
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files
#' @param fastq.read1, a character string indicating the filename of read 1 fastq or fastq.gz containing barcodes sequencing
#' @param DIs, integer, a minimum of reads associated to a clone to show it in plots of clones, for visualization purposes.
#' @param clones, a character string indicating the filename of txt file containing a newline separated list of clones of interest in the format of barcode_UCI
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return plots and stats
#'
#' @examples
#'\dontrun{
#'
#' catcheR_step2QC(
#'   group=("docker"),
#'   folder = "/20tb/ratto/catcheR/napoli_final/", 
#'   fastq.read1 = "V350180591_L04_SPIKEIN_1.fq",
#'   DIs = 1000, 
#'   clones = "clones.txt")
#'
#' @export

catcheR_step2QC <- function(
  group=c("docker","sudo"),
  folder, fastq.read1, DIs = 1000, clones = NULL){ 
  
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
      "/home/plasmid_final.sh",
      "/data/scratch",
      fastq.read1,
      DIs,
      clones
    )
  )
  #executing the docker job
  #docker run --platform linux/amd64 -v /20tb/ratto/catcheR/test_CM5/:/data/scratch repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_slicing.sh /data/scratch 1st2nd_hiPSC_CM_S5_R1_001.fastq 1st2nd_hiPSC_CM_S5_R2_001.fastq y12.csv GGCGCGTTCATCTGGGGGAGCCG 6 12
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data/scratch -d docker.io/repbioinfo/catcher_barcode_pipeline_update /home/plasmid_final.sh /data/scratch ", fastq.read1, " ", DIs, " ",clones, sep="")
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data -d docker.io/repbioinfo/desc.2018.01 Rscript /bin/top.R ", matrixName," ",format," ",separator, " ", logged, " ", threshold," ",type, sep="")
  setwd(home)
} 
