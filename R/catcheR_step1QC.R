#' @title catcheR_step1QC
#' @description For the analysis of the plasmids in their final form, i.e. where only the barcode is present
#' @param group, a character string. Two options: sudo or docker, depending to which group the user belongs
#' @param folder, a character string indicating the path of the working folder containing the input files
#' @param fastq.read1, a character string indicating the filename of read 1 fastq or fastq.gz containing barcodes sequencing
#' @param DIs, integer of the minimum number of diversity indexes (DIs, pseudo-unique reads) of the most represented shRNA matching to a given UCI-BC; in combination with "ratio", it selects UCI-BCs for which it is possible to reliably assign an shRNA. Default is 100
#' @param ratio, integer of the minimum ratio between the number of DIs of the most represented and second most represented shRNAs matching to a given UCI-BC. in combination with DI, it is employed to filter UCI-BCs for which it is possible to reliably assign an shRNA. Default is 10
#' @param plot.threshold, integer of the minimum number of DIs per UCI-BC for output plot
#' @param clones, a character string indicating the filename of txt file containing a newline separated list of clones of interest in the format of barcode_UCI
#' 
#' @author Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
#'
#' @return plots and stats
#'
#' @examples
#'\dontrun{
#'
#' catcheR_step1QC(
#'   group=("docker"),
#'   folder = "/20tb/ratto/catcheR/napoli_inter/", 
#'   fastq.read1 = "V350180591_L04_SPIKEIN_1.fq", 
#'   DIs = 100,
#'   ratio = 10,
#'   plot.threshold = 2000,
#'   clones = "clones.txt")
#'
#' @export


catcheR_step1QC <- function(
  group=c("docker","sudo"),
  folder, fastq.read1, DIs = 100, ratio = 10, plot.threshold = 2000, clones = NULL){ 
  
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
      "/home/plasmid_inter.sh",
      "/data/scratch",
      fastq.read1,
      plot.threshold,
      DIs,
      ratio,
      clones
    )
  )
  #executing the docker job
  #docker run --platform linux/amd64 -v /20tb/ratto/catcheR/test_CM5/:/data/scratch repbioinfo/catcher_barcode_pipeline /home/barcode_silencing_slicing.sh /data/scratch 1st2nd_hiPSC_CM_S5_R1_001.fastq 1st2nd_hiPSC_CM_S5_R2_001.fastq y12.csv GGCGCGTTCATCTGGGGGAGCCG 6 12
  #params <- paste("--cidfile ",folder,"/dockerID -v ",folder, ":/data/scratch -d docker.io/repbioinfo/catcher_barcode_pipeline_update /home/plasmid_inter.sh /data/scratch ", fastq.read1, " ", plot.threshold, " ", DIs, " ", ratio, " ",clones,sep="")
  setwd(home)
} 
