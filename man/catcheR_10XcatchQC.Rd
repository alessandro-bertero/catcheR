% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/catcheR_10XQC.R
\name{catcheR_10XcatchQC}
\alias{catcheR_10XcatchQC}
\title{catcheR_10XcatchQC}
\usage{
catcheR_10XcatchQC(
  group = c("docker", "sudo"),
  folder,
  reference = "GGCGCGTTCATCTGGGGGAGCCG",
  mode = "bimodal"
)
}
\arguments{
\item{group, }{a character string. Two options: sudo or docker, depending to which group the user belongs}

\item{folder, }{a character string indicating the path of the working folder containing the input files}

\item{reference, }{a character string indicating the sequence to identify reads containing barcodes. Should be found at at the beginning of read2. Use reverse complement!}

\item{mode, }{a character string. Two options: "bimodal" or "noise". To evaluate a threshold number of UMIs to consider a UCI valid there are 2 options: "bimodal" (default) which sets the threshold at the valley of the UMIxUCI distribution, or "noise", which sets the threshold at 1.35 * number of UCI supported by a single UMI.}
}
\value{
QC plots and intermediate data
}
\description{
Explorative step of the data analysis pipeline for the iPS2-seq methods from HEDGe lab.
}
\author{
Maria Luisa Ratto, marialuisa.ratto [at] unito [dot] it, UNITO
}
