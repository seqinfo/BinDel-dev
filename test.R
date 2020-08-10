
bam_location1 <- "C:/Users/Priit/Dropbox/Informaatika/Helsingi Ülikool/Käsikiri 2/CNV/bams/A749N2.bam"
bed_location1 <- "coordinates/chr15.bed"
out_location1 <- "test/"


commandArgs <- function(...){
  c(bam_location1, bed_location1, out_location1)
}

source("count_reads.R")


bam_location2 <- "C:/Users/Priit/Dropbox/Informaatika/Helsingi Ülikool/Käsikiri 2/CNV/bams/A749N2.bam"
bed_location2 <- "coordinates/chr15.bed" 
reference_location <- "reference.tsv" 

commandArgs <- function(...){
  c(bam_location2, bed_location2, reference_location)
}

source("analyse_bam.R")



results_location <- "results.A749N2.bam.tsv"
bed_location3 <- "coordinates/as_pws.bed" 

commandArgs <- function(...){
  c(results_location, bed_location3)
}

source("visualize.R")
