library(tidyverse)
library(fs)

args = commandArgs(trailingOnly=TRUE)

if (length(args) != 2) {
  stop("Please provide (only) the location of the folder where reference counts are stored.", call.=FALSE)
}

folder_location <- args[1]
ref_name <- args[2]


ref_count_files <- fs::dir_ls(path = folder_location, regexp = "count\\..*\\.tsv$")


reference <- ref_count_files %>% 
  map_dfr(read_tsv)

write_tsv(reference, ref_name)

