library(tidyverse)
library(fs)

if (length(args) != 1) {
  stop("Please provide (only) the location of the folder where reference counts are stored.", call.=FALSE)
}

folder_location <- args[1]


ref_count_files <- fs::dir_ls(path = folder_location, regexp = "count\\..*\\.tsv$")


reference <- ref_count_files %>% 
  map_dfr(read_tsv)

write_tsv(reference, "reference.tsv")

