source(here::here("R/util.R"))

args = commandArgs(trailingOnly = TRUE)


if (length(args) != 2) {
  stop(
    "Please provide (1) binnable BAM file (.bam) and (2) regions to bin (.bed)",
    call. = FALSE
  )
}

bam_location <- args[1]
bed_location <- args[2]

reads_per_bin <- bin_counts(bam_location, bed_location)


write_tsv(path = paste0(basename(bam_location), ".tsv"), x = reads_per_bin)
