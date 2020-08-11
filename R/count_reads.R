source(here::here("util.R"))

args = commandArgs(trailingOnly = TRUE)


if (length(args) != 3) {
  stop(
    "Please provide (1) binnable BAM file (.bam), (2) regions to bin (.bed) and output folder location.",
    call. = FALSE
  )
}

bam_location <- args[1]
bed_location <- args[2]
out_location <- args[3]

reads_per_bin <- bin_counts(bam_location, bed_location)


write_tsv(path = paste0(out_location, "count.", basename(bam_location), ".tsv"), x = reads_per_bin)
