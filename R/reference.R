# Copyright: Priit Paluoja and Priit Palta; University of Tartu and University of Helsinki.
# All rights reserved, unauthorised usage and distribution are prohibited.
# Contact: priit.palta@gmail.com


#' Create a BinDel reference
#' 
#' @importFrom magrittr %>%
#' @param bam_locations The vector of .bam paths to write into the reference.
#' @param coordinates_file A location to the coordinates file  with columns: \emph{chr}, \emph{start}, \emph{end}, \emph{focus}, \emph{length}.
#' @param output_name The name of the reference file. If .gz is appended to the output name, the reference file is compressed.
#' @param col_names output column names?
#' @param anonymise include original sample name in the reference?
#' @param prefix if anonymise is T, set prefix to the sample name
#'
#' @export
#' @examples
#' write_reference(c("sample.bam"), "coordinates.tsv", "reference.gz")
write_reference <-
  function(bam_locations,
           coordinates_file,
           output_name,
           col_names = T,
           anonymise = T,
           prefix = "S") {
    message("Creating BinDel reference.")
    
    if (length(bam_locations) == 0) {
      stop("No .bam files found in '", bam_locations, "'.")
    }
    
    if (!file.exists(coordinates_file)) {
      stop("Coordinates file not found in '", coordinates_file, "'.")
    }
    
    if (file.exists(output_name)) {
      stop("Reference already exists, aborting.")
    }
    
    bed <- divide_bins(coordinates_file)
    
    if (col_names) {
      df_cols <- c("chr", "start", "end", "focus", "reads", "sample")
      empty_df <-
        data.frame(matrix(ncol = length(df_cols), nrow = 0))
      colnames(empty_df) <- df_cols
      readr::write_tsv(x = empty_df, output_name)
    }
    
    
    i <- 1
    lapply(bam_locations, function(x) {
      message("Processing '", x, "'.")
      
      binned <- bin_bam(x, bed)
      
      if (anonymise) {
        binned <- binned %>%
          dplyr::mutate(sample = paste0(prefix, i))
        i <<- i + 1
      }
      
      readr::write_tsv(
        file = output_name,
        x = binned,
        append = TRUE,
        col_names = F
      )
    })
  }


#' Create a reference file.
#'
#' Takes a location to the folder with .bam files and .bed file and bins the
#' .bam files with GRCh38. 
#' 
#' @param bam_locations The path to the folder, where the reference files exists.
#' @param bed_location A location to the .bed file with columns: \emph{chr}, \emph{start}, \emph{end}, \emph{focus}.
#' @param reference_name The name of the output file. File must not exist.
#'
#' @export
create_reference <-
  function(bam_locations,
           bed_location,
           reference_name) {
    files <- list.files(
      path = bam_locations,
      pattern = "*.bam",
      full.names = TRUE,
      recursive = FALSE
    )
    
    .Deprecated("write_reference")
    
    if (length(files) == 0) {
      stop("No .bam files found in '", bam_locations, "'.")
    }
    
    if (!file.exists(bed_location)) {
      stop(".bed file not found in '", bed_location, "'.")
    }
    
    if (file.exists(reference_name)) {
      stop("Output '", reference_name, "' already exists.")
    }
    
    bed <- readr::read_tsv(bed_location)
    
    
    df_cols <- c("chr", "start", "end", "focus", "reads", "sample")
    empty_df <- data.frame(matrix(ncol = length(df_cols), nrow = 0))
    colnames(empty_df) <- df_cols
    readr::write_tsv(x = empty_df, reference_name)
    
    
    lapply(files, function(x) {
      message("Processing '", x, "'.")
      readr::write_tsv(
        file = reference_name,
        x = bin_bam(x, bed),
        append = TRUE,
        col_names = FALSE
      )
    })
  }

