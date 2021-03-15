#' Infer sample high-risk probability.
#'
#' This function does the following:
#' \enumerate{
#' \item Apply bin-based GC% correct.
#' \item Normalizing by read count.
#' \item Normalize by bin length.
#' \item Apply PCA-based normalization
#' \item Calculate reference group statistics per focus region.
#' \item Calculate sample Mahalanobis distance from the reference group.
#' \item Calculate -log10 chi-squared distribution probabilities.
#' }
#'
#' Outputs results to different files.
#'
#' @importFrom magrittr %>%
#' @param bam_location Path to the input file
#' @param reference_location Path to the reference file
#' @param use_pca Use PCA based normalization?
#' @param nComp How many components to use in PCA-based normalization?
#' @param bin_plot Create and save detailed bin plots?
#' @param result_plot Create and save detailed result?
#' @param save_bins Save bins?
#' @param pretty_output_columns Output 'nice' human readable column names?
#' @export
#' @examples
#'
#' bam <- "sample.bam"
#' reference <- "reference.tsv"
#' infer_normality(bam, reference)
#' head("sample.bam.tsv")
infer_normality <- function(bam_location,
                            reference_location,
                            use_pca = TRUE,
                            nComp = 80,
                            bin_plot = TRUE,
                            result_plot = TRUE,
                            save_bins = FALSE,
                            pretty_output_columns = TRUE)  {
  message_package_version()
  
  
  message("Reading reference set from: ", reference_location)
  reference <-
    readr::read_tsv(reference_location) %>%
    dplyr::mutate(reference = TRUE) %>%
    dplyr::group_by(sample) %>%
    # de-identify reference.
    dplyr::mutate(sample = as.character(dplyr::cur_group_id())) %>%
    dplyr::ungroup()
  
  check_reference(reference, use_pca)
  
  message("Reading and binning: ", bam_location)
  binned_reads <- bin_bam(
    bam_location,
    reference %>%
      dplyr::select(chr, start, end, focus) %>%
      dplyr::distinct(chr, start, end, focus)
  ) %>%
    dplyr::mutate(reference = FALSE)
  
  check_uniq(binned_reads, reference)
  
  
  message("Merging BAM and reference for calculations.")
  samples <- reference %>%
    dplyr::bind_rows(binned_reads)
  
  rm(binned_reads, reference)
  
  
  samples <- samples %>%
    gc_correct() %>%
    normalize_reads() %>%
    # Optimize memory
    dplyr::select(chr, focus, start, sample, reference, gc_corrected)
  
  if (use_pca) {
    samples <- pca_correct(samples)
  }
  
  
  samples <- calculate_bin_stat(samples)
  
  sample_name <- basename(bam_location)
  
  
  if (bin_plot) {
    message("Creating and saving region plots.")
    save_bin_plot(samples, sample_name)
  }
  
  
  if (save_bins) {
    message("Writing bins to file.")
    readr::write_tsv(
      samples %>%
        dplyr::ungroup() %>%
        dplyr::filter(sample == sample_name) %>%
        dplyr::select(chr, focus, start, PPDX_norm) %>%
        dplyr::rename(`Normalized Z-Score` = PPDX_norm),
      paste0(sample_name, ".bins", ".tsv")
    )
  }
  
  
  samples <- calculate_summary(samples)
  
  if (result_plot) {
    message("Creating and saving sample specific plot.")
    save_result_plot(samples, sample_name)
  }
  
  
  message("Saving metrics.")
  if (pretty_output_columns) {
    samples %>%
      dplyr::filter(!reference) %>%
      dplyr::select(-reference) %>%
      dplyr::rename(Sample = sample) %>%
      dplyr::rename(Chromosome = chr) %>%
      dplyr::rename(Subregion = focus) %>%
      dplyr::rename(`% of bins over reference bin mean median` = affected_over) %>%
      dplyr::rename(`High risk probability` = p) %>%
      dplyr::rename(`Z-Score` = PPDX) %>%
      dplyr::rename(`Normalized Z-Score` = PPDX_norm) %>%
      dplyr::mutate_if(is.numeric, round, 3) %>%
      readr::write_tsv(paste0(sample_name, ".tsv"))
  } else{
    samples %>%
      dplyr::filter(!reference) %>%
      dplyr::select(-reference) %>%
      dplyr::rename(subregion = focus) %>%
      dplyr::rename(`percentage_bins_over_ref_bin_mean_median` = affected_over) %>%
      dplyr::rename(`high_risk_probability` = p) %>%
      dplyr::rename(`z_score` = PPDX) %>%
      dplyr::rename(`normalized_z_score` = PPDX_norm) %>%
      dplyr::mutate_if(is.numeric, round, 3) %>%
      readr::write_tsv(paste0(sample_name, ".tsv"))
  }
}
