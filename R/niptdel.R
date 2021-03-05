#' Find GC% for GRCh38 based .bed data frame.
#'
#' @importFrom magrittr %>%
#' @param bed A data frame in .bed format with columns: 'chr', 'start', 'end', 'focus'.
#' @return A data frame in bed format with GC%.
find_gc <- function(bed) {
  reads <- bed %>%
    dplyr::mutate(gc = Biostrings::letterFrequency(
      Biostrings::getSeq(
        BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38,
        GenomicRanges::GRanges(chr, IRanges::IRanges(start = start, end = end))
      ),
      "GC",
      as.prob = T
    )) %>%
    dplyr::mutate(gc = round(gc, 1))
}


#' Bin aligned sequences (from .bam) into genomic bins based on the .bed file.
#'
#'
#' @importFrom magrittr %>%
#' @param bam_location A location to the BAM-file to bin.
#' @param bed A data frame in .bed format with columns: 'chr', 'start', 'end', 'focus'.
#' @return A data frame in bed format with GC%.
#' @export
bin_bam <- function(bam_location, bed) {
  bam <-
    GenomicAlignments::readGAlignments(bam_location,
                                       param = Rsamtools::ScanBamParam(
                                         flag = Rsamtools::scanBamFlag(isDuplicate = FALSE, isSecondaryAlignment = FALSE),
                                         what = c("pos")
                                       ))
  
  binned_counts <- bed %>%
    dplyr::mutate(reads = SummarizedExperiment::assay(
      GenomicAlignments::summarizeOverlaps(GenomicRanges::makeGRangesFromDataFrame(.), bam, mode = "IntersectionStrict")
    )) %>%
    dplyr::mutate(sample = basename(bam_location))
  
  return(binned_counts)
}


#' Infer sample normality probability.
#'
#' This function does the following:
#' 1. Apply bin-based GC% correct.
#' 2. Normalizing by read count.
#' 3. Normalize by bin length.
#' 4. Apply PCA-based normalization
#' 5. Calculate reference group statistics per focus region.
#' 6. Calculate sample Mahalanobis distance from the reference group.
#' 7. Calculate -log10 chi-squared distribution probabilities.
#'
#' @importFrom magrittr %>%
#' @param bam_location Path to the input file
#' @param reference_location Path to the reference file
#' @param use_pca Use PCA based normalization?
#' @param nComp How many components to use in PCA-based normalization?
#' @param plot_results Create and save detailed plots?
#' @return A data frame with scores for the provided BAM.
#' @export
#' @examples
#'
#' bam <- "sample.bam"
#' reference <- "reference.tsv"
#' metrics <- infer_normality(bam, reference)
infer_normality <- function(bam_location,
                            reference_location,
                            use_pca = TRUE,
                            nComp = 80,
                            plot_results = TRUE)  {
  sample_name <- basename(bam_location)
  
  
  message("Reading reference file from:", reference_location)
  # Reference samples to be used to calculate z-scores
  reference <-
    readr::read_tsv(reference_location) %>%
    dplyr::mutate(reference = TRUE)
  
  # Check if reference file has all the required columns.
  ref_expected_cols <-
    c("chr", "start", "end", "focus", "reads", "sample")
  for (col in ref_expected_cols) {
    if (!col %in% colnames(reference)) {
      stop(paste0("Reference file has a missing column: '", col, "'."))
    }
  }
  
  number_of_reference_samples <- nrow(reference %>%
                                        dplyr::select(sample) %>%
                                        dplyr::distinct(sample))
  
  if (nrow(reference) == 0) {
    stop("Reference file cannot be empty")
  }
  
  if (number_of_reference_samples < 10) {
    warning("Reference group has less than 10 samples.")
  } else {
    message("Reference group has ",
            number_of_reference_samples,
            " samples.")
  }
  
  
  if (use_pca) {
    if (number_of_reference_samples < nComp) {
      stop(
        paste0(
          "nComp (",
          nComp,
          ") must be lower than number of reference samples(",
          number_of_reference_samples,
          ")."
        )
      )
    }
  }
  
  message("Reading and binning: ", bam_location)
  # Bin BAM under investigation
  binned_reads <- bin_bam(
    bam_location,
    reference %>%
      dplyr::select(chr, start, end, focus) %>%
      dplyr::distinct(chr, start, end, focus)
  ) %>%
    dplyr::mutate(reference = FALSE)
  
  
  message("Checking if bam name is unique in terms of reference.")
  # Note, reference samples names must not overlap with the analyzable sample.
  if (nrow(
    binned_reads %>%
    dplyr::select(sample) %>%
    dplyr::distinct(sample) %>%
    dplyr::inner_join(reference %>% dplyr::select(sample))
  ) != 0) {
    stop("BAM name is present in the reference group. Stopping.")
  }
  
  message("Merging BAM and reference for calculations.")
  # Merge: BAM + reference
  samples <- reference %>%
    dplyr::bind_rows(binned_reads)
  
  
  rm(binned_reads)
  
  
  # GC-correct (sample wise) (PMID: 28500333 and PMID: 20454671)
  message("Applying GC% correct.")
  # Find GC% of the genome (HG38)
  samples <- samples %>%
    dplyr::left_join(find_gc(
      dplyr::select(.data = ., chr, start, end, focus) %>%
        dplyr::distinct(chr, start, end, focus)
    )) %>%
    # Do GC correct
    dplyr::group_by(sample, gc) %>%
    dplyr::mutate(avg_reads_gc_interval = mean(reads)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(weights = mean(reads) / avg_reads_gc_interval) %>%
    dplyr::mutate(gc_corrected = reads * weights) %>%
    dplyr::filter(!is.na(gc_corrected)) %>%
    dplyr::ungroup()
  
  
  message("Normalizing by read count and bin length.")
  samples <- samples %>%
    # Sample read count correct
    dplyr::group_by(sample) %>%
    dplyr::mutate(gc_corrected = gc_corrected / sum(gc_corrected)) %>%
    dplyr::ungroup() %>%
    # Sample bin length correct
    dplyr::mutate(gc_corrected = gc_corrected / (end - start)) %>%
    # Optimize memory
    dplyr::select(chr, focus, start, sample, reference, gc_corrected)
  
  
  if (use_pca) {
    message("Applying PCA normalization with ", nComp, " components.")
    # For PCA sort ()
    samples <- samples %>%
      dplyr::arrange(reference)
    
    # Pivot wide for PCA normalization
    wider <- samples %>%
      dplyr::select(focus, start, sample, reference, gc_corrected) %>%
      tidyr::pivot_wider(
        names_from = c(focus, start),
        id_cols = c(sample, reference),
        values_from = gc_corrected,
        names_sep = ":"
      )
    
    # https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com
    # Train PCA
    ref <- wider %>%
      dplyr::filter(reference) %>%
      dplyr::select(-reference, -sample)
    
    mu <- colMeans(ref, na.rm = T)
    refPca <- stats::prcomp(ref)
    
    
    Xhat <- refPca$x[, 1:nComp] %*% t(refPca$rotation[, 1:nComp])
    Xhat <- scale(Xhat, center = -mu, scale = FALSE)
    
    # Use trained PCA on other samples
    pred <- wider %>%
      dplyr::filter(!reference) %>%
      dplyr::select(-reference, -sample)
    
    Yhat <-
      stats::predict(refPca, pred)[, 1:nComp] %*% t(refPca$rotation[, 1:nComp])
    Yhat <- scale(Yhat, center = -mu, scale = FALSE)
    
    # Actual PCA normalization and conversion back to long:
    normalized <-
      dplyr::bind_rows(as.data.frame(as.matrix(pred) / as.matrix(Yhat)),
                       as.data.frame(as.matrix(ref) / as.matrix(Xhat))) %>%
      tidyr::pivot_longer(
        names_sep = ":",
        names_to = c("focus", "start"),
        cols = dplyr::everything(),
        values_to = "gc_corrected"
      )
    
    normalized$sample <- samples$sample
    normalized$reference <- samples$reference
    normalized$chr <- samples$chr
    samples <- normalized
    
    
    rm(wider)
    rm(pred)
    rm(ref)
    rm(Yhat)
    rm(Xhat)
    rm(normalized)
    
    
  }
  
  message("Calulating reference group statistics.")
  # Calculate each reference bin i mean
  reference <- samples %>%
    dplyr::filter(reference) %>%
    dplyr::group_by(chr, start) %>%
    dplyr::summarise(mean_ref_bin = mean(gc_corrected),
                     mean_ref_sd = sd(gc_corrected)) %>%
    dplyr::ungroup()
  
  
  message("Calculating metrics.")
  samples <- samples %>%
    dplyr::right_join(reference) %>%
    dplyr::mutate(
      z_score = (gc_corrected - mean_ref_bin) / mean_ref_sd,
      over_median = as.integer(gc_corrected >= mean_ref_bin)
    ) %>%
    dplyr::filter(!is.na(z_score)) %>%
    dplyr::group_by(sample, chr, focus, reference) %>%
    dplyr::mutate(
      z_score_PPDX = z_score / sqrt(dplyr::n()),
      # Normalize z_score_PPDX with over_median + Laplace smoothing
      z_score_PPDX_norm = ((z_score / sqrt(dplyr::n(
        
      ))) + 1 / n()) / (sum(over_median) + 2 / n())
    )
  
  # Create detailed plots of the regions
  if (plot_results) {
    message("Creating and saving region plots.")
    
    pdf(paste0(sample_name, ".details", ".pdf"), title = sample_name)
    samples %>%
      dplyr::mutate(chr = as.numeric(stringr::str_remove(chr, "chr"))) %>%
      dplyr::arrange(chr) %>%
      dplyr::group_by(focus) %>%
      dplyr::do({
        print(
          ggplot2::ggplot(
            .,
            ggplot2::aes(
              as.numeric(start),
              z_score_PPDX_norm,
              group = sample,
              color = reference
            )
          ) +
            ggrastr::rasterise(ggplot2::geom_line(), dpi = 300) +
            ggplot2::geom_hline(yintercept = 0) +
            ggplot2::scale_x_continuous(
              labels = function(x)
                format(
                  x,
                  big.mark = " ",
                  decimal.mark = ",",
                  scientific = FALSE
                )
            ) +
            ggplot2::scale_color_manual(values = c("red", "grey")) +
            ggplot2::ylab("Normalized Z-score") +
            ggplot2::ggtitle(paste(sample_name, unique(.$focus))) +
            ggplot2::theme_bw() +
            ggplot2::theme(
              panel.border = ggplot2::element_blank(),
              axis.title.x = ggplot2::element_blank(),
              axis.line = ggplot2::element_line(),
              strip.background = ggplot2::element_blank(),
              axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)
            )
        )
        invisible(.)
      })
    dev.off()
  }
  
  samples <- samples %>%
    dplyr::summarise(
      z_score_PPDX = sum(z_score_PPDX),
      z_score_PPDX_norm = sum(z_score_PPDX_norm),
      affected_over = sum(over_median) / dplyr::n() * 100,
    ) %>%
    dplyr::group_by(reference, chr, focus) %>%
    dplyr::mutate(mean_x = mean(z_score_PPDX_norm),
                  mean_y = mean(z_score_PPDX)) %>%
    dplyr::ungroup()
  
  
  
  
  samples <- samples %>%
    dplyr::group_by(chr, focus) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(~ {
      cov <-
        stats::cov(
          .x %>%
            dplyr::filter(reference) %>%
            dplyr::select(z_score_PPDX_norm, z_score_PPDX)
        )
      
      center <- .x %>%
        dplyr::filter(reference) %>%
        dplyr::select(mean_x, mean_y) %>%
        dplyr::distinct()
      
      distances <- stats::mahalanobis(
        .x %>% dplyr::select(z_score_PPDX_norm, z_score_PPDX),
        c(center$mean_x, center$mean_y),
        cov
      )
      
      dplyr::bind_cols(
        sample = .x$sample,
        chr = .x$chr,
        focus = .x$focus,
        reference = .x$reference,
        z_score_PPDX_norm = .x$z_score_PPDX_norm,
        z_score_PPDX = .x$z_score_PPDX,
        affected_over = .x$affected_over,
        affected_under = 100 - .x$affected_over,
        p = -log10(stats::pchisq(
          distances, df = 2, lower.tail = FALSE
        ) + 1e-100)
      )
    })
  
  
  if (plot_results) {
    ordered <- samples %>%
      dplyr::mutate(chr = as.numeric(stringr::str_remove(chr, "chr"))) %>%
      dplyr::mutate(sign = sign(z_score_PPDX)) %>%
      dplyr::mutate(shape = ifelse(sign > 0, 24L, ifelse(sign < 0, 25L, 18L))) %>%
      dplyr::mutate(shape = ifelse(reference, 21L, shape)) %>%
      dplyr::mutate(color = ifelse(reference, "grey", "black")) %>%
      dplyr::mutate(alpha = ifelse(reference, 0.5, 1))
    
    
    # Plot results
    ggplot2::ggsave(
      paste0(sample_name, ".summary", ".png"),
      ggplot2::ggplot(ordered, ggplot2::aes(
        x = stats::reorder(focus, chr), y = p
      )) +
        ggplot2::geom_point(
          alpha = ordered$alpha,
          shape = ordered$shape,
          color = ordered$color,
          fill = ordered$color,
          size = 1.6
        ) +
        ggplot2::ylab("High risk probability") +
        ggplot2::ggtitle(basename(bam_path)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          panel.border = ggplot2::element_blank(),
          axis.line = ggplot2::element_line(),
          strip.background = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          legend.position = "none",
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
          )
        )
    )
  }
  
  return(samples %>% dplyr::filter(!reference) %>% dplyr::select(-reference))
}