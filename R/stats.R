#' Calculate bin statistics. Can be only called right after normalization.
#'
#' Required columns "chr", "start", "reference", "gc_corrected".
#'
#'
#' 1. Calculates reference group each (by chromosome and start coordinate) bin
#' mean and sd.
#' 2. Calculates regular Z-score for each sample (inc reference sample) and
#' if this samples bin is over ref mean bin.
#' 3. Group by "sample", "chr", "focus", "reference" and calculate
#' * PPDX (Z-score divided by sqrt())
#' * PPDX_norm ((PPDX + 1) / (n_bins_over_mean + 2))
#
#' @param samples A normalized data frame to be used in bin calculation.
#' @return A data frame that has added columns PPDX, PPDX_norm-
calculate_bin_stat <- function(samples) {
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
      over_mean = as.integer(gc_corrected >= mean_ref_bin)
    ) %>%
    dplyr::filter(!is.na(z_score)) %>%
    dplyr::group_by(sample, chr, focus, reference) %>%
    dplyr::mutate(
      PPDX = z_score / sqrt(dplyr::n()),
      # Normalize PPDX with over_mean + Laplace smoothing
      PPDX_norm = ((z_score / sqrt(dplyr::n(
        
      ))) + 1 / dplyr::n()) / (sum(over_mean) + 2 / dplyr::n())
    )
  
  return(samples)
}


#' Calculate summary statistics (probabilities)
#'
#' 1. Sums sample "PPDX", "PPDX_norm", "over_mean"
#' 2. Calculates reference mean of "PPDX" and "PPDX_norm"
#' 3. For each subregion (focus) calculates Mahalanobis distance from the 
#' reference group (PPDX, PPDX_norm) and transforms it with chi-square(df = 2)
#' to -log10 probability.
#'
#' @param samples A binned statistics data frame.
#' @return A data frame with probabilities
calculate_summary <- function(samples) {
  samples <- samples %>%
    dplyr::summarise(
      PPDX = sum(PPDX),
      PPDX_norm = sum(PPDX_norm),
      affected_over = sum(over_mean) / dplyr::n() * 100
    ) %>%
    dplyr::group_by(reference, chr, focus) %>%
    dplyr::mutate(mean_x = mean(PPDX_norm),
                  mean_y = mean(PPDX)) %>%
    dplyr::ungroup()
  
  
  samples <- samples %>%
    dplyr::group_by(chr, focus) %>%
    dplyr::group_split() %>%
    purrr::map_dfr( ~ {
      cov <-
        stats::cov(.x %>%
                     dplyr::filter(reference) %>%
                     dplyr::select(PPDX_norm, PPDX))
      
      center <- .x %>%
        dplyr::filter(reference) %>%
        dplyr::select(mean_x, mean_y) %>%
        dplyr::distinct()
      
      distances <-
        stats::mahalanobis(.x %>% dplyr::select(PPDX_norm, PPDX),
                           c(center$mean_x, center$mean_y),
                           cov)
      
      dplyr::bind_cols(
        sample = .x$sample,
        chr = .x$chr,
        focus = .x$focus,
        reference = .x$reference,
        PPDX_norm = .x$PPDX_norm,
        PPDX = .x$PPDX,
        affected_over = .x$affected_over,
        p = -log10(stats::pchisq(
          distances, df = 2, lower.tail = FALSE
        ) + 1e-100)
      )
    })
}
