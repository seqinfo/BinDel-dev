#' Calculate bin statistics.
#'
#' Calculate each genomic bin statistics compared to the reference group. Can
#'  be only called right after normalization.
#'
#' Required columns \emph{chr}, \emph{start}, \emph{reference}, \emph{gc_corrected}.
#'
#'
#' \enumerate{
#'  \item Calculates reference group each bin (by chromosome and start coordinate)
#' mean and SD.
#'  \item Calculates regular Z-score for each sample (including reference sample) and
#' if the bin is over the reference mean bin.
#' \item Group by \emph{sample}, \emph{chr}, \emph{focus}, \emph{reference} and calculate
#' * \eqn{PPDX = \frac{Z-score}{\sqrt{n}}}
#' * \eqn{PPDX_norm = \frac{PPDX + 1}{n_{bins_over_mean} + 2}} # Laplace smoothing
#' }
#
#' @param samples A normalized data frame to be used in bin calculation.
#' @return A data frame that has added columns \emph{PPDX}, \emph{PPDX_norm}.
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
    ) %>%
    dplyr::ungroup()
  
  return(samples)
}


#' Calculate summary statistics (probabilities)
#'
#' Aggregates bin statistics to region high probability risks. Requires columns
#' \emph{PPDX}, \emph{PPDX_norm} and \emph{over_mean}.
#'
#'\enumerate{
#' \item Sums sample \emph{PPDX}, \emph{PPDX_norm}, \emph{over_mean}.
#' \item Calculates reference mean of \emph{PPDX} and \emph{PPDX_norm}.
#' \item For each subregion (focus) calculates Mahalanobis distance from the
#' reference group (\emph{PPDX}, \emph{PPDX_norm}) and transforms it with chi-square(df = 2)
#' to -log10 probability.
#'}
#' @param samples A binned statistics data frame.
#' @param amplify_sensitivity used both the \emph{PPDX} and \emph{PPDX_norm}?
#' @return A data frame with probabilities
calculate_summary <- function(samples, amplify_sensitivity) {
  samples <- samples %>%
    dplyr::group_by(sample, chr, focus, reference) %>%
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
      if (amplify_sensitivity) {
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
        
      } else{
        cov <-
          stats::cov(.x %>%
                       dplyr::filter(reference) %>%
                       dplyr::select(PPDX_norm))
        
        center <- .x %>%
          dplyr::filter(reference) %>%
          dplyr::select(mean_x) %>%
          dplyr::distinct()
        
        distances <-
          stats::mahalanobis(.x %>% dplyr::select(PPDX_norm),
                             c(center$mean_x),
                             cov)
      }
      
      
      dplyr::bind_cols(
        sample = .x$sample,
        chr = .x$chr,
        focus = .x$focus,
        reference = .x$reference,
        PPDX_norm = .x$PPDX_norm,
        PPDX = .x$PPDX,
        affected_over = .x$affected_over,
        p = -log10(
          stats::pchisq(distances,
                        df = if (amplify_sensitivity)
                          2
                        else
                          1,
                        lower.tail = FALSE) + 1e-100
        )
      )
    })
}
