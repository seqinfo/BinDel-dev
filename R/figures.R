#' Creates and saves binned normalized Z-Score figure per region.
#'
#' @importFrom magrittr %>%
#' @param data A binned Z-scored data frame with columns "focus"," start",
#'  "PPDX_norm", "sample" and reference.
#' @param sample_name sample name of interest
save_bin_plot <- function(samples, sample_name) {
  capt <-
    paste(
      "Supplementary bin figure complementing the final result",
      "figure. Each subregion depicts the bins from the region of interest",
      "and their deviation from the reference group."
    )
  
  ggplot2::ggsave(
    paste0(sample_name, ".bins", ".png"),
    ggplot2::ggplot(
      samples,
      ggplot2::aes(
        as.numeric(start),
        PPDX_norm,
        group = sample,
        color = reference
      )
    ) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap( ~ focus, scales = "free", ncol = 2) +
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
      ggplot2::labs(
        title = paste(
          basename(bam_path),
          "scientific",
          utils::packageName(),
          "bin report"
        ),
        subtitle = paste("Version: ", utils::packageVersion(utils::packageName())),
        caption = stringr::str_wrap(capt, width = 100)
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(),
        strip.background = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5),
        plot.caption = ggplot2::element_text(
          hjust = 0,
          vjust = 0,
          face = "italic",
          color = "black"
        )
      )
    ,
    width = 10,
    height = max(1, round(length(
      unique(samples$focus)
    ) / 2)) * 2
  )
}


#' Create and save combined probability result figure.
#'
#' @importFrom magrittr %>%
#' @param data A probability data frame with columns "chr", "PPDX",
#' and "focus".
#' @param sample_name sample name of interest
save_result_plot <- function(samples, sample_name) {
  ordered <- samples %>%
    dplyr::mutate(chr = as.numeric(stringr::str_remove(chr, "chr"))) %>%
    dplyr::mutate(sign = sign(PPDX)) %>%
    dplyr::mutate(shape = ifelse(sign > 0, 24L, ifelse(sign < 0, 25L, 18L))) %>%
    dplyr::mutate(shape = ifelse(reference, 21L, shape)) %>%
    dplyr::mutate(color = ifelse(reference, "grey", "black")) %>%
    dplyr::mutate(alpha = ifelse(reference, 0.5, 1))
  
  
  capt <- paste(
    "High-risk probability of the regions of interest.",
    "The higher the region's risk probability, the more probable the",
    "deviation from the healthy reference group. The upwards triangle hints",
    "for duplication, downwards triangle hints for the deletion. High-risk",
    "candidates need to be verified by outputted genomic bin figure, which",
    "depicts a sample of interest compared to the healthy reference group."
  )
  
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
      ggplot2::ylim(0, 105) +
      ggplot2::ylab("High risk probability") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste(
          basename(bam_path),
          "scientific",
          utils::packageName(),
          "report"
        ),
        subtitle = paste("Version: ", utils::packageVersion(utils::packageName())),
        caption = stringr::str_wrap(capt, width = 100)
      ) +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(),
        strip.background = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "none",
        plot.caption = ggplot2::element_text(
          hjust = 0,
          vjust = 0,
          face = "italic",
          color = "black"
        ),
        axis.title.x = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 90,
          hjust = 1,
          vjust = 0.5
        )
      )
  )
}
