#' Plot the transposed knee plot and inflection point
#' 
#' This function plots a transposed knee plot, showing the inflection point and
#' the number of remaining cells after inflection point filtering. It's
#' transposed since it's more generalizable to multi-modal data.
#' 
#' @param df The data frame from \code{\link{get_knee_df}}.
#' @param inflection Output of \code{\link{get_inflection}}.
#' @return A \code{ggplot2} object.
#' @export
#' @importFrom ggplot2 ggplot aes geom_path geom_vline geom_hline 
#' scale_x_log10 scale_y_log10 labs annotation_logticks geom_text
#' @examples 
#' # Download dataset already in BUS format
#' library(TENxBUSData)
#' TENxBUSData(".", dataset = "retina")
#' tr2g <- transcript2gene("Mus musculus", type = "vertebrate",
#'   ensembl_version = 99, kallisto_out_path = "./out_retina")
#' m <- make_sparse_matrix("./out_retina/output.sorted.txt",
#' tr2g = tr2g, est_ncells = 1e5,
#' est_ngenes = nrow(tr2g))
#' df <- get_knee_df(m)
#' infl <- get_inflection(df)
#' knee_plot(df, infl)
knee_plot <- function(df, inflection) {
  annot <- tibble(inflection = inflection,
                  rank_cutoff = max(df$rank[df$total > inflection]))
  ggplot(df, aes(total, rank)) +
    geom_path() +
    geom_vline(aes(xintercept = inflection), data = annot, linetype = 2, color = "gray40") +
    geom_hline(aes(yintercept = rank_cutoff), data = annot, linetype = 2, color = "gray40") +
    geom_text(aes(inflection, rank_cutoff, 
                  label = paste(rank_cutoff, "'cells'")),
              data = annot, vjust = 1) +
    scale_x_log10() +
    scale_y_log10() +
    labs(y = "Rank", x = "Total UMIs") +
    annotation_logticks()
}

#' @rdname knee_plot
#' @param mat Gene count matrix, a dgCMatrix.
#' @return A tibble with two columns: \code{total} for total UMI counts for 
#' each barcode, and \code{rank} for rank of the total counts, with number 1
#' for the barcode with the most counts.
#' @export
#' 
get_knee_df <- function(mat) {
  tibble(total = colSums(mat),
         rank = row_number(desc(total))) %>%
    distinct() %>%
    filter(total > 0) %>% 
    arrange(rank)
}

#' Compute inflection point of knee plot
#'
#'
#' @param df The data frame from \code{\link{get_knee_df}}.
#' @param lower Minimum total UMI counts for barcode for it to be considered
#' when calculating the inflection point; this helps to avoid the noisy part of
#' the curve for barcodes with very few counts.
#' @return A \code{numeric(1)} for the total UMI count at the inflection point.
#' @note Code in part adapted from \code{barcodeRanks} from \code{DropetUtils}.
#' @export
#' @importFrom dplyr transmute
#' 
get_inflection <- function(df, lower = 100) {
  df_fit <- df %>% 
    filter(total > lower) %>% 
    transmute(log_total = log10(total),
              log_rank = log10(rank))
  d1n <- diff(df_fit$log_total)/diff(df_fit$log_rank)
  right.edge <- which.min(d1n)
  10^(df_fit$log_total[right.edge])
}
