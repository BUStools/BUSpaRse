#' Knee plot for filtering empty droplets
#' 
#' Visualizes the inflection point to filter empty droplets. This function plots 
#' different datasets with a different color. Facets can be added after calling
#' this function with `facet_*` functions. Will be added to the next release
#' version of BUSpaRse.
#' 
#' @param bc_rank A `DataFrame` output from `DropletUtil::barcodeRanks` or a
#' named list of such `DataFrame`s.
#' @return A ggplot2 object.
#' @export
setGeneric("knee_plot", function(bc_rank) standardGeneric("knee_plot"),
           signature = "bc_rank")

#' @rdname knee_plot
#' @importFrom tidyr unnest
#' @importFrom ggplot2 ggplot aes geom_line geom_hline geom_vline scale_x_log10
#' scale_y_log10 labs
#' @export
#' @examples 
#' 
setMethod("knee_plot", "DataFrame",
         function(bc_rank) {
  knee_plt <- tibble(rank = bc_rank[["rank"]],
                     total = bc_rank[["total"]]) %>% 
    distinct() %>% 
    dplyr::filter(total > 0)
  annot <- tibble(inflection = metadata(bc_rank)[["inflection"]],
                  rank_cutoff = max(bc_rank$rank[bc_rank$total > metadata(bc_rank)[["inflection"]]]))
  rank <- total <- inflection <- rank_cutoff <- NULL
  p <- ggplot(knee_plt, aes(rank, total)) +
    geom_line() +
    geom_hline(aes(yintercept = inflection), data = annot, linetype = 2) +
    geom_vline(aes(xintercept = rank_cutoff), data = annot, linetype = 2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = "Rank", y = "Total UMIs")
  return(p)
})

#' @rdname knee_plot
#' @export
setMethod("knee_plot", "list",
          function(bc_rank) {
            knee_plt <- tibble(rank = lapply(bc_rank, function(x) x[["rank"]]),
                               total = lapply(bc_rank, function(x) x[["total"]]),
                               dataset = names(bc_rank)) %>% 
              unnest(cols = c(rank, total)) %>% 
              distinct() %>% 
              dplyr::filter(total > 0)
            annot <- tibble(inflection = vapply(bc_rank, function(x) metadata(x)[["inflection"]], 
                                                FUN.VALUE = numeric(1L)),
                            rank_cutoff = vapply(bc_rank, 
                                                 function(x) max(x$rank[x$total >
                                                                          metadata(x)[["inflection"]]]),
                                                 FUN.VALUE = numeric(1L)),
                            dataset = names(bc_rank))
            rank <- total <- inflection <- rank_cutoff <- dataset <- NULL
            p <- ggplot(knee_plt, aes(rank, total, color = dataset)) +
              geom_line() +
              geom_hline(aes(yintercept = inflection, color = dataset), 
                         data = annot, linetype = 2) +
              geom_vline(aes(xintercept = rank_cutoff, color = dataset),
                         data = annot, linetype = 2) +
              scale_x_log10() +
              scale_y_log10() +
              labs(x = "Rank", y = "Total UMIs")
            return(p)
          }
)
