#' @useDynLib BUStoolsR
#' @importFrom Rcpp sourceCpp
NULL

#' Convert the Output of \code{kallisto bus} into Gene by Gell Matrix
#' 
#' This function takes the output file of \code{kallisto bus}, after being 
#' sorted and converted into text witth \code{bustools}. See vignettes on the
#' website of this package for a tutorial. The \code{bustools} output has 4 
#' columns: barcode, UMI, equivalence class, and counts. This function converts
#' that file into a sparse matrix that can be used in downstream analyses.
#' 
#' For 10x data sets, you can find a barcode whitelist file that comes with
#' CellRanger installation. You don't need to run CellRanger to get that. An 
#' example path to get the whitelist file is
#' \code{cellranger-2.1.0/cellranger-cs/2.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt}
#' for v2 chemistry.
#' 
#' @param fn File name of output of \code{bustools}, namely \code{output.bus}
#' from \code{kallisto bus} after sorted and converted to text by
#'  \code{bustools}.
#' @param genes A list with each element a string vector of genes that an 
#' equivalence class maps to, generated earlier in R.
#' @param whitelist A character vector with valid cell barcodes.
#' @param est_ncells Estimated number of cells; providing this argument will
#' speed up computation as it minimizes memory reallocation as vectors grow.
#' @param est_ngenes Estimated number of genes.
#' @param display_progress Whether to display progress.
#' @param progress_unit How many iteration to print one progress update when
#' reading in the \code{kallisto bus} file.
#' @return A sparse matrix with genes in rows and cells in columns.
#' @importFrom zeallot %<-%
#' @importClassesFrom Matrix dgCMatrix
#' @export
make_sparse_matrix <- function(fn, genes, est_ncells, est_ngenes, 
                               whitelist = NULL, 
                               display_progress = TRUE,
                               progress_unit = 5e6) {
  fn <- normalizePath(fn, mustWork = TRUE)
  if (is.null(whitelist)) {
    whitelist <- ""
  }
  c(res_mat, barcodes, genes) %<-% fill_cell_gene(fn, genes, 
                                                  est_ncells, est_ngenes, 
                                                  whitelist, 
                                                  display_progress,
                                                  progress_unit)
  rownames(res_mat) <- genes
  colnames(res_mat) <- barcodes
  res_mat
}
