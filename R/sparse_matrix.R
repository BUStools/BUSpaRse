#' @useDynLib BUSpaRse
#' @importFrom Rcpp sourceCpp
NULL

#' Convert the Output of \code{kallisto bus} into Gene by Gell Matrix
#'
#' This function takes the output file of \code{kallisto bus}, after being
#' sorted and converted into text with \code{bustools}. See vignettes on the
#' [website of this package](https://bustools.github.io/BUS_notebooks_R/) for a
#' tutorial. The \code{bustools} output has 4 columns: barcode, UMI, equivalence
#' class, and counts. This function converts that file into a sparse matrix that
#' can be used in downstream analyses.
#'
#' This function can generate both the gene count matrix and the transcript
#' compatibility count (TCC) matrix. The TCC matrix has barcodes in the columns
#' and equivalence classes in the rows. See
#' [Ntranos et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0970-8)
#' for more information about the RCC matrix.
#'
#' For 10x data sets, you can find a barcode whitelist file that comes with
#' CellRanger installation. You don't need to run CellRanger to get that. An
#' example path to get the whitelist file is
#' \code{cellranger-2.1.0/cellranger-cs/2.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt}
#' for v2 chemistry.
#'
#' @param bus_path Path to the sorted text `bus` output file.
#' @param tr2g A Data frame with columns \code{gene} and \code{transcript}, in
#' the same order as in the transcriptome index for \code{kallisto}. This
#' argument can be missing or is ignored if only the TCC matrix, not the gene
#' count matrix, is made.
#' @param whitelist A character vector with valid cell barcodes. This is an
#' optional argument, that defaults to \code{NULL}. When it is \code{NULL},
#' all cell barcodes present that have some UMI assignable to genes or ECs will
#' be included in the sparse matrix whether they are known to be valid or not.
#' Barcodes with only UMIs that are not assignable to genes or ECs will still be
#' excluded.
#' @param gene_count Logical, whether the gene count matrix should be returned.
#' @param TCC Logical, whether the TCC matrix should be returned.
#' @param single_gene Logical, whether to use single gene mode. In single gene
#' mode, only UMIs that can be uniquely mapped to one gene are kept. Without
#' single gene mode, UMIs mapped to multiple genes will be evenly distributed to
#' those genes.
#' @param est_ncells Estimated number of cells; providing this argument will
#' speed up computation as it minimizes memory reallocation as vectors grow.
#' @param est_ngenes Estimated number of genes or equivalence classes.
#' @param verbose Whether to display progress.
#' @param progress_unit How many iteration to print one progress update when
#' reading in the \code{kallisto bus} file.
#' @return If both gene count and TCC matrices are returned, then this function
#' returns a list with two matrices, each with genes/equivalence classes in the
#' rows and barcodes in the columns. If only one of gene count and TCC matrices
#' is returned, then a \code{dgCMatrix} with genes/equivalence classes in the
#' rows and barcodes in the columns. These matrices are unfiltered. Please filter
#' the empty droplets before downstream analysis.
#' @seealso \code{\link{EC2gene}}
#' @importFrom zeallot %<-%
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' # Load toy example for testing
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' load(paste(toy_path, "toy_example.RData", sep = "/"))
#' out_fn <- paste0(toy_path, "/output.sorted.txt")
#' # With whitelist
#' m <- make_sparse_matrix(out_fn, tr2g_toy, 10, 3, whitelist = whitelist,
#'   gene_count = TRUE, TCC = FALSE, single_gene = TRUE,
#'   verbose = FALSE)
make_sparse_matrix <- function(bus_path, tr2g, est_ncells,
                               est_ngenes, whitelist = NULL, gene_count = TRUE,
                               TCC = TRUE, single_gene = TRUE, 
                               verbose = TRUE, progress_unit = 5e6) {
  bus_path <- normalizePath(bus_path, mustWork = TRUE)
  kallisto_out_path <- dirname(bus_path)
  bus_fn <- basename(bus_path)
  if (!grepl(".txt$", bus_fn)) {
    stop("Argument fn must point to a text file. Please run bustools text.")
  }
  if (TCC && !gene_count) {
    tr2g <- data.frame(gene = character(0),
      transcript = character(0))
  }
  if (is.null(whitelist)) {
    whitelist <- ""
  }
  # Prevent the no visible binding of global variable note in R CMD check
  barcodes_gc <- geneIDs <- barcodes_tcc <- ec_inds <- NULL
  if (gene_count && TCC) {
    c(c(gc_mat, barcodes_gc, geneIDs),
      c(tcc_mat, barcodes_tcc, ec_inds)) %<-% fill_cell_gene(bus_fn, kallisto_out_path,
      tr2g,
      est_ncells, est_ngenes,
      whitelist, gene_count, TCC,
      single_gene, 
      verbose, progress_unit)
    rownames(gc_mat) <- geneIDs
    colnames(gc_mat) <- barcodes_gc
    rownames(tcc_mat) <- ec_inds
    colnames(tcc_mat) <- barcodes_tcc
    return(list(gene_count = gc_mat, TCC = tcc_mat))
  } else if (gene_count) {
    c(gc_mat, barcodes_gc, geneIDs) %<-% fill_cell_gene(bus_fn, kallisto_out_path,
      tr2g,
      est_ncells, est_ngenes,
      whitelist, gene_count,
      FALSE, single_gene, 
      verbose, progress_unit)
    rownames(gc_mat) <- geneIDs
    colnames(gc_mat) <- barcodes_gc
    return(gc_mat)
  } else {
    c(tcc_mat, barcodes_tcc, ec_inds) %<-% fill_cell_gene(bus_fn, kallisto_out_path,
      tr2g,
      est_ncells, est_ngenes,
      whitelist, FALSE, TCC,
      FALSE, 
      verbose, progress_unit)
    rownames(tcc_mat) <- ec_inds
    colnames(tcc_mat) <- barcodes_tcc
    return(tcc_mat)
  }
}
