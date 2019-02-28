#' @useDynLib BUSpaRse
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
#' @inheritParams EC2gene
#' @param bus_path Path to the sorted text `bus` output file.
#' @param whitelist A character vector with valid cell barcodes. This is an
#' optional argument, that defaults to \code{NULL}. When it is \code{NULL},
#' all cell barcodes present will be included in the sparse matrix whether they
#' are known to be valid or not.
#' @param gene_count Logical, whether the gene count matrix should be returned.
#' @param TCC Logical, whether the TCC matrix should be returned.
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
#' rows and barcodes in the columns.
#' @seealso \code{\link{EC2gene}}
#' @family functions to generate sparse matrix from outputs of other 
#' \code{BUSpaRse} functions
#' @importFrom zeallot %<-%
#' @importClassesFrom Matrix dgCMatrix
#' @export
#' @examples
#' \dontrun{
#' # Download dataset already in BUS format
#' library(TENxhgmmBUS)
#' library(Matrix)
#' download_hgmm(".", "hgmm100")
#' tr2g <- transcript2gene(c("Homo sapiens", "Mus musculus"),
#'                         kallisto_out_path = "./out_hgmm100")
#' res_mat <- make_sparse_matrix("./out_hgmm100/output.sorted.txt",
#'                               tr2g = tr2g, est_ncells = 3e5, 
#'                               est_ngenes = nrow(tr2g), gene_count = TRUE,
#'                               TCC = FALSE)
#' # Remove empty droplets
#' tot_counts <- colSums(res_mat)
#' res_mat <- res_mat[,tot_counts > 500]
#' }

make_sparse_matrix <- function(bus_path, tr2g, est_ncells, 
                               est_ngenes, whitelist = NULL, gene_count = TRUE,
                               TCC = TRUE, ncores = 0, verbose = TRUE,
                               progress_unit = 5e6) {
  bus_path <- normalizePath(bus_path, mustWork = TRUE)
  kallisto_out_path <- dirname(bus_path)
  bus_fn <- basename(bus_path)
  if (!grepl(".txt$", bus_fn)) {
    stop("Argument fn must point to a text file. Please run bustools text.")
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
                                                             ncores,
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
                                                        FALSE, ncores,
                                                        verbose, progress_unit)
    rownames(gc_mat) <- geneIDs
    colnames(gc_mat) <- barcodes_gc
    return(gc_mat)
  } else {
    c(tcc_mat, barcodes_tcc, ec_inds) %<-% fill_cell_gene(bus_fn, kallisto_out_path,
                                                          tr2g, 
                                                          est_ncells, est_ngenes, 
                                                          whitelist, FALSE, TCC,
                                                          ncores,
                                                          verbose, progress_unit)
    rownames(tcc_mat) <- ec_inds
    colnames(tcc_mat) <- barcodes_tcc
    return(tcc_mat)
  }
}

#' Get gene count or TCC matrix in one step
#' 
#' The \code{bustools} output has 4 columns: barcode, UMI, equivalence class, 
#' and counts. This function directly converts that file into a sparse matrix 
#' that can be used in downstream analyses in one step for species that are in
#' the Ensembl database. This function condenses a multi-step workflow 
#' implemented in this package into one step. For non-model organisms absent 
#' from Ensembl, please run the individual steps separately, as this function
#' either queries Ensembl through biomart or parses Ensembl FASTA sequence names
#' for transcript and gene information required to aggregate read counts mapped 
#' to transcripts into counts for genes. The vignette has a tutorial of running 
#' the individual steps in this workflow.
#' 
#' For 10x data sets, you can find a barcode whitelist file that comes with
#' CellRanger installation. You don't need to run CellRanger to get that. An 
#' example path to get the whitelist file is
#' \code{cellranger-2.1.0/cellranger-cs/2.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt}
#' for v2 chemistry.
#' 
#' Passing FASTA files to \code{fasta_file} is faster than passing \code{species}
#' since with the latter, this function will query the Ensembl biomart database,
#' which is really slow.
#' 
#' By default, this function does not save the data frame that maps
#' transcripts to genes. As this information, along with human readable gene
#' names corresponding to gene IDs and other information passed to \code{other_attrs}
#' that can be retrieved from Ensembl might be useful later in the analysis,
#' this function has an option to save the data frame. Set \code{save_tr2g = TRUE}
#' to save this to disk.
#' 
#' @inheritParams make_sparse_matrix
#' @inheritParams transcript2gene
#' @inheritParams EC2gene
#' @param save_tr2g Logical, whether to save the data frame that maps transcripts
#' to genes to disk.
#' @return If both gene count and TCC matrices are returned, then this function
#' returns a list with two matrices, each with genes/equivalence classes in the
#' rows and barcodes in the columns. If only one of gene count and TCC matrices
#' is returned, then a \code{dgCMatrix} with genes/equivalence classes in the 
#' rows and barcodes in the columns.
#' @family functions to generate sparse matrix in one step
#' @importFrom devtools install_github
#' @export

busparse_matrix <- function(species, fasta_file, bus_path, 
                            est_ncells, est_ngenes, whitelist = NULL, 
                            gene_count = TRUE, TCC = TRUE,
                            ensembl_version = NULL, other_attrs = NULL,
                            save_tr2g = FALSE, 
                            file_save = "./tr2g_sorted.csv",
                            verbose = TRUE, ncores = 0,
                            progress_unit = 5e6, ...) {
  bus_path <- normalizePath(bus_path, mustWork = TRUE)
  kallisto_out_path <- dirname(bus_path)
  #bus_fn <- basename(bus_path)
  # Get transcript and gene information
  tr2g <- transcript2gene(species, fasta_file, kallisto_out_path, 
                          other_attrs, ensembl_version,
                          save_tr2g, file_save, verbose, ...)
  genes <- EC2gene(tr2g, kallisto_out_path, ncores = ncores, verbose = verbose)
  make_sparse_matrix(bus_path, tr2g, est_ncells = est_ncells, 
                     est_ngenes = est_ngenes, whitelist = whitelist,
                     gene_count = gene_count, TCC = TCC, ncores = ncores,
                     verbose = verbose, progress_unit = progress_unit)
}

# To do: 
# Unit test one step functions, write examples for all exported functions.
# Unit test tr2g_ensembl related functions.
# Unit test file saving
