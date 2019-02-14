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
#' For 10x data sets, you can find a barcode whitelist file that comes with
#' CellRanger installation. You don't need to run CellRanger to get that. An 
#' example path to get the whitelist file is
#' \code{cellranger-2.1.0/cellranger-cs/2.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt}
#' for v2 chemistry.
#' 
#' @param bus_fn File name of output of \code{bustools}, namely \code{output.bus}
#' from \code{kallisto bus} after sorted and converted to text by
#'  \code{bustools}.
#' @param genes A list with each element a string vector of genes that an 
#' equivalence class maps to, generated earlier in \code{EC2gene}.
#' @param whitelist A character vector with valid cell barcodes. This is an
#' optional argument, that defaults to \code{NULL}. When it is \code{NULL},
#' all cell barcodes present will be included in the sparse matrix, not only
#' the barcodes known to be valid.
#' @param est_ncells Estimated number of cells; providing this argument will
#' speed up computation as it minimizes memory reallocation as vectors grow.
#' @param est_ngenes Estimated number of genes.
#' @param verbose Whether to display progress.
#' @param progress_unit How many iteration to print one progress update when
#' reading in the \code{kallisto bus} file.
#' @return A sparse matrix with genes in rows and cells in columns.
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
#' genes <- EC2gene(tr2g, "./out_hgmm100", ncores = 1, verbose = FALSE)
#' res_mat <- make_sparse_matrix("./out_hgmm100/output.sorted.txt",
#'                               genes = genes, est_ncells = 3e5, 
#'                               est_ngenes = nrow(tr2g))
#' # Remove empty droplets
#' tot_counts <- colSums(res_mat)
#' res_mat <- res_mat[,tot_counts > 500]
#' }

make_sparse_matrix <- function(bus_fn, genes, est_ncells, est_ngenes, 
                               whitelist = NULL, 
                               verbose = TRUE,
                               progress_unit = 5e6) {
  bus_fn <- normalizePath(bus_fn, mustWork = TRUE)
  if (!grepl(".txt$", bus_fn)) {
    stop("Argument fn must point to a text file. Please run bustools text.")
  }
  if (is.null(whitelist)) {
    whitelist <- ""
  }
  # Prevent the no visible binding of global variable note in R CMD check
  barcodes <- geneIDs <- NULL
  c(res_mat, barcodes, geneIDs) %<-% fill_cell_gene(bus_fn, genes, 
                                                  est_ncells, est_ngenes, 
                                                  whitelist, 
                                                  verbose, progress_unit)
  rownames(res_mat) <- geneIDs
  colnames(res_mat) <- barcodes
  res_mat
}

#' Get gene count matrix in one step
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
#' @inheritParams make_sparse_matrix
#' @inheritParams transcript2gene
#' @inheritParams EC2gene
#' @param fasta_file Character vector of paths to the transcriptome FASTA files
#' used to build the kallisto index. If multiple FASTA files were used for the
#' kallisto index, then the order of paths here must match that of the FASTA
#' files when building the kallisto index. Exactly one of \code{species} and
#' \code{fasta_file} can be missing.
#' @param save_tr2g Logical, whether to save the data frame that maps transcripts
#' to genes to disk.
#' @note By default, this function does not save the data frame that maps
#' transcripts to genes. As this information, along with human readable gene
#' names corresponding to gene IDs and other information passed to \code{other_attrs}
#' that can be retrieved from Ensembl might be useful later in the analysis,
#' this function has an option to save the data frame. Set \code{save_tr2g = TRUE}
#' to save this to disk.
#' @return A sparse matrix with genes in rows and cells in columns.
#' @family functions to generate sparse matrix in one step
#' @export

busparse_gene_count <- function(species, fasta_file, kallisto_out_path, bus_fn,
                                est_ncells, est_ngenes, whitelist = NULL, 
                                ensembl_version = NULL, other_attrs = NULL,
                                save_tr2g = FALSE, 
                                file_save = "./tr2g_sorted.csv",
                                verbose = TRUE, ncores = 1,
                                progress_unit = 5e6, ...) {
  if (!xor(missing(species), missing(fasta_file))) {
    stop("Exactly one of species and fasta_file can be missing.\n")
  }
  # Get transcript and gene information
  if (missing(fasta_file)) {
    tr2g <- transcript2gene(species, kallisto_out_path,
                            other_attrs = other_attrs,
                            ensembl_version = ensembl_version,
                            save = save_tr2g, file_save = file_save, 
                            verbose = verbose, ...)
  } else {
    fls <- lapply(fasta_file, tr2g_fasta, verbose = verbose)
    tr2g <- rbindlist(fls)
  }
  genes <- EC2gene(tr2g, kallisto_out_path, ncores = ncores, verbose = verbose)
  make_sparse_matrix(bus_fn, genes, est_ncells = est_ncells, 
                     est_ngenes = est_ngenes, whitelist = whitelist,
                     verbose = verbose, progress_unit = progress_unit)
}

# To do: Make TCC matrix, parallelize read count processing
# Unit test one step functions, write examples for all exported functions.
# Unit test tr2g_ensembl related functions.
# Unit test file saving
