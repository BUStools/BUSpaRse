#' Map Transcript ID to Gene ID
#' 
#' In \code{kallisto bus}, the Ensembl transcript ID (with version number) is
#' used, but usually users are interested in number of UMI per gene rather than
#' transcript. This function queries biomart to get the gene IDs (with version
#' number) corresponding to the transcripts in the transcriptome assembly. 
#' 
#' This function applies only to species that have gene and transcript IDs on
#' Ensembl.
#' 
#' @param species A character vector of Latin names of species present in this
#' scRNA-seq dataset.
#' @param kallisto_out_path Path to the \code{kallisto bus} output directory.
#' @return A data frame with two columns: \code{gene} and \code{transcript},
#' with Ensembl gene and transcript IDs (with version number), in the same order
#' as in the transcriptome index used in \code{kallisto}.
#' @importFrom biomaRt useMart getBM
#' @importFrom data.table fread :=
#' @export
transcript2gene <- function(species, kallisto_out_path) {
  marts <- lapply(species, function(x) {
    mart_name <- species2dataset(x)
    useMart("ensembl", mart_name)
  })
  tr2g_list <- lapply(marts, function(mart) {
    getBM(attributes = c("ensembl_gene_id_version", 
                         "ensembl_transcript_id_version"), 
          mart = mart)
  })
  tr2g <- Reduce(rbind, tr2g_list)
  names(tr2g) <- c("gene", "transcript")
  # Get transcripts used in kallisto
  path_use <- normalizePath(kallisto_out_path)
  trs <- fread(paste(path_use, "transcripts.txt", sep = "/"),
               col.names = "transcript")
  # Subset and sort tr2g
  merge(trs, tr2g, all = FALSE, sort = FALSE)
}

#' Map EC Index to Genes Compatible with the EC
#' 
#' In the output file \code{output.bus}, equivalence classes (EC) are denoted by
#' an index, which is related to the set of transcripts the EC is compatible to 
#' in the output file \code{matrix.ec}. This function further relates the set of
#' transcripts to the set of genes the EC is compatible to. This function first
#' reads in \code{matrix.ec}, and then translates the transcripts into genes.
#' 
#' The data frame passed to \code{tr2g} can be generated from function 
#' \code{transcript2gene} in this package for any organism that has gene and
#' transcript ID on Ensembl. The output of this function should be passed to
#' \code{make_sparse_matrix} in the next step in the workflow, which will 
#' produce the sparse matrix that can be used in \code{Seurat}.
#' 
#' @inheritParams transcript2gene
#' @param tr2g A Data frame with columns \code{gene} and \code{transcript}, in
#' the same order as in the transcriptome index for \code{kallisto}.
#' @param ncores Number of cores to use, defaults to 1.
#' @param verbose Whether to display progress, though this function usually
#' does not take very long. Defaults to \code{TRUE}.
#' @return A list each element of whom is the set of genes the corresponding EC
#' is compatible to. The genes are in Ensembl ID with version number. The 
#' elements of this list are in the same order as the ECs listed in the
#' \code{kallisto bus} output file \code{matrix.ec}. 
#' @seealso \code{\link{transcript2gene}}
#' @importFrom parallel mclapply
#' @export
#' 
EC2gene <- function(tr2g, kallisto_out_path, ncores = 1, verbose = TRUE) {
  genes <- tr2g$gene
  # Read in matrix.ec
  if (verbose) cat("Reading matrix.ec\n")
  path_use <- normalizePath(kallisto_out_path, mustWork = TRUE)
  ECs <- fread(paste(path_use, "matrix.ec", sep = "/"), 
               col.names = c("EC_index", "EC"),
               data.table = TRUE, showProgress = verbose)
  if (verbose) cat("Processing genes\n")
  # Prevent R CMD check note no visible binding for global variable
  EC_index <- EC <- NULL
  ECs[, c("EC_index", "EC") := list(EC_index, 
                                 strsplit(EC, ","))]
  ECs[, genes := mclapply(EC, 
                          function(x) {
                            inds <- as.integer(x) + 1
                            unique(genes[inds])
                          }, mc.cores = ncores)]
  ECs$genes
}
