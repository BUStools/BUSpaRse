#' @include sparse_matrix.R
NULL

#' Get transcript and gene info from Ensembl
#' 
#' The memoised version of this function is exported, since biomart queries are
#' really slow.
#' 
#' @param species Character vector of length 1, Latin name of the species of
#' interest.
#' @param ensembl_version Integer version number of Ensembl (e.g. 94 for the
#' October 2018 release). This argument defaults to \code{NULL}, which will use
#' the current release of Ensembl. Use 
#' \code{\link[biomaRt]{listEnsemblArchives}} to see the version number corresponding
#' to the Ensembl release of a particular date.
#' @param other_attrs Character vector. Other attributes to get from Ensembl, 
#' such as gene symbol and position on the genome.
#' @param verbose Whether to display progress.
#' @param \dots Othe arguments to be passed to \code{\link[biomaRt]{useEnsembl}},
#' such as host and mirror.
#' Use \code{\link[biomaRt]{listAttributes}} to see which attributes are available.
#' @importFrom biomaRt useMart getBM
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. If \code{other_attrs}
#' has been specified, then those will also be columns in the data frame returned.
#' @seealso \code{\link{tr2g_ensembl}}
#' 
.tr2g_ensembl <- function(species, other_attrs = NULL, ensembl_version = NULL, 
                          verbose = TRUE, ...) {
  # Validate arguments
  if (!is.null(ensembl_version) && !is.numeric(ensembl_version)) {
    stop("ensembl_version must be integer.\n")
  }
  if (!is.null(other_attrs) && 
      (!is.atomic(other_attrs) || !is.character(other_attrs))) {
    stop("other_attrs must be an atomic character vector.\n")
  }
  mart_name <- species2dataset(species)
  if (verbose) {
    message("Querying biomart\n")
  }
  mart <- useEnsembl(biomart = "ensembl", dataset = mart_name, 
                     version = ensembl_version, ...)
  out <- getBM(c("ensembl_gene_id_version", 
                 "ensembl_transcript_id_version",
                 "external_gene_name",
                 other_attrs), mart = mart)
  names(out)[1:3] <- c("gene", "transcript", "gene_name")
  out
}

#' Get transcript and gene info from Ensembl
#' 
#' This function queries Ensembl biomart to convert transcript IDs to gene IDs.
#' 
#' The output from this function can be cached, so subsequent calls with exactly
#' the same parameters will retrieve result from the cache, unless the cache is 
#' cleared. To clear cache, see \code{\link[R.cache]{clearCache}}.
#' 
#' @inheritParams .tr2g_ensembl
#' @param cache Logical, whether to cache results, defaults to \code{TRUE}.
#' @importFrom biomaRt useMart getBM
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. The gene and transcript
#' IDs will have the version number. If \code{other_attrs} has been specified, 
#' then those will also be columns in the data frame returned.
#' @importFrom R.cache addMemoization
#' @family functions to retrieve transcript and gene info
#' @seealso \code{\link{tr2g_gtf}} \code{\link{tr2g_fasta}} \code{\link{transcript2gene}}
#' @export
tr2g_ensembl <- function(species, other_attrs = NULL, ensembl_version = NULL, 
                         verbose = TRUE, cache = TRUE, ...) {
  if (cache) {
    f <- addMemoization(.tr2g_ensembl)
  } else {
    f <- .tr2g_ensembl
  }
  f(species, other_attrs, ensembl_version, verbose, ...)
}

#' Get transcript and gene info from GTF file
#' 
#' This function reads a GTF file and extracts the transcript ID and
#' corresponding gene ID. The memoised version of this function is exported 
#' since reading large gtf files can take a while.
#' 
#' Transcript and gene versions may not be present in all GTF files, so these
#' arguments are optional. This function has arguments for transcript and gene 
#' version numbers because Ensembl IDs have version numbers. For Ensembl IDs, we
#' recommend including the version number, since a change in version number 
#' signals a change in the entity referred to by the ID after reannotation. If a
#' version is used, then it will be appended to the ID, separated by 
#' \code{version_sep}.
#' 
#' @param file Path to a GTF file to be read. The file can remain gzipped.
#' @param type_use Character vector, the values taken by the \code{type} field 
#' in the GTF file that denote the desired transcripts. This can be "exon",
#' "transcript", "mRNA", and etc.
#' @param transcript_id Character vector of length 1. Tag in \code{attribute} 
#' field corresponding to transcript IDs. This argument must be supplied and
#' cannot be \code{NA} or \code{NULL}. Will throw error if tag indicated in this
#' argument does not exist.
#' @param gene_id Character vector of length 1. Tag in \code{attribute} 
#' field corresponding to gene IDs. This argument must be supplied and
#' cannot be \code{NA} or \code{NULL}. Note that this is different from gene
#' symbols, which do not have to be unique. This can be Ensembl or Entrez IDs.
#' However, if the gene symbols are in fact unique for each gene, you may
#' supply the tag for human readable gene symbols to this argument. Will throw 
#' error if tag indicated in this argument does not exist.
#' @param gene_name Character vector of length 1. Tag in \code{attribute} 
#' field corresponding to gene symbols. This argument can be \code{NA} or
#' \code{NULL} if you are fine with non-human readable gene IDs. 
#' @param transcript_version Character vector of length 1. Tag in \code{attribute} 
#' field corresponding to _transcript_ version number. If your GTF file does not
#' include transcript version numbers, or if you do not wish to include the
#' version number, then use \code{NULL} for this argument. 
#' @param gene_version Character vector of length 1. Tag in \code{attribute} 
#' field corresponding to _gene_ version number. If your GTF file does not
#' include gene version numbers, or if you do not wish to include the
#' version number, then use \code{NULL} for this argument. 
#' @param version_sep Character to separate bewteen the main ID and the version
#' number. Defaults to ".", as in Ensembl.
#' @inheritParams .tr2g_ensembl
#' @note The transcript and gene IDs are The \code{attribute} field (the last 
#' field) of GTF files can be complicated and inconsistent across different 
#' sources. Please check the \code{attribute} tags in your GTF file and consider
#' the arguments of this function carefully. The defaults are set according to 
#' Ensembl GTF files; defaults may not work for files from other sources. Due to
#' the general lack of standards for the \code{attribute} field, you may need to
#' further clean up the output of this function.
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. 
#' @importFrom plyranges read_gff filter
#' @importFrom magrittr %>% 
#' @importFrom stringr str_detect
#' @importFrom dplyr distinct
#' @importFrom purrr walk2
#' @importFrom S4Vectors mcols
#' @seealso \code{\link{tr2g_gtf}}
#' 
.tr2g_gtf <- function(file, type_use = "exon", transcript_id = "transcript_id",
                      gene_id = "gene_id", gene_name = "gene_name",
                      transcript_version = "transcript_version",
                      gene_version = "gene_version", version_sep = ".",
                      verbose = TRUE) {
  # Check that some arguments are length 1 character vectors
  args_used <- get_args()
  args_check <- c("file", "transcript_id", "gene_id", "gene_name",
                  "transcript_version", "gene_version", "version_sep")
  check_char1(unlist(args_used[args_check]))
  
  if (!str_detect(file, "\\.gtf")) {
    stop("file must be a GTF file.\n")
  }
  file <- normalizePath(file, mustWork = TRUE)
  if (verbose) {
    message("Reading GTF file\n")
  }
  gr <- read_gff(file)
  tags <- names(mcols(gr))
  if (is.null(transcript_id)) stop("transcript_id cannot be NULL.\n")
  if (is.null(gene_id)) stop("gene_id cannot be NULL.\n")
  check_tag_present(c(transcript_id, gene_id), tags, error = TRUE)
  # Will do nothing if all are NULL
  check_tag_present(c(gene_name, transcript_version, gene_version), 
                    tags, error = FALSE)
  
  gr <- gr %>% 
    filter(!is.na(transcript_id), type %in% type_use)
  if (length(gr) == 0) {
    stop(paste("No transcript has types", paste(type_use, collapse = ", "), 
               ".\n"))
  }
  
  out <- data.frame(gene = mcols(gr)[[gene_id]],
                    transcript = mcols(gr)[[transcript_id]],
                    stringsAsFactors = FALSE)
  if (!is.null(gene_name)) {
    out$gene_name <- mcols(gr)[[gene_name]]
  }
  if (!is.null(transcript_version)) {
    tv <- mcols(gr)[[transcript_version]]
    out$transcript <- paste(out$transcript, tv, sep = version_sep)
  }
  if (!is.null(gene_version)) {
    gv <- mcols(gr)[[gene_version]]
    out$gene <- paste(out$gene, gv, sep = version_sep)
  }
  distinct(out)
}

#' Get transcript and gene info from GTF file
#' 
#' This function reads a GTF file and extracts the transcript ID and
#' corresponding gene ID. The memoised version of this function is exported 
#' since reading large gtf files can take a while.
#' 
#' Transcript and gene versions may not be present in all GTF files, so these
#' arguments are optional. This function has arguments for transcript and gene 
#' version numbers because Ensembl IDs have version numbers. For Ensembl IDs, we
#' recommend including the version number, since a change in version number 
#' signals a change in the entity referred to by the ID after reannotation.
#' 
#' @inheritParams .tr2g_gtf
#' @inheritParams tr2g_ensembl
#' @note The transcript and gene IDs are The \code{attribute} field (the last 
#' field) of GTF files can be complicated and inconsistent across different 
#' sources. Please check the \code{attribute} tags in your GTF file and consider
#' the arguments of this function carefully. The defaults are set according to 
#' Ensembl GTF files; defaults may not work for files from other sources.
#' @seealso \code{\link{tr2g_ensembl}} \code{\link{tr2g_fasta}} \code{\link{transcript2gene}}
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. 
#' @family functions to retrieve transcript and gene info
#' @export
tr2g_gtf <- function(file, type_use = "exon", transcript_id = "transcript_id",
                     gene_id = "gene_id", gene_name = "gene_name",
                     transcript_version = "transcript_version",
                     gene_version = "gene_version", version_sep = ".",
                     verbose = TRUE, cache = TRUE) {
  if (cache) {
    f <- addMemoization(.tr2g_gtf)
  } else {
    f <- .tr2g_gtf
  }
  f(file, type_use, transcript_id, gene_id, gene_name, transcript_version,
    gene_version, version_sep, verbose)
}

# To do: tr2g_gff

#' Get transcript and gene info from names in FASTA files
#' 
#' FASTA files, such as those for cDNA and ncRNA from Ensembl, have genome
#' annotation information in the name of each sequence entry. This function
#' extracts the transcript and gene IDs from such FASTA files. 
#' 
#' At present, this function only works with FASTA files from Ensembl, and uses
#' regex to extract Ensembl IDs. Sequence names should be formatted as follows:
#' 
#' ```
#' ENST00000632684.1 cdna chromosome:GRCh38:7:142786213:142786224:1 
#' gene:ENSG00000282431.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene 
#' gene_symbol:TRBD1 description:T cell receptor beta diversity 1 
#' [Source:HGNC Symbol;Acc:HGNC:12158]
#' ```
#' 
#' If your FASTA file sequence names are formatted differently, then you must
#' extract the transcript and gene IDs by some other means. The Bioconductor
#' package \code{Biostrings} is recommended; after reading the FASTA file into 
#' R, the sequence names can be accessed by the \code{names} function.
#' 
#' @inheritParams tr2g_ensembl
#' @param file Path to the FASTA file to be read. The file can remain gzipped.
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. 
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_extract
#' @seealso \code{\link{tr2g_fasta}}
#' 
.tr2g_fasta <- function(file, verbose = TRUE) {
  if (!str_detect(file, "(\\.fasta)|(\\.fa)|(\\.fna)")) {
    stop("file must be a FASTA file.\n")
  }
  file <- normalizePath(file, mustWork = TRUE)
  s <- readDNAStringSet(file)
  out <- data.frame(gene = str_extract(names(s), "ENS[A-Z]*G\\d+\\.\\d+"),
                    transcript = str_extract(names(s), "ENS[A-Z]*T\\d+\\.\\d+"),
                    gene_name = str_extract(names(s), "gene_symbol:[a-zA-Z\\d-\\.]+"),
                    stringsAsFactors = FALSE) %>% 
    separate(gene_name, into = c("g", "gene_name"), sep = ":") %>% 
    dplyr::select(-g) %>% 
    distinct()
  out
}

#' Get transcript and gene info from names in FASTA files
#' 
#' FASTA files, such as those for cDNA and ncRNA from Ensembl, have genome
#' annotation information in the name of each sequence entry. This function
#' extracts the transcript and gene IDs from such FASTA files. 
#' 
#' At present, this function only works with FASTA files from Ensembl, and uses
#' regex to extract Ensembl IDs. Sequence names should be formatted as follows:
#' 
#' ```
#' ENST00000632684.1 cdna chromosome:GRCh38:7:142786213:142786224:1 
#' gene:ENSG00000282431.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene 
#' gene_symbol:TRBD1 description:T cell receptor beta diversity 1 
#' [Source:HGNC Symbol;Acc:HGNC:12158]
#' ```
#' 
#' If your FASTA file sequence names are formatted differently, then you must
#' extract the transcript and gene IDs by some other means. The Bioconductor
#' package \code{Biostrings} is recommended; after reading the FASTA file into 
#' R, the sequence names can be accessed by the \code{names} function.
#' 
#' @inheritParams tr2g_ensembl
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. 
#' @seealso \code{\link{tr2g_ensembl}} \code{\link{tr2g_gtf}} \code{\link{transcript2gene}}
#' @family functions to retrieve transcript and gene info
tr2g_fasta <- function(file, verbose = TRUE, cache = TRUE) {
  if (cache) {
    f <- addMemoization(.tr2g_fasta)
  } else {
    f <- .tr2g_fasta
  }
  f(file, verbose)
}

#' Map Transcript ID to Gene ID
#' 
#' This function first retrieves gene and transcript ID information with the
#' \code{tr2g_*} family of functions, and then sorts the transcripts so they are
#' in the same order as in the \code{kallisto} index.
#' 
#' You can supply both species names and file paths. Eventually, all the
#' transcript and gene information will be combined regardless of source and
#' sorted.
#' 
#' The \code{tr2g_*} family of functions called by this function
#' have the option to cache results, determined by the \code{cache} argument
#' of this function. To clear cache, see \code{\link[R.cache]{clearCache}}.
#' 
#' @inheritParams tr2g_ensembl
#' @param species A character vector of Latin names of species present in this
#' scRNA-seq dataset. This is used to retrieve Ensembl information from biomart.
#' @param files A character vector of paths to GTF/GFF/FASTA files. The files 
#' can remain gzipped.
#' @param kallisto_out_path Path to the \code{kallisto bus} output directory.
#' @param verbose Whether to display progress. Defaults to \code{TRUE}.
#' @return A data frame with two columns: \code{gene} and \code{transcript},
#' with Ensembl gene and transcript IDs (with version number), in the same order
#' as in the transcriptome index used in \code{kallisto}.
#' @importFrom data.table fread := rbindlist
#' @export
#' @family functions to retrieve transcript and gene info
#' @seealso \code{\link{tr2g_ensembl}} \code{\link{tr2g_gtf}} \code{\link{tr2g_fasta}}
#' @examples
#' \dontrun{
#' # Download dataset already in BUS format
#' library(TENxhgmmBUS)
#' download_hgmm(".", "hgmm100")
#' tr2g <- transcript2gene(c("Homo sapiens", "Mus musculus"),
#'                         kallisto_out_path = "./out_hgmm100")
#' }
#' 
transcript2gene <- function(species, files, kallisto_out_path, 
                            other_attrs = NULL, ensembl_version = NULL,
                            verbose = TRUE, cache = TRUE, ...) {
  kallisto_out_path <- normalizePath(kallisto_out_path, mustWork = TRUE)
  ensembl_data <- file_data <- NULL
  if (!missing(species)) {
    ensembl_dfs <- lapply(species, tr2g_ensembl, other_attrs = other_attrs, 
                  ensembl_version = ensembl_version, 
                  verbose = verbose, cache = cache, ...)
    ensembl_data <- rbindlist(ensembl_dfs) %>% unique()
  }
  if (!missing(files)) {
    if (any(!str_detect(files, "(\\.gff)|(\\.gtf)|(\\.fa)|(\\.fasta)|(\\.fna)"))) {
      stop("Files must be GTF/GFF/FASTA.\n")
    }
    file_dfs <- lapply(files, function(x) {
      if (str_detect(x, "(\\.gff)|(\\.gtf)")) {
        tr2g_gtf(x, verbose = verbose, cache = cache)
      } else {
        tr2g_fasta(x, verbose = verbose, cache = cache)
      }
    })
    file_data <- rbindlist(file_dfs) %>% unique()
  }
  if (!is.null(ensembl_data) && !is.null(file_data)) {
    out <- merge(ensembl_data, file_data, 
                 by = c("gene", "transcript", "gene_name"),
                 all = TRUE) %>% unique(by = c("gene", "transcript", "gene_name"))
  }
  
  # Sort according to the order in kallisto index
  if (verbose) {
    message("Sorting\n")
  }
  trs <- fread(paste(kallisto_out_path, "transcripts.txt", sep = "/"),
               col.names = "transcript")
  merge(trs, out, by = c("gene", "transcript"), sort = FALSE)
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
#' produce the sparse matrix that can be used in 
#' \code{\href{https://satijalab.org/seurat/}{Seurat}}.
#' 
#' @inheritParams transcript2gene
#' @param tr2g A Data frame with columns \code{gene} and \code{transcript}, in
#' the same order as in the transcriptome index for \code{kallisto}.
#' @param ncores Number of cores to use, defaults to 1.
#' @return A list each element of whom is the set of genes the corresponding EC
#' is compatible to. The genes are in Ensembl ID with version number. The 
#' elements of this list are in the same order as the ECs listed in the
#' \code{kallisto bus} output file \code{matrix.ec}. 
#' @seealso \code{\link{transcript2gene}}
#' @importFrom parallel mclapply
#' @export
#' @examples
#' \dontrun{
#' # Download dataset already in BUS format
#' library(TENxhgmmBUS)
#' download_hgmm(".", "hgmm100")
#' tr2g <- transcript2gene(c("Homo sapiens", "Mus musculus"),
#'                         kallisto_out_path = "./out_hgmm100")
#' genes <- EC2gene(tr2g, "./out_hgmm100", ncores = 1, verbose = FALSE)
#' }
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
