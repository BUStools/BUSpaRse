#' @include tr2g.R
NULL

#' Convert Latin species name to dataset name
#'
#' This function converts Latin species name to a dataset name in biomart to
#' query gene and transcript ID.
#'
#' @inheritParams transcript2gene
#' @return The appropriate dataset name for biomart.
#' @export
#' @examples
#' species2dataset(species = "Homo sapiens")
species2dataset <- function(species, type = c("vertebrate", "metazoa", "plant",
                              "fungus", "protist")) {
  type <- match.arg(type)
  species <- strsplit(species, " ")[[1]]
  if (length(species) != 2) {
    stop("Please use the Latin binomial convention for species rather than the colloquial name.\n")
  }
  species[1] <- tolower(substr(species[1], 1, 1))
  species <- paste(species, collapse = "")
  if (type == "vertebrate") {
    return(paste(species, "gene_ensembl", sep = "_"))
  } else {
    return(paste(species, "eg_gene", sep = "_"))
  }
}

#' Check that a tag is present in attribute field of GTF/GFF
#'
#' The attribute field of GTF/GFF files are very complicated and is very
#' inconsistent between sources. This function is to make sure that transcript
#' and gene IDs can be extracted properly.
#'
#' @param tags_use The tags to be checked.
#' @param tags The tags present in attribute field.
#' @param error Whether to throw an error for absent tags. If \code{FALSE}, then
#' a warning will be given.
#' @return Error or warning if tag is absent.
check_tag_present <- function(tags_use, tags, error = TRUE) {
  present <- tags_use %in% tags
  if (!all(present)) {
    s <- paste("Tags", paste(tags_use[!present], collapse = ", "),
      "are absent from the attribute field.\n")
    if (error) {
      stop(s)
    } else {
      warning(paste(s, "These tags are ignored.\n"))
    }
  }
}

#' Check that an object is a character vector of length 1
#'
#' Just in case the user passes something with length more than 1 and messes up
#' everything thanks to vectorization.
#'
#' @param x Named vector of arguments to be checked.
#' @return Error if \code{x} is not a character vector with length 1.
check_char1 <- function(x) {
  arg_names <- names(x)
  inds <- vapply(x, function(x) {
    !is.character(x) | length(x) > 1
  },
  FUN.VALUE = logical(1))
  if (any(inds)) {
    stop(paste(paste(arg_names[inds], sep = ", "),
      "must be a character vector with length 1.\n"))
  }
}

#' Check inputs to tr2g_gtf and tr2g_gff3
#'
#' This function validates inputs to tr2g_gtf and tr2g_gff3 and throws error
#' early if some inputs are wrong.
#'
#' @inheritParams tr2g_gtf
#' @param format Whether it's gtf or gff3.
#' @return Nothing, will throw error if there's a problem.

check_gff <- function(format, file, transcript_id, gene_id) {
  if (is.null(transcript_id)) stop("transcript_id cannot be NULL.\n")
  if (is.null(gene_id)) stop("gene_id cannot be NULL.\n")

  if (!str_detect(file, paste0("\\.", format))) {
    stop(paste("file must be a", toupper(format), "file.\n"))
  }
}

#' Read matrix along with barcode and gene names
#'
#' This function takes in a directory and name and reads the mtx file, genes,
#' and barcodes from the output of `bustools` to return a sparse matrix with
#' column names and row names.
#'
#' @param dir Directory with the bustools count outputs.
#' @param name The files in the output directory should be <name>.mtx, <name>.genes.txt,
#' and <name>.barcodes.txt.
#' @param tcc Logical, whether the matrix of interest is a TCC matrix. Defaults
#' to \code{FALSE}.
#' @return A dgCMatrix with barcodes as column names and genes as row names.
#' @importFrom Matrix readMM
#' @importFrom methods as
#' @export
#' @examples
#' # Internal toy data used for unit testing
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' m <- read_count_output(toy_path, name = "genes", tcc = FALSE)
read_count_output <- function(dir, name, tcc = TRUE) {
  dir <- normalizePath(dir, mustWork = TRUE)
  m <- readMM(paste0(dir, "/", name, ".mtx"))
  m <- Matrix::t(m)
  m <- as(m, "dgCMatrix")
  # The matrix read has cells in rows
  ge <- if (tcc) ".ec.txt" else ".genes.txt"
  genes <- fread(paste0(dir, "/", name, ge), header = FALSE)$V1
  barcodes <- fread(paste0(dir, "/", name, ".barcodes.txt"), header = FALSE)$V1
  colnames(m) <- barcodes
  rownames(m) <- genes
  return(m)
}

#' Read intronic and exonic matrices into R
#'
#' @param spliced_dir Directory with `mtx` file for UMI counts of spliced
#' transcripts.
#' @param unspliced_dir Directory with `mtx` file for UMI counts of unspliced
#' transcripts.
#' @param spliced_name The files in the splicedd directory should be
#' <spliced_name>.mtx, <spliced_name>.genes.txt, and
#' <spliced_name>.barcodes.txt.
#' @param unspliced_name The files in the unsplicedd directory should be
#' <unspliced_name>.mtx, <unspliced_name>.genes.txt, and
#' <unspliced_name>.barcodes.txt.
#' @return A list of two dgCMatrix with barcodes as column names and genes as
#' row names. The elements of the list will be `spliced` and `unspliced`.
#' @export
#' @examples
#' # Internal toy data used for unit testing
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' m <- read_velocity_output(toy_path, toy_path,
#'   spliced_name = "genes",
#'   unspliced_name = "genes")
read_velocity_output <- function(spliced_dir, unspliced_dir, spliced_name,
                                 unspliced_name) {
  spliced <- read_count_output(spliced_dir, spliced_name, FALSE)
  unspliced <- read_count_output(unspliced_dir, unspliced_name, FALSE)
  list(spliced = spliced, unspliced = unspliced)
}
