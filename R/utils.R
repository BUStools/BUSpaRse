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
species2dataset <- function(species) { 
  species <- strsplit(species, " ")[[1]]
  if (length(species) != 2) {
    stop("Please use the Latin binomial convention for species rather than the colloquial name.\n")
  }
  species[1] <- tolower(substr(species[1], 1, 1))
  species <- paste(species, collapse = "")
  paste(species, "gene", "ensembl", sep = "_")
}

#' Check that a tag is present in attribute field of GTF/GFF
#' 
#' The attribute field of GTF/GFF files are very complicated and is very
#' inconsistent between sources. This function is to make sure that transcript 
#' and gene IDs can be extracted properly.
#' 
#' @param tags_use The tags to be checked.
#' @param tags The tags present in attribute field.
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

#' Get values passed to arguments in a function
#' 
#' This function is for giving informative error messages about which argument
#' is wrong. This is from the StackOverflow question
#' \url{https://stackoverflow.com/questions/17256834/getting-the-arguments-of-a-parent-function-in-r-with-names}
#' @return A named list mapping arguments to values passed to them.
get_args <- function() {
  as.list(match.call(
    definition = sys.function(-1),
    call = sys.call(-1)))[-1]
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
  inds <- !is.character(x) | length(x) > 1
  if (any(inds)) {
    stop(paste(paste(arg_names[inds], sep = ", "), 
               "must be a character vector with length 1.\n"))
  }
}
