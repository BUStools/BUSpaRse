#' @include tr2g.R
NULL

#' Convert Latin species name to dataset name
#' 
#' This function converts Latin species name to a dataset name in biomart to
#' query gene and transcript ID.
#' 
#' @inheritParams transcript2gene
species2dataset <- function(species) { 
  species <- strsplit(species, " ")[[1]]
  if (length(species) != 2) {
    stop("Please use the Latin binomial convention for species rather than the colloquial name.")
  }
  species[1] <- tolower(substr(species[1], 1, 1))
  species <- paste(species, collapse = "")
  paste(species, "gene", "ensembl", sep = "_")
}