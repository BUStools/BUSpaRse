err0 <- function(n) {
  if (n == 0) {
    stop("Chromosone names in annotation and genome do not overlap. ",
         "Consider checking that formatting of chromosome names of ",
         "the annotation and the genome match.")
  }
}

msg1 <- function(n) {
  message("There are ", n, " chromosomes in the annotation that are absent ",
          "from the genome. Consider checking that formatting of chromosome ",
          "names of the annotation and the genome match.")
}

msg2 <- function(n) {
  message("There are ", n, " chromosomes in the genome are absent from the ",
          "annotation. Consider checking that formatting of chromosome names ",
          "of the annotation and the genome match.")
}

#' Remove chromosomes in anotation absent from genome
#' 
#' @param chrs_use Character vector of names of chromosomes present in both the
#' annotation and the genome.
#' @param annot Either a `GRanges` object or a `TxDb` object for gene 
#' annotation.
#' @return A subsetted genome annotation of the same type ofo the input genome
#' annotation.
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicFeatures seqlevels<-
sub_annot <- function(chrs_use, annot) {
  err0(length(chrs_use))
  n <- length(seqlevels(annot)) - length(chrs_use)
  if (n > 0) {
    msg1(n)
    if (is(annot, "GRanges")) {
      annot <- annot[seqnames(annot) %in% chrs_use]
    } else if (is(annot, "TxDb")) {
      seqlevels(annot) <- chrs_use
    }
  }
  return(annot)
}

#' Check for chromosomes in genome but not annotation
#' 
#' @inheritParams get_velocity_files
#' @inheritParams sub_annot
#' @return Nothing. Will emit message if the genome contains chromosomes absent
#' from the annotation.
check_genome <- function(chrs_use, genome) {
  if (is(genome, "DNAStringSet")) {
    n <- length(genome)
  } else if (is(genome, "BSgenome")) {
    n <- length(seqlevels(genome))
  }
  n2 <- n - length(chrs_use)
  if (n2 > 0) {
    msg2(n2)
  }
}

#' Subset genome annotation
#' 
#' Exclude chromosomes present in the annotation but absent from the genonme.
#' 
#' @inheritParams get_velocity_files
#' @inheritParams sub_annot
#' @return A subsetted genome annotation of the same type ofo the input genome
#' annotation.
setGeneric("subset_annot", 
           function(genome, annot)
             standardGeneric("subset_annot"))

#' @rdname subset_annot
setMethod("subset_annot", signature = c("DNAStringSet", "ANY"),
          function(genome, annot) {
            chrs_use <- intersect(seqlevels(annot), names(genome))
            annot <- sub_annot(chrs_use, annot)
            check_genome(chrs_use, genome)
            return(annot)
          })

#' @rdname subset_annot
setMethod("subset_annot", signature = c("BSgenome", "ANY"),
          function(genome, annot) {
            chrs_use <- intersect(seqlevels(annot), seqlevels(genome))
            annot <- sub_annot(chrs_use, annot)
            check_genome(chrs_use, genome)
            return(annot)
          })

#' Extract transcriptomic sequences from genome
#' 
#' This is a thin wrapper around \code{\link{extractTranscriptSeqs}}. This 
#' wrapper takes care of object type of the genome.
#' 
#' @inheritParams get_velocity_files
#' @inheritParams subset_annot
#' @param \dots Extra argument for methods.
#' @return A `DNAStringSet` object with transcript sequences, whose names are
#' the transcript IDs. 
setGeneric("extract_tx", 
           function(genome, annot, ...) 
             standardGeneric("extract_tx"))

#' @rdname extract_tx
#' @importFrom GenomicRanges sort
#' @importFrom GenomicFeatures extractTranscriptSeqs
setMethod("extract_tx", signature = c("ANY", "GRanges"),
          function(genome, annot, transcript_id = "transcript_id") {
            annot <- subset_annot(genome, annot)
            gr <- sort(annot[annot$type == "exon"])
            gr <- split(gr, mcols(gr)[[transcript_id]])
            extractTranscriptSeqs(genome, gr)
          })

#' @rdname extract_tx
#' @importFrom GenomicFeatures exonsBy
setMethod("extract_tx", signature = c("ANY", "TxDb"),
          function(genome, annot) {
            annot <- subset_annot(genome, annot)
            gr <- exonsBy(annot, by = "tx", use.names = TRUE)
            extractTranscriptSeqs(genome, gr)
          })
