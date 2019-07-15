err0 <- function() {
  stop("Chromosone names in annotation and genome do not overlap. ",
       "Consider checking that formatting of chromosome names of ",
       "the annotation and the genome match.")
}

msg1 <- function() {
  message("Some chromosomes in the annotation are absent from the genome. ",
          "Consider checking that formatting of chromosome names of ",
          "the annotation and the genome match.")
}

msg2 <- function() {
  message("Some chromosomes in the genome are absent from the annotation. ",
          "Consider checking that formatting of chromosome names of ",
          "the annotation and the genome match.")
}

#' Subset genome annotation
#' 
#' Exclude chromosomes present in the annotation but absent from the genonme.
#' 
#' @inheritParams get_velocity_files
#' @param annot Either a `GRanges` object or a `TxDb` object for gene 
#' annotation.
#' @return A subsetted genome annotation of the same type ofo the input genome
#' annotation.
#' @importFrom GenomeInfoDb seqlevels
setGeneric("subset_annot", 
           function(genome, annot)
             standardGeneric("subset_annot"))

#' @rdname subset_annot
setMethod("subset_annot", signature = c("DNAStringSet", "GRanges"),
          function(genome, annot) {
            chrs_use <- intersect(seqlevels(annot), names(genome))
            if (length(chrs_use) == 0) {
              err0()
            } else if (length(chrs_use) < length(seqlevels(annot))) {
              msg1()
              annot <- annot[seqnames(annot) %in% chrs_use]
            } 
            if (length(chrs_use) < length(genome)) {
              msg2()
            }
            return(annot)
          })

#' @rdname subset_annot
#' @importFrom GenomicFeatures seqlevels<-
setMethod("subset_annot", signature = c("DNAStringSet", "TxDb"),
          function(genome, annot) {
            chrs_use <- intersect(seqlevels(annot), names(genome))
            if (length(chrs_use) == 0) {
              err0()
            } else if (length(chrs_use) < length(seqlevels(annot))) {
              msg1()
              seqlevels(annot) <- chrs_use
            } 
            if (length(chrs_use) < length(genome)) {
              msg2()
            }
            return(annot)
          })

#' @rdname subset_annot
setMethod("subset_annot", signature = c("BSgenome", "GRanges"),
          function(genome, annot) {
            chrs_use <- intersect(seqlevels(annot), seqlevels(genome))
            if (length(chrs_use) == 0) {
              err0()
            } else if (length(chrs_use) < length(seqlevels(annot))) {
              msg1()
              annot <- annot[seqnames(annot) %in% chrs_use]
            } 
            if (length(chrs_use) < length(seqlevels(genome))) {
              msg2()
            }
            return(annot)
          })

#' @rdname subset_annot
setMethod("subset_annot", signature = c("BSgenome", "TxDb"),
          function(genome, annot) {
            chrs_use <- intersect(seqlevels(annot), seqlevels(genome))
            if (length(chrs_use) == 0) {
              err0()
            } else if (length(chrs_use) < length(seqlevels(annot))) {
              msg1()
              seqlevels(annot) <- chrs_use
            } 
            if (length(chrs_use) < length(seqlevels(genome))) {
              msg2()
            }
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
