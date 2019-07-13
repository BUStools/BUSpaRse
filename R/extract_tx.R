#' Extract transcriptomic sequences from genome
#' 
#' This is a thin wrapper around \code{\link{extractTranscriptSeqs}}. This 
#' wrapper takes care of object type of the genome.
#' 
#' @inheritParams get_velocity_files
#' @param annot Either a `GRanges` object or a `TxDb` object for gene 
#' annotation.
#' @param \dots Extra argument for methods.
#' @return A `DNAStringSet` object with transcript sequences, whose names are
#' the transcript IDs. 
setGeneric("extract_tx", 
           function(genome, annot, ...) 
             standardGeneric("extract_tx"))

#' @rdname extract_tx
#' @importFrom GenomicRanges sort
#' @importFrom GenomicFeatures extractTranscriptSeqs
setMethod("extract_tx", signature = c("DNAStringSet", "GRanges"),
          function(genome, annot, transcript_id = "transcript_id") {
            gr <- sort(annot[annot$type == "exon"])
            gr <- split(gr, mcols(gr)[[transcript_id]])
            extractTranscriptSeqs(genome, gr)
          })

#' @rdname extract_tx
#' @importFrom GenomicFeatures exonsBy
setMethod("extract_tx", signature = c("DNAStringSet", "TxDb"),
          function(genome, annot) {
            gr <- exonsBy(annot, by = "tx", use.names = TRUE)
            extractTranscriptSeqs(genome, gr)
          })

#' @rdname extract_tx
setMethod("extract_tx", signature = c("BSgenome", "GRanges"),
          function(genome, annot, transcript_id = "transcript_id") {
            gr <- sort(annot[annot$type == "exon"])
            gr <- split(gr, mcols(gr)[[transcript_id]])
            extractTranscriptSeqs(genome, gr)
          })

#' @rdname extract_tx
setMethod("extract_tx", signature = c("BSgenome", "TxDb"),
          function(genome, annot) {
            gr <- exonsBy(annot, by = "tx", use.names = TRUE)
            extractTranscriptSeqs(genome, gr)
          })
