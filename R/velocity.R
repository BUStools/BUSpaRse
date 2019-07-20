err0 <- function(n) {
  if (n == 0) {
    stop("Chromosone names in annotation and genome do not overlap. ",
         "Consider checking that formatting of chromosome names of ",
         "the annotation and the genome match.")
  }
}

msg1 <- function(n) {
  message(n, " sequences in the annotation absent from the genome were dropped.")
}

msg2 <- function(n) {
  message(n, " sequences in the genome are absent from the annotation.")
}

#' Match chromosome naming styles of annotation and genome
#' 
#' Internal use. This function matches chromosome naming styles. It will also 
#' give the genome and the annotation the same `genome` slot. This function 
#' assumes that the annotation and the genome refer to the same version of 
#' genome. If more than one style, then the first element will be used.
#' 
#' @inheritParams get_velocity_files
#' @param annot Genome annotation, an object of a class with a 
#' \code{\link{seqlevels}} method, such as `GRanges`, `TxDb`, and `EnsDb`.
#' @return A list of two. The first element is the genome with the proper style,
#' and the second element is the annotation with the proper style.
#' @importFrom GenomeInfoDb seqlevelsStyle genome isCircular 
#' @importFrom GenomeInfoDb seqlevelsStyle<- genome<- isCircular<-
match_style <- function(Genome, annot, style) {
  # Match styles
  gr_style <- tryCatch(seqlevelsStyle(annot),
                       error = function(e) return("other"))
  gn_style <- tryCatch(seqlevelsStyle(Genome),
                       error = function(e) return("other"))
  if (all(c(gr_style, gn_style) != "other")) {
    if (style == "annotation") {
      seqlevelsStyle(Genome) <- gr_style[1]
    } else if (style == "genome") {
      seqlevelsStyle(annot) <- gn_style[1]
    } else {
      seqlevelsStyle(Genome) <- style
      seqlevelsStyle(annot) <- style
    }
  }
  list(genome = Genome, annot = annot)
}

#' Remove chromosomes in anotation absent from genome
#' 
#' @param chrs_use Character vector of names of chromosomes present in both the
#' annotation and the genome.
#' @param annot Either a `GRanges` object or a `TxDb` object for gene 
#' annotation.
#' @return A subsetted genome annotation of the same type ofo the input genome
#' annotation.
#' @importFrom GenomeInfoDb seqlevels seqlevelsInUse
#' @importFrom GenomicFeatures seqlevels<-
sub_annot <- function(chrs_use, annot) {
  err0(length(chrs_use))
  n <- length(seqlevels(annot)) - length(chrs_use)
  if (n > 0) {
    msg1(n)
    if (is(annot, "GRanges")) {
      annot <- annot[as.vector(seqnames(annot)) %in% chrs_use]
      seqlevels(annot) <- seqlevelsInUse(annot)
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
check_genome <- function(chrs_use, Genome) {
  if (is(Genome, "DNAStringSet")) {
    n <- length(Genome)
  } else if (is(Genome, "BSgenome")) {
    n <- length(seqlevels(Genome))
  }
  n2 <- n - length(chrs_use)
  if (n2 > 0) {
    msg2(n2)
  }
}

#' Subset genome annotation
#' 
#' Exclude chromosomes present in the annotation but absent from the genome and
#' add information about circular chromosomes.
#' 
#' @inheritParams get_velocity_files
#' @inheritParams sub_annot
#' @return A subsetted genome annotation of the same type ofo the input genome
#' annotation.
setGeneric("subset_annot", 
           function(Genome, annot)
             standardGeneric("subset_annot"))

#' @rdname subset_annot
setMethod("subset_annot", signature = c("DNAStringSet", "ANY"),
          function(Genome, annot) {
            chrs_use <- intersect(seqlevels(annot), names(Genome))
            annot <- sub_annot(chrs_use, annot)
            check_genome(chrs_use, Genome)
            return(annot)
          })

#' @rdname subset_annot
setMethod("subset_annot", signature = c("BSgenome", "ANY"),
          function(Genome, annot) {
            chrs_use <- intersect(seqlevels(annot), seqlevels(Genome))
            annot <- sub_annot(chrs_use, annot)
            check_genome(chrs_use, Genome)
            return(annot)
          })

#' Extract transcriptomic sequences from genome
#' 
#' This is a thin wrapper around \code{\link{extractTranscriptSeqs}}. This 
#' wrapper takes care of object type of the genome.
#' 
#' @inheritParams subset_annot
#' @param annot A `GRanges` object.
#' @return A `DNAStringSet` object with transcript sequences, whose names are
#' the transcript IDs. 
extract_tx <- function(Genome, annot) {
  gr <- sort(annot[annot$type == "exon"])
  gr <- GenomicRanges::split(gr, gr$transcript_id)
  extractTranscriptSeqs(Genome, gr)
}

#' Transfer information about circular chromosomes between genome and annotation
#' 
#' Internal use, called after calling \code{\link{subset_annot}}. 
#' 
#' @inheritParams match_style
#' @return If neither genome nor annotation indicates which chromosome is 
#' circular, then the input will be returned unchanged. If only one of genome 
#' and annotation has such information, then it will be transferred to the one
#' that does not. If both do have such information, the information from the 
#' genome will be transferred to the annotation if they're different.
annot_circular <- function(Genome, annot) {
  if (all(is.na(isCircular(Genome))) && all(is.na(isCircular(annot)))) {
    message("isCircular information is absent.")
  } else if (all(is.na(isCircular(Genome)))) {
    isCircular(Genome) <- isCircular(annot)[seqlevels(Genome)]
  } else {
    isCircular(annot) <- isCircular(Genome)[seqlevels(annot)]
  }
  return(list(genome = Genome, annot = annot))
}

#' Check if transcript ID in transcriptome and annotation match
#' 
#' This function throws an error if transcript IDs in transcriptome and
#' annotation do not overlap. If they do overlap, this function will give a 
#' message about transcript IDs that do not agree in the transcriptome and the
#' annotation
#' 
#' @param tx_annot Character vector of transcript IDs from the annotation.
#' @param tx Character vector of transcript IDs from the transcriptome.
#' @return Character vector of the overlapping transcript IDs.
check_tx <- function(tx_annot, tx) {
  tx_overlap <- intersect(tx_annot, tx)
  if (length(tx_overlap) == 0) {
    stop("Transcripts in gene annotation do not overlap with those in the ",
         "transcriptome. Consider checking the type of transcript ID used ",
         "or whether version number is included.")
  }
  len_diff <- length(unique(tx_annot)) - length(tx_overlap)
  if (len_diff > 0) {
    message("There are ", len_diff, " transcripts in the gene annotation is ",
            "absent from the transcriptome. These transcripts are removed.")
  }
  return(tx_overlap)
}

#' Standardize GRanges field names
#' 
#' To avoid introducing rlang as another dependency for tidyeval. This function
#' will also convert exon numbers to integer.
#' 
#' @param gr A `GRanges` object.
#' @param gene_id Name of the metadata field for gene ID.
#' @param transcript_id Name of the metadata field for transcript ID.
#' @return A `GRanges` object with standardized names: gene ID as `gene_id`,
#' and transcript ID as `transcript_id`.
#' @importFrom GenomicRanges mcols<-
standardize_tags <- function(gr, gene_id, transcript_id) {
  fields <- names(mcols(gr))
  if (gene_id != "gene_id") {
    names(mcols(gr))[fields == gene_id] <- "gene_id"
  }
  if (transcript_id != "transcript_id") {
    names(mcols(gr))[fields == transcript_id] <- "transcript_id"
  }
  return(gr)
}

#' Get flanked intronic ranges
#' 
#' Adapted from Julien Roux.
#' 
#' @param gr Exon `GRangesList`, each element of which is for a transcript or a
#' gene.
#' @param L Read length.
#' @return A `GRangesList` object with ranges for flanked intronic regions.
#' Each element of the list is for a transcript or a gene.
#' @importFrom GenomicRanges resize setdiff
#' @importFrom GenomicFeatures intronsByTranscript
#' @importFrom S4Vectors elementNROWS
get_flanked_introns <- function(gr, L) {
  exonExt <- 2 * (L - 1)
  newWidths <- width(gr) - exonExt
  # Exclude flanking regions from exons
  newWidths <- lapply(newWidths, function(x) {x[x < 0] <- 0; x}) 
  gr <- resize(gr, width = newWidths, fix = "center")
  # Get introns
  introns <- setdiff(range(gr), gr)
  introns <- introns[!elementNROWS(introns) == '0']
  return(introns)
}

#' Write the files for RNA velocity to disk
#' 
#' Write the files for RNA velocity to disk, in the specified output directory.
#' 
#' @inheritParams get_flanked_introns
#' @param out_path Directory to save the outputs written to disk. If this
#' directory does not exist, then it will be created.
#' @param introns Intronic ranges plus flanking region, returned by 
#' \code{\link{get_flanked_introns}}.
#' @param tr2g_cdna A data frame with columns `transcript` and `gene` that maps
#' transcripts to genes for spliced transcripts.
#' @return Nothing into the R session. The files are written to disk.
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings writeXStringSet
write_velocity_output <- function(out_path, introns, Genome, Transcriptome,
                                  isoform_action, tr2g_cdna, compress_fa,
                                  width) {
  message("Extracting flanked intronic sequences")
  intron_seqs <- unlist(getSeq(Genome, introns))
  ns <- names(intron_seqs)
  # Not to confuse with version number
  names(intron_seqs) <- make.unique(paste0(ns, "-I"), sep = "")
  if (isoform_action == "collapse") {
    tr2g_intron <- data.frame(transcript = names(intron_seqs),
                              gene = ns,
                              stringsAsFactors = FALSE)
  } else {
    # ns will be transcript ID in this case
    tr2g_intron <- data.frame(transcript = names(intron_seqs),
                              gene = tr2g_cdna$gene[match(ns, tr2g_cdna$transcript)],
                              stringsAsFactors = FALSE)
  }
  out_fa <- paste(out_path, "cDNA_introns.fa", sep = "/")
  if (compress_fa) {
    out_fa <- paste0(out_fa, ".gz")
  }
  message("Writing outputs")
  if (is.character(Transcriptome)) {
    file.copy(Transcriptome, out_fa)
    writeXStringSet(intron_seqs, out_fa, append = TRUE, compress = compress_fa,
                    width = width)
  } else {
    writeXStringSet(c(Transcriptome, intron_seqs), out_fa, compress = compress_fa,
                    width = width)
  }
  writeLines(names(intron_seqs), paste(out_path, "introns_tx_to_capture.txt",
                                   sep = "/"))
  writeLines(unique(tr2g_cdna$transcript),
             paste(out_path, "cDNA_tx_to_capture.txt", sep = "/"))
  fwrite(rbind(tr2g_cdna, tr2g_intron), paste(out_path, "tr2g.tsv", sep = "/"),
         quote = FALSE, sep = "\t", col.names = FALSE)
}

#' Validate input to get_velocity_files
#' 
#' @inheritParams get_velocity_files
#' @return Will throw error if validation fails. Returns a named list whose
#' first element is the normalized path to output directory, and whose second
#' element is the normalized path to the transcriptome file if specified.
#' @importFrom methods is
validate_velocity_input <- function(L, Genome, Transcriptome, out_path,
                                    compress_fa, width) {
  if (length(L) > 1 || L %% 1 > sqrt(.Machine$double.eps) || !is.atomic(L)) {
    stop("L must be an integer vector with length 1.")
  }
  if (length(width) > 1 || width %% 1 > sqrt(.Machine$double.eps) || 
      !is.atomic(width)) {
    stop("width must be an integer vector with length 1.")
  }
  if (!is.logical(compress_fa) || length(compress_fa) > 1) {
    stop("compress_fa must be logical with length 1.")
  }
  out_path <- normalizePath(out_path, mustWork = FALSE)
  if (!dir.exists(out_path)) {
    dir.create(out_path)
  }
  if (!is(Genome, "BSgenome") && !is(Genome, "DNAStringSet")) {
    stop("Genome must be either a BSgenome object or a DNAStringSet.")
  }
  if (is.character(Transcriptome)) {
    Transcriptome <- normalizePath(Transcriptome, mustWork = TRUE)
  } else if (!is(Transcriptome, "DNAStringSet") && !is.null(Transcriptome)) {
    stop("Transcriptome must be either a valid file path or a DNAStringSet.")
  } else {
    Transcriptome <- NULL
  }
  return(list(out_path = out_path, tx_path = Transcriptome))
}
