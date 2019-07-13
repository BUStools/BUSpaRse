#' Standardize GRanges field names
#' 
#' To avoid introducing rlang as another dependency for tidyeval. This function
#' will also convert exon numbers to integer.
#' 
#' @param gr A `GRanges` object.
#' @param gene_id Name of the metadata field for gene ID.
#' @param transcript_id Name of the metadata field for transcript ID.
#' @param exon_number Name of the metadata field for exon number.
#' @return A `GRanges` object with standardized names: gene ID as `gene_id`,
#' transcript ID as `transcript_id`, and exon number as `exon_number`.
#' @importFrom GenomicRanges mcols<-
standardize_tags <- function(gr, gene_id, transcript_id, exon_number) {
  fields <- names(mcols(gr))
  if (gene_id != "gene_id") {
    names(mcols(gr))[fields == gene_id] <- "gene_id"
  }
  if (transcript_id != "transcript_id") {
    names(mcols(gr))[fields == transcript_id] <- "transcript_id"
  }
  if (exon_number != "exon_number") {
    names(mcols(gr))[fields == exon_number] <- "exon_number"
  }
  if (!is.integer(gr$exon_number)) {
    gr$exon_number <- as.integer(gr$exon_number)
  }
  return(gr)
}

#' Collapse different isoforms of a gene
#' 
#' This function takes the union of all exon ranges of different isoform of each
#' gene.
#' 
#' @inheritParams standardize_tags
#' @param gr A `GRanges` object with only exon ranges. As this function is for 
#' internal use only and the fields are checked elsewhere, fields will not be 
#' checked in this function.
#' @return A `GRanges` object with only exon ranges and different isoforms
#' collapsed. Now the "transcript ID" in the metadata columns will be the same
#' as the gene ID.
#' @importFrom GenomicRanges split reduce strand mcols start
#' @importFrom dplyr mutate group_by min_rank
collapse_isoforms <- function(gr, gene_id = "gene_id", 
                              transcript_id = "transcript_id",
                              exon_number = "exon_number") {
  gr <- standardize_tags(gr, gene_id, transcript_id, exon_number)
  gr <- GenomicRanges::split(gr, gr$gene_id)
  gr <- GenomicRanges::reduce(gr)
  gr <- unlist(gr)
  gr$type <- "exon"
  gr$gene_id <- gr$transcript_id <- names(gr)
  names(gr) <- NULL
  # Address R CMD check note
  start_arr <- NULL
  metas <- as.data.frame(mcols(gr)) %>% 
    mutate(start_arr = as.integer(paste0(as.vector(strand(gr)), 
                                         start(gr)))) %>% 
    group_by(transcript_id) %>% 
    mutate(exon_number = min_rank(start_arr))
  gr$exon_number <- metas$exon_number
  return(gr)
}

#' Exclude exons that are too short
#' 
#' Exons that are too short are excluded from the `GRanges` object, except the
#' first and the last exons of a transcript.
#' 
#' @inheritParams collapse_isoforms
#' @inheritParams get_velocity_files
#' @return A `GRanges` object without exons that are too short.
#' @importFrom GenomicRanges width
#' @importFrom dplyr desc
exclude_short_exons <- function(gr, L, gene_id = "gene_id", 
                                transcript_id = "transcript_id", 
                                exon_number = "exon_number",
                                exon_id = "exon_id") {
  gr <- standardize_tags(gr, gene_id, transcript_id, exon_number)
  if (!exon_id %in% names(mcols(gr))) {
    gr$exon_id <- paste(gr$transcript_id, gr$exon_number, sep = "-")
  }
  exon_number_rev <- NULL
  # Remove exons that are too short except the first and the last exon of
  # the transcript
  metas <- mcols(gr) %>% 
    as.data.frame() %>% 
    mutate(width = width(gr)) %>% 
    group_by(transcript_id) %>% 
    mutate(exon_number_rev = min_rank(dplyr::desc(exon_number)),
           include = (width > L - 1) | 
             (exon_number == 1 | exon_number_rev == 1))
  exons_use <- unique(metas[[exon_id]][metas$include])
  gr <- gr[mcols(gr)[[exon_id]] %in% exons_use]
  return(gr)
}

#' Get lengths of flanking regions
#' 
#' Get lengths of flanking regions at both 5' and 3' ends of each intron, and 
#' sort the flanking lengths in the order the introns will be sorted.
#' 
#' @inheritParams get_velocity_files
#' @param df A data frame with the following columns:
#' \describe{
#' \item{chr}{Which chromosome the gene is from.}
#' \item{start}{Start of the exonic range. For genes on the minus strand, this
#' is the 5' end on the plus strand, and thus actually the end of the exonic
#' range.}
#' \item{end}{End of the exonic range.}
#' \item{strand}{Strand of the transcript of interest.}
#' \item{gene_id}{Gene ID.}
#' \item{transcript_id}{Transcript ID.}
#' \item{exon_number}{Order of the exon from 5' to 3' on the transcript. Should
#' be integer.}
#' }
#' The column names must be as described above.
#' @return A data frame with the following columns:
#' \describe{
#' \item{gene_id}{Gene ID.}
#' \item{intron_id}{Intron ID, constructed from transcript ID and the order
#' of introns from 5' annd 3' of the transcript.}
#' \item{flank5}{Length of the flanking region on the 5' end of the intron.}
#' \item{flank3}{Length of the flanking region on the 3' end of the intron.}
#' }
#' @importFrom dplyr lead case_when arrange ungroup
get_flank_lengths <- function(df, L) {
  next_same <- too_short <- chr <- min_start <- exon_number_sort <- 
    intron_id <- flank5 <- flank3 <- exon_number <- 
    transcript_id <- gene_id <- NULL
  if (!is.integer(df$exon_number)) {
    df$exon_number <- as.integer(df$exon_number)
  }
  # Calculate the appropriate flanking lengths
  metas <- df %>% 
    mutate(too_short = width < L - 1,
           exon_number_sort = as.integer(paste0(strand, exon_number))) %>% 
    arrange(transcript_id, exon_number) %>% 
    mutate(next_same = lead(exon_number) != 1, # FALSE for the last exon of a transcript
           next_same = ifelse(is.na(next_same), FALSE, next_same),
           # Whether to truncate the 5' flanking region
           trunc5 = too_short & next_same,
           # Whether to truncate the 3' flanking region
           trunc3 = lead(too_short) & next_same,
           # Length of flanking region by 5' end of the intron
           flank5 = case_when(
             trunc5 ~ width,
             TRUE ~ L - 1L
           ),
           # Length of flannking region by 3' end of the intron
           flank3 = case_when(
             trunc3 ~ lead(width),
             TRUE ~ L - 1L
           ))
  metas <- metas[metas$next_same,]
  # Sort to be the same order as the intronic ranges.
  metas <- metas %>% 
    group_by(transcript_id) %>% 
    mutate(min_start = min(start)) %>% 
    ungroup() %>% 
    arrange(chr, desc(strand), min_start, transcript_id, exon_number_sort) %>% 
    mutate(intron_id = paste0(transcript_id, "-I", exon_number)) %>% 
    dplyr::select(gene_id, intron_id, flank5, flank3)
  return(metas)
}

#' Get flanked intronic ranges
#' 
#' @param txdb A `TxDb` object with exonic annotations.
#' @param metas The data frame returned by \code{\link{get_flank_lengths}}.
#' @return A `GRanges` object with ranges for flanked intronic regions.
#' @importFrom GenomicRanges start<- end<- strand end
#' @importFrom GenomicFeatures intronsByTranscript
get_flanked_introns <- function(txdb, metas) {
  introns <- intronsByTranscript(txdb)
  introns <- unlist(introns)
  start(introns) <- ifelse(strand(introns) == "+", start(introns) - metas$flank5,
                           start(introns) - metas$flank3)
  end(introns) <- ifelse(strand(introns) == "+", end(introns) + metas$flank3,
                         end(introns) + metas$flank5)
  names(introns) <- metas$intron_id
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
#' @param genome Either a \code{\link{BSgenome}} or a \code{\link{XStringSet}}
#' object of genomic sequences, where the intronic sequences will be extracted
#' from. Naming convention of chromosomes in the genome must match that in the
#' GTF file, e.g. consistently use 1 for chromosome 1 or consistently use chr1.
#' @param transcriptome Either a \code{\link{XStringSet}} or a path to a fasta
#' file (can be gzipped) of the transcriptome, which contains sequences of
#' spliced transcripts. This will be concatenated with the intronic sequences
#' to give one fasta file.
#' @param tr2g_cdna A data frame with columns `transcript` and `gene` that maps
#' transcripts to genes for spliced transcripts.
#' @return Nothing into the R session. The files are written to disk.
#' @importFrom BSgenome getSeq
#' @importFrom Biostrings writeXStringSet
write_velocity_output <- function(out_path, introns, genome, transcriptome,
                                  tr2g_cdna, metas) {
  intron_seqs <- getSeq(genome, introns)
  out_fa <- paste(out_path, "cDNA_introns.fa", sep = "/")
  if (is.character(transcriptome)) {
    file.copy(transcriptome, out_fa)
    writeXStringSet(intron_seqs, out_fa, append = TRUE)
  } else {
    writeXStringSet(c(transcriptome, intron_seqs), out_fa)
  }
  writeLines(names(introns), paste(out_path, "introns_tx_to_capture.txt",
                                   sep = "/"))
  writeLines(unique(tr2g_cdna$transcript),
             paste(out_path, "cDNA_tx_to_capture.txt", sep = "/"))
  tr2g_intron <- setNames(metas[,c("intron_id", "gene_id")], 
                          c("transcript", "gene"))
  fwrite(rbind(tr2g_cdna, tr2g_intron), paste(out_path, "tr2g.tsv", sep = "/"),
         quote = FALSE, sep = "\t", col.names = FALSE)
}

#' Validate input to get_velocity_files
#' 
#' @inheritParams get_velocity_files
#' @return Nothing. Will throw error if validation fails.
#' @importFrom methods is
validate_velocity_input <- function(L, genome, transcriptome, out_path) {
  if (length(L) > 1 || L %% 1 > sqrt(.Machine$double.eps) || !is.atomic(L)) {
    stop("L must be an integer vector with length 1.")
  }
  out_path <- normalizePath(out_path, mustWork = FALSE)
  if (!dir.exists(out_path)) {
    dir.create(out_path)
  }
  if (!is(genome, "BSgenome") && !is(genome, "DNAStringSet")) {
    stop("genome must be either a BSgenome object or a DNAStringSet.")
  }
  if (is.character(transcriptome)) {
    transcriptome <- normalizePath(transcriptome, mustWork = TRUE)
  } else if (!is(transcriptome, "DNAStringSet")) {
    stop("transcriptome must be either a valid file path or a DNAStringSet.")
  }
}
