#' Get files required for RNA velocity with bustools
#' 
#' Computation of RNA velocity requires the number of unspliced transcripts,
#' which can be quantified with reads containing intronic sequences. This 
#' function extracts intronic sequences flanked by L-1 bases of exonic sequences
#' where L is the biological read length of the single cell technology of 
#' interest. The flanking exonic sequences are included for reads partially
#' mapping to an intron and an exon.
#' 
#' @inheritParams tr2g_gtf
#' @param file Path to a GTF file with annotation of exon coordinates of 
#' transcripts, preferably from Ensembl. The file must be formatted as in Ensembl.
#' In the metadata, the following fields are required: type (e.g. whether the 
#' range of interest is a gene or transcript or exon or CDS), gene ID, 
#' transcript ID, and exon number. These fields need not to have standard names,
#' as long as their names are specified in arguments of this function.
#' @param L Length of the biological read. For instance, 98 nt for 10x v2 and 91
#' nt for 10x v3.
#' @param genome Either a \code{\link{BSgenome}} or a \code{\link{XStringSet}}
#' object of genomic sequences, where the intronic sequences will be extracted
#' from. Naming convention of chromosomes in the genome must match that in the
#' GTF file, e.g. consistently use 1 for chromosome 1 or consistently use chr1.
#' @param transcriptome Either a \code{\link{XStringSet}} or a path to a fasta
#' file (can be gzipped) of the transcriptome, which contains sequences of
#' spliced transcripts. This will be concatenated with the intronic sequences
#' to give one fasta file.
#' @param out_path Directory to save the outputs written to disk. If this
#' directory does not exist, then it will be created.
#' @param short_exon_action Character, indicating action to take with exons
#' that are shorter than L-1. Must be one of the following:
#' \describe{
#' \item{truncate}{The flanking region involving this exon will be truncated;
#' the entire exon will be included in the flanking region, but not the intron
#' on the other side of this exon.}
#' \item{include}{The short exon will be included as part of the intronic 
#' sequence between the surrounding exons. Note that if the short exon is the
#' first or the last exon of the transcript, then the flanking region in this
#' exon will always be truncated to not include intergenic regions.}
#' }
#' @param isoform_action Character, indicating action to take with different
#' transcripts of the same gene. Must be one of the following:
#' \describe{
#' \item{separate}{Introns from different transcripts will be kept separate.}
#' \item{collapse}{First, the union of all exons of different transcripts of a
#' gene will be taken. Then the introns will be inferred from this union.}
#' }
#' @param exon_number Character vector of length 1. Tag in \code{attribute} 
#' field corresponding to exon numbers, namely the order of exons from 5' to 3'
#' on the transcript. Must be coercable to integer.
#' @param exon_id Character vector of length 1. Tag in \code{attribute} 
#' field corresponding to exon ID. If the specified tag does not exist in the
#' GTF file, then a new metadata column that can uniquely identify exons will
#' be created.
#' @return The following files will be written to disk in the directory 
#' `out_path`:
#' \describe{
#' \item{cDNA_introns.fa}{A fasta file containing both the spliced transcripts
#' and the flanked intronic sequences. This will be used to build the `kallisto`
#' index.}
#' \item{cDNA_tx_to_capture.txt}{A text file of transcript IDs of spliced
#' transcripts.}
#' \item{introns_tx_to_capture.txt}{A text file of IDs of introns. The names 
#' will have the pattern <transcript ID>-Ix, where x is a number differentiating
#' between introns of the same transcript. If all transcripts of the same gene
#' are collapsed before inferring intronic sequences, gene ID will be used in
#' place of transcript ID.}
#' \item{tr2g.txt}{A text file with two columns matching transcripts and introns
#' to genes. The first column is transcript or intron ID, and the second column
#' is the corresponding gene ID. The part for transcripts are generated from
#' the GTF file.}
#' }
#' Nothing is returned into the R session.
#' @export
#' @importFrom BSgenome getSeq
#' @importFrom GenomicFeatures intronsByTranscript makeTxDbFromGRanges
#' @importFrom dplyr group_by arrange pull min_rank desc lead case_when ungroup
#' @importFrom BiocGenerics width start end start<- end<- 
#' @importFrom GenomicRanges strand reduce split seqnames mcols mcols<-
#' @importFrom methods is
#' @importFrom Biostrings writeXStringSet
#' @examples 
#' # Use toy example
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file <- paste0(toy_path, "/velocity_annot.gtf")
#' genome <- Biostrings::readDNAStringSet(paste0(toy_path, "/velocity_genome.fa"))
#' transcriptome <- paste0(toy_path, "/velocity_tx.fa")
#' get_velocity_files(file, 11, genome, transcriptome, ".",
#'                    gene_version = NULL, transcript_version = NULL)
get_velocity_files <- function(file, L, genome, transcriptome, out_path,
                               short_exon_action = c("truncate", "include"),
                               isoform_action = c("separate", "collapse"),
                               transcript_id = "transcript_id",
                               gene_id = "gene_id", 
                               transcript_version = "transcript_version",
                               gene_version = "gene_version", 
                               version_sep = ".", exon_number = "exon_number",
                               exon_id = "exon_id") {
  # Validate arguments------------------------
  file <- normalizePath(file, mustWork = TRUE)
  if (length(L) > 1 || L %% 1 > sqrt(.Machine$double.eps) || !is.atomic(L)) {
    stop("L must be an integer vector with length 1.")
  }
  L <- as.integer(L)
  if (!is(genome, "BSgenome") && !is(genome, "DNAStringSet")) {
    stop("genome must be either a BSgenome object or a DNAStringSet.")
  }
  if (is.character(transcriptome)) {
    transcriptome <- normalizePath(transcriptome, mustWork = TRUE)
  } else if (!is(transcriptome, "DNAStringSet")) {
    stop("transcriptome must be either a valid file path or a DNAStringSet.")
  }
  out_path <- normalizePath(out_path, mustWork = FALSE)
  if (!dir.exists(out_path)) {
    dir.create(out_path)
  }
  short_exon_action <- match.arg(short_exon_action)
  isoform_action <- match.arg(isoform_action)
  
  # Process GTF file to get intronic ranges-----------------
  gr <- plyranges::read_gff(file)
  check_gff("gtf", file, transcript_id, gene_id)
  fields <- names(mcols(gr))
  check_tag_present(exon_number, fields, error = TRUE)
  if (isoform_action == "collapse") {
    check_tag_present(exon_id, fields, error = TRUE)
  }
  # I need this anyway and this will validate the GTF file.
  tr2g_cdna <- tr2g_GRanges(gr, gene_name = NULL, 
                            transcript_version = transcript_version,
                            gene_version = gene_version,
                            version_sep = version_sep)
  gr <- gr[gr$type == "exon" & (as.vector(strand(gr)) %in% c("+", "-"))]
  # Switch to standard naming
  names(mcols(gr))[fields == gene_id] <- "gene_id"
  names(mcols(gr))[fields == transcript_id] <- "transcript_id"
  names(mcols(gr))[fields == exon_number] <- "exon_number"
  gr$exon_number <- as.integer(gr$exon_number)
  if (!is.null(gene_version)) {
    names(mcols(gr))[fields == gene_version] <- "gene_version"
    gr$gene_id <- paste(gr$gene_id, gr$gene_version, sep = version_sep)
  }
  if (!is.null(transcript_version)) {
    names(mcols(gr))[fields == transcript_version] <- "transcript_version"
    gr$transcript_id <- paste(gr$transcript_id, gr$transcript_version,
                              sep = version_sep)
  }
  if (isoform_action == "collapse") {
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
  }
  if (short_exon_action == "include") {
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
    if (!exon_id %in% names(metas)) {
      metas$exon_id <- paste(metas$transcript_id, metas$exon_number, sep = "-")
    }
    exons_use <- unique(metas$exon_id[metas$include])
    gr <- gr[gr$exon_id %in% exons_use]
  }
  next_same <- too_short <- chr <- min_start <- exon_number_sort <- 
    intron_id <- flank5 <- flank3 <- NULL
  # Calculate the appropriate flanking lengths
  metas <- mcols(gr) %>% 
    as.data.frame() %>% 
    mutate(width = width(gr),
           start = start(gr),
           strand = as.vector(strand(gr)),
           chr = as.vector(seqnames(gr)),
           too_short = width < L - 1,
           exon_number = as.integer(exon_number),
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
    
  # Get intronic ranges
  txdb <- makeTxDbFromGRanges(gr)
  introns <- intronsByTranscript(txdb, use.names = TRUE)
  introns <- unlist(introns)
  start(introns) <- ifelse(strand(introns) == "+", start(introns) - metas$flank5,
                           start(introns) - metas$flank3)
  end(introns) <- ifelse(strand(introns) == "+", end(introns) + metas$flank3,
                         end(introns) + metas$flank5)
  names(introns) <- metas$intron_id
  
  # Prepare output---------------------------
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

# To do: Write method for TxDb as that's a more convenient way to get gene
# annotation than GTF files.