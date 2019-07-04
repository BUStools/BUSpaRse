#' Get files required for RNA velocity with bustools
#' 
#' Computation of RNA velocity requires the number of unspliced transcripts,
#' which can be quantified with reads containing intronic sequences. This 
#' function extracts intronic sequences flanked by L-1 bases of exonic sequences
#' where L is the biological read length of the single cell technology of 
#' interest. The flanking exonic sequences are included for reads partially
#' mapping to an intron and an exon.
#' 
#' @inheritParams tr2g_fasta
#' @param file Path to a GTF file with annotation of exon coordinates of 
#' transcripts, preferably from Ensembl. The file must be formatted as in Ensembl.
#' In the metadata, gene ID must be named `gene_id`, transcript ID must be
#' `transcript_id`, gene version must be `gene_version`, and transcript version must
#' be `transcript_version`. There also must be a metadata field called "type",
#' which indicates whether the range is for gene or transcript or exon.
#' @param L Length of the biological read. 
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
#' is the corresponding gene ID.}
#' }
#' Nothing is returned into the R session.
#' @export
#' @importFrom BSgenome getSeq
#' @importFrom GenomicFeatures intronsByTranscript makeTxDbFromGRanges
#' @importFrom dplyr group_by arrange pull min_rank desc lead lag case_when
#' @importFrom BiocGenerics width start end
#' @importFrom GenomicRanges strand
#' 
get_velocity_files <- function(file, L, genome, transcriptome, out_path,
                               short_exon_action = c("truncate", "include"),
                               isoform_action = c("separate", "collapse"),
                               use_transcript_version = TRUE, 
                               use_gene_version = TRUE) {
  # Validate arguments------------------------
  file <- normalizePath(file, mustWork = TRUE)
  if (!is(genonme, "BSgenome") && !is(gennome, "DNAStringSet")) {
    stop("genonme must be either a BSgenome object or a DNAStringSet.")
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
  gr <- read_gff(file)
  gr <- gr[gr$type == "exon" & strand(gr) %in% c("+", "-")]
  if (use_transcript_version) {
    gr$transcript_id <- paste(gr$transcript_id, gr$transcript_version,
                              sep = ".")
  }
  if (use_gene_version) {
    gr$gene_id <- paste(gr$gene_id, gr$gene_version, sep = ".")
  }
  if (isoform_action == "collapse") {
    gr <- split(gr, gr$gene_id)
    gr <- reduce(gr)
    gr <- unlist(gr)
    gr$type <- "exon"
    gr$gene_id <- gr$transcript_id <- names(gr)
    names(gr) <- NULL
    metas <- as.data.frame(mcols(gr)) %>% 
      mutate(start_arr = as.integer(paste0(strand(gr), start(gr)))) %>% 
      group_by(transcript_id) %>% 
      mutate(exon_number = min_rank(start_arr))
    gr$exon_number <- metas$exon_number
  }
  if (short_exon_action == "include") {
    # Remove exons that are too short except the first and the last exon of
    # the transcript
    metas <- mcols(gr) %>% 
      as.data.frame() %>% 
      group_by(transcript_id) %>% 
      mutate(exon_number = as.integer(exon_number),
             exon_number_rev = min_rank(dplyr::desc(exon_number)),
             width = width(gr),
             include = (width > L - 1) | 
               (exon_number == 1 | exon_number_rev == 1))
    exons_use <- metas %>% 
      filter(include) %>% 
      pull(exon_id) %>% unique()
    gr <- gr[gr$exon_id %in% exons_use]
  }
  metas <- mcols(gr) %>% 
    as.data.frame() %>% 
    mutate(width = width(gr),
           start = start(gr),
           strand = strand(gr),
           too_short = width < L - 1,
           exon_number = as.integer(exon_number)) %>% 
    arrange(transcript_id, exon_number) %>% 
    mutate(next_same = lead(exon_number) != 1, # FALSE for the last exon of a transcript
           next_same = ifelse(is.na(next_same), FALSE, next_same),
           # Whether to truncate by the 5' end of the intron + flank
           trunc5 = too_short & next_same,
           # Whether to truncate by the 3' end of the intron + flank
           trunc3 = lead(too_short) & next_same,
           # Length of flanking region by 5' end of the intron
           flank5 = case_when(
             trunc5 ~ L - 1 - width,
             !next_same ~ NA_real_,
             TRUE ~ L - 1
           ),
           # Lenngth of flannking region by 3' end of the intron
           flank_end = case_when(
             trunc3 ~ L - 1 - lead(width),
             !next_same ~ NA_real_,
             TRUE ~ L - 1
           )) %>% 
    group_by(transcript_id) %>% 
    mutate(min_start = min(start)) %>% 
    # Sort to be the same order as the intronic ranges.
    arrange(desc(strand), min_start, transcript_id, desc(exon_number))
  
  # Get intronic rannges
  txdb <- makeTxDbFromGRanges(gr)
  introns <- intronsByTranscript(txdb, use.names = TRUE)
  fl5 <- metas$flank5[!is.na(metas$flank5)]
  fl3 <- metas$flank3[!is.na(metas$flank3)]
  start(introns) <- ifelse(strand(gr) == "+", start(introns) - fl5,
                           start(introns) - fl3)
  end(introns) <- ifelse(strand(gr) == "+", end(introns) + fl3,
                         end(introns) + fl5)
  metas <- data.frame(gene_id = names(introns)) %>% 
    group_by(gene_id) %>% 
    mutate(intron_number = row_number(gene_id),
           intron_id = paste0(gene_id, "-I", intron_number))
  names(introns) <- metas$intron_id
  
  # Prepare output---------------------------
  # Get transcriptome
  if (is.character(transcriptome)) {
    transcriptome <- readDNAStringSet(transcriptome)
  }
  intron_seqs <- getSeq(genome, introns)
  writeXStringSet(c(transcriptome, intron_seqs), 
                  paste(out_path, "cDNA_introns.fa"))
  writeLines(names(introns), paste(out_path, "introns_tx_to_capture.txt",
                                   sep = "/"))
  if (use_transcript_version) {
    tx_regex <- "^[a-zA-Z\\d-\\.]+"
  } else {
    tx_regex <- "^[a-zA-Z\\d-]+"
  }
  writeLines(str_extract(names(transcriptome), tx_regex),
             paste(out_path, "cDNA_tx_to_capture.txt", sep = "/"))
  tr2g_cdna <- tr2g_DNAStringSet(transcriptome)
  tr2g_intron <- setNames(metas[,c("intron_id", "gene_id")], 
                          c("transcript", "gene"))
  fwrite(rbind(tr2g_cdna, tr2g_intron), paste(out_path, "tr2g.txt", sep = "/"),
         quote = FALSE, sep = "\t", col.names = FALSE)
}
