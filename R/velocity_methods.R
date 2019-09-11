#' Generate RNA velocity files for GRanges
#'
#' @inheritParams tr2g_gtf
#' @param gr A `GRanges` object for gene annotation.
#' @param L Length of the biological read. For instance, 10xv1: 98 nt,
#' 10xv2: 98 nt, 10xv3: 91 nt, Drop-seq: 50 nt. If in doubt check read length
#' in a fastq file for biological reads with the `bash` commands:
#' If the fastq file is gzipped, then do `zcat your_file.fastq.gz | head` on
#' Linux. If on Mac, then `zcat < your_file.fastq.gz | head`. Then you will see
#' lines with nucleotide bases. Copy one of those lines and determine its length
#' with \code{\link{str_length}} in R or `echo -n <the sequence> | wc -c` in
#' `bash`. Which file corresponds to biological reads depends on the particular
#' technology.
#' @param Genome Either a \code{\link{BSgenome}} or a \code{\link{XStringSet}}
#' object of genomic sequences, where the intronic sequences will be extracted
#' from. Use \code{\link{genomeStyles}} to check which styles are supported for
#' your organism of interest; supported styles can be interconverted. If the
#' style in your genome or annotation is not supported, then the style of
#' chromosome names in the genome and annotation should be manually set to be
#' consistent.
#' @param Transcriptome A \code{\link{XStringSet}}, a path to a fasta
#' file (can be gzipped) of the transcriptome which contains sequences of
#' spliced transcripts, or `NULL`. The transcriptome here will be concatenated
#' with the intronic sequences to give one fasta file. When `NULL`, the
#' transriptome sequences will be extracted from the genome
#' given the gene annotation, so it will be guaranteed that transcript IDs in
#' the transcriptome and in the annotation match. Otherwise, the type of
#' transcript ID in the transcriptome must match that in the gene annotation
#' supplied via argument `X`.
#' @param out_path Directory to save the outputs written to disk. If this
#' directory does not exist, then it will be created. Defaults to the current
#' working directory.
#' @param style Formatting of chromosome names. Use
#' \code{\link{genomeStyles}} to check which styles are supported for your
#' organism of interest and what those styles look like. This can also be a
#' style supported for your organism different from the style used by the
#' annotation and the genome. Then this style will be used for both the
#' annotation and the genome. Can take the following values:
#' \describe{
#' \item{annotation}{If style of the annnotation is different from that of the
#' genome, then the style of the annotation will be used.}
#' \item{genome}{If style of the annnotation is different from that of the
#' genome, then the style of the genome will be used.}
#' \item{other}{Custom style, need to manually ensure that the style in
#' annotation matches that of the genome.}
#' \item{Ensembl}{Or `UCSC` or `NCBI`, whichever is supported by your species
#' of interest.}
#' }
#' @param isoform_action Character, indicating action to take with different
#' transcripts of the same gene. Must be one of the following:
#' \describe{
#' \item{collapse}{First, the union of all exons of different transcripts of a
#' gene will be taken. Then the introns will be inferred from this union. Only
#' the flanked intronic sequences are affected; isoforms will always be taken
#' into account for spliced sequences or exon-exon junctions.}
#' \item{separate}{Introns from different transcripts will be kept separate.}
#' }
#' @param exon_option Character, indicating how exonic sequences should be
#' included in the kallisto index. Must be one of the following:
#' \describe{
#' \item{full}{The full cDNA sequences, which include the full exonic sequences,
#' will be used. This is the default.}
#' \item{junction}{Only the exon-exon junctions, with L-1 bases on each side
#' of the junctions, will be used.}
#' }
#' @param compress_fa Logical, whether to compress the output fasta file of
#' transcriptome and flanked intronic sequenncess. If `TRUE`, then the fasta
#' file will be gzipped.
#' @param width Maximum number of letters per line of sequence in the output
#' fasta file. Must be an integer.
#' @return See \code{\link{get_velocity_files}}
#' @importFrom GenomicRanges seqnames strand
#' @importFrom GenomicFeatures extractTranscriptSeqs
#' @importFrom S4Vectors split
.get_velocity_files <- function(gr, L, Genome, Transcriptome = NULL,
                                out_path = ".", style = c("annotation",
                                  "genome", "Ensembl",
                                  "UCSC", "NCBI",
                                  "other"),
                                isoform_action = c("collapse", "separate"),
                                exon_option = c("full", "junction"),
                                transcript_id = "transcript_id",
                                gene_id = "gene_id",
                                transcript_version = "transcript_version",
                                gene_version = "gene_version",
                                version_sep = ".",
                                compress_fa = FALSE, width = 80L) {
  tx_path <- NULL
  L <- as.integer(L)
  width <- as.integer(width)
  isoform_action <- match.arg(isoform_action)
  exon_option <- match.arg(exon_option)
  c(out_path, tx_path) %<-% validate_velocity_input(L, Genome, Transcriptome,
    out_path, compress_fa, width, exon_option)
  style <- match.arg(style)
  fields <- names(mcols(gr))
  gr <- standardize_tags(gr, gene_id, transcript_id)
  gr$exon_rank <- NULL # not to cause trouble with extractTranscriptSeqs
  gr <- gr[gr$type == "exon" & (as.vector(strand(gr)) %in% c("+", "-"))]
  c(Genome, gr) %<-% match_style(Genome, gr, style)
  gr <- subset_annot(Genome, gr)
  c(Genome, gr) %<-% annot_circular(Genome, gr)
  genome(gr) <- genome(Genome)[seqlevels(gr)]
  if (!is.null(gene_version)) {
    names(mcols(gr))[fields == gene_version] <- "gene_version"
    gr$gene_id <- paste(gr$gene_id, gr$gene_version, sep = version_sep)
  }
  if (!is.null(transcript_version)) {
    names(mcols(gr))[fields == transcript_version] <- "transcript_version"
    gr$transcript_id <- paste(gr$transcript_id, gr$transcript_version,
      sep = version_sep)
  }
  if (is(Transcriptome, "DNAStringSet")) {
    # Check that transcript IDs in GRanges match those in the transcriptome
    tx_overlap <- check_tx(gr$transcript_id, names(Transcriptome))
    gr <- gr[gr$transcript_id %in% tx_overlap]
  } else {
    Transcriptome <- tx_path
  }
  tr2g_cdna <- tr2g_GRanges(gr, gene_name = NULL,
    transcript_version = NULL,
    gene_version = NULL) # version already added
  if (isoform_action == "collapse") {
    message("Collapsing gene isoforms")
    grl <- GenomicRanges::split(gr, gr$gene_id)
    grl <- GenomicRanges::reduce(grl)
    # Get intronic ranges
    introns <- get_intron_flanks(grl, L, FALSE)
    if (exon_option == "junction") {
      message("Extracting exon-exon junctions")
      grt <- GenomicRanges::split(gr, gr$transcript_id)
      exons <- get_intron_flanks(grt, L, TRUE)[[2]]
      Transcriptome <- extractTranscriptSeqs(Genome, exons)
      tr2g_cdna <- tr2g_junction(tr2g_cdna, names(Transcriptome))
    }
  } else {
    grl <- GenomicRanges::split(gr, gr$transcript_id)
    grl <- revElements(grl, any(strand(grl) == "-"))
    if (exon_option == "junction") {
      message("Extracting exon-exon junctions")
      c(introns, exons) %<-% get_intron_flanks(grl, L, TRUE)
      Transcriptome <- extractTranscriptSeqs(Genome, exons)
      tr2g_cdna <- tr2g_junction(tr2g_cdna, names(Transcriptome))
    } else {
      introns <- get_intron_flanks(grl, L, FALSE)
    }
  }
  if (is.null(Transcriptome) && exon_option == "full") {
    message("Extracting transcriptome from genome")
    if (isoform_action != "collapse") {
      Transcriptome <- extractTranscriptSeqs(Genome, grl)
    } else {
      # To distinguish from grl, which is split at gene level
      grt <- GenomicRanges::split(gr, gr$transcript_id)
      grt <- revElements(grt, any(strand(grt) == "-"))
      Transcriptome <- extractTranscriptSeqs(Genome, grt)
    }
    # Again, make sure that all transcripts in tr2g are in the transcriptome
    tr2g_cdna <- tr2g_cdna[tr2g_cdna$transcript %in% names(Transcriptome), ]
  }
  write_velocity_output(out_path, introns, Genome, Transcriptome,
    isoform_action, exon_option, tr2g_cdna, compress_fa, width)
}

#' Get files required for RNA velocity with bustools
#'
#' Computation of RNA velocity requires the number of unspliced transcripts,
#' which can be quantified with reads containing intronic sequences. This
#' function extracts intronic sequences flanked by L-1 bases of exonic sequences
#' where L is the biological read length of the single cell technology of
#' interest. The flanking exonic sequences are included for reads partially
#' mapping to an intron and an exon. 
#' 
#' @inheritParams .get_velocity_files
#' @inheritParams tr2g_EnsDb
#' @param X Gene annotation with transcript and exon information. It can be a
#' path to a GTF file with annotation of exon coordinates of
#' transcripts, preferably from Ensembl. In the metadata, the following fields
#' are required: type (e.g. whether the range of interest is a gene or
#' transcript or exon or CDS), gene ID, and transcript ID. These
#' fields need not to have standard names, as long as their names are specified
#' in arguments of this function. It can also be a \code{\link{TxDb}} object,
#' such as from the Bioconductor package
#' \code{TxDb.Hsapiens.UCSC.hg38.knownGene}. It can also be a
#' \code{\link{EnsDb}} object.
#' @param use_transcript_version Logical, whether to include version number in
#' the Ensembl transcript ID.
#' @param \dots Extra arguments for methods.
#' @return The following files will be written to disk in the directory
#' `out_path`:
#' \describe{
#' \item{cDNA_introns.fa}{A fasta file containing both the spliced transcripts
#' and the flanked intronic sequences. The intronic sequences are flanked by L-1
#' nt of exonic sequences to capture reads from nascent transcript partially
#' mapping to exons. If the exon is shorter than 2*(L-1) nt, then the entire
#' exon will be included in the intronic sequence. This will be used to build
#' the `kallisto` index.}
#' \item{cDNA_tx_to_capture.txt}{A text file of transcript IDs of spliced
#' transcripts. If `exon_option == "junction"`, then IDs of the exon-exon
#' junctions. These IDs will have the pattern <transcript ID>-Jx, where x is a
#' number differentiating between different junctions of the same transcript.
#' Here x will always be ordered from 5' to 3' as on the plus strand.}
#' \item{introns_tx_to_capture.txt}{A text file of IDs of introns. The names
#' will have the pattern <transcript ID>-Ix, where x is a number differentiating
#' between introns of the same transcript. If all transcripts of the same gene
#' are collapsed before inferring intronic sequences, gene ID will be used in
#' place of transcript ID. Here x will always be ordered from 5' to 3' as on the
#' plus strand.}
#' \item{tr2g.txt}{A text file with two columns matching transcripts and introns
#' to genes. The first column is transcript or intron ID, and the second column
#' is the corresponding gene ID. The part for transcripts are generated from
#' the gene annotation supplied.}
#' }
#' Nothing is returned into the R session.
#' @export
#' @examples
#' # Use toy example
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file <- paste0(toy_path, "/velocity_annot.gtf")
#' genome <- Biostrings::readDNAStringSet(paste0(toy_path, "/velocity_genome.fa"))
#' transcriptome <- paste0(toy_path, "/velocity_tx.fa")
#' get_velocity_files(file, 11, genome, transcriptome, ".",
#'   gene_version = NULL, transcript_version = NULL)
setGeneric("get_velocity_files",
  function(X, L, Genome, Transcriptome = NULL, out_path = ".",
             style = c("annotation", "genome", "Ensembl", "UCSC", "NCBI",
               "other"),
             isoform_action = c("collapse", "separate"),
             exon_option = c("full", "junction"),
             compress_fa = FALSE, width = 80L, ...)
    standardGeneric("get_velocity_files"),
  signature = "X")

#' @rdname get_velocity_files
#' @export
setMethod("get_velocity_files", "GRanges",
  function(X, L, Genome, Transcriptome = NULL, out_path = ".",
             style = c("annotation", "genome", "Ensembl", "UCSC", "NCBI",
               "other"),
             isoform_action = c("collapse", "separate"),
             exon_option = c("full", "junction"),
             compress_fa = FALSE, width = 80L,
             transcript_id = "transcript_id", gene_id = "gene_id",
             transcript_version = "transcript_version",
             gene_version = "gene_version", version_sep = ".") {
    .get_velocity_files(X, L, Genome, Transcriptome, out_path, style,
      isoform_action, exon_option, transcript_id, gene_id,
      transcript_version, gene_version, version_sep,
      compress_fa, width)
  }
)

#' @rdname get_velocity_files
#' @param is_circular Logical vector of the same length as the number of
#' sequences in the annotation and with the same names as the sequences,
#' indicating whether the sequence is circular. If `NULL`, then all sequences
#' will be assumed to be linear.
#' @export
setMethod("get_velocity_files", "character",
  function(X, L, Genome, Transcriptome = NULL, out_path = ".",
             style = c("annotation", "genome", "Ensembl", "UCSC", "NCBI",
               "other"),
             isoform_action = c("collapse", "separate"),
             exon_option = c("full", "junction"),
             compress_fa = FALSE, width = 80L,
             is_circular = NULL,
             transcript_id = "transcript_id", gene_id = "gene_id",
             transcript_version = "transcript_version",
             gene_version = "gene_version", version_sep = ".") {
    file <- normalizePath(X, mustWork = TRUE)
    check_gff("gtf", file, transcript_id, gene_id)
    gr <- plyranges::read_gff(file)
    if (is.null(is_circular)) {
      message("Assuming that all chromosomes are linenar.")
      isCircular(gr) <- setNames(rep(FALSE, length(seqlevels(gr))),
        seqlevels(gr))
    }
    .get_velocity_files(gr, L, Genome, Transcriptome, out_path, style,
      isoform_action, exon_option, transcript_id, gene_id,
      transcript_version, gene_version, version_sep,
      compress_fa, width)
  }
)

#' @rdname get_velocity_files
#' @importFrom GenomicFeatures exonsBy
#' @export
setMethod("get_velocity_files", "TxDb",
  function(X, L, Genome, Transcriptome, out_path,
             style = c("annotation", "genome", "Ensembl", "UCSC", "NCBI",
               "other"),
             isoform_action = c("collapse", "separate"),
             exon_option = c("full", "junction"),
             compress_fa = FALSE, width = 80L) {
    exons_by_tx <- function(X, tx_id, tx) {
      gr <- exonsBy(X, by = "tx") # Will be numbers
      # Remove transcripts that aren't mapped to genes
      gr <- gr[names(gr) %in% tx_id]
      names(gr) <- tx[match(names(gr), tx_id)]
      gr
    }
    tx_path <- NULL
    L <- as.integer(L)
    width <- as.integer(width)
    isoform_action <- match.arg(isoform_action)
    exon_option <- match.arg(exon_option)
    c(out_path, tx_path) %<-%
      validate_velocity_input(L, Genome, Transcriptome, out_path,
        compress_fa, width, exon_option)
    style <- match.arg(style)
    c(Genome, X) %<-% match_style(Genome, X, style)
    X <- subset_annot(Genome, X)
    tr2g_cdna <- tr2g_TxDb(X)
    if (is(Transcriptome, "DNAStringSet")) {
      tx_overlap <- check_tx(tr2g_cdna$transcript, names(Transcriptome))
      tr2g_cdna <- tr2g_cdna[tr2g_cdna$transcript %in% tx_overlap, ]
    } else if (is.character(Transcriptome)) {
      Transcriptome <- tx_path
    }
    if (isoform_action == "collapse") {
      message("Collapsing gene isoforms")
      grl <- exonsBy(X, by = "gene")
      grl <- GenomicRanges::reduce(grl)
      c(Genome, grl) %<-% annot_circular(Genome, grl)
      genome(grl) <- genome(Genome)[seqlevels(grl)]
      # Get intronic ranges
      introns <- get_intron_flanks(grl, L, FALSE)
      if (exon_option == "junction") {
        message("Extracting exon-exon junctions")
        grt <- exons_by_tx(X, tr2g_cdna$tx_id, tr2g_cdna$transcript)
        c(Genome, grt) %<-% annot_circular(Genome, grt)
        genome(grt) <- genome(Genome)[seqlevels(grt)]
        exons <- get_intron_flanks(grt, L, TRUE)[[2]]
        Transcriptome <- extractTranscriptSeqs(Genome, exons)
        tr2g_cdna <- tr2g_junction(tr2g_cdna, names(Transcriptome))
      }
    } else {
      # tr2g_cdna has already been filtered if it needs to be filtered
      grl <- exons_by_tx(X, tr2g_cdna$tx_id, tr2g_cdna$transcript)
      c(Genome, grl) %<-% annot_circular(Genome, grl)
      genome(grl) <- genome(Genome)[seqlevels(grl)]
      if (exon_option == "junction") {
        message("Extracting exon-exon junctions")
        c(introns, exons) %<-% get_intron_flanks(grl, L, TRUE)
        Transcriptome <- extractTranscriptSeqs(Genome, exons)
        tr2g_cdna <- tr2g_junction(tr2g_cdna, names(Transcriptome))
      } else {
        introns <- get_intron_flanks(grl, L, FALSE)
      }
    }
    if (is.null(Transcriptome) && exon_option == "full") {
      message("Extracting transcriptome from genome")
      if (isoform_action != "collapse") {
        Transcriptome <- extractTranscriptSeqs(Genome, grl)
      } else {
        grt <- exons_by_tx(X, tr2g_cdna$tx_id, tr2g_cdna$transcript)
        c(Genome, grt) %<-% annot_circular(Genome, grt)
        genome(grt) <- genome(Genome)[seqlevels(grt)]
        Transcriptome <- extractTranscriptSeqs(Genome, grt)
      }
      # Again, make sure that all transcripts in tr2g are in the transcriptome
      tr2g_cdna <- tr2g_cdna[tr2g_cdna$transcript %in% names(Transcriptome), ]
    }
    tr2g_cdna <- tr2g_cdna[, c("transcript", "gene")]
    write_velocity_output(out_path, introns, Genome, Transcriptome,
      isoform_action, exon_option, tr2g_cdna, compress_fa, width)
  }
)

#' @rdname get_velocity_files
#' @importFrom AnnotationFilter AnnotationFilterList SeqNameFilter TxIdFilter
#' @importFrom ensembldb exonsBy
#' @export
setMethod("get_velocity_files", "EnsDb",
  function(X, L, Genome, Transcriptome, out_path,
             style = c("annotation", "genome", "Ensembl", "UCSC", "NCBI",
               "other"),
             isoform_action = c("collapse", "separate"),
             exon_option = c("full", "junction"),
             compress_fa = FALSE, width = 80L,
             use_transcript_version = TRUE, use_gene_version = TRUE) {
    exons_by_tx <- function(X, tx, chrs_use, use_transcript_version) {
      if (use_transcript_version) {
        # tr2g_cdna has already been filtered if it needs to be filtered
        tx_nv <- str_remove(tx, "\\.\\d+$")
        filter_use <- AnnotationFilterList(SeqNameFilter(chrs_use),
          TxIdFilter(tx_nv))
      } else {
        filter_use <- AnnotationFilterList(SeqNameFilter(chrs_use),
          TxIdFilter(tx))
      }
      gr <- exonsBy(X, by = "tx", filter = filter_use)
      seqlevels(gr) <- seqlevelsInUse(gr)
      if (use_transcript_version) {
        # Add transcript version
        names(gr) <- tx[match(names(gr), tx_nv)]
      }
      gr
    }
    tx_path <- NULL
    L <- as.integer(L)
    width <- as.integer(width)
    isoform_action <- match.arg(isoform_action)
    exon_option <- match.arg(exon_option)
    c(out_path, tx_path) %<-% validate_velocity_input(L, Genome,
      Transcriptome, out_path, compress_fa, width, exon_option)
    style <- match.arg(style)
    c(Genome, X) %<-% match_style(Genome, X, style)
    chrs_use <- intersect(seqlevels(X), seqlevels(Genome))
    tr2g_cdna <- tr2g_EnsDb(X, use_gene_name = FALSE,
      use_transcript_version = use_transcript_version,
      use_gene_version = use_gene_version)
    if (is(Transcriptome, "DNAStringSet")) {
      tx_overlap <- check_tx(tr2g_cdna$transcript, names(Transcriptome))
      tr2g_cdna <- tr2g_cdna[tr2g_cdna$transcript %in% tx_overlap, ]
    } else if (is.character(Transcriptome)) {
      Transcriptome <- tx_path
    }
    # extractTranscriptSeqs does not use transcript version numbers
    if (isoform_action == "collapse") {
      message("Collapsing gene isoforms")
      filter_use <- AnnotationFilterList(SeqNameFilter(chrs_use))
      grl <- exonsBy(X, by = "gene", filter = filter_use)
      # Add version number if used
      if (use_gene_version) {
        g_nv <- str_remove(tr2g_cdna$gene, "\\.\\d+$")
        names(grl) <- tr2g_cdna$gene[match(names(grl), g_nv)]
      }
      grl <- GenomicRanges::reduce(grl)
      seqlevels(grl) <- seqlevelsInUse(grl)
      c(Genome, grl) %<-% annot_circular(Genome, grl)
      genome(grl) <- genome(Genome)[seqlevels(grl)]
      # Get intronic ranges
      introns <- get_intron_flanks(grl, L, FALSE)
      if (exon_option == "junction") {
        message("Extracting exon-exon junctions")
        grt <- exons_by_tx(X, tr2g_cdna$transcript, chrs_use,
          use_transcript_version)
        c(Genome, grt) %<-% annot_circular(Genome, grt)
        genome(grt) <- genome(Genome)[seqlevels(grt)]
        exons <- get_intron_flanks(grt, L, TRUE)[[2]]
        Transcriptome <- extractTranscriptSeqs(Genome, exons)
        tr2g_cdna <- tr2g_junction(tr2g_cdna, names(Transcriptome))
      }
    } else {
      grl <- exons_by_tx(X, tr2g_cdna$transcript, chrs_use,
        use_transcript_version)
      c(Genome, grl) %<-% annot_circular(Genome, grl)
      genome(grl) <- genome(Genome)[seqlevels(grl)]
      if (exon_option == "junction") {
        message("Extracting exon-exon junctions")
        c(introns, exons) %<-% get_intron_flanks(grl, L, TRUE)
        Transcriptome <- extractTranscriptSeqs(Genome, exons)
        tr2g_cdna <- tr2g_junction(tr2g_cdna, names(Transcriptome))
      } else {
        introns <- get_intron_flanks(grl, L, FALSE)
      }
    }
    if (is.null(Transcriptome) && exon_option == "full") {
      message("Extracting transcriptome from genome")
      if (isoform_action != "collapse") {
        Transcriptome <- extractTranscriptSeqs(Genome, grl)
      } else {
        grt <- exons_by_tx(X, tr2g_cdna$transcript, chrs_use,
          use_transcript_version)
        c(Genome, grt) %<-% annot_circular(Genome, grt)
        genome(grt) <- genome(Genome)[seqlevels(grt)]
        Transcriptome <- extractTranscriptSeqs(Genome, grt)
      }
      tr2g_cdna <- tr2g_cdna[tr2g_cdna$transcript %in% names(Transcriptome), ]
    }
    write_velocity_output(out_path, introns, Genome, Transcriptome,
      isoform_action, exon_option, tr2g_cdna, compress_fa, width)
  })
