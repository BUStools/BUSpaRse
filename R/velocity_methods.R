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
#' @param genome Either a \code{\link{BSgenome}} or a \code{\link{XStringSet}}
#' object of genomic sequences, where the intronic sequences will be extracted
#' from. Naming convention of chromosomes in the genome must match that in the
#' annotation, e.g. consistently use 1 for chromosome 1 or consistently use chr1.
#' @param transcriptome Either a \code{\link{XStringSet}} or a path to a fasta
#' file (can be gzipped) of the transcriptome, which contains sequences of
#' spliced transcripts. This will be concatenated with the intronic sequences
#' to give one fasta file. For the exported function 
#' \code{\link{get_velocity_files}}, this argument can be missing. If it is
#' missing, then the transriptome sequences will be extracted from the genome
#' given the gene annotation. 
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
#' @return See \code{\link{get_velocity_files}}
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicFeatures makeTxDbFromGRanges

.get_velocity_files <- function(gr, L, genome, transcriptome, out_path,
                                short_exon_action = c("truncate", "include"),
                                isoform_action = c("separate", "collapse"),
                                transcript_id = "transcript_id",
                                gene_id = "gene_id", 
                                transcript_version = "transcript_version",
                                gene_version = "gene_version", 
                                version_sep = ".", exon_number = "exon_number",
                                exon_id = "exon_id") {
  L <- as.integer(L)
  validate_velocity_input(L, genome, transcriptome, out_path)
  short_exon_action <- match.arg(short_exon_action)
  isoform_action <- match.arg(isoform_action)
  fields <- names(mcols(gr))
  check_tag_present(exon_number, fields, error = TRUE)
  # I need this anyway and this will validate the GTF file.
  tr2g_cdna <- tr2g_GRanges(gr, gene_name = NULL, gene_id = gene_id,
                            transcript_id = transcript_id,
                            transcript_version = transcript_version,
                            gene_version = gene_version,
                            version_sep = version_sep)
  gr <- gr[gr$type == "exon" & (as.vector(strand(gr)) %in% c("+", "-"))]
  gr <- standardize_tags(gr, gene_id, transcript_id, exon_number)
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
    gr <- collapse_isoforms(gr)
  }
  if (short_exon_action == "include") {
    gr <- exclude_short_exons(gr, L, exon_id = exon_id)
  }
  # Calculate the appropriate flanking lengths
  metas <- mcols(gr) %>% 
    as.data.frame() %>% 
    mutate(width = width(gr),
           start = start(gr),
           strand = as.vector(strand(gr)),
           chr = as.vector(seqnames(gr)))
  metas <- get_flank_lengths(metas, L)
  
  # Get intronic ranges
  txdb <- makeTxDbFromGRanges(gr)
  introns <- get_flanked_introns(txdb, metas)
  write_velocity_output(out_path, introns, genome, transcriptome, tr2g_cdna,
                        metas)
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
#' @param X Gene annotation with transcript and exon information. It can be a 
#' path to a GTF file with annotation of exon coordinates of 
#' transcripts, preferably from Ensembl. The file must be formatted as in Ensembl.
#' In the metadata, the following fields are required: type (e.g. whether the 
#' range of interest is a gene or transcript or exon or CDS), gene ID, 
#' transcript ID, and exon number. These fields need not to have standard names,
#' as long as their names are specified in arguments of this function. It can
#' also be a \code{\link{TxDb}} object, such as from the Bioconductor package
#' \code{TxDb.Hsapiens.UCSC.hg38.knownGene}.
#' @param \dots Extra parameters for methods.
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
#' @examples 
#' # Use toy example
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file <- paste0(toy_path, "/velocity_annot.gtf")
#' genome <- Biostrings::readDNAStringSet(paste0(toy_path, "/velocity_genome.fa"))
#' transcriptome <- paste0(toy_path, "/velocity_tx.fa")
#' get_velocity_files(file, 11, genome, transcriptome, ".",
#'                    gene_version = NULL, transcript_version = NULL)
setGeneric("get_velocity_files", 
           function(X, L, genome, transcriptome, out_path,
                    short_exon_action = c("truncate", "include"),
                    isoform_action = c("separate", "collapse"), ...) 
             standardGeneric("get_velocity_files"),
           signature = "X")

#' @rdname get_velocity_files
#' @export
setMethod("get_velocity_files", "GRanges",
          function(X, L, genome, transcriptome, out_path,
                   short_exon_action = c("truncate", "include"),
                   isoform_action = c("separate", "collapse"),
                   transcript_id = "transcript_id",
                   gene_id = "gene_id", 
                   transcript_version = "transcript_version",
                   gene_version = "gene_version", 
                   version_sep = ".", exon_number = "exon_number",
                   exon_id = "exon_id") {
            X <- subset_annot(genome, X)
            if (missing(transcriptome)) {
              transcriptome <- extract_tx(genome, X)
            }
            .get_velocity_files(X, L, genome, transcriptome, out_path,
                                short_exon_action, isoform_action, 
                                transcript_id, gene_id, 
                                transcript_version, gene_version, version_sep, 
                                exon_number, exon_id)
          }
)

#' @rdname get_velocity_files
#' @export
setMethod("get_velocity_files", "character",
          function(X, L, genome, transcriptome, out_path,
                   short_exon_action = c("truncate", "include"),
                   isoform_action = c("separate", "collapse"),
                   transcript_id = "transcript_id",
                   gene_id = "gene_id", 
                   transcript_version = "transcript_version",
                   gene_version = "gene_version", 
                   version_sep = ".", exon_number = "exon_number",
                   exon_id = "exon_id") {
            file <- normalizePath(X, mustWork = TRUE)
            check_gff("gtf", file, transcript_id, gene_id)
            gr <- plyranges::read_gff(file)
            gr <- subset_annot(genome,gr)
            if (missing(transcriptome)) {
              transcriptome <- extract_tx(genome, gr)
            }
            .get_velocity_files(gr, L, genome, transcriptome, out_path,
                                short_exon_action, isoform_action, 
                                transcript_id, gene_id, 
                                transcript_version, gene_version, version_sep, 
                                exon_number, exon_id)
          }
)

#' @rdname get_velocity_files
#' @importFrom AnnotationDbi keys
#' @importFrom GenomicFeatures exons
#' @importFrom dplyr inner_join
#' @export
setMethod("get_velocity_files", "TxDb",
          function(X, L, genome, transcriptome, out_path,
                   short_exon_action = c("truncate", "include"),
                   isoform_action = c("separate", "collapse")) {
            short_exon_action <- match.arg(short_exon_action)
            isoform_action <- match.arg(isoform_action)
            X <- subset_annot(genome, X)
            if (missing(transcriptome)) {
              transcriptome <- extract_tx(genome, X)
            }
            if (short_exon_action == "include" || isoform_action == "collapse") {
              gr <- exons(X, columns = c("TXNAME", "GENEID", "EXONRANK", 
                                         "EXONID"))
              gr$type <- "exon"
              gr$TXNAME <- unlist(gr$TXNAME)
              gr$GENEID <- unlist(gr$GENEID)
              gr$EXONRANK <- unlist(gr$EXONRANK)
              .get_velocity_files(gr, L, genome, transcriptome, out_path,
                                  short_exon_action, isoform_action, 
                                  transcript_id = "TXNAME", gene_id = "GENEID", 
                                  transcript_version = NULL, 
                                  gene_version = NULL,
                                  exon_number = "EXONRANK",
                                  exon_id = "EXONID")
            } else {
              if (length(L) > 1 || L %% 1 > sqrt(.Machine$double.eps) || !is.atomic(L)) {
                stop("L must be an integer vector with length 1.")
              }
              L <- as.integer(L)
              out_path <- normalizePath(out_path, mustWork = FALSE)
              if (!dir.exists(out_path)) {
                dir.create(out_path)
              }
              tr2g_cdna <- tr2g_TxDb(X)
              metas <- AnnotationDbi::select(X, keys = keys(X, keytype = "EXONID"),
                                             keytype = "EXONID",
                                             columns = c("EXONCHROM", "EXONSTART",
                                                         "EXONEND", "EXONSTRAND",
                                                         "EXONRANK",
                                                         "GENEID", "TXNAME"))
              metas <- metas[complete.cases(metas),]
              names_use <- data.frame(txdb = c("EXONID", "EXONCHROM", 
                                               "EXONSTART", "EXONEND", 
                                               "EXONSTRAND", "EXONRANK",
                                               "GENEID", "TXNAME"),
                                      wanted = c("exon_id", "chr", "start", 
                                                 "end", "strand", "exon_number",
                                                 "gene_id", "transcript_id"))
              names(metas) <- names_use$wanted[match(names(metas), names_use$txdb)]
              metas$width <- metas$end - metas$start + 1L
              metas <- get_flank_lengths(metas, L)
              introns <- get_flanked_introns(X, metas)
              write_velocity_output(out_path, introns, genome, transcriptome, 
                                    tr2g_cdna, metas)
            }
          }
)
