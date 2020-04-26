.annots_from_fa <- function(fa) {
  # Avoid R CMD check note
  transcript_id <- type <- cr <- gene_biotype <- transcript_biotype <- 
    gene_id <- gene_symbol <- description <- source <- accession <- 
    start <- end <- strand <- NULL
  out <- tibble(transcript_id = str_extract(names(fa), "^[a-zA-Z\\d-\\.]+"),
                type = str_extract(names(fa), paste0("(?<=", transcript_id, " ).*?(?=\\s)")),
                cr = str_extract(names(fa),
                                 "(?<=((chromosome)|(scaffold)):).*?(?=\\s)"),
                gene_biotype = str_extract(names(fa), "(?<=gene_biotype:).*?(?=\\s)"),
                transcript_biotype = str_extract(names(fa), "(?<=transcript_biotype:).*?(?=\\s)"),
                gene_id = str_extract(names(fa), "(?<=gene:).*?(?=\\s)"),
                gene_symbol = str_extract(names(fa), "(?<=gene_symbol:).*?(?=\\s)"),
                description = str_extract(names(fa), "(?<=description:).*?(?= \\[)"),
                source = str_extract(names(fa), "(?<=Source:).*?(?=;)"),
                accession = str_extract(names(fa), "(?<=Acc:).*?(?=\\])")) %>% 
    tidyr::separate(cr, into = c("genome", "seqnames", "start", "end", "strand"), sep = ":") %>% 
    mutate(start = as.integer(start),
           end = as.integer(end),
           strand = dplyr::case_when(
             strand == "1" ~ "+",
             strand == "-1" ~ "-",
             TRUE ~ "*"
           ))
}

#' Get genome annotation from Ensembl FASTA file
#' 
#' Ensembl FASTA files for RNA contain much of the information contained in GTF
#' files, such as chromosome, genome assembly version, coordinates, strand,
#' gene ID, gene symbol, and gene description. Given such a FASTA file, this
#' function can extract all the genome annotation information and return a
#' data frame or a `GRanges` object.
#' 
#' @param file Path to the FASTA file to be read. The file can remain gzipped.
#' @return `annots_from_fa_df` returns a data frame and 
#' `annots_from_fa_GRanges` returns `GRanges`.
#' @importFrom tidyr separate
#' @importFrom dplyr case_when
#' @export
#' @examples
#' fn <- system.file("testdata/fasta_test.fasta", package = "BUSpaRse")
#' annots_from_fa_df(fn)
#' annots_from_fa_GRanges(fn)
#' 
#' @return A data frame with genome annotations.
#' @export

annots_from_fa_df <- function(file) {
  fa <- readDNAStringSet(file)
  .annots_from_fa(fa)
}

#' @rdname annots_from_fa_df
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
annots_from_fa_GRanges <- function(file) {
  df <- annots_from_fa_df(file)
  gn <- unique(df$genome)
  df <- dplyr::select(df, -genome)
  DF <- as(df[, c("type", "transcript_id", "gene_id", "gene_symbol", 
                  "transcript_biotype", "gene_biotype", "description",
                  "source", "accession")], "DataFrame")
  ranges <- IRanges(start = df$start, end = df$end)
  gr <- GRanges(seqnames = as.factor(as(df$seqnames, "Rle")),
                ranges = ranges, strand = as.factor(as(df$strand, "Rle")),
                mcols = DF)
  names(mcols(gr)) <- str_remove(names(mcols(gr)), "^mcols\\.")
  genome(gr) <- gn
  gr
}
