#' @include sparse_matrix.R
NULL

#' Get transcript and gene info from Ensembl
#'
#' This function queries Ensembl biomart to convert transcript IDs to gene IDs.
#'
#' @param species Character vector of length 1, Latin name of the species of
#' interest.
#' @param type Character, must be one of "vertebrate", "metazoa", "plant",
#' "fungus" and "protist". Passing "vertebrate" will use the default
#' www.ensembl.org host. Gene annotation of some common invertebrate model
#' organisms, such as _Drosophila melanogaster_, are available on www.ensembl.org
#' so for these invertebrate model organisms, "vertebrate" can be used for this
#' argument. Passing values other than "vertebrate" will use other Ensembl hosts.
#' For animals absent from www.ensembl.org, try "metazoa".
#' @param ensembl_version Integer version number of Ensembl (e.g. 94 for the
#' October 2018 release). This argument defaults to \code{NULL}, which will use
#' the current release of Ensembl. Use
#' \code{\link{listEnsemblArchives}} to see the version number corresponding
#' to the Ensembl release of a particular date. The version specified here must
#' match the version of Ensembl where the transcriptome used to build the
#' kallisto index was downloaded.
#' @param other_attrs Character vector. Other attributes to get from Ensembl,
#' such as gene symbol and position on the genome.
#' Use \code{\link{listAttributes}} to see which attributes are available.
#' @param use_gene_name Logical, whether to get gene names.
#' @param use_transcript_version Logical, whether to include version number in
#' the Ensembl transcript ID. To decide whether to
#' include transcript version number, check whether version numbers are included
#' in the `transcripts.txt` in the `kallisto` output directory. If that file
#' includes version numbers, then trannscript version numbers must be included
#' here as well. If that file does not include version numbers, then transcript
#' version numbers must not be included here.
#' @param use_gene_version Logical, whether to include version number in the
#' Ensembl gene ID. Unlike transcript
#' version number, it's up to you whether to include gene version number.
#' @param verbose Whether to display progress.
#' @param \dots Othe arguments to be passed to \code{\link{useEnsembl}},
#' such as mirror. Note that setting mirrors other than the default, e.g. uswest,
#' does not work for archived versions.
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom stats setNames
#' @return A data frame with at least 2 columns: \code{gene} for gene ID,
#' \code{transcript} for transcript ID, and optionally \code{gene_name}
#' for gene names. If \code{other_attrs} has been specified, then those will
#' also be columns in the data frame returned.
#' @family functions to retrieve transcript and gene info
#' @export
#' @examples
#' tr2g <- tr2g_ensembl(species = "Felis catus", other_attrs = "description")
#' # This will use plants.ensembl.org as host instead of www.ensembl.org
#' tr2g <- tr2g_ensembl(species = "Arabidopsis thaliana", type = "plant")
tr2g_ensembl <- function(species, type = c("vertebrate", "metazoa", "plant",
                           "fungus", "protist"),
                         other_attrs = NULL,
                         use_gene_name = TRUE,
                         use_transcript_version = TRUE,
                         use_gene_version = TRUE,
                         ensembl_version = NULL,
                         verbose = TRUE, ...) {
  # Validate arguments
  check_char1(setNames(c(species, type), c("species", "type")))
  type <- match.arg(type)
  if (!is.null(ensembl_version) && !is.numeric(ensembl_version)) {
    stop("ensembl_version must be integer.")
  }
  if (!is.null(other_attrs) &&
    (!is.atomic(other_attrs) || !is.character(other_attrs))) {
    stop("other_attrs must be an atomic character vector.")
  }
  if (type != "vertebrate" && (use_transcript_version || use_gene_version)) {
    message("Version is only available to vertebrates.")
    use_transcript_version <- use_gene_version <- FALSE
  }
  ds_name <- species2dataset(species, type)
  host_pre <- switch(type,
    vertebrate = "www",
    metazoa = "metazoa",
    plant = "plants",
    fungus = "fungi",
    protist = "protists")
  mart_use <- paste(host_pre, "mart", sep = "_")
  host_use <- paste0(host_pre, ".ensembl.org")
  if (type == "vertebrate") mart_use <- "ensembl"
  if (verbose) {
    message(paste("Querying biomart for transcript and gene IDs of",
      species))
  }
  mart <- useEnsembl(biomart = mart_use, dataset = ds_name, host = host_use,
    version = ensembl_version, ...)
  attrs_use <- c("ensembl_transcript_id", "ensembl_gene_id", other_attrs)
  if (use_gene_name) {
    attrs_use <- c(attrs_use, "external_gene_name")
  }
  if (use_transcript_version) {
    attrs_use[1] <- paste(attrs_use[1], "version", sep = "_")
  }
  if (use_gene_version) {
    attrs_use[2] <- paste(attrs_use[2], "version", sep = "_")
  }
  out <- getBM(attrs_use, mart = mart)
  names(out)[seq_len(2)] <- c("transcript", "gene")
  names(out)[names(out) == "external_gene_name"] <- "gene_name"
  out
}

#' Get transcript and gene info from GRanges
#'
#' Internal use, for GRanges from GTF files
#'
#' @param gr A \code{\link{GRanges}} object. The metadata columns should be
#' atomic vectors, not lists.
#' @param type_use Character vector, the values taken by the \code{type} field
#' in the GTF file that denote the desired transcripts. This can be "exon",
#' "transcript", "mRNA", and etc.
#' @param transcript_id Character vector of length 1. Tag in \code{attribute}
#' field corresponding to transcript IDs. This argument must be supplied and
#' cannot be \code{NA} or \code{NULL}. Will throw error if tag indicated in this
#' argument does not exist.
#' @param gene_id Character vector of length 1. Tag in \code{attribute}
#' field corresponding to gene IDs. This argument must be supplied and
#' cannot be \code{NA} or \code{NULL}. Note that this is different from gene
#' symbols, which do not have to be unique. This can be Ensembl or Entrez IDs.
#' However, if the gene symbols are in fact unique for each gene, you may
#' supply the tag for human readable gene symbols to this argument. Will throw
#' error if tag indicated in this argument does not exist.
#' @param gene_name Character vector of length 1. Tag in \code{attribute}
#' field corresponding to gene symbols. This argument can be \code{NA} or
#' \code{NULL} if you are fine with non-human readable gene IDs and do not wish
#' to extract human readable gene symbols.
#' @param transcript_version Character vector of length 1. Tag in \code{attribute}
#' field corresponding to _transcript_ version number. If your GTF file does not
#' include transcript version numbers, or if you do not wish to include the
#' version number, then use \code{NULL} for this argument. To decide whether to
#' include transcript version number, check whether version numbers are included
#' in the `transcripts.txt` in the `kallisto` output directory. If that file
#' includes version numbers, then trannscript version numbers must be included
#' here as well. If that file does not include version numbers, then transcript
#' version numbers must not be included here.
#' @param gene_version Character vector of length 1. Tag in \code{attribute}
#' field corresponding to _gene_ version number. If your GTF file does not
#' include gene version numbers, or if you do not wish to include the
#' version number, then use \code{NULL} for this argument. Unlike transcript
#' version number, it's up to you whether to include gene version number.
#' @param version_sep Character to separate bewteen the main ID and the version
#' number. Defaults to ".", as in Ensembl.
#' @return A data frame at least 2 columns: \code{gene} for gene ID,
#' \code{transcript} for transcript ID, and optionally, \code{gene_name} for
#' gene names.
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect
#' @importFrom dplyr distinct
#' @importFrom S4Vectors mcols
tr2g_GRanges <- function(gr, type_use = "exon", transcript_id = "transcript_id",
                         gene_id = "gene_id", gene_name = "gene_name",
                         transcript_version = "transcript_version",
                         gene_version = "gene_version", version_sep = ".") {
  tags <- names(mcols(gr))
  check_tag_present(c(transcript_id, gene_id), tags, error = TRUE)
  # Will do nothing if all are NULL
  check_tag_present(c(gene_name, transcript_version, gene_version),
    tags, error = FALSE)
  gr <- gr[!is.na(mcols(gr)[[transcript_id]])]
  gr <- gr[gr$type %in% type_use]
  if (length(gr) == 0) {
    stop(paste("No entry has types", paste(type_use, collapse = ", ")))
  }
  out <- data.frame(transcript = mcols(gr)[[transcript_id]],
    gene = mcols(gr)[[gene_id]],
    stringsAsFactors = FALSE)
  if (!is.null(gene_name) && gene_name %in% tags) {
    out$gene_name <- mcols(gr)[[gene_name]]
  }
  if (!is.null(transcript_version) && transcript_version %in% tags) {
    tv <- mcols(gr)[[transcript_version]]
    out$transcript <- paste(out$transcript, tv, sep = version_sep)
  }
  if (!is.null(gene_version) && gene_version %in% tags) {
    gv <- mcols(gr)[[gene_version]]
    out$gene <- paste(out$gene, gv, sep = version_sep)
  }
  distinct(out)
}

#' Get transcript and gene info from GTF file
#'
#' This function reads a GTF file and extracts the transcript ID and
#' corresponding gene ID. This function assumes that the GTF file is properly
#' formatted. See \url{http://mblab.wustl.edu/GTF2.html} for a detailed
#' description of proper GTF format. Note that GFF3 files have a somewhat
#' different and more complicated format in the attribute field, which this
#' function does not support. See \url{http://gmod.org/wiki/GFF3} for a detailed
#' description of proper GFF3 format. To extract transcript and gene information
#' from GFF3 files, see the function \code{\link{tr2g_gff3}} in this package.
#'
#' Transcript and gene versions may not be present in all GTF files, so these
#' arguments are optional. This function has arguments for transcript and gene
#' version numbers because Ensembl IDs have version numbers. For Ensembl IDs, we
#' recommend including the version number, since a change in version number
#' signals a change in the entity referred to by the ID after reannotation. If a
#' version is used, then it will be appended to the ID, separated by
#' \code{version_sep}.
#'
#' The transcript and gene IDs are The \code{attribute} field (the last
#' field) of GTF files can be complicated and inconsistent across different
#' sources. Please check the \code{attribute} tags in your GTF file and consider
#' the arguments of this function carefully. The defaults are set according to
#' Ensembl GTF files; defaults may not work for files from other sources. Due to
#' the general lack of standards for the \code{attribute} field, you may need to
#' further clean up the output of this function.
#'
#' @param file Path to a GTF file to be read. The file can remain gzipped.
#' @inheritParams tr2g_GRanges
#' @inheritParams tr2g_ensembl
#' @return A data frame at least 2 columns: \code{gene} for gene ID,
#' \code{transcript} for transcript ID, and optionally, \code{gene_name} for
#' gene names.
#' @importFrom plyranges read_gff
#' @family functions to retrieve transcript and gene info
#' @export
#' @examples
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file_use <- paste(toy_path, "gtf_test.gtf", sep = "/")
#' # Default
#' tr2g <- tr2g_gtf(file = file_use, verbose = FALSE)
#' # Excluding version numbers
#' tr2g <- tr2g_gtf(file = file_use, transcript_version = NULL,
#'   gene_version = NULL)
tr2g_gtf <- function(file, type_use = "exon", transcript_id = "transcript_id",
                     gene_id = "gene_id", gene_name = "gene_name",
                     transcript_version = "transcript_version",
                     gene_version = "gene_version", version_sep = ".",
                     verbose = TRUE) {
  # Validate arguments
  check_char1(setNames(file, "file"))
  file <- normalizePath(file, mustWork = TRUE)
  check_gff("gtf", file, transcript_id, gene_id)
  if (verbose) {
    message(paste("Reading GTF file."))
  }
  gr <- read_gff(file)
  tr2g_GRanges(gr, type_use, transcript_id, gene_id, gene_name,
    transcript_version, gene_version, version_sep)
}

#' Get transcript and gene info from GFF3 file
#'
#' This function reads a GFF3 file and extracts the transcript ID and
#' corresponding gene ID. This function assumes that the GFF3 file is properly
#' formatted. See \url{http://gmod.org/wiki/GFF3} for a detailed
#' description of proper GFF3 format. Note that GTF files have a somewhat
#' different and simpler format in the attribute field, which this function does
#' not support. See \url{http://mblab.wustl.edu/GTF2.html} for a detailed
#' description of proper GTF format. To extract transcript and gene information
#' from GTF files, see the function \code{\link{tr2g_gtf}} in this package.
#' Some files bearing the \code{.gff3} are in fact more like the GTF format. If
#' this is so, then change the extension to \code{.gtf} and use the function
#' \code{\link{tr2g_gtf}} in this package instead.
#'
#' Transcript and gene versions may not be present in all GTF files, so these
#' arguments are optional. This function has arguments for transcript and gene
#' version numbers because Ensembl IDs have version numbers. For Ensembl IDs, we
#' recommend including the version number, since a change in version number
#' signals a change in the entity referred to by the ID after reannotation. If a
#' version is used, then it will be appended to the ID, separated by
#' \code{version_sep}.
#'
#' The transcript and gene IDs are The \code{attribute} field (the last
#' field) of GTF files can be complicated and inconsistent across different
#' sources. Please check the \code{attribute} tags in your GTF file and consider
#' the arguments of this function carefully. The defaults are set according to
#' Ensembl GTF files; defaults may not work for files from other sources. Due to
#' the general lack of standards for the \code{attribute} field, you may need to
#' further clean up the output of this function.
#'
#' @inheritParams tr2g_gtf
#' @return A data frame at least 2 columns: \code{gene} for gene ID,
#' \code{transcript} for transcript ID, and optionally, \code{gene_name} for
#' gene names.
#' @family functions to retrieve transcript and gene info
#' @importFrom plyranges read_gff3
#' @importFrom stringr str_split
#' @importFrom dplyr left_join distinct
#' @importFrom tidyr unite
#' @export
#' @examples
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file_use <- paste(toy_path, "gff3_test.gff3", sep = "/")
#' # Default
#' tr2g <- tr2g_gff3(file = file_use, verbose = FALSE)
#' # Excluding version numbers
#' tr2g <- tr2g_gff3(file = file_use, transcript_version = NULL,
#'   gene_version = NULL)
tr2g_gff3 <- function(file, type_use = "mRNA", transcript_id = "transcript_id",
                      gene_id = "gene_id", gene_name = "Name",
                      transcript_version = "version",
                      gene_version = "version", version_sep = ".",
                      verbose = TRUE) {
  # Validate arguments
  check_char1(setNames(file, "file"))
  file <- normalizePath(file, mustWork = TRUE)
  check_gff("gff3", file, transcript_id, gene_id)
  if (verbose) {
    message(paste("Reading GFF3 file."))
  }
  gr <- read_gff3(file)
  tags <- names(mcols(gr))
  check_tag_present(c(transcript_id, gene_id), tags, error = TRUE)
  # Will do nothing if all are NULL
  check_tag_present(c(gene_name, transcript_version, gene_version),
    tags, error = FALSE)
  # Get transcript ID
  gr_tx <- gr[!is.na(mcols(gr)[[transcript_id]])]
  gr_tx <- gr_tx[gr_tx$type %in% type_use]
  if (length(gr_tx) == 0) {
    stop(paste("No entry has types", paste(type_use, collapse = ", ")))
  }
  genes <- str_split(gr_tx$Parent, ":", simplify = TRUE)[, 2]
  out <- data.frame(transcript = mcols(gr_tx)[[transcript_id]],
    gene = genes,
    stringsAsFactors = FALSE)
  if (!is.null(transcript_version) && transcript_version %in% tags) {
    tv <- mcols(gr_tx)[[transcript_version]]
    out$transcript <- paste(out$transcript, tv, sep = version_sep)
  }
  # Get gene name and version
  get_gene_name <- !is.null(gene_name) && gene_name %in% tags
  get_gene_version <- !is.null(gene_version) && gene_version %in% tags
  if (get_gene_name || get_gene_version) {
    gr_g <- gr[!is.na(mcols(gr)[[gene_id]])]
    gs <- data.frame(gene = mcols(gr_g)[[gene_id]],
      stringsAsFactors = FALSE)
    if (get_gene_name) {
      gs$gene_name <- mcols(gr_g)[[gene_name]]
    }
    # Add gene names to output
    out <- out %>%
      left_join(gs, by = "gene")
    if (get_gene_version) {
      gs$gv <- mcols(gr_g)[[gene_version]]
      # Add gene version to output
      # Avoid R CMD check note
      gene <- gv <- NULL
      out <- out %>%
        left_join(gs, by = c("gene", "gene_name")) %>%
        unite("gene", gene, gv, sep = version_sep)
    }
  }
  distinct(out)
}

#' Get transcript and gene info from names in FASTA files
#'
#' FASTA files, such as those for cDNA and ncRNA from Ensembl, might have genome
#' annotation information in the name of each sequence entry. This function
#' extracts the transcript and gene IDs from such FASTA files.
#'
#' At present, this function only works with FASTA files from Ensembl, and uses
#' regex to extract vertebrate Ensembl IDs. Sequence names should be formatted
#' as follows:
#'
#' ```
#' ENST00000632684.1 cdna chromosome:GRCh38:7:142786213:142786224:1
#' gene:ENSG00000282431.1 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene
#' gene_symbol:TRBD1 description:T cell receptor beta diversity 1
#' [Source:HGNC Symbol;Acc:HGNC:12158]
#' ```
#'
#' If your FASTA file sequence names are formatted differently, then you must
#' extract the transcript and gene IDs by some other means. The Bioconductor
#' package \code{Biostrings} is recommended; after reading the FASTA file into
#' R, the sequence names can be accessed by the \code{names} function.
#'
#' While normally, you should call \code{\link{sort_tr2g}} to sort the
#' transcript IDs from the output of the \code{tr2g_*} family of functions, If
#' the FASTA file supplied here is the same as the one used to build the
#' kallisto index, then the transcript IDs in the output of this function are in
#' the same order as in the kallisto index, so you can skip \code{\link{sort_tr2g}}
#' and proceed directly to \code{\link{EC2gene}} with the output of this
#' function.
#'
#' @inheritParams tr2g_ensembl
#' @param file Path to the FASTA file to be read. The file can remain gzipped.
#' @return A data frame with at least 2 columns: \code{gene} for gene ID,
#' \code{transcript} for transcript ID, and optionally \code{gene_name} for gene
#' names.
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_extract str_remove str_replace
#' @importFrom dplyr select mutate
#' @family functions to retrieve transcript and gene info
#' @export
#' @examples
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file_use <- paste(toy_path, "fasta_test.fasta", sep = "/")
#' tr2g <- tr2g_fasta(file = file_use, verbose = FALSE)
tr2g_fasta <- function(file, use_gene_name = TRUE,
                       use_transcript_version = TRUE,
                       use_gene_version = TRUE,
                       verbose = TRUE) {
  check_char1(setNames(file, "file"))
  file <- normalizePath(file, mustWork = TRUE)
  if (!str_detect(file, "(\\.fasta)|(\\.fa)|(\\.fna)")) {
    stop("file must be a FASTA file.")
  }
  file <- normalizePath(file, mustWork = TRUE)
  if (verbose) {
    message("Reading FASTA file.")
  }
  s <- readDNAStringSet(file)
  is_ens <- all(str_detect(names(s), "^ENS[A-Z]*T\\d+"))
  if (!is_ens && (use_transcript_version || use_gene_version)) {
    message("Version is not applicable to IDs not of the form ENS[species prefix][feature type prefix][a unique eleven digit number].")
    use_transcript_version <- use_gene_version <- FALSE
  }
  # Avoid R CMD check note
  g <- gene_name <- NULL
  out <- tibble(transcript = str_extract(names(s), "^[a-zA-Z\\d-\\.]+"),
    gene = str_replace(names(s), "^.*gene:", "") %>%
      str_replace("\\s+.*$", ""))
  if (use_gene_name) {
    out$gene_name <- str_replace(names(s), "^.*gene_symbol:", "") %>%
      str_replace("\\s+.*$", "")
  }
  out <- distinct(out)
  # Remove version number
  if (is_ens) {
    # Prevent R CMD check note of no visible binding for global variable
    transcript <- gene <- NULL
    if (!use_transcript_version) {
      out <- out %>%
        mutate(transcript = str_remove(transcript, "\\.\\d+$"))
    }
    if (!use_gene_version) {
      out <- out %>%
        mutate(gene = str_remove(gene, "\\.\\d+$"))
    }
  }
  out
}

#' Get transcript and gene info from TxDb objects
#'
#' The genome and gene annotations of some species can be conveniently obtained
#' from Bioconductor packages. This is more convenient than downloading GTF
#' files from Ensembl and reading it into R. In these packages, the gene
#' annotation is stored in a \code{\link{TxDb}} object, which has standardized
#' names for gene IDs, transcript IDs, exon IDs, and so on, which are stored in
#' the metadata fields in GTF and GFF3 files, which are not standardized.
#' This function extracts transcript and corresponding gene information from
#' gene annotation stored in a \code{\link{TxDb}} object.
#'
#' @param txdb A \code{\link{TxDb}} object with gene annotation.
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{tx_id} for internal transcript IDs used to avoid
#' duplicate transcript names. For TxDb packages from Bioconductor, gene ID is
#' Entrez ID, while transcript IDs are Ensembl IDs with version numbers for
#' `TxDb.Hsapiens.UCSC.hg38.knownGene`. In some cases, the transcript ID
#' have duplicates, and this is resolved by adding numbers to make the IDs
#' unique.
#' @importFrom AnnotationDbi columns keys keytypes
#' @importFrom stats complete.cases
#' @family functions to retrieve transcript and gene info
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. If \code{other_attrs}
#' has been specified, then those will also be columns in the data frame returned.
#' @family functions to retrieve transcript and gene info
#' @export
#' @examples
#' library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#' tr2g_TxDb(TxDb.Hsapiens.UCSC.hg38.knownGene)
tr2g_TxDb <- function(txdb) {
  df <- AnnotationDbi::select(txdb, AnnotationDbi::keys(txdb, keytype = "TXID"),
    keytype = "TXID",
    columns = c("TXNAME", "GENEID", "TXID"))
  if (anyDuplicated(df$TXNAME)) {
    df$TXNAME <- make.unique(df$TXNAME, sep = "_")
  }
  df <- df[complete.cases(df), c("TXNAME", "GENEID", "TXID")]
  names(df) <- c("transcript", "gene", "tx_id")
  df
}

#' Get transcript and gene info from EnsDb objects
#'
#' Bioconductor provides Ensembl genome annotation in `AnnotationHub`; older
#' versions of Ensembl annotation can be obtained from packages like
#' `EnsDb.Hsapiens.v86`. This is an alternative to querying Ensembl with
#' biomart; Ensembl's server seems to be less stable than that of Bioconductor.
#' However, more information and species are available on Ensembl biomart than
#' on `AnnotationHub`.
#'
#' @inheritParams tr2g_ensembl
#' @param ensdb Ann `EnsDb` object, such as from `AnnotationHub` or
#' `EnsDb.Hsapiens.v86`.
#' @param other_attrs Character vector. Other attributes to get from the `EnsDb`
#' object, such as gene symbol and position on the genome.
#' Use \code{\link{columns}} to see which attributes are available.
#' @return A data frame with at least 2 columns: \code{gene} for gene ID,
#' \code{transcript} for transcript ID, and optionally \code{gene_name}
#' for gene names. If \code{other_attrs} has been specified, then those will
#' also be columns in the data frame returned.
#' @family functions to retrieve transcript and gene info
#' @export
#' @examples
#' library(EnsDb.Hsapiens.v86)
#' tr2g_EnsDb(EnsDb.Hsapiens.v86, use_transcript_version = FALSE,
#'   use_gene_version = FALSE)
tr2g_EnsDb <- function(ensdb, other_attrs = NULL, use_gene_name = TRUE,
                       use_transcript_version = TRUE,
                       use_gene_version = TRUE) {
  attrs_use <- c("TXID", "GENEID", other_attrs)
  if (use_gene_name) {
    attrs_use <- c(attrs_use, "GENENAME")
  }
  if (use_transcript_version) {
    attrs_use[1] <- "TXIDVERSION"
  }
  if (use_gene_version) {
    attrs_use[2] <- "GENEIDVERSION"
  }
  df <- AnnotationDbi::select(ensdb, AnnotationDbi::keys(ensdb, keytype = "TXID"),
    keytype = "TXID",
    columns = attrs_use)
  if (use_transcript_version) {
    df$TXID <- NULL
  }
  names(df)[str_detect(names(df), "^TXID")] <- "transcript"
  names(df)[str_detect(names(df), "^GENEID")] <- "gene"
  names(df)[names(df) == "GENENAME"] <- "gene_name"
  df
}

#' Sort transcripts to the same order as in kallisto index
#'
#' This function takes the data frame output from the \code{tr2g_*} family of
#' functions in this package as the input, and sorts it so the transcripts are
#' in the same order as in the kallisto index used to generate the \code{bus}
#' file. Sorting is vital to obtain the correct sparse matrix from the \code{bus}
#' file as equivalence class notations are based on the index of transcripts
#' in the kallisto index.
#'
#' Since the attribute field of GTF and GFF3 files varies across sources, output
#' from \code{\link{tr2g_gtf}} and \code{\link{tr2g_gff3}} may need further
#' clean up. You may also supply gene and transcript IDs from other sources.
#' This function should be used after the clean up, when the transcript IDs in
#' the cleaned up data frame have the same format as those in \code{transcript}
#'
#' @param tr2g The data frame output from the \code{tr2g_*} family of functions.
#' @param file Character vector of length 1, path to a tsv file with
#' transcript IDs and the corresponding gene IDs, in the format required for
#' `bustools`, or written by \code{\link{save_tr2g_bustools}}.
#' @param kallisto_out_path Character vector of length 1, path to the directory
#' for the outputs of kallisto bus.
#' @param verbose Whether to display progress.
#' @return A data frame with columns \code{transcript} and \code{gene} and the
#' other columns present in \code{tr2g} or the data frame in \code{file}, with
#' the transcript IDs sorted to be in the same order as in the kallisto index.
#' @importFrom data.table fread fwrite
#' @export
#' @family functions to retrieve transcript and gene info
#' @examples
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file_use <- paste(toy_path, "gtf_test.gtf", sep = "/")
#' tr2g <- tr2g_gtf(file = file_use, verbose = FALSE,
#'   transcript_version = NULL)
#' tr2g <- sort_tr2g(tr2g, kallisto_out_path = toy_path, verbose = FALSE)
sort_tr2g <- function(tr2g, file, kallisto_out_path, verbose = TRUE) {
  if (!xor(missing(tr2g), missing(file))) {
    stop("Exactly one of tr2g and file should be missing.")
  }
  kallisto_out_path <- normalizePath(kallisto_out_path, mustWork = TRUE)
  trs_path <- paste(kallisto_out_path, "transcripts.txt", sep = "/")
  if (!file.exists(trs_path)) {
    stop("The file transcripts.txt does not exist in",
      kallisto_out_path, "")
  }
  if (missing(tr2g)) {
    tr2g <- fread(file, header = FALSE, col.names = c("transcript", "gene"))
  }
  trs <- fread(trs_path, header = FALSE, col.names = "transcript")
  if (verbose) {
    message("Sorting transcripts")
  }
  out <- merge(trs, tr2g, by = "transcript", sort = FALSE)
  if (nrow(trs) != nrow(out)) {
    stop("Some transcripts in the kallisto index are absent from tr2g.")
  }
  out
}

#' Save transcript to gene file for use in `bustools`
#'
#' This function saves the transcript to gene data frame generated by this package
#' in whatever means in a format required by `bustools`. In order to use
#' `bustools` to generate the gene count or TCC matrix, a file
#' that maps transcripts to genes is required. This should be a tsv file with 2
#' columns: the first column for transcript ID and the second for gene ID. The
#' order of transcripts in this file must be the same as the order in the
#' kallisto index, and this ordering can be ensured by the function
#' \code{\link{sort_tr2g}}. There must also be no headers. All columns other than
#' `transcript` and `gene` will be discarded. To save a file with those columns,
#' directly save the transcript to gene data frame with function like
#' \code{\link{write.table}}, \code{readr::write_delim}, and
#' \code{\link{fwrite}}.
#'
#' @inheritParams sort_tr2g
#' @param \dots Other arguments passed to \code{\link{fwrite}}, such
#' as \code{sep}, \code{quote}, and \code{col.names}.
#' @param file_save File name of the file to be saved. The directory in which
#' the file is to be saved must exist.
#' @return Nothing is returned into the R session. A tsv file of the format
#' required by `bustools` with the name and directory specified will be written
#' to disk.
#' @export
#' @examples
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file_use <- paste(toy_path, "gtf_test.gtf", sep = "/")
#' tr2g <- tr2g_gtf(file = file_use, verbose = FALSE)
#' save_tr2g_bustools(tr2g, file_save = "./tr2g.tsv")
save_tr2g_bustools <- function(tr2g, file_save = "./tr2g.tsv", ...) {
  file_save <- normalizePath(file_save, mustWork = FALSE)
  fwrite(tr2g[, c("transcript", "gene")], file = file_save, sep = "\t",
    col.names = FALSE)
}

#' Map Ensembl transcript ID to gene ID
#'
#' This function is a shortcut to get the correctly sorted data frame with
#' transcript IDs and the corresponding gene IDs from Ensembl biomart or Ensembl
#' transcriptome FASTA files. For biomart query, it calls
#' \code{\link{tr2g_ensembl}} and then \code{\link{sort_tr2g}}. For FASTA files,
#' it calls \code{\link{tr2g_fasta}} and then \code{\link{sort_tr2g}}. Unlike in
#' \code{\link{tr2g_ensembl}} and \code{\link{tr2g_fasta}}, multiple species can
#' be supplied if cells from different species were sequenced together. This
#' function should only be used if the kallisto inidex was built with
#' transcriptomes from Ensembl. Also, if querying biomart, please make sure to set
#' \code{ensembl_version} to match the version where the transcriptomes were
#' downloaded.
#'
#' @inheritParams tr2g_ensembl
#' @inheritParams sort_tr2g
#' @param species A character vector of Latin names of species present in this
#' scRNA-seq dataset. This is used to retrieve Ensembl information from biomart.
#' @param type A character vector indicating the type of each species. Each
#' element must be one of "vertebrate", "metazoa", "plant", "fungus", and
#' "protist". If length is 1, then this type will be used for all species specified
#' here. Can be missing if `fasta_file` is specified.
#' @param fasta_file Character vector of paths to the transcriptome FASTA files
#' used to build the kallisto index. Exactly one of \code{species} and
#' \code{fasta_file} can be missing.
#' @param kallisto_out_path Path to the \code{kallisto bus} output directory.
#' @param verbose Whether to display progress. Defaults to \code{TRUE}.
#' @return A data frame with two columns: \code{gene} and \code{transcript},
#' with Ensembl gene and transcript IDs (with version number), in the same order
#' as in the transcriptome index used in \code{kallisto}.
#' @param \dots Other arguments passed to `tr2g_ensembl` such as `other_attrs`,
#' `ensembl_version`, and arguments passed to \code{\link{useEnsembl}}. If
#' `fasta_files` is supplied instead of `species`, then this will be extra
#' argumennts to \code{\link{tr2g_fasta}}, such as `use_transcript_version` and
#' `use_gene_version`.
#' @importFrom data.table rbindlist
#' @export
#' @family functions to retrieve transcript and gene info
#' @examples
#' # Download dataset already in BUS format
#' library(TENxBUSData)
#' TENxBUSData(".", dataset = "retina")
#' tr2g <- transcript2gene("Mus musculus", type = "vertebrate",
#'   ensembl_version = 94, kallisto_out_path = "./out_retina")
transcript2gene <- function(species, fasta_file, kallisto_out_path,
                            type = "vertebrate",
                            verbose = TRUE, ...) {
  if (!xor(missing(species), missing(fasta_file))) {
    stop("Exactly one of species and fasta_file can be missing.")
  }
  if (missing(fasta_file)) {
    if (length(type) != 1 && length(species) != length(type)) {
      stop("species and type must have the same length.")
    }
    if (length(type) == 1) {
      type <- rep(type, length(species))
    }
    kallisto_out_path <- normalizePath(kallisto_out_path, mustWork = TRUE)
    MoreArgs <- list(...)
    fls <- mapply(tr2g_ensembl, species, type,
      verbose = verbose,
      MoreArgs = MoreArgs,
      SIMPLIFY = FALSE)
    tr2g <- rbindlist(fls)
    return(sort_tr2g(tr2g, kallisto_out_path = kallisto_out_path, verbose = verbose))
  } else {
    fls <- lapply(fasta_file, tr2g_fasta, verbose = verbose, ...)
    tr2g <- rbindlist(fls)
    # Just to be safe, to make sure that the transcripts are in the right order
    return(sort_tr2g(tr2g, kallisto_out_path = kallisto_out_path,
      verbose = verbose))
  }
}

#' Map EC Index to Genes Compatible with the EC
#'
#' In the output file \code{output.bus}, equivalence classes (EC) are denoted by
#' an index, which is related to the set of transcripts the EC is compatible to
#' in the output file \code{matrix.ec}. This function further relates the set of
#' transcripts to the set of genes the EC is compatible to. This function first
#' reads in \code{matrix.ec}, and then translates the transcripts into genes.
#'
#' The data frame passed to \code{tr2g} can be generated from function
#' \code{\link{transcript2gene}} in this package for any organism that has gene and
#' transcript ID on Ensembl, or from the \code{tr2g_*} family of function.
#' You no longer need to use this function before running \code{make_sparse_matrix};
#' the purpose of this function is to query which genes equivalence classes map
#' to.
#'
#' Calling this function is unnessary when working with gene count matrices.
#' However, this function is useful for finding genes the ECs map to in TCC
#' matrices, such as when finding species-specific ECs in mixed species datasets
#' and identifying ECs mapped to known marker genes of cell types.
#'
#' @inheritParams transcript2gene
#' @param tr2g A Data frame with columns \code{gene} and \code{transcript}, in
#' the same order as in the transcriptome index for \code{kallisto}.
#' @param ncores Number of cores to use, defaults to 0, which means the system
#' will automatically determine the number of cores as it sees fit. Negative
#' numbers are interpreted as 0. Positive numbers will limit the number of cores
#' used. This might not speed up `EC2gene` very much unless there are many genes
#' or ECs detected.
#' @return A data frame with 3 columns:
#' \describe{
#' \item{EC_ind}{Index of the EC as appearing in the `matrix.ec` file.}
#' \item{EC}{A list column each element of which is a numeric vector of the
#' transcripts in the EC corresponding to the EC index. To learn more about list
#' columns, see the [relevant section in the R for Data Science book](https://r4ds.had.co.nz/many-models.html#list-columns-1).}
#' \item{gene}{A list column each element of which is a character vector of genes
#' the EC maps to.}
#' }
#' @seealso \code{\link{transcript2gene}}
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom tibble tibble
#' @export
#' @examples
#' # Load toy example for testing
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' load(paste(toy_path, "toy_example.RData", sep = "/"))
#' EC2gene(tr2g_toy, toy_path, verbose = FALSE, ncores = 1)
EC2gene <- function(tr2g, kallisto_out_path, ncores = 0, verbose = TRUE) {
  kallisto_out_path <- normalizePath(kallisto_out_path, mustWork = TRUE)
  c(ec_vec, genes) %<-% EC2gene_export(tr2g, kallisto_out_path, ncores, verbose)
  # Sort according to indices
  EC_inds <- 0:(length(genes) - 1)
  genes <- genes[as.character(EC_inds)]
  names(genes) <- NULL
  ec_vec <- ec_vec[as.character(EC_inds)]
  names(ec_vec) <- NULL
  ec_vec <- lapply(ec_vec, as.numeric)
  tibble(EC_ind = EC_inds,
    EC = ec_vec,
    gene = genes)
}
