#' @include sparse_matrix.R
NULL

write_tr2g_fun <- function(out, out_path, overwrite) {
  fn <- paste(out_path, "tr2g.tsv", sep = "/")
  if (file.exists(fn) && !overwrite) {
    message("File ", fn, " already exists.")
  } else {
    write.table(out, file = fn, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
}

check_out_path <- function(out_path) {
  out_path <- normalizePath(out_path, mustWork = FALSE)
  if (!dir.exists(out_path)) dir.create(out_path)
  out_path
}

#' Get transcript and gene info from Ensembl
#'
#' This function queries Ensembl biomart to convert transcript IDs to gene IDs.
#' 
#' @inheritParams tr2g_GRanges
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
#' kallisto index was downloaded. This only works for vertebrates and the most
#' common invertebrate model organisms like _Drosophila melanogaster_ and 
#' _C. elegans_ (i.e. www.ensembl.org and its mirrors), not the other Ensembl
#' sites for plants, protists, fungi, and metazoa. 
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
#' @seealso dl_transcriptome
#' @importFrom BiocGenerics %in%
#' @importFrom GenomeInfoDb genomeStyles
#' @export
#' @examples
#' tr2g <- tr2g_ensembl(species = "Felis catus", other_attrs = "description",
#'   write_tr2g = FALSE)
#' # This will use plants.ensembl.org as host instead of www.ensembl.org
#' tr2g <- tr2g_ensembl(species = "Arabidopsis thaliana", type = "plant",
#'   write_tr2g = FALSE)
tr2g_ensembl <- function(species, type = c("vertebrate", "metazoa", "plant",
                           "fungus", "protist"), out_path = ".",
                         write_tr2g = TRUE,
                         other_attrs = NULL,
                         use_gene_name = TRUE,
                         use_transcript_version = TRUE,
                         use_gene_version = TRUE,
                         transcript_biotype_col = "transcript_biotype",
                         gene_biotype_col = "gene_biotype", 
                         transcript_biotype_use = "all",
                         gene_biotype_use = "all", 
                         chrs_only = TRUE,
                         ensembl_version = NULL, overwrite = FALSE,
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
  if (type != "vertebrate" & !is.null(ensembl_version)) {
    warning("Archive only works for vertebrates. Using current version instead.")
    ensembl_version <- NULL
  }
  if (write_tr2g) {
    out_path <- check_out_path(out_path)
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
  if (transcript_biotype_use != "all") {
    attrs_use <- c(attrs_use, transcript_biotype_col)
  }
  if (gene_biotype_use != "all") {
    attrs_use <- c(attrs_use, gene_biotype_col)
  }
  if (chrs_only) {
    chrs_use <- try(genomeStyles(species)$Ensembl)
    if (is(chrs_use, "try-error")) {
      chrs_only <- FALSE
    } else {
      attrs_use <- c(attrs_use, "chromosome_name")
    }
  }
  attrs_use <- unique(attrs_use)
  out <- getBM(attrs_use, mart = mart)
  
  if (gene_biotype_use != "all") {
    gbt_use <- which_biotypes(gene_biotype_use, out[[gene_biotype_col]])
    out <- out[out[[gene_biotype_col]] %in% gbt_use,]
  }
  if (transcript_biotype_use != "all") {
    tbt_use <- which_biotypes(transcript_biotype_use, 
                              out[[transcript_biotype_col]])
    out <- out[out[[transcript_biotype_col]] %in% tbt_use,]
  }
  if (chrs_only) {
    out <- out[out$chromosome_name %in% chrs_use,]
  }
  names(out)[seq_len(2)] <- c("transcript", "gene")
  names(out)[names(out) == "external_gene_name"] <- "gene_name"
  if (write_tr2g) {
    write_tr2g_fun(out, out_path, overwrite)
  }
  out
}

which_biotypes <- function(bt, bt_col) {
  u_bt_col <- unique(bt_col)
  u_bt_col <- u_bt_col[!is.na(u_bt_col)]
  if (length(bt) == 1L) {
    if (bt == "all") {
      return(u_bt_col)
    } else if (bt == "cellranger") {
      bt <- c("protein_coding", "lincRNA", "antisense", "IG_LV_gene", 
             "IG_V_gene", "IG_V_pseudogene", "IG_D_gene", "IG_J_gene",
             "IG_J_pseudogene", "IG_C_gene", "IG_C_pseudogene", "TR_V_gene", 
             "TR_V_pseudogene", "TR_D_gene", "TR_J_gene", "TR_J_pseudogene", 
             "TR_C_gene")
    }
  }
  out <- intersect(bt, u_bt_col)
  if (length(out) == 0L) {
    stop("Biotype must be one of ", paste(u_bt_col, collapse = ", "))
  }
  out
}

filter_chr <- function(gr, chrs_only) {
  if (chrs_only) {
    gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  }
  gr
}

filter_biotype <- function(gr, transcript_biotype_col,
                           gene_biotype_col, transcript_biotype_use,
                           gene_biotype_use, gff3 = FALSE) {
  if (gene_biotype_use == "all" && transcript_biotype_use == "all") {
    return(gr)
  }
  if (gene_biotype_use != "all") {
    gbt_use <- which_biotypes(gene_biotype_use, mcols(gr)[[gene_biotype_col]])
    gr <- gr[mcols(gr)[[gene_biotype_col]] %in% gbt_use]
  }
  if (transcript_biotype_use != "all") {
    if (transcript_biotype_use == "cellranger") {
      warning("Transcript biotypes are not the same as gene biotypes. ",
              "Cell Ranger reference filters by gene biotypes.")
    }
    tbt_use <- which_biotypes(transcript_biotype_use, mcols(gr)[[transcript_biotype_col]])
    gr <- gr[mcols(gr)[[transcript_biotype_col]] %in% tbt_use]
  }
  gr
}

filter_biotype_gff3 <- function(gr, transcript_id, gene_id, transcript_biotype_col,
                                gene_biotype_col, transcript_biotype_use,
                                gene_biotype_use, source) {
  grt <- gr[!is.na(mcols(gr)[[transcript_id]])]
  if (gene_biotype_use == "all" && transcript_biotype_use == "all") {
    grg <- gr[gr$type == "gene"]
    return(list(gr_tx = grt, gr_g = grg))
  }
  sep <- if (source == "ensembl") ":" else "-"
  if (gene_biotype_use != "all") {
    grg <- gr[gr$type == "gene"]
    gbt_use <- which_biotypes(gene_biotype_use, mcols(grg)[[gene_biotype_col]])
    genes_use <- unique(mcols(grg)[[gene_id]][mcols(grg)[[gene_biotype_col]] %in% gbt_use])
    genes_use <- genes_use[!is.na(genes_use)]
    if (transcript_biotype_use == "all") {
      genes <- paste("gene", genes_use, sep = sep)
      return(list(gr_tx = grt[grt$Parent %in% genes], gr_g = grg))
    }
  }
  if (transcript_biotype_use != "all") {
    if (transcript_biotype_use == "cellranger") {
      warning("Transcript biotypes are not the same as gene biotypes. ",
              "Cell Ranger reference filters by gene biotypes.")
    }
    tbt_use <- which_biotypes(transcript_biotype_use, mcols(grt)[[transcript_biotype_col]])
    grt <- grt[mcols(grt)[[transcript_biotype_col]] %in% tbt_use]
    if (gene_biotype_use == "all") {
      grg <- gr[gr$type == "gene"]
      return(list(gr_tx = grt, gr_g = grg))
    } else {
      genes <- paste("gene", genes_use, sep = sep)
      return(list(gr_tx = grt[grt$Parent %in% genes], gr_g = grg))
    }
  }
}

check_genome_present <- function(Genome, get_transcriptome) {
  if (get_transcriptome) {
    if (is.null(Genome)) {
      stop("Genome is required for transcriptome extraction.")
    } else if (!is(Genome, "BSgenome") && !is(Genome, "DNAStringSet")) {
      stop("Genome must be either a BSgenome or a DNAStringSet.")
    }
  }
}

#' Get transcript and gene info from GRanges
#'
#' Internal use, for GRanges from GTF files
#'
#' @param gr A \code{\link{GRanges}} object. The metadata columns should be
#' atomic vectors, not lists.
#' @param get_transcriptome Logical, whether to extract transcriptome from
#' genome with the GTF file. If filtering biotypes or chromosomes, the filtered
#' `GRanges` will be used to extract transcriptome.
#' @param out_path Directory to save the outputs written to disk. If this
#' directory does not exist, then it will be created. Defaults to the current
#' working directory.
#' @param write_tr2g Logical, whether to write tr2g to disk. If `TRUE`, then
#' a file `tr2g.tsv` will be written into `out_path`.
#' @param Genome Either a \code{\link{BSgenome}} or a \code{\link{XStringSet}}
#' object of genomic sequences, where the intronic sequences will be extracted
#' from. Use \code{\link{genomeStyles}} to check which styles are supported for
#' your organism of interest; supported styles can be interconverted. If the
#' style in your genome or annotation is not supported, then the style of
#' chromosome names in the genome and annotation should be manually set to be
#' consistent.
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
#' error if tag indicated in this argument does not exist. This is typically
#' "gene_id" for annotations from Ensembl and "gene" for refseq.
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
#' @param transcript_biotype_col Character vector of length 1. Tag in 
#' \code{attribute} field corresponding to _transcript_ biotype. 
#' @param gene_biotype_col Character vector of length 1. Tag in \code{attribute}
#' field corresponding to _gene_ biotype. 
#' @param transcript_biotype_use Character, can be "all" or
#' a vector of _transcript_ biotypes to be used. Transcript biotypes aren't
#' entirely the same as gene biotypes. For instance, in Ensembl annotation,
#' `retained_intron` is a transcript biotype, but not a gene biotype. If 
#' "cellranger", then a warning will be given. See `data("ensembl_tx_biotypes")`
#' for all available transcript biotypes from Ensembl.
#' @param gene_biotype_use Character, can be "all", "cellranger", or
#' a vector of _gene_ biotypes to be used. If "cellranger", then the biotypes 
#' used by Cell Ranger's reference are used. See `data("cellranger_biotypes")` 
#' for gene biotypes the Cell Ranger reference uses. See 
#' `data("ensembl_gene_biotypes")` for all available gene biotypes from Ensembl.
#' Note that gene biotypes and transcript biotypes are not always the same.
#' @param chrs_only Logical, whether to include chromosomes only, for GTF and
#' GFF files can contain annotations for scaffolds, which are not incorporated
#' into chromosomes. This will also exclude haplotypes. Defaults to `TRUE`. 
#' Only applicable to species found in `genomeStyles()`.
#' @param compress_fa Logical, whether to compress the output fasta file. If 
#' `TRUE`, then the fasta file will be gzipped.
#' @param save_filtered_gtf Logical. If filtering type, biotypes, and/or 
#' chromosomes, whether to save the filtered `GRanges` as a GTF file.
#' @param overwrite Logical, whether to overwrite if files with names of outputs
#' written to disk already exist.
#' @return A data frame at least 2 columns: \code{gene} for gene ID,
#' \code{transcript} for transcript ID, and optionally, \code{gene_name} for
#' gene names.
#' @importFrom magrittr %>%
#' @importFrom stringr str_detect
#' @importFrom dplyr distinct
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb mapSeqlevels
#' @importFrom plyranges write_gff
#' @importFrom GenomeInfoDb keepStandardChromosomes
tr2g_GRanges <- function(gr, Genome = NULL, get_transcriptome = TRUE, 
                         out_path = ".", write_tr2g = TRUE, 
                         transcript_id = "transcript_id",
                         gene_id = "gene_id", gene_name = "gene_name",
                         transcript_version = "transcript_version",
                         gene_version = "gene_version", version_sep = ".",
                         transcript_biotype_col = "transcript_biotype",
                         gene_biotype_col = "gene_biotype", 
                         transcript_biotype_use = "all",
                         gene_biotype_use = "all", 
                         chrs_only = TRUE, compress_fa = FALSE,
                         save_filtered_gtf = TRUE, overwrite = FALSE) {
  tags <- names(mcols(gr))
  check_tag_present(c(transcript_id, gene_id), tags, error = TRUE)
  if (transcript_biotype_use != "all") {
    check_tag_present(transcript_biotype_col, tags, error = TRUE)
  }
  if (gene_biotype_use != "all") {
    check_tag_present(gene_biotype_col, tags, error = TRUE)
  }
  # Will do nothing if all are NULL
  check_tag_present(c(gene_name, transcript_version, gene_version),
    tags, error = FALSE)
  if (any(get_transcriptome, write_tr2g, save_filtered_gtf)) {
    out_path <- check_out_path(out_path)
  }
  check_genome_present(Genome = Genome, get_transcriptome = get_transcriptome)
  gr <- gr[!is.na(mcols(gr)[[transcript_id]])]
  gr <- gr[gr$type == "exon"]
  if (length(gr) == 0) {
    stop(paste("No entry has type exon."))
  }
  # Filter chromosomes
  gr <- filter_chr(gr, chrs_only)
  # Filter by biotypes
  gr <- filter_biotype(gr, transcript_biotype_col, gene_biotype_col, 
                       transcript_biotype_use, gene_biotype_use)
  if (get_transcriptome) {
    gr <- sort(gr)
  }
  # Prepare output
  out <- tibble(transcript = mcols(gr)[[transcript_id]],
                gene = mcols(gr)[[gene_id]])
  if (!is.null(gene_name) && gene_name %in% tags) {
    out$gene_name <- mcols(gr)[[gene_name]]
  }
  if (!is.null(transcript_version) && transcript_version %in% tags) {
    tv <- mcols(gr)[[transcript_version]]
    mcols(gr)[[transcript_id]] <- out$transcript <- 
      paste(out$transcript, tv, sep = version_sep)
    mcols(gr)[[transcript_version]] <- NULL
  }
  if (!is.null(gene_version) && gene_version %in% tags) {
    gv <- mcols(gr)[[gene_version]]
    mcols(gr)[[gene_id]] <- out$gene <- 
      paste(out$gene, gv, sep = version_sep)
    mcols(gr)[[gene_version]] <- NULL
  }
  out <- distinct(out)
  do_filter <- gene_biotype_use != "all" | transcript_biotype_col != "all" |
    chrs_only
  if (save_filtered_gtf && do_filter) {
    gtf_save <- paste(out_path, "gtf_filtered.gtf", sep = "/")
    if (file.exists(gtf_save) && !overwrite) {
      message("File ", gtf_save, " already exists.")
    } else write_gff(gr, gtf_save)
  }
  if (get_transcriptome) {
    tx_save <- paste(out_path, "transcriptome.fa", sep = "/")
    if (compress_fa) tx_save <- paste0(tx_save, ".gz")
    if (file.exists(tx_save) && !overwrite) {
      message("File ", tx_save, " already exists.")
    } else {
      c(Genome, gr) %<-% match_style(Genome, gr, style = "annotation")
      gr <- subset_annot(Genome, gr)
      c(Genome, gr) %<-% annot_circular(Genome, gr)
      genome(gr) <- genome(Genome)[seqlevels(gr)]
      grl <- GenomicRanges::split(gr, gr$transcript_id)
      grl <- revElements(grl, any(strand(grl) == "-"))
      tx <- extractTranscriptSeqs(Genome, grl)
      out <- out[out$transcript %in% names(tx),]
      tx <- tx[out$transcript]
      writeXStringSet(tx, tx_save, compress = compress_fa)
    }
  }
  if (write_tr2g) {
    write_tr2g_fun(out, out_path, overwrite)
  }
  out
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
#' @param file Path to a GTF file to be read. The file can remain gzipped. Use
#' \code{getGTF} from the \code{biomartr} package to download GTF files
#' from Ensembl, and use \code{getGFF} from \code{biomartr} to download
#' GFF3 files from Ensembl and RefSeq.
#' @inheritParams tr2g_GRanges
#' @inheritParams tr2g_ensembl
#' @return A data frame at least 2 columns: \code{gene} for gene ID,
#' \code{transcript} for transcript ID, and optionally, \code{gene_name} for
#' gene names.
#' @importFrom plyranges read_gff
#' @family functions to retrieve transcript and gene info
#' @seealso ensembl_gene_biotypes ensembl_tx_biotypes cellranger_biotypes
#' @export
#' @examples
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file_use <- paste(toy_path, "gtf_test.gtf", sep = "/")
#' # Default
#' tr2g <- tr2g_gtf(file = file_use, get_transcriptome = FALSE,
#'   write_tr2g = FALSE, save_filtered_gtf = FALSE)
#' # Excluding version numbers
#' tr2g <- tr2g_gtf(file = file_use, transcript_version = NULL,
#'   gene_version = NULL, get_transcriptome = FALSE,
#'   write_tr2g = FALSE, save_filtered_gtf = FALSE)
tr2g_gtf <- function(file, Genome = NULL, get_transcriptome = TRUE, 
                     out_path = ".", write_tr2g = TRUE, 
                     transcript_id = "transcript_id",
                     gene_id = "gene_id", gene_name = "gene_name",
                     transcript_version = "transcript_version",
                     gene_version = "gene_version", version_sep = ".",
                     transcript_biotype_col = "transcript_biotype",
                     gene_biotype_col = "gene_biotype", 
                     transcript_biotype_use = "all",
                     gene_biotype_use = "all", 
                     chrs_only = TRUE, compress_fa = FALSE,
                     save_filtered_gtf = TRUE, overwrite = FALSE) {
  # Validate arguments
  check_char1(setNames(file, "file"))
  file <- normalizePath(file, mustWork = TRUE)
  check_gff("gtf", file, transcript_id, gene_id)
  gr <- read_gff(file)
  tr2g_GRanges(gr = gr, Genome = Genome, get_transcriptome = get_transcriptome, 
               out_path = out_path, write_tr2g = write_tr2g, 
               transcript_id = transcript_id,
               gene_id = gene_id, gene_name = gene_name,
               transcript_version = transcript_version,
               gene_version = gene_version, version_sep = version_sep,
               transcript_biotype_col = transcript_biotype_col,
               gene_biotype_col = gene_biotype_col, 
               transcript_biotype_use = transcript_biotype_use,
               gene_biotype_use = gene_biotype_use, 
               chrs_only = chrs_only, compress_fa = compress_fa,
               save_filtered_gtf = save_filtered_gtf,
               overwrite = overwrite)
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
#' @param source Name of the database where this GFF3 file was downloaded. Must
#' be either "ensembl" or "refseq".
#' @param save_filtered_gff Logical. If filtering type, biotypes, and/or 
#' chromosomes, whether to save the filtered `GRanges` as a GFF3 file.
#' @return A data frame at least 2 columns: \code{gene} for gene ID,
#' \code{transcript} for transcript ID, and optionally, \code{gene_name} for
#' gene names.
#' @family functions to retrieve transcript and gene info
#' @importFrom plyranges read_gff3
#' @importFrom stringr str_split
#' @importFrom dplyr left_join distinct
#' @importFrom tidyr unite
#' @importFrom plyranges write_gff3
#' @export
#' @seealso ensembl_gene_biotypes ensembl_tx_biotypes cellranger_biotypes
#' ensembl_gtf_mcols ensembl_gff_mcols refseq_gff_mcols
#' @note The defaults here are for Ensembl GFF3 files. To see all attribute
#' tags for Ensembl and RefSeq GFF3 files, see `data("ensembl_gff_mcols")` and
#' `data("refseq_gff_mcols")`.
#' @examples
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file_use <- paste(toy_path, "gff3_test.gff3", sep = "/")
#' # Default
#' tr2g <- tr2g_gff3(file = file_use, write_tr2g = FALSE, 
#' get_transcriptome = FALSE, save_filtered_gff = FALSE)
#' # Excluding version numbers
#' tr2g <- tr2g_gff3(file = file_use, transcript_version = NULL,
#'   gene_version = NULL, write_tr2g = FALSE, get_transcriptome = FALSE,
#'   save_filtered_gff = FALSE)
tr2g_gff3 <- function(file, Genome = NULL, get_transcriptome = TRUE,
                      out_path = ".", write_tr2g = TRUE, 
                      transcript_id = "transcript_id",
                      gene_id = "gene_id", gene_name = "Name",
                      transcript_version = "version",
                      gene_version = "version", version_sep = ".",
                      transcript_biotype_col = "biotype",
                      gene_biotype_col = "biotype", 
                      transcript_biotype_use = "all",
                      gene_biotype_use = "all", 
                      chrs_only = TRUE, compress_fa = FALSE,
                      save_filtered_gff = TRUE, overwrite = FALSE,
                      source = c("ensembl", "refseq")) {
  # Validate arguments
  check_char1(setNames(file, "file"))
  file <- normalizePath(file, mustWork = TRUE)
  check_gff("gff3", file, transcript_id, gene_id)
  source <- match.arg(source)
  check_genome_present(Genome = Genome, get_transcriptome = get_transcriptome)
  gr <- read_gff3(file)
  if (source == "refseq") {
    # refseq has its own seqnames; use chromosome column instead
    chr <- sn <- NULL
    sn2chr <- tibble(sn = as.vector(seqnames(gr)), chr = gr$chromosome) %>% 
      dplyr::filter(!is.na(chr)) %>% 
      distinct() %>% 
      dplyr::filter(str_detect(sn, "^NC"))
    gr <- gr[seqnames(gr) %in% sn2chr$sn]
    seqlevels(gr, pruning.mode = "coarse") <- seqlevelsInUse(gr)
    seqlevels(gr) <- sn2chr$chr[as.vector(match(seqlevels(gr), sn2chr$sn))]
    transcript_version <- gene_version <- NULL
  }
  tags <- names(mcols(gr))
  check_tag_present(c(transcript_id, gene_id), tags, error = TRUE)
  # Will do nothing if all are NULL
  check_tag_present(c(gene_name, transcript_version, gene_version),
    tags, error = FALSE)
  if (get_transcriptome && !"exon" %in% unique(gr$type)) {
    stop("The type 'exon' must be present to extract transcriptome.")
  }
  if (write_tr2g || save_filtered_gff || get_transcriptome) {
    out_path <- check_out_path(out_path)
  }
  # Filter by chromosome
  gr <- filter_chr(gr, chrs_only)
  # Filter by biotype
  if (is(gr$Parent, "List")) {
    if (!any(lengths(gr$Parent) > 1)) {
      ind3 <- lengths(gr$Parent) < 1
      gr$Parent[ind3] <- NA
      gr$Parent <- unlist(gr$Parent)
    } else {
      if (get_transcriptome) {
        ind2 <- gr$type == "exon"
        gr$Parent[ind2] <- 
          lapply(gr$Parent[ind2], 
                 function(x) x[str_detect(x, "^(transcript)|(rna)")])
      }
      ind1 <- !is.na(mcols(gr)[[transcript_id]])
      gr$Parent[ind1] <- lapply(gr$Parent[ind1], 
                          function(x) x[str_detect(x, "^gene")])
      ind3 <- lengths(gr$Parent) < 1
      gr$Parent[ind3] <- NA
      gr$Parent <- unlist(gr$Parent)
    }
  }
  gr_tx <- gr_g <- NULL
  c(gr_tx, gr_g) %<-% filter_biotype_gff3(gr, transcript_id, gene_id, 
                                          transcript_biotype_col, 
                                          gene_biotype_col, transcript_biotype_use, 
                                          gene_biotype_use, source)
  # Get transcript ID
  genes <- str_split(gr_tx$Parent, ":|-", simplify = TRUE)[, 2]
  out <- tibble(transcript = mcols(gr_tx)[[transcript_id]],
                gene = genes)
  if (!is.null(transcript_version) && transcript_version %in% tags) {
    tv <- mcols(gr_tx)[[transcript_version]]
    out$transcript_version <- paste(out$transcript, tv, sep = version_sep)
  }
  # Get gene name and version
  get_gene_name <- !is.null(gene_name) && gene_name %in% tags
  get_gene_version <- !is.null(gene_version) && gene_version %in% tags
  if (get_gene_name || get_gene_version) {
    gs <- tibble(gene = mcols(gr_g)[[gene_id]])
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
        left_join(gs, by = c("gene", "gene_name"))
      out$gene_version <- paste(out$gene, out$gv, sep = version_sep)
      out$gv <- NULL
    }
  }
  out <- distinct(out)
  do_filter <- gene_biotype_use != "all" | transcript_biotype_col != "all" |
    chrs_only
  if ((save_filtered_gff || get_transcriptome)) {
    # The exons
    gre <- gr[gr$type == "exon"]
    txs <- if (source == "ensembl") {
      paste0("transcript:", out$transcript)
    } else {
      paste0("rna-", out$transcript)
    }
    gre <- gre[gre$Parent %in% txs]
  }
  if (save_filtered_gff && do_filter) {
    file_save <- paste(out_path, "gff_filtered.gff3", sep = "/")
    if (file.exists(file_save) && !overwrite) {
      message("File ", file_save, " already exists.")
    } else {
      gr_g <- gr_g[mcols(gr_g)[[gene_id]] %in% out$gene]
      gr_out <- c(gr_g, gr_tx, gre)
      write_gff3(gr_out, file_save)
    }
  }
  if (get_transcriptome) {
    tx_save <- paste(out_path, "transcriptome.fa", sep = "/")
    if (compress_fa) tx_save <- paste0(tx_save, ".gz")
    if (file.exists(tx_save) && !overwrite) {
      message("File ", tx_save, " already exists.")
    } else {
    gre <- sort(gre)
    c(Genome, gre) %<-% match_style(Genome, gre, style = "annotation")
    gre <- subset_annot(Genome, gre)
    c(Genome, gre) %<-% annot_circular(Genome, gre)
    genome(gre) <- genome(Genome)[seqlevels(gre)]
    grl <- split(gre, gre$Parent)
    names(grl) <- str_remove(names(grl), "(transcript:)|(rna-)")
    grl <- revElements(grl, any(strand(grl) == "-"))
    # Transcript version number
    if (!is.null(transcript_version)) {
      names(grl) <- out$transcript_version[match(names(grl), out$transcript)]
      out$transcript <- NULL
      names(out)[names(out) == "transcript_version"] <- "transcript"
    }
    tx <- extractTranscriptSeqs(Genome, grl)
    out <- out[out$transcript %in% names(tx),]
    tx <- tx[out$transcript]
    writeXStringSet(tx, tx_save, compress = compress_fa)
    }
  }
  if ("transcript_version" %in% names(out)) {
    out$transcript <- NULL
    names(out)[names(out) == "transcript_version"] <- "transcript"
  }
  if ("gene_version" %in% names(out)) {
    out$gene <- NULL
    names(out)[names(out) == "gene_version"] <- "gene"
  }
  out <- out[, c("transcript", "gene", setdiff(names(out), 
                                               c("transcript", "gene")))]
  if (write_tr2g) { 
    write_tr2g_fun(out, out_path, overwrite)
  }
  out
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
#' @inheritParams tr2g_GRanges
#' @inheritParams annots_from_fa_df
#' @param save_filtered If filtering for biotype and chromosomes, whether to
#' save the filtered fasta file. If `TRUE`, the file will be `tx_filtered.fa` in
#' `out_path`.
#' @return A data frame with at least 2 columns: \code{gene} for gene ID,
#' \code{transcript} for transcript ID, and optionally \code{gene_name} for gene
#' names.
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_extract str_remove str_replace
#' @importFrom dplyr select mutate
#' @family functions to retrieve transcript and gene info
#' @export
#' @seealso ensembl_gene_biotypes ensembl_tx_biotypes cellranger_biotypes
#' @examples
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file_use <- paste(toy_path, "fasta_test.fasta", sep = "/")
#' tr2g <- tr2g_fasta(file = file_use, save_filtered = FALSE, write_tr2g = FALSE)
tr2g_fasta <- function(file, out_path = ".", write_tr2g = TRUE,
                       use_gene_name = TRUE, use_transcript_version = TRUE,
                       use_gene_version = TRUE,
                       transcript_biotype_use = "all",
                       gene_biotype_use = "all", 
                       chrs_only = TRUE, save_filtered = TRUE, 
                       compress_fa = FALSE, overwrite = FALSE) {
  check_char1(setNames(file, "file"))
  file <- normalizePath(file, mustWork = TRUE)
  if (!str_detect(file, "(\\.fasta)|(\\.fa)|(\\.fna)")) {
    stop("file must be a FASTA file.")
  }
  if (write_tr2g || save_filtered) {
    out_path <- check_out_path(out_path)
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
                gene = str_extract(names(s), "(?<=gene:).*?(?=\\s)"))
  if (use_gene_name) {
    out$gene_name <- str_extract(names(s), "(?<=gene_symbol:).*?(?=(\\s|$))")
  }
  
  inds <- TRUE
  do_filter <- transcript_biotype_use != "all" | gene_biotype_use != "all" |
    chrs_only
  if (chrs_only) {
    sns <- str_extract(names(s), "(?<=((chromosome)|(scaffold)):).*?(?=\\s)") %>% 
      str_split(pattern = ":", simplify = TRUE)
    sns <- sns[,2]
    inds <- !is.na(mapSeqlevels(sns, style = "Ensembl"))
  }
  if (transcript_biotype_use != "all") {
    tbts <- str_extract(names(s), "(?<=transcript_biotype:).*?(?=\\s)")
    tbts_use <- which_biotypes(transcript_biotype_use, tbts)
    inds <- inds & (tbts %in% tbts_use)
  }
  if (gene_biotype_use != "all") {
    gbts <- str_extract(names(s), "(?<=gene_biotype:).*?(?=\\s)")
    gbts_use <- which_biotypes(gene_biotype_use, gbts)
    inds <- inds & (gbts %in% gbts_use)
  }
  if (do_filter) out <- distinct(out[inds, ])
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
  if (save_filtered & do_filter) {
    file_save <- paste(out_path, "tx_filtered.fa", sep = "/")
    fa_out <- s[inds]
    if (compress_fa) file_save <- paste0(file_save, ".gz")
    if (file.exists(file_save) && !overwrite) {
      message("File ", file_save, " already exists.")
    } else writeXStringSet(fa_out, file_save, compress = compress_fa)
  }
  if (write_tr2g) {
    write_tr2g_fun(out, out_path, overwrite)
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
#' @inheritParams tr2g_GRanges
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
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' tr2g_TxDb(TxDb.Hsapiens.UCSC.hg38.knownGene, BSgenome.Hsapiens.UCSC.hg38)
tr2g_TxDb <- function(txdb, Genome = NULL, get_transcriptome = TRUE, 
                      out_path = ".", write_tr2g = TRUE, chrs_only = TRUE, 
                      compress_fa = FALSE, overwrite = FALSE) {
  check_genome_present(Genome = Genome, get_transcriptome = get_transcriptome)
  if (get_transcriptome || write_tr2g) {
    out_path <- check_out_path(out_path)
  }
  attrs_use <- c("TXNAME", "GENEID", "TXID")
  if (chrs_only) {
    chrs_use <- mapSeqlevels(seqlevels(txdb), seqlevelsStyle(txdb)[1])
    chrs_use <- unname(chrs_use[!is.na(chrs_use)])
    attrs_use <- c(attrs_use, "TXCHROM")
    txdb <- keepStandardChromosomes(txdb, pruning.mode = "coarse")
  }
  df <- AnnotationDbi::select(txdb, AnnotationDbi::keys(txdb, keytype = "TXID"),
    keytype = "TXID",
    columns = attrs_use)
  if (anyDuplicated(df$TXNAME)) {
    df$TXNAME <- make.unique(df$TXNAME, sep = "_")
  }
  df <- df[complete.cases(df), ]
  if (chrs_only) {
    df <- df[df$TXCHROM %in% chrs_use, ]
  }
  names(df)[1:3] <- c("tx_id", "gene", "transcript")
  if (chrs_only) {
    names(df)[4] <- "seqnames"
  }
  if (get_transcriptome) {
    tx_save <- paste0(out_path, "/transcriptome.fa")
    if (compress_fa) tx_save <- paste0(tx_save, ".gz")
    if (file.exists(tx_save) && !overwrite) {
      message("File ", tx_save, " already exists.")
    } else {
      c(Genome, txdb) %<-% match_style(Genome, txdb, "annotation")
      txdb <- subset_annot(Genome, txdb)
      grl <- exonsBy(txdb, by = "tx") # Will be numbers
      # Remove transcripts that aren't mapped to genes
      grl <- grl[names(grl) %in% df$tx_id]
      names(grl) <- df$transcript[match(names(grl), df$tx_id)]
      c(Genome, grl) %<-% annot_circular(Genome, grl)
      genome(grl) <- genome(Genome)[seqlevels(grl)]
      tx <- extractTranscriptSeqs(Genome, grl)
      df <- df[df$transcript %in% names(tx),]
      tx <- tx[df$transcript]
      writeXStringSet(tx, tx_save, compress = compress_fa)
    }
  }
  df <- df[, c("transcript", "gene", setdiff(names(df), c("transcript", "gene")))]
  if (write_tr2g) {
    write_tr2g_fun(df, out_path, overwrite)
  }
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
#' @inheritParams tr2g_GRanges
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
#' @seealso ensembl_gene_biotypes ensembl_tx_biotypes cellranger_biotypes
#' @examples
#' library(EnsDb.Hsapiens.v86)
#' tr2g_EnsDb(EnsDb.Hsapiens.v86, get_transcriptome = FALSE, write_tr2g = FALSE,
#'  use_transcript_version = FALSE,
#'  use_gene_version = FALSE)
tr2g_EnsDb <- function(ensdb, Genome = NULL, get_transcriptome = TRUE, 
                       out_path = ".", write_tr2g = TRUE,
                       other_attrs = NULL, use_gene_name = TRUE,
                       use_transcript_version = TRUE,
                       use_gene_version = TRUE,
                       transcript_biotype_col = "TXBIOTYPE",
                       gene_biotype_col = "GENEBIOTYPE", 
                       transcript_biotype_use = "all",
                       gene_biotype_use = "all", chrs_only = TRUE, 
                       compress_fa = FALSE, overwrite = FALSE) {
  check_genome_present(Genome, get_transcriptome)
  if (get_transcriptome || write_tr2g) {
    out_path <- check_out_path(out_path)
  }
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
  if (transcript_biotype_use != "all") {
    attrs_use <- c(attrs_use, transcript_biotype_col)
  }
  if (gene_biotype_use != "all") {
    attrs_use <- c(attrs_use, gene_biotype_col)
  }
  if (chrs_only) {
    chrs_use <- mapSeqlevels(seqlevels(ensdb), seqlevelsStyle(ensdb)[1])
    chrs_use <- unname(chrs_use[!is.na(chrs_use)])
    attrs_use <- c(attrs_use, "SEQNAME")
  }
  attrs_use <- unique(attrs_use)
  cls <- columns(ensdb)
  if (any(!attrs_use %in% cls)) {
    stop("Attributes must be one of ", paste(cls, collapse = ", "),
         "; attributes ", paste(setdiff(attrs_use, cls), collapse = " ,"),
         " are absent.")
  }
  df <- AnnotationDbi::select(ensdb, AnnotationDbi::keys(ensdb, keytype = "TXID"),
    keytype = "TXID",
    columns = attrs_use)
  if (use_transcript_version) {
    df$TXID <- NULL
  }
  if (gene_biotype_use != "all") {
    gbt_use <- which_biotypes(gene_biotype_use, df[[gene_biotype_col]])
    df <- df[df[[gene_biotype_col]] %in% gbt_use,]
  }
  if (transcript_biotype_use != "all") {
    tbt_use <- which_biotypes(transcript_biotype_use, 
                              df[[transcript_biotype_col]])
    df <- df[df[[transcript_biotype_col]] %in% tbt_use,]
  }
  if (chrs_only) {
    df <- df[df$SEQNAME %in% chrs_use,]
  }
  names(df)[str_detect(names(df), "^TXID")] <- "transcript"
  names(df)[str_detect(names(df), "^GENEID")] <- "gene"
  names(df)[names(df) == "GENENAME"] <- "gene_name"
  names(df)[names(df) == transcript_biotype_col] <- "transcript_biotype"
  names(df)[names(df) == gene_biotype_col] <- "gene_biotype"
  names(df)[names(df) == "SEQNAME"] <- "seqnames"
  if (get_transcriptome) {
    tx_save <- paste0(out_path, "/transcriptome.fa")
    if (compress_fa) tx_save <- paste0(tx_save, ".gz")
    if (file.exists(tx_save) && !overwrite) {
      message("File ", tx_save, " already exists.")
    } else {
      c(Genome, ensdb) %<-% match_style(Genome, ensdb, "annotation")
      chrs_use <- intersect(unique(df$seqnames), seqlevels(Genome))
      if (use_transcript_version) {
        # tr2g_cdna has already been filtered if it needs to be filtered
        tx_nv <- str_remove(df$transcript, "\\.\\d+$")
        filter_use <- AnnotationFilterList(SeqNameFilter(chrs_use),
                                           TxIdFilter(tx_nv))
      } else {
        filter_use <- AnnotationFilterList(SeqNameFilter(chrs_use),
                                           TxIdFilter(df$transcript))
      }
      grl <- exonsBy(ensdb, by = "tx", filter = filter_use)
      seqlevels(grl) <- seqlevelsInUse(grl)
      if (use_transcript_version) {
        # Add transcript version
        names(grl) <- df$transcript[match(names(grl), tx_nv)]
      }
      c(Genome, grl) %<-% annot_circular(Genome, grl)
      genome(grl) <- genome(Genome)[seqlevels(grl)]
      tx <- extractTranscriptSeqs(Genome, grl)
      writeXStringSet(tx, tx_save, compress = compress_fa)
    }
  }
  df <- df[, c("transcript", "gene", setdiff(names(df), c("transcript", "gene")))]
  if (write_tr2g) {
    write_tr2g_fun(df, out_path, overwrite)
  }
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
#' @return A data frame with columns \code{transcript} and \code{gene} and the
#' other columns present in \code{tr2g} or the data frame in \code{file}, with
#' the transcript IDs sorted to be in the same order as in the kallisto index.
#' @export
#' @importFrom utils read.table
#' @family functions to retrieve transcript and gene info
#' @note This function has been superseded by the new version of tr2g_* 
#' functions that can extract transcriptome for only the biotypes specified and
#' with only the standard chromosomes. The new version of tr2g_* functions also
#' sorts the transcriptome so the tr2g and the transcriptome have transcripts in
#' the same order.
#' @examples
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file_use <- paste(toy_path, "gtf_test.gtf", sep = "/")
#' tr2g <- tr2g_gtf(file = file_use, get_transcriptome = FALSE,
#'   write_tr2g = FALSE, save_filtered_gtf = FALSE, transcript_version = NULL)
#' tr2g <- sort_tr2g(tr2g, kallisto_out_path = toy_path)
sort_tr2g <- function(tr2g, file, kallisto_out_path) {
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
    tr2g <- read.table(file, header = FALSE, col.names = c("transcript", "gene"),
                       stringsAsFactors = FALSE)
  }
  trs <- read.table(trs_path, header = FALSE, col.names = "transcript",
                    stringsAsFactors = FALSE)
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
#' \code{\link{write.table}}, \code{readr::write_delim}.
#'
#' @inheritParams sort_tr2g
#' @param file_save File name of the file to be saved. The directory in which
#' the file is to be saved must exist.
#' @return Nothing is returned into the R session. A tsv file of the format
#' required by `bustools` with the name and directory specified will be written
#' to disk.
#' @export
#' @importFrom utils write.table
#' @note This function has been superseded by the new version of tr2g_* 
#' functions that can extract transcriptome for only the biotypes specified and
#' with only the standard chromosomes. The new version of tr2g_* functions also
#' sorts the transcriptome so the tr2g and the transcriptome have transcripts in
#' the same order, and write the tr2g.tsv file in the bustools format.
#' @examples
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' file_use <- paste(toy_path, "gtf_test.gtf", sep = "/")
#' tr2g <- tr2g_gtf(file = file_use, get_transcriptome = FALSE, 
#'   write_tr2g = FALSE, save_filtered_gtf = FALSE)
#' save_tr2g_bustools(tr2g, file_save = "./tr2g.tsv")
save_tr2g_bustools <- function(tr2g, file_save = "./tr2g.tsv") {
  file_save <- normalizePath(file_save, mustWork = FALSE)
  write.table(tr2g[, c("transcript", "gene")], file = file_save, sep = "\t",
    col.names = FALSE, quote = FALSE, row.names = FALSE)
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
#' @return A data frame with two columns: \code{gene} and \code{transcript},
#' with Ensembl gene and transcript IDs (with version number), in the same order
#' as in the transcriptome index used in \code{kallisto}.
#' @param \dots Other arguments passed to `tr2g_ensembl` such as `other_attrs`,
#' `ensembl_version`, and arguments passed to \code{\link{useEnsembl}}. If
#' `fasta_files` is supplied instead of `species`, then this will be extra
#' argumennts to \code{\link{tr2g_fasta}}, such as `use_transcript_version` and
#' `use_gene_version`.
#' @importFrom dplyr bind_rows
#' @export
#' @family functions to retrieve transcript and gene info
#' @note This function has been superseded by the new version of tr2g_* 
#' functions that can extract transcriptome for only the biotypes specified and
#' with only the standard chromosomes. The new version of tr2g_* functions also
#' sorts the transcriptome so the tr2g and the transcriptome have transcripts in
#' the same order.
#' @examples
#' # Download dataset already in BUS format
#' library(TENxBUSData)
#' TENxBUSData(".", dataset = "hgmm100")
#' tr2g <- transcript2gene(c("Homo sapiens", "Mus musculus"), 
#'   type = "vertebrate",
#'   ensembl_version = 99, kallisto_out_path = "./out_hgmm100")
transcript2gene <- function(species, fasta_file, kallisto_out_path,
                            type = "vertebrate", ...) {
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
      MoreArgs = MoreArgs,
      SIMPLIFY = FALSE)
    tr2g <- bind_rows(fls)
    return(sort_tr2g(tr2g, kallisto_out_path = kallisto_out_path))
  } else {
    fls <- lapply(fasta_file, tr2g_fasta, ...)
    tr2g <- bind_rows(fls)
    # Just to be safe, to make sure that the transcripts are in the right order
    return(sort_tr2g(tr2g, kallisto_out_path = kallisto_out_path))
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
#' @param verbose Logical, whether to display progress.
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
#' @importFrom tibble tibble
#' @export
#' @examples
#' # Load toy example for testing
#' toy_path <- system.file("testdata", package = "BUSpaRse")
#' load(paste(toy_path, "toy_example.RData", sep = "/"))
#' EC2gene(tr2g_toy, toy_path, verbose = FALSE)
EC2gene <- function(tr2g, kallisto_out_path, verbose = TRUE) {
  kallisto_out_path <- normalizePath(kallisto_out_path, mustWork = TRUE)
  c(ec_vec, genes) %<-% EC2gene_export(tr2g, kallisto_out_path, verbose)
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
