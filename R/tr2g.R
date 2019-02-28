#' @include sparse_matrix.R
NULL

#' Get transcript and gene info from Ensembl
#' 
#' This function queries Ensembl biomart to convert transcript IDs to gene IDs.
#' 
#' @param species Character vector of length 1, Latin name of the species of
#' interest.
#' @param ensembl_version Integer version number of Ensembl (e.g. 94 for the
#' October 2018 release). This argument defaults to \code{NULL}, which will use
#' the current release of Ensembl. Use 
#' \code{\link[biomaRt]{listEnsemblArchives}} to see the version number corresponding
#' to the Ensembl release of a particular date. The version specified here must
#' match the version of Ensembl where the transcriptome used to build the
#' kallisto index was downloaded.
#' @param other_attrs Character vector. Other attributes to get from Ensembl, 
#' such as gene symbol and position on the genome.
#' @param verbose Whether to display progress.
#' @param \dots Othe arguments to be passed to \code{\link[biomaRt]{useEnsembl}},
#' such as host and mirror.
#' Use \code{\link[biomaRt]{listAttributes}} to see which attributes are available.
#' @importFrom biomaRt useEnsembl getBM
#' @importFrom stats setNames
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. If \code{other_attrs}
#' has been specified, then those will also be columns in the data frame returned.
#' @family functions to retrieve transcript and gene info
#' @export
tr2g_ensembl <- function(species, other_attrs = NULL, ensembl_version = NULL, 
                          verbose = TRUE, ...) {
  # Validate arguments
  check_char1(setNames(species, "species"))
  if (!is.null(ensembl_version) && !is.numeric(ensembl_version)) {
    stop("ensembl_version must be integer.\n")
  }
  if (!is.null(other_attrs) && 
      (!is.atomic(other_attrs) || !is.character(other_attrs))) {
    stop("other_attrs must be an atomic character vector.\n")
  }
  mart_name <- species2dataset(species)
  if (verbose) {
    message(paste("Querying biomart for transcript and gene IDs of",
                  species))
  }
  mart <- useEnsembl(biomart = "ensembl", dataset = mart_name, 
                     version = ensembl_version, ...)
  out <- getBM(c("ensembl_gene_id_version", 
                 "ensembl_transcript_id_version",
                 "external_gene_name",
                 other_attrs), mart = mart)
  names(out)[1:3] <- c("gene", "transcript", "gene_name")
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
#' @param file Path to a GTF file to be read. The file can remain gzipped.
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
#' \code{NULL} if you are fine with non-human readable gene IDs. 
#' @param transcript_version Character vector of length 1. Tag in \code{attribute} 
#' field corresponding to _transcript_ version number. If your GTF file does not
#' include transcript version numbers, or if you do not wish to include the
#' version number, then use \code{NULL} for this argument. 
#' @param gene_version Character vector of length 1. Tag in \code{attribute} 
#' field corresponding to _gene_ version number. If your GTF file does not
#' include gene version numbers, or if you do not wish to include the
#' version number, then use \code{NULL} for this argument. 
#' @param version_sep Character to separate bewteen the main ID and the version
#' number. Defaults to ".", as in Ensembl.
#' @inheritParams tr2g_ensembl
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. 
#' @importFrom plyranges read_gff
#' @importFrom magrittr %>% 
#' @importFrom stringr str_detect
#' @importFrom dplyr distinct
#' @importFrom S4Vectors mcols
#' @family functions to retrieve transcript and gene info
#' @export
tr2g_gtf <- function(file, type_use = "exon", transcript_id = "transcript_id",
                      gene_id = "gene_id", gene_name = "gene_name",
                      transcript_version = "transcript_version",
                      gene_version = "gene_version", version_sep = ".",
                      verbose = TRUE) {
  # Validate arguments
  check_char1(setNames(file, "file"))
  file <- normalizePath(file, mustWork = TRUE)
  check_gff("gtf", file, transcript_id, gene_id, gene_name,
            transcript_version, gene_version, version_sep)
  if (verbose) {
    message(paste("Reading GTF file."))
  }
  gr <- read_gff(file)
  tags <- names(mcols(gr))
  check_tag_present(c(transcript_id, gene_id), tags, error = TRUE)
  # Will do nothing if all are NULL
  check_tag_present(c(gene_name, transcript_version, gene_version), 
                    tags, error = FALSE)
  gr <- gr[!is.na(mcols(gr)[[transcript_id]])]
  gr <- gr[gr$type %in% type_use]
  if (length(gr) == 0) {
    stop(paste("No entry has types", paste(type_use, collapse = ", "), 
               "\n"))
  }
  out <- data.frame(gene = mcols(gr)[[gene_id]],
                    transcript = mcols(gr)[[transcript_id]],
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
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. 
#' @family functions to retrieve transcript and gene info
#' @importFrom plyranges read_gff3
#' @importFrom stringr str_split
#' @importFrom dplyr left_join distinct
#' @importFrom tidyr unite
#' @export
tr2g_gff3 <- function(file, type_use = "mRNA", transcript_id = "transcript_id",
                      gene_id = "gene_id", gene_name = "Name",
                      transcript_version = "version",
                      gene_version = "version", version_sep = ".",
                      verbose = TRUE) {
  # Validate arguments
  check_char1(setNames(file, "file"))
  file <- normalizePath(file, mustWork = TRUE)
  check_gff("gff3", file, transcript_id, gene_id, gene_name,
            transcript_version, gene_version, version_sep)
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
    stop(paste("No entry has types", paste(type_use, collapse = ", "), 
               "\n"))
  }
  genes <- str_split(gr_tx$Parent, ":", simplify = TRUE)[,2]
  out <- data.frame(gene = genes,
                    transcript = mcols(gr_tx)[[transcript_id]],
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
#' regex to extract Ensembl IDs. Sequence names should be formatted as follows:
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
#' @return A data frame with 3 columns: \code{gene} for gene ID, \code{transcript}
#' for transcript ID, and \code{gene_name} for gene names. 
#' @importFrom Biostrings readDNAStringSet
#' @importFrom stringr str_extract
#' @importFrom tidyr separate
#' @importFrom dplyr select
#' @family functions to retrieve transcript and gene info
#' 
tr2g_fasta <- function(file, verbose = TRUE) {
  check_char1(setNames(file, "file"))
  file <- normalizePath(file, mustWork = TRUE)
  if (!str_detect(file, "(\\.fasta)|(\\.fa)|(\\.fna)")) {
    stop("file must be a FASTA file.\n")
  }
  file <- normalizePath(file, mustWork = TRUE)
  if (verbose) {
    message("Reading FASTA file.")
  }
  s <- readDNAStringSet(file)
  # Avoid R CMD check note
  g <- gene_name <- NULL
  out <- data.frame(gene = str_extract(names(s), "ENS[A-Z]*G\\d+\\.\\d+"),
                    transcript = str_extract(names(s), "ENS[A-Z]*T\\d+\\.\\d+"),
                    gene_name = str_extract(names(s), "gene_symbol:[a-zA-Z\\d-\\.]+"),
                    stringsAsFactors = FALSE) %>% 
    tidyr::separate(gene_name, into = c("g", "gene_name"), sep = ":") %>% 
    dplyr::select(-g) %>% 
    distinct()
  out
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
#' Exactly one of \code{tr2g} and \code{file} should be missing.
#' @param file Character vector of length 1, path to a csv or tsv file with
#' transcript IDs and the corresponding gene IDs. Headers \code{transcript} and
#' \code{gene} must be present in the file.
#' @param kallisto_out_path Character vector of length 1, path to the directory
#' for the outputs of kallisto bus.
#' @param save Whether to save the output.
#' @param verbose Whether to display progress.
#' @param \dots Other arguments passed to \code{\link[data.table]{fwrite}}, such
#' as \code{sep}, \code{quote}, and \code{col.names}.
#' @param file_save File name of the file to be saved. If the directory in which
#' the file is to be saved does not exist, then the directory will be created.
#' @return A data frame with columns \code{transcript} and \code{gene} and the
#' other columns present in \code{tr2g} or the data frame in \code{file}, with
#' the transcript IDs sorted to be in the same order as in the kallisto index.
#' When \code{save = TRUE}, the data frame is not only saved on disk but also
#' returned in the R session.
#' @importFrom data.table fread fwrite
#' @export
#' @family functions to retrieve transcript and gene info
#' 
sort_tr2g <- function(tr2g, file, kallisto_out_path, save = FALSE,
                      file_save = "./tr2g_sorted.csv", verbose = TRUE, ...) {
  if (!xor(missing(tr2g), missing(file))) {
    stop("Exactly one of tr2g and file should be missing.\n")
  }
  kallisto_out_path <- normalizePath(kallisto_out_path, mustWork = TRUE)
  trs_path <- paste(kallisto_out_path, "transcripts.txt", sep = "/")
  if (!file.exists(trs_path)) {
    stop("The file transcripts.txt does not exist in",
               kallisto_out_path, "\n")
  }
  if (missing(tr2g)) {
    tr2g <- fread(file)
  }
  trs <- fread(trs_path, header = FALSE, col.names = "transcript")
  if (verbose) {
    message("Sorting transcripts")
  }
  out <- merge(trs, tr2g, by = "transcript", sort = FALSE)
  if (nrow(trs) != nrow(out)) {
    stop("Some transcripts in the kallisto index absent from tr2g.\n")
  }
  if (save) {
    file_save <- normalizePath(file_save, mustWork = FALSE)
    dn <- dirname(file_save)
    if (!dir.exists(dn)) dir.create(dn)
    fwrite(out, file_save, ...)
  } else {
    extra_args <- list(...)
    if (length(extra_args) > 0) {
      arg_names <- paste(names(extra_args), collapse = ", ")
      message(paste("Arguments", arg_names, "are ignored."))
    }
  }
  out
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
#' Note that here, the arguments passed to \dots will be passed to 
#' \code{\link[biomaRt]{useEnsembl}} rather than \code{\link[data.table]{fwrite}},
#' so default settings of \code{\link[data.table]{fwrite}} will be used if 
#' \code{save = TRUE}, which will save a csv file that includes the column names
#' to \code{file_save}.
#' 
#' @inheritParams tr2g_ensembl
#' @inheritParams sort_tr2g
#' @param species A character vector of Latin names of species present in this
#' scRNA-seq dataset. This is used to retrieve Ensembl information from biomart.
#' @param fasta_file Character vector of paths to the transcriptome FASTA files
#' used to build the kallisto index. Exactly one of \code{species} and
#' \code{fasta_file} can be missing.
#' @param kallisto_out_path Path to the \code{kallisto bus} output directory.
#' @param verbose Whether to display progress. Defaults to \code{TRUE}.
#' @return A data frame with two columns: \code{gene} and \code{transcript},
#' with Ensembl gene and transcript IDs (with version number), in the same order
#' as in the transcriptome index used in \code{kallisto}.
#' @importFrom data.table rbindlist
#' @export
#' @family functions to retrieve transcript and gene info
#' @examples
#' \dontrun{
#' # Download dataset already in BUS format
#' library(TENxhgmmBUS)
#' download_hgmm(".", "hgmm100")
#' tr2g <- transcript2gene(c("Homo sapiens", "Mus musculus"),
#'                         kallisto_out_path = "./out_hgmm100")
#' }
#' 
transcript2gene <- function(species, fasta_file, kallisto_out_path, 
                            other_attrs = NULL, ensembl_version = NULL,
                            save = FALSE, file_save = "./tr2g_sorted.csv",
                            verbose = TRUE, ...) {
  if (!xor(missing(species), missing(fasta_file))) {
    stop("Exactly one of species and fasta_file can be missing.\n")
  }
  if (missing(fasta_file)) {
    kallisto_out_path <- normalizePath(kallisto_out_path, mustWork = TRUE)
    fls <- lapply(species, tr2g_ensembl, other_attrs = other_attrs,
                  ensembl_version = ensembl_version, verbose = verbose, ...)
    tr2g <- rbindlist(fls)
    return(sort_tr2g(tr2g, kallisto_out_path = kallisto_out_path, save = save,
              file_save = file_save, verbose = verbose))
  } else {
    if (!is.null(ensembl_version) || length(list(...)) > 0) {
      message("Arguments related to Ensembl biomart queries are ignored.")
    }
    fls <- lapply(fasta_file, tr2g_fasta, verbose = verbose)
    tr2g <- rbindlist(fls)
    # Just to be safe, to make sure that the transcripts are in the right order
    return(sort_tr2g(tr2g, kallisto_out_path = kallisto_out_path, 
                      save = save, file_save = file_save, 
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
#' used.
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
#' \dontrun{
#' # Download dataset already in BUS format
#' library(TENxhgmmBUS)
#' download_hgmm(".", "hgmm100")
#' tr2g <- transcript2gene(c("Homo sapiens", "Mus musculus"),
#'                         kallisto_out_path = "./out_hgmm100")
#' genes <- EC2gene(tr2g, "./out_hgmm100", ncores = 1, verbose = FALSE)
#' }
#' 
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
