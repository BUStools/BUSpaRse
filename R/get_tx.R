#' Download transcriptome from Ensembl
#' 
#' This function downloads the cDNA fasta file from specific version of Ensembl.
#' It can also filter the fasta file by gene and transcript biotype and remove
#' scaffolds and haplotypes.
#' 
#' @inheritParams tr2g_ensembl
#' @inheritParams .get_velocity_files
#' @param ... Other arguments passed to `tr2g_fasta`.
#' @return Invisibly returns the path to the fasta file. The following files are
#' written to disk, in the `out_path` directory:
#' \describe{
#' \item{species.genome.cdna.all.fa.gz}{The cDNA fasta file from Ensembl, from
#' the specified version.}
#' \item{cdna_filtered.fa}{The filtered cDNA fasta file, only keeping the 
#' desired biotypes and without scaffolds and haplotypes (if 
#' `chrs_only = TRUE`). This file will not be written if all gene and transcript
#' biotypes are used and scaffolds and haplotypes are not removed.}
#' \item{tr2g.tsv}{The transcript to gene file, without headers so can be
#' directly used for `bustools`.}
#' }
#' @export
#' @importFrom biomaRt searchDatasets listMarts listEnsemblArchives
#' @importFrom utils download.file
#' @examples 
#' dl_transcriptome("Drosophila melanogaster", gene_biotype_use = "cellranger")
dl_transcriptome <- function(species, out_path = ".",
                             type = c("vertebrate", "metazoa", "plant",
                                      "fungus", "protist"),
                             transcript_biotype_use = "all",
                             gene_biotype_use = "all", 
                             chrs_only = TRUE,
                             ensembl_version = NULL,
                             verbose = TRUE, ...) {
  type <- match.arg(type)
  out_path <- normalizePath(out_path, mustWork = FALSE)
  if (!dir.exists(out_path)) dir.create(out_path)
  if (length(strsplit(species, " ")[[1]]) != 2L) {
    stop("Please use the Latin binomial convention for species rather than the colloquial name.")
  }
  if (!is.null(ensembl_version) && !is.numeric(ensembl_version)) {
    stop("ensembl_version must be integer.")
  }
  styles <- try(genomeStyles(species))
  if (is(styles, "try-error")) {
    warning("All seqnames are used.")
    chrs_only <- FALSE
  }
  # Get current ensembl version
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
  if (is.null(ensembl_version)) {
    ds <- listMarts(host = host_use)
    ev <- str_extract(ds$version[1], "\\d+") %>% as.integer()
  } else ev <- ensembl_version
  if (type != "vertebrate" & !is.null(ensembl_version)) {
    warning("Cannot retrieve genome version of non-vertebrates from archive. ",
            "Using genome of current version instead.")
    ensembl_version <- NULL
  }
  mart <- my_useMart(ensembl_version, mart_use, ds_name, host_use)
  # Get genome version
  ds2<- searchDatasets(mart = mart, ds_name)
  gv <- ds2$version
  if (str_detect(gv, "^GRC")) {
    gv <- str_extract(gv, "^GRC[a-z\\d]+")
  }
  
  # Build URL for download
  species_url <- str_replace(species, " ", "_") %>% tolower()
  if (type == "vertebrate") {
    part1 <- "ftp://ftp.ensembl.org/pub/release-"
  }
  else {
    part1 <- paste0("ftp://ftp.ensemblgenomes.org/pub/", host_pre, "/release-")
  }
  fn <- paste0(str_replace(species, " ", "_"), ".", gv, ".cdna.all.fa.gz")
  url <- paste0(part1, ev, "/fasta/", species_url, "/cdna/", fn)
  destfile <- paste(out_path, fn, sep = "/")
  # Download fasta file
  if (file.exists(destfile)) {
    message("File ", destfile, " already exists, skipping download.")
  } else {
    download.file(url, destfile = destfile, quiet = !verbose)
  }
  
  # Filter biotype
  do_filter <- transcript_biotype_use != "all" | gene_biotype_use != "all" |
    chrs_only
  if (do_filter) {
    file_filtered <- paste(out_path, "tx_filtered.fa", sep = "/")
    tr2g <- tr2g_fasta(destfile, out_path = out_path, write_tr2g = TRUE,
                       gene_biotype_use = gene_biotype_use,
                       transcript_biotype_use = transcript_biotype_use,
                       chrs_only = chrs_only, save_filtered = TRUE, ...)
    file_return <- file_filtered
  } else {
    tr2g <- tr2g_fasta(destfile, out_path = out_path, gene_biotype_use = "all",
                       transcript_biotype_use = "all", write_tr2g = TRUE,
                       chrs_only = FALSE, save_filtered = FALSE, ...)
    file_return <- destfile
  }
  invisible(file_return)
}
