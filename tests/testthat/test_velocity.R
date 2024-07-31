context("Test for correct files output for RNA velocity")

library(Biostrings)
library(GenomicFeatures)
library(GenomicRanges)
library(plyranges)

# Load toy examples
toy_path <- system.file("testdata", package = "BUSpaRse")
toy_genome <- readDNAStringSet(file.path(toy_path, "velocity_genome.fa"))
L <- 16
txdb <- txdbmaker::makeTxDbFromGFF(paste(toy_path, "velocity_annot.gtf", sep = "/"))

# Some functions just for these tests
test_fasta <- function(toy_path, out_path, option = "normal") {
  fa_fn <- switch(option,
    normal = "velocity_introns.fa",
    coll = "velocity_introns_coll.fa",
    junc = "junctions.fa",
    junc_coll = "junction_coll.fa")
  expected_fa <- readDNAStringSet(file.path(toy_path, fa_fn))
  actual_fa <- readDNAStringSet(file.path(out_path, "cDNA_introns.fa"))
  # sort
  actual_fa <- actual_fa[names(expected_fa)]
  expect_true(all(abs(pcompare(actual_fa, expected_fa)) < 1e-8))
}

test_tr2g <- function(toy_path, out_path, option = "normal") {
  tr2g_fn <- switch(option,
    normal = "/velocity_tr2g.tsv",
    coll = "/velocity_tr2g_coll.tsv",
    junc = "/junction_tr2g.tsv",
    junc_coll = "/junction_tr2g_coll.tsv")
  expected_tr2g <- read.table(paste0(toy_path, tr2g_fn),
    header = FALSE, sep = "\t",
    col.names = c("transcript", "gene"),
    stringsAsFactors = FALSE)
  actual_tr2g <- read.table(file.path(out_path, "tr2g.tsv"),
    header = FALSE, sep = "\t",
    col.names = c("transcript", "gene"),
    stringsAsFactors = FALSE)
  actual_tr2g <- actual_tr2g[match(expected_tr2g$transcript, actual_tr2g$transcript), ]
  expect_equivalent(actual_tr2g, expected_tr2g)
}

test_tx_capture <- function(toy_path, out_path, option = "normal") {
  cap_fn <- switch(option,
    normal = "velocity_cdna_tx.txt",
    junc = "junction_capture.txt")
  expected_tx_capture <- readLines(file.path(toy_path, cap_fn))
  actual_tx_capture <- readLines(file.path(out_path, "cDNA_tx_to_capture.txt"))
  actual_tx_capture <- actual_tx_capture[match(expected_tx_capture, actual_tx_capture)]
  expect_equal(actual_tx_capture, expected_tx_capture)
}

test_intron_capture <- function(toy_path, out_path, option = "normal") {
  ins_fn <- switch(option,
    normal = "/velocity_introns_capture.txt",
    coll = "/velocity_introns_coll_capture.txt")
  expected_intron_capture <- readLines(paste0(toy_path, ins_fn))
  actual_intron_capture <- readLines(file.path(out_path,
                                               "introns_tx_to_capture.txt"))
  actual_intron_capture <- actual_intron_capture[match(expected_intron_capture, actual_intron_capture)]
  expect_equal(actual_intron_capture, expected_intron_capture)
}

test_that("GRanges, keep isoforms separate", {
  out_path <- "./trunc_sep"
  get_velocity_files(file.path(toy_path, "velocity_annot.gtf"), L = L,
    Genome = toy_genome,
    Transcriptome = file.path(toy_path, "velocity_tx.fa"),
    out_path = out_path,
    isoform_action = "separate",
    exon_option = "full",
    transcript_version = NULL,
    gene_version = NULL)
  # fasta file
  test_fasta(toy_path, out_path)
  # tr2g
  test_tr2g(toy_path, out_path)
  # transcript to capture
  test_tx_capture(toy_path, out_path)
  # introns to capture
  test_intron_capture(toy_path, out_path)
  unlink(out_path, recursive = TRUE)
})

test_that("GRanges, collapse isoforms", {
  out_path <- "./trunc_coll"
  get_velocity_files(file.path(toy_path, "velocity_annot.gtf"), L = L,
    Genome = toy_genome,
    Transcriptome = file.path(toy_path, "velocity_tx.fa"),
    out_path = out_path,
    isoform_action = "collapse",
    transcript_version = NULL,
    gene_version = NULL)
  # fasta file
  test_fasta(toy_path, out_path, option = "coll")
  # tr2g
  test_tr2g(toy_path, out_path, option = "coll")
  # transcript to capture
  test_tx_capture(toy_path, out_path)
  # introns to capture
  test_intron_capture(toy_path, out_path, option = "coll")
  unlink(out_path, recursive = TRUE)
})

test_that("TxDb, keep isoforms separate", {
  out_path <- "./trunc_sep"
  expect_warning(get_velocity_files(txdb, L = L,
    Genome = toy_genome,
    Transcriptome = file.path(toy_path, "velocity_tx.fa"),
    out_path = out_path, isoform_action = "separate", chrs_only = FALSE),
    regexp = "isCircular information")
  # fasta file
  test_fasta(toy_path, out_path)
  # tr2g
  test_tr2g(toy_path, out_path)
  # transcript to capture
  test_tx_capture(toy_path, out_path)
  # introns to capture
  test_intron_capture(toy_path, out_path)
  unlink(out_path, recursive = TRUE)
})

test_that("TxDb, collapse isoforms", {
  out_path <- "./trunc_coll"
  expect_warning(get_velocity_files(txdb, L = L,
    Genome = toy_genome,
    Transcriptome = file.path(toy_path, "velocity_tx.fa"),
    out_path = out_path, isoform_action = "collapse", chrs_only = FALSE),
    regexp = "isCircular information")
  # fasta file
  test_fasta(toy_path, out_path, "coll")
  # tr2g
  test_tr2g(toy_path, out_path, "coll")
  # transcript to capture
  test_tx_capture(toy_path, out_path)
  # introns to capture
  test_intron_capture(toy_path, out_path, "coll")
  unlink(out_path, recursive = TRUE)
})

test_that("When transcriptome is missing, GRanges", {
  out_path <- "./trunc_sep"
  get_velocity_files(file.path(toy_path, "velocity_annot.gtf"), L = L,
    Genome = toy_genome,
    out_path = out_path,
    isoform_action = "separate",
    transcript_version = NULL,
    gene_version = NULL)
  # fasta file
  test_fasta(toy_path, out_path, "normal")
  # tr2g
  test_tr2g(toy_path, out_path, "normal")
  # transcript to capture
  test_tx_capture(toy_path, out_path)
  # introns to capture
  test_intron_capture(toy_path, out_path, "normal")
  unlink(out_path, recursive = TRUE)
})

test_that("When transcriptome is missing, TxDb", {
  out_path <- "./trunc_coll"
  expect_warning(get_velocity_files(txdb, L = L,
    Genome = toy_genome,
    out_path = out_path,
    isoform_action = "collapse"),
  regexp = "isCircular information")
  # fasta file
  test_fasta(toy_path, out_path, "coll")
  # tr2g
  test_tr2g(toy_path, out_path, "coll")
  # transcript to capture
  test_tx_capture(toy_path, out_path)
  # introns to capture
  test_intron_capture(toy_path, out_path, "coll")
  unlink(out_path, recursive = TRUE)
})

test_that("Get exon-exon junctions, GRanges", {
  out_path <- "./trunc_sep"
  get_velocity_files(file.path(toy_path, "velocity_annot.gtf"), L = L,
    Genome = toy_genome,
    out_path = out_path,
    isoform_action = "separate",
    exon_option = "junction",
    transcript_version = NULL,
    gene_version = NULL)
  # fasta file
  test_fasta(toy_path, out_path, "junc")
  # tr2g
  test_tr2g(toy_path, out_path, "junc")
  # transcript to capture
  test_tx_capture(toy_path, out_path, "junc")
  # introns to capture
  test_intron_capture(toy_path, out_path)
  unlink(out_path, recursive = TRUE)
})

test_that("Get exon-exon junctions, GRanges, collapse isoform", {
  out_path <- "./trunc_coll"
  get_velocity_files(file.path(toy_path, "velocity_annot.gtf"), L = L,
    Genome = toy_genome,
    out_path = out_path,
    isoform_action = "collapse",
    exon_option = "junction",
    transcript_version = NULL,
    gene_version = NULL)
  # fasta file
  test_fasta(toy_path, out_path, "junc_coll")
  # tr2g
  test_tr2g(toy_path, out_path, "junc_coll")
  # transcript to capture
  test_tx_capture(toy_path, out_path, "junc")
  # introns to capture
  test_intron_capture(toy_path, out_path, "coll")
  unlink(out_path, recursive = TRUE)
})

test_that("Get exon-exon junctions, TxDb", {
  out_path <- "./trunc_sep"
  expect_warning(get_velocity_files(txdb, L = L,
    Genome = toy_genome,
    out_path = out_path,
    isoform_action = "separate",
    exon_option = "junction"),
  regexp = "isCircular information")
  # fasta file
  test_fasta(toy_path, out_path, "junc")
  # tr2g
  test_tr2g(toy_path, out_path, "junc")
  # transcript to capture
  test_tx_capture(toy_path, out_path, "junc")
  # introns to capture
  test_intron_capture(toy_path, out_path)
  unlink(out_path, recursive = TRUE)
})

test_that("Get exon-exon junctions, TxDb, collapse isoforms", {
  out_path <- "./trunc_coll"
  expect_warning(get_velocity_files(txdb, L = L,
    Genome = toy_genome,
    out_path = out_path,
    isoform_action = "collapse",
    exon_option = "junction"),
  regexp = "isCircular information")
  # fasta file
  test_fasta(toy_path, out_path, "junc_coll")
  # tr2g
  test_tr2g(toy_path, out_path, "junc_coll")
  # transcript to capture
  test_tx_capture(toy_path, out_path, "junc")
  # introns to capture
  test_intron_capture(toy_path, out_path, "coll")
  unlink(out_path, recursive = TRUE)
})
