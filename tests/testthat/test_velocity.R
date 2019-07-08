context("Test for correct files output for RNA velocity")

library(Biostrings)

# Load toy examples
toy_path <- system.file("testdata", package = "BUSpaRse")
toy_genome <- readDNAStringSet(paste(toy_path, "velocity_genome.fa",
                                     sep = "/"))
L <- 11

# Some functions just for these tests
test_fasta <- function(toy_path, out_path, 
                       option = c("trunc", "inc", "coll")) {
  option <- match.arg(option)
  expected_fa <- readDNAStringSet(paste0(toy_path, "/velocity_introns_",
                                         option, ".fa"))
  actual_fa <- readDNAStringSet(paste(out_path, "cDNA_introns.fa", sep = "/"))
  # sort
  actual_fa <- actual_fa[names(expected_fa)]
  expect_equivalent(actual_fa, expected_fa)
}

test_tr2g <- function(toy_path, out_path,
                      option = c("trunc", "inc", "coll")) {
  option <- match.arg(option)
  expected_tr2g <- read.table(paste0(toy_path, "/velocity_tr2g_", option, ".tsv"),
                              header = FALSE, sep = "\t",
                              col.names = c("transcript", "gene"),
                              stringsAsFactors = FALSE)
  actual_tr2g <- read.table(paste(out_path, "tr2g.tsv", sep = "/"), 
                            header = FALSE, sep = "\t",
                            col.names = c("transcript", "gene"),
                            stringsAsFactors = FALSE)
  actual_tr2g <- actual_tr2g[match(expected_tr2g$transcript, actual_tr2g$transcript),]
  expect_equivalent(actual_tr2g, expected_tr2g)
}

test_tx_capture <- function(toy_path, out_path) {
  expected_tx_capture <- readLines(paste(toy_path, "velocity_cdna_tx.txt", 
                                         sep = "/"))
  actual_tx_capture <- readLines(paste(out_path, "cDNA_tx_to_capture.txt",
                                       sep = "/"))
  actual_tx_capture <- actual_tx_capture[match(expected_tx_capture, actual_tx_capture)]
  expect_equal(actual_tx_capture, expected_tx_capture)
}

test_intron_capture <- function(toy_path, out_path,
                                option = c("trunc", "inc", "coll")) {
  option <- match.arg(option)
  expected_intron_capture <- readLines(paste0(toy_path, 
                                             "/velocity_introns_", option,
                                             "_capture.txt"))
  actual_intron_capture <- readLines(paste(out_path, 
                                           "introns_tx_to_capture.txt",
                                           sep = "/"))
  actual_intron_capture <- actual_intron_capture[match(expected_intron_capture, actual_intron_capture)]
  expect_equal(actual_intron_capture, expected_intron_capture)
}

test_that("Truncate short exons, keep isoforms separate", {
  out_path <- "./trunc_sep"
  get_velocity_files(paste(toy_path, "velocity_annot.gtf", sep = "/"), L = L,
                     genome = toy_genome, 
                     transcriptome = paste(toy_path, "velocity_tx.fa", sep = "/"),
                     out_path = out_path,
                     short_exon_action = "truncate",
                     isoform_action = "separate",
                     transcript_version = NULL,
                     gene_version = NULL)
  # fasta file
  test_fasta(toy_path, out_path, "trunc")
  # tr2g
  test_tr2g(toy_path, out_path, "trunc")
  # transcript to capture
  test_tx_capture(toy_path, out_path)
  # introns to capture
  test_intron_capture(toy_path, out_path, "trunc")
  unlink(out_path, recursive = TRUE)
})

test_that("Bypass short exons, keep isoforms separate", {
  out_path <- "./inc_sep"
  get_velocity_files(paste(toy_path, "velocity_annot.gtf", sep = "/"), L = L,
                     genome = toy_genome, 
                     transcriptome = paste(toy_path, "velocity_tx.fa", sep = "/"),
                     out_path = out_path,
                     short_exon_action = "include",
                     isoform_action = "separate",
                     transcript_version = NULL,
                     gene_version = NULL)
  # fasta file
  test_fasta(toy_path, out_path, "inc")
  # tr2g
  test_tr2g(toy_path, out_path, "inc")
  # transcript to capture
  test_tx_capture(toy_path, out_path)
  # introns to capture
  test_intron_capture(toy_path, out_path, "inc")
  unlink(out_path, recursive = TRUE)
})

test_that("Truncate short exons, collapse isoforms", {
  out_path <- "./trunc_coll"
  get_velocity_files(paste(toy_path, "velocity_annot.gtf", sep = "/"), L = L,
                     genome = toy_genome, 
                     transcriptome = paste(toy_path, "velocity_tx.fa", sep = "/"),
                     out_path = out_path,
                     short_exon_action = "truncate",
                     isoform_action = "collapse",
                     transcript_version = NULL,
                     gene_version = NULL)
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
