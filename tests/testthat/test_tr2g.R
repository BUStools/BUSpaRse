context("Test function to convert transcript to gene")
library(testthat)
library(dplyr)

# Load toy example for testing
toy_path <- system.file("testdata", package = "BUSpaRse")
load(paste(toy_path, "toy_example.RData", sep = "/"))
tr2g_expected <- read.csv(paste(toy_path, "tr2g_expected.csv", sep = "/"),
                          header = TRUE, stringsAsFactors = FALSE)
tr2g_expected_version <- read.csv(paste(toy_path, "tr2g_expected_version.csv", sep = "/"),
                                  header = TRUE, stringsAsFactors = FALSE)
fa_tr2g_expected <- read.csv(paste(toy_path, "fa_tr2g_expected.csv", sep = "/"),
                             header = TRUE, stringsAsFactors = FALSE)

test_that("Correct data set name for biomart", {
  error_use <- "Please use the Latin binomial convention"
  expect_equal(species2dataset("Homo sapiens"), "hsapiens_gene_ensembl")
  expect_equal(species2dataset("Felis catus"), "fcatus_gene_ensembl")
  expect_error(species2dataset("cats are so cute"), error_use)
  expect_error(species2dataset("mouse"), error_use)
})

test_that("Extract transcript and gene ID from GTF file", {
  # No version number
  fn <- paste(toy_path, "gtf_test.gtf", sep = "/")
  tr2g_no_vn <- tr2g_gtf(fn, gene_version = NULL, transcript_version = NULL) %>% 
    arrange(gene)
  # No gene name
  expect_equal(tr2g_no_vn, tr2g_expected)
  tr2g_no_gn <- tr2g_gtf(fn, gene_name = NULL, gene_version = NULL,
                         transcript_version = NULL) %>% 
    arrange(gene)
  expect_equal(tr2g_no_gn, tr2g_expected[,-3])
  # With version number
  expect_equal(tr2g_gtf(fn) %>% arrange(gene), tr2g_expected_version)
  
  # Test error messages
  # transcript_id or gene_id is wrong
  expect_error(tr2g_gtf(fn, transcript_id = "foo"), "Tags foo are absent")
  expect_error(tr2g_gtf(fn, gene_id = "bar"), "Tags bar are absent")
  # transcript_id or gene_id is NULL
  expect_error(tr2g_gtf(fn, transcript_id = NULL), "transcript_id cannot be NULL.")
  expect_error(tr2g_gtf(fn, gene_id = NULL), "gene_id cannot be NULL.")
  # gene_name, gene_version, or transcript_version is wrong
  expect_warning(tr2g_gtf(fn, gene_name = "foo", gene_version = "bar"),
                 "Tags foo, bar are absent")
  # Error when input file is not a gtf file
  expect_error(tr2g_gtf(paste(toy_path, "fasta_test.fasta", sep = "/")),
               "file must be a GTF file.")
  expect_error(tr2g_gtf(fn, type_use = "foo"),
               "No entry has types foo")
})

test_that("Extract transcript and gene ID from GFF3 file", {
  # No version number
  fn <- paste(toy_path, "gff3_test.gff3", sep = "/")
  tr2g_no_vn <- tr2g_gff3(fn, type_use = c("mRNA", "lnc_RNA"),
                          gene_version = NULL, transcript_version = NULL) %>% 
    arrange(gene)
  # No gene name
  expect_equal(tr2g_no_vn, tr2g_expected)
  tr2g_no_gn <- tr2g_gff3(fn, type_use = c("mRNA", "lnc_RNA"),
                          gene_name = NULL, gene_version = NULL,
                         transcript_version = NULL) %>% 
    arrange(gene)
  expect_equal(tr2g_no_gn, tr2g_expected[,-3])
  # With version number
  expect_equal(tr2g_gff3(fn, type_use = c("mRNA", "lnc_RNA")) %>% arrange(gene),
               tr2g_expected_version)
  
  # Test error messages
  # transcript_id or gene_id is wrong
  expect_error(tr2g_gff3(fn, transcript_id = "foo"), "Tags foo are absent")
  expect_error(tr2g_gff3(fn, gene_id = "bar"), "Tags bar are absent")
  # transcript_id or gene_id is NULL
  expect_error(tr2g_gff3(fn, transcript_id = NULL), "transcript_id cannot be NULL.")
  expect_error(tr2g_gff3(fn, gene_id = NULL), "gene_id cannot be NULL.")
  # gene_name, gene_version, or transcript_version is wrong
  expect_warning(tr2g_gff3(fn, gene_name = "foo", gene_version = "bar"),
                 "Tags foo, bar are absent")
  # Error when input file is not a gff3 file
  expect_error(tr2g_gff3(paste(toy_path, "fasta_test.fasta", sep = "/")),
               "file must be a GFF3 file.")
  expect_error(tr2g_gff3(fn, type_use = "foo"),
               "No entry has types foo")
})

test_that("Extract transcript and gene ID from FASTA file", {
  fn <- paste(toy_path, "fasta_test.fasta", sep = "/")
  expect_equal(tr2g_fasta(fn), fa_tr2g_expected)
  expect_error(tr2g_fasta(file = 3), 
               "file must be a character vector with length 1")
  expect_error(tr2g_fasta(paste(toy_path, "gtf_test.gtf", sep = "/")),
               "file must be a FASTA file.")
})

test_that("Correct gene list output", {
  expect_equal(EC2gene(tr2g_toy, toy_path, verbose = FALSE),
               genes_toy)
  # Test multicore
  expect_equal(EC2gene(tr2g_toy, toy_path, ncores = 2, verbose = FALSE),
               genes_toy)
})
