context("Test function to convert transcript to gene")
library(testthat)
library(dplyr)
library(stringr)
library(tibble)

# Load toy example for testing
toy_path <- system.file("testdata", package = "BUSpaRse")
load(paste(toy_path, "toy_example.RData", sep = "/"))
tr2g_expected <- read.csv(paste(toy_path, "tr2g_expected.csv", sep = "/"),
  header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()
tr2g_expected_version <- read.csv(paste(toy_path, "tr2g_expected_version.csv", sep = "/"),
  header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()
fa_tr2g_expected <- read.csv(paste(toy_path, "fa_tr2g_expected.csv", sep = "/"),
  header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()
fa_tr2g_no_version <- read.csv(paste(toy_path, "fa_tr2g_no_version.csv", sep = "/"),
  header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()
fa_tr2g_dm <- read.csv(paste(toy_path, "fa_tr2g_dm.csv", sep = "/"),
  header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()

test_that("Correct data set name for biomart", {
  error_use <- "Please use the Latin binomial convention"
  expect_equal(species2dataset("Homo sapiens"), "hsapiens_gene_ensembl")
  expect_equal(species2dataset("Felis catus"), "fcatus_gene_ensembl")
  expect_equal(species2dataset("Arabidopsis thaliana", "plant"), "athaliana_eg_gene")
  expect_equal(species2dataset("Aspergillus nidulans", "fungus"), "anidulans_eg_gene")
  expect_error(species2dataset("Mus musculus", "foo"))
  expect_error(species2dataset("cats are so cute"), error_use)
  expect_error(species2dataset("mouse"), error_use)
})

test_that("Sensible results from biomart query", {
  # Specify version of Ensembl
  tr2g <- tr2g_ensembl(species = "Danio rerio",
    use_gene_version = FALSE, use_transcript_version = TRUE,
    ensembl_version = 94, chrs_only = FALSE)
  expect_equal(names(tr2g), c("transcript", "gene", "gene_name"))
  expect_gt(nrow(tr2g), 1)
  # version numbers
  expect_false(any(str_detect(tr2g$gene, "\\.\\d+$")))
  expect_true(all(str_detect(tr2g$transcript, "\\.\\d+$")))
  expect_message({
    tr2g <- tr2g_ensembl(species = "Arabidopsis thaliana", type = "plant",
      use_transcript_version = TRUE)
  },
  "Version is only available to vertebrates.")
  expect_equal(names(tr2g), c("transcript", "gene", "gene_name",
                              "chromosome_name"))
  expect_gt(nrow(tr2g), 1)
  rm(tr2g)
})

test_that("Extract transcript and gene ID from GTF file", {
  # No version number
  fn <- paste(toy_path, "gtf_test.gtf", sep = "/")
  tr2g_no_vn <- tr2g_gtf(fn, get_transcriptome = FALSE, write_tr2g = FALSE,
                         save_filtered_gtf = FALSE,
                         gene_version = NULL, transcript_version = NULL) %>%
    arrange(gene)
  # No gene name
  expect_equal(tr2g_no_vn, tr2g_expected)
  tr2g_no_gn <- tr2g_gtf(fn, get_transcriptome = FALSE, write_tr2g = FALSE,
                         save_filtered_gtf = FALSE,
                         gene_name = NULL, gene_version = NULL,
                         transcript_version = NULL) %>%
    arrange(gene)
  expect_equal(tr2g_no_gn, tr2g_expected[, -3])
  # With version number
  expect_equal(tr2g_gtf(fn, get_transcriptome = FALSE, write_tr2g = FALSE,
                        save_filtered_gtf = FALSE) %>% 
                 arrange(gene), tr2g_expected_version)

  # Test error messages
  # transcript_id or gene_id is wrong
  expect_error(tr2g_gtf(fn, transcript_id = "foo", 
                        get_transcriptome = FALSE, write_tr2g = FALSE), 
               "Tags foo are absent")
  expect_error(tr2g_gtf(fn, gene_id = "bar", write_tr2g = FALSE, 
                        get_transcriptome = FALSE), 
               "Tags bar are absent")
  # transcript_id or gene_id is NULL
  expect_error(tr2g_gtf(fn, transcript_id = NULL, write_tr2g = FALSE, 
                        get_transcriptome = FALSE), 
               "transcript_id cannot be NULL.")
  expect_error(tr2g_gtf(fn, gene_id = NULL, write_tr2g = FALSE, 
                        get_transcriptome = FALSE), 
               "gene_id cannot be NULL.")
  # gene_name, gene_version, or transcript_version is wrong
  expect_warning(tr2g_gtf(fn, gene_name = "foo", gene_version = "bar",
                          write_tr2g = FALSE, get_transcriptome = FALSE,
                          save_filtered_gtf = FALSE),
    "Tags foo, bar are absent")
  # Error when input file is not a gtf file
  expect_error(tr2g_gtf(paste(toy_path, "fasta_test.fasta", sep = "/"),
                        write_tr2g = FALSE, get_transcriptome = FALSE,
                        save_filtered_gtf = FALSE),
    "file must be a GTF file.")
  rm(list = c("fn", "tr2g_no_vn", "tr2g_no_gn"))
})

test_that("Extract transcript and gene ID from GFF3 file", {
  # No version number
  fn <- paste(toy_path, "gff3_test.gff3", sep = "/")
  tr2g_no_vn <- tr2g_gff3(fn, get_transcriptome = FALSE, save_filtered_gff = FALSE,
    gene_version = NULL, transcript_version = NULL, write_tr2g = FALSE) %>%
    arrange(gene)
  # No gene name
  expect_equal(tr2g_no_vn, tr2g_expected)
  tr2g_no_gn <- tr2g_gff3(fn, get_transcriptome = FALSE, save_filtered_gff = FALSE,
    gene_name = NULL, gene_version = NULL, write_tr2g = FALSE,
    transcript_version = NULL) %>%
    arrange(gene)
  expect_equal(tr2g_no_gn, tr2g_expected[, -3])
  # With version number
  expect_equal(tr2g_gff3(fn, get_transcriptome = FALSE, save_filtered_gff = FALSE,
                         write_tr2g = FALSE) %>% arrange(gene),
    tr2g_expected_version)

  # Test error messages
  # transcript_id or gene_id is wrong
  expect_error(tr2g_gff3(fn, transcript_id = "foo", get_transcriptome = FALSE), 
               "Tags foo are absent")
  expect_error(tr2g_gff3(fn, gene_id = "bar", get_transcriptome = FALSE), 
               "Tags bar are absent")
  # transcript_id or gene_id is NULL
  expect_error(tr2g_gff3(fn, transcript_id = NULL, get_transcriptome = FALSE), 
               "transcript_id cannot be NULL.")
  expect_error(tr2g_gff3(fn, gene_id = NULL, get_transcriptome = FALSE), 
               "gene_id cannot be NULL.")
  # gene_name, gene_version, or transcript_version is wrong
  expect_warning(tr2g_gff3(fn, get_transcriptome = FALSE, save_filtered_gff = FALSE,
                           gene_name = "foo", gene_version = "bar",
                           write_tr2g = FALSE),
    "Tags foo, bar are absent")
  # Error when input file is not a gff3 file
  expect_error(tr2g_gff3(paste(toy_path, "fasta_test.fasta", sep = "/"),
                         get_transcriptome = FALSE),
    "file must be a GFF3 file.")
  rm(list = c("fn", "tr2g_no_vn", "tr2g_no_gn"))
})

test_that("Extract transcript and gene ID from FASTA file", {
  fn <- paste(toy_path, "fasta_test.fasta", sep = "/")
  expect_equal(tr2g_fasta(fn, write_tr2g = FALSE, save_filtered = FALSE), 
               fa_tr2g_expected)
  expect_error(tr2g_fasta(file = 3),
    "file must be a character vector with length 1")
  expect_error(tr2g_fasta(paste(toy_path, "gtf_test.gtf", sep = "/")),
    "file must be a FASTA file.")
  # Version number
  expect_equal(tr2g_fasta(fn, use_transcript_version = FALSE,
    use_gene_version = FALSE, write_tr2g = FALSE, save_filtered = FALSE),
  fa_tr2g_no_version)
  # Non-ENS IDs
  fn_dm <- paste(toy_path, "fasta_dm_test.fasta", sep = "/")
  m <- capture_messages(tr2g_fasta(fn_dm, use_transcript_version = TRUE, 
                                   write_tr2g = FALSE,
                                   save_filtered = FALSE))
  expect_match(m, "Version is not applicable.*", all = FALSE)
  tr2g_dm <- read.csv(paste(toy_path, "fa_tr2g_dm.csv", sep = "/"),
    header = TRUE, stringsAsFactors = FALSE) %>% as_tibble()
  expect_equal(tr2g_fasta(fn_dm, use_transcript_version = TRUE, 
                          write_tr2g = FALSE, save_filtered = FALSE), tr2g_dm)
  expect_equal(tr2g_fasta(fn_dm, use_transcript_version = FALSE, 
                          write_tr2g = FALSE, save_filtered = FALSE), tr2g_dm)
  rm(list = c("fn", "fn_dm", "tr2g_dm"))
})

test_that("Correct gene list output", {
  foo <- EC2gene(tr2g_toy, toy_path, verbose = FALSE)
  expect_equal(foo$EC_ind, EC2g_toy$EC_ind)
  expect_equal(foo$EC, EC2g_toy$EC)
  expect_equal(foo$gene, EC2g_toy$gene)
  rm(foo)
})

test_that("Correct bustools tr2g tsv file", {
  save_tr2g_bustools(tr2g_expected, file_save = "./tr2g.tsv")
  tr2g1 <- read.table(paste(toy_path, "tr2g_bustools.tsv", sep = "/"),
    header = FALSE, stringsAsFactors = FALSE)
  tr2g2 <- read.table("./tr2g.tsv", header = FALSE, stringsAsFactors = FALSE)
  expect_equal(tr2g1, tr2g2)
  file.remove("./tr2g.tsv")
  rm(list = c("tr2g1", "tr2g2"))
})

test_that("Correct sorting of transcripts", {
  file_use <- paste(toy_path, "gtf_test.gtf", sep = "/")
  tr2g <- tr2g_gtf(file = file_use, save_filtered_gtf = FALSE, 
                   get_transcriptome = FALSE, write_tr2g = FALSE,
                   transcript_version = NULL, gene_version = NULL)
  s <- sort_tr2g(tr2g, kallisto_out_path = toy_path)
  tx <- readLines(paste(toy_path, "transcripts.txt", sep = "/"))
  expect_equal(s$transcript, tx)
  rm(list = c("tr2g", "s", "tx"))
})
