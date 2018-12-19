context("Test function to convert transcript to gene")
library(testthat)

# Load toy example for testing
toy_path <- system.file("testdata", package = "BUSpaRse")
load(paste(toy_path, "toy_example.RData", sep = "/"))

test_that("Correct data set name for biomart", {
  error_use <- "Please use the Latin binomial convention for species rather than the colloquial name."
  expect_equal(species2dataset("Homo sapiens"), "hsapiens_gene_ensembl")
  expect_equal(species2dataset("Felis catus"), "fcatus_gene_ensembl")
  expect_error(species2dataset("cats are so cute", error_use))
  expect_error(species2dataset("mouse", error_use))
})

test_that("Correct gene list output", {
  expect_equal(EC2gene(tr2g_toy, toy_path, verbose = FALSE),
               genes_toy)
})
