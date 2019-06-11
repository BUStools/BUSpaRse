context("Test make_sparse_matrix function")
library(testthat)

# Load toy example for testing
toy_path <- system.file("testdata", package = "BUSpaRse")
load(paste(toy_path, "toy_example.RData", sep = "/"))

test_that("Check file name input", {
  expect_error(make_sparse_matrix(paste0(toy_path, "/matrix.ec"),
                                  genes_toy, 11, 3, verbose = FALSE),
            "Argument fn must point to a text file. Please run bustools text.")
  expect_error(make_sparse_matrix("foo.txt", genes_toy, 11, 3, 
                                  verbose = FALSE),
               "(No such file or directory)|(cannot find)")
})

test_that("Check for correct gene count matrix", {
  out_fn <- paste0(toy_path, "/output.sorted.txt")
  # With whitelist
  m <- make_sparse_matrix(out_fn, tr2g_toy, 10, 3, whitelist = whitelist,
                          gene_count = TRUE, TCC = FALSE, verbose = FALSE,
                          ncores = 2)
  expect_equal(dim(m), dim(expected_mat))
  # Reoroder rows and columns
  m <- m[rownames(expected_mat), colnames(expected_mat)]
  expect_equal(m, expected_mat)
  # Without whitelist
  m2 <- make_sparse_matrix(out_fn, tr2g_toy, 11, 3, 
                           gene_count = TRUE, TCC = FALSE, verbose = FALSE,
                           ncores = 2)
  expect_equal(dim(m2), dim(expected_mat_full))
  # Reorder
  m2 <- m2[rownames(expected_mat_full), colnames(expected_mat_full)]
  expect_equal(m2, expected_mat_full)
  rm(list = c("m", "m2"))
})

test_that("Check for correct TCC matrix", {
  out_fn <- paste0(toy_path, "/output.sorted.txt")
  # With whitelist
  m <- make_sparse_matrix(out_fn, est_ncells = 10, est_ngenes = 9, 
                          whitelist = whitelist,
                          gene_count = FALSE, TCC = TRUE, verbose = FALSE,
                          ncores = 2)
  expect_equal(dim(m), dim(expected_tcc))
  # Reoroder rows and columns
  m <- m[rownames(expected_tcc), colnames(expected_tcc)]
  expect_equal(m, expected_tcc)
  # Without whitelist
  m2 <- make_sparse_matrix(out_fn, 
                           est_ncells = 11, est_ngenes = 9, 
                           gene_count = FALSE, TCC = TRUE, verbose = FALSE,
                           ncores = 2)
  expect_equal(dim(m2), dim(expected_tcc_full))
  # Reorder
  m2 <- m2[rownames(expected_tcc_full), colnames(expected_tcc_full)]
  expect_equal(m2, expected_tcc_full)
  rm(list = c("m", "m2"))
})
