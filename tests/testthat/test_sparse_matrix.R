context("Test make_sparse_matrix function")
library(testthat)

# Load toy example for testing
toy_path <- system.file("testdata", package = "BUSpaRse")
load(paste(toy_path, "toy_example.RData", sep = "/"))

test_that("Check file name input", {
  expect_error(make_sparse_matrix(paste0(toy_path, "/matrix.ec"),
                                  genes_toy, 11, 3, display_progress = FALSE),
            "Argument fn must point to a text file. Please run bustools text.")
  expect_error(make_sparse_matrix("foo.txt", genes_toy, 11, 3, 
                                  display_progress = FALSE),
               "(No such file or directory)|(cannot find)")
})

test_that("Check for correct matrix", {
  out_fn <- paste0(toy_path, "/output.sorted.txt")
  # With whitelist
  m <- make_sparse_matrix(out_fn, genes_toy, 10, 3, whitelist = whitelist,
                          display_progress = FALSE)
  expect_equal(dim(m), dim(expected_mat))
  # Reoroder rows and columns
  m <- m[rownames(expected_mat), colnames(expected_mat)]
  expect_equal(m, expected_mat)
  # Without whitelist
  m2 <- make_sparse_matrix(out_fn, genes_toy, 11, 3, display_progress = FALSE)
  expect_equal(dim(m2), dim(expected_mat_full))
  # Reorder
  m2 <- m2[rownames(expected_mat_full), colnames(expected_mat_full)]
  expect_equal(m2, expected_mat_full)
})