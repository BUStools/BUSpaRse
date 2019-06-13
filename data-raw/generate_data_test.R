# Generate fake data for unit testing
# Cases to be tested for:
# 1. Multiple UMIs for one EC, which has only one transcript
# 2. One UMI for one EC, which has one transcript
# 3. One UMI for one EC, which has multiple transcripts of the same gene
# 4. Multiple UMIs for one EC, which has multiple transcripts of the same gene
# 5. One UMI for one EC, which has transcripts of different genes
# 6. Multiple UMIs for one EC, which has transcripts of different genes
# 7. Some UMIs exclusively for one transcript, and some other UMIs exclusively
# for another transcript of the same gene
# 8. Same UMI maps to different ECs that don't have intersecting genes
# 9. Same UMI maps to different ECs with intersection that has multiple genes
# 10. Same UMI maps to different ECs with intersection with 1 gene
library(dplyr)
library(tibble)

tr2g_toy <- data.frame(
  gene = c("meow", "kitty", "kitty", "kitty", "purr"),
  transcript = c("meow", "kitty", "kedi", "qitah", "purr"),
  stringsAsFactors = FALSE
)

ECs_toy <- tibble(
  EC_ind = 0:9,
  EC = c(as.character(0:4), "1,2", "2,3", "0,1", "0,4", "0,2,4")
)
EC2g_toy <- ECs_toy %>% 
  mutate(EC = lapply(EC, function(x) as.numeric(strsplit(x, ",")[[1]])),
         gene = list("meow", "kitty", "kitty", "kitty", "purr",
                     "kitty", "kitty", c("kitty", "meow"), c("meow", "purr"),
                     c("kitty", "meow", "purr")))

output_sorted_toy <- data.frame(
  barcode = c("A", "A", "B", "C", "D", "D", "E", "F", "F", "G", "G", "H", "H", 
              "I", "I", "J", "J", "K"),
  UMI = c("a1", "a2", "b", "c", "d1", "d2", "e", "f1", "f2", "g1", "g2", "h", "h", 
          "i", "i", "j", "j", "k"),
  EC_index = c("0", "0", "0", "5", "5", "6", "7", "7", "7", "2", "3", "1", "8",
               "9", "8", "4", "8", "1"),
  counts = 1,
  stringsAsFactors = FALSE
)
whitelist <- LETTERS[1:10]

# The expected gene_count matrix
rowinds <- c(1, 1, 2, 2, 1, 2, 1, 2, 2, 1, 3, 3)
colinds <- c(1, 2, 3, 4, 5, 5, 6, 6, 7, 8, 8, 9)
values <- c(2, 1, 1, 2, 0.5, 0.5, 1, 1, 2, 0.5, 0.5, 1)
expected_mat <- Matrix::sparseMatrix(i = rowinds, j = colinds, x = values,
                                     dimnames = list(c("meow", "kitty", "purr"),
                                                     LETTERS[c(1:7, 9:10)]))
expected_mat_full <- Matrix::sparseMatrix(i = c(rowinds, 2), j = c(colinds, 10),
                                          x = c(values, 1),
                                          dimnames = list(c("meow", "kitty", "purr"),
                                                          LETTERS[c(1:7, 9:11)]))
# The expected gene_count matrix in single gene wide
rowinds <- c(1,1,2,2,2,3)
colinds <- 1:6
values <- c(2,1,1,2,2,1)
expected_single <- Matrix::sparseMatrix(i = rowinds, j = colinds, x = values,
                                        dimnames = list(c("meow", "kitty", "purr"),
                                                        LETTERS[c(1:4, 7, 10)]))

# The expected TCC matrix
rowinds <- c(1,1,5,5,6,7,7,2,3,8,4)
rowinds_full <- c(1,1,6,6,7,8,8,3,4,9,5,2)
colinds <- c(1,2,3,4,4,5,6,7,7,8,9)
values <- c(2,1,1,1,1,1,2,1,1,1,1)
expected_tcc <- Matrix::sparseMatrix(i = rowinds, j = colinds, x = values,
                                     dimnames = list(ECs_toy$EC[c(1, 3:9)],
                                                     LETTERS[c(1:7, 9:10)]))
# Without whitelist
expected_tcc_full <- Matrix::sparseMatrix(i = rowinds_full, j = c(colinds, 10), 
                                          x = c(values, 1),
                                          dimnames = list(ECs_toy$EC[1:9],
                                                          LETTERS[c(1:7, 9:11)]))

# Save data
write.table(ECs_toy, file = "./inst/testdata/matrix.ec", quote = FALSE, 
            row.names = FALSE, col.names = FALSE, sep = "\t")
write.table(output_sorted_toy, file = "./inst/testdata/output.sorted.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
save(tr2g_toy, whitelist, expected_mat, expected_mat_full, expected_single,
     expected_tcc, expected_tcc_full, EC2g_toy,
     file = "./inst/testdata/toy_example.RData")
