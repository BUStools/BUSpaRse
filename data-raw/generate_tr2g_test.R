# Generate toy gtf, gff3, and fasta files for testing
library(biomartr)
library(Biostrings)
library(plyranges)
library(dplyr)
library(tidyr)
library(stringr)
getGTF(organism = "Homo sapiens", path = ".", db = "ensembl")
getGFF(organism = "Homo sapiens", path = ".", db = "ensembl")
# fasta file
getRNA(db = "ensembl", organism = "Homo sapiens", path = ".", reference = 95)
# For non-ENS gene symbols
getRNA(db = "ensembl", organism = "Drosophila melanogaster", path = ".", reference = 95)

gtf <- read_gff("./Homo_sapiens.GRCh38.96_ensembl.gtf.gz")
gff3 <- plyranges::read_gff3("./Homo_sapiens.GRCh38.96_ensembl.gff3.gz")
fa <- readDNAStringSet("./Homo_sapiens.GRCh38.ncrna.fa.gz")
fa_dm <- readDNAStringSet("./Drosophila_melanogaster.BDGP6.22.ncrna.fa.gz")

# Get a subset of the files for testing
# Randomly select 5 transcripts
genes <- unique(gff3$gene_id[gff3$type == "gene"])
set.seed(1)
genes_use <- sample(genes, 5)
gtf_sub <- gtf %>%
  filter(gene_id %in% genes_use)
tx_use <- unique(gtf_sub$transcript_id)
tx_use <- tx_use[!is.na(tx_use)]
write_gff(gtf_sub, "./inst/testdata/gtf_test.gtf")

gff3_sub <- gff3 %>%
  filter(gene_id %in% genes_use | transcript_id %in% tx_use)
write_gff3(gff3_sub, "./inst/testdata/gff3_test.gff3")

# Expected tr2g for GTF and GFF3 files
tr2g_expected <- mcols(gtf_sub) %>%
  as.data.frame() %>%
  select(transcript = transcript_id, gene = gene_id, gene_name) %>%
  filter(complete.cases(.)) %>%
  distinct() %>%
  arrange(gene)
write.csv(tr2g_expected, "./inst/testdata/tr2g_expected.csv", row.names = FALSE,
  quote = FALSE)
writeLines(tr2g_expected$transcript, "./inst/testdata/transcripts.txt", sep = "\n")

# Save in the bustools format
write.table(tr2g_expected[, c("transcript", "gene")],
  "./inst/testdata/tr2g_bustools.tsv", sep = "\t",
  row.names = FALSE, col.names = FALSE, quote = FALSE)

# With version number
tr2g_expected_version <- mcols(gtf_sub) %>%
  as.data.frame() %>%
  select(transcript = transcript_id, gene = gene_id, gene_name,
    gene_version, transcript_version) %>%
  filter(complete.cases(.)) %>%
  distinct() %>%
  unite("gene", gene, gene_version, sep = ".") %>%
  unite("transcript", transcript, transcript_version, sep = ".") %>%
  arrange(gene)
write.csv(tr2g_expected_version, "./inst/testdata/tr2g_expected_version.csv",
  row.names = FALSE, quote = FALSE)

# Use a subset of fasta file to test tr2g_fasta
fa_sub <- fa[1:5]
writeXStringSet(fa_sub, "./inst/testdata/fasta_test.fasta")
fa_tr2g_expected <- data.frame(transcript = str_extract(names(fa_sub), "^[a-zA-Z\\d-\\.]+"),
  gene = str_replace(names(fa_sub), "^.*gene:", "") %>%
    str_replace("\\s+.*$", ""),
  gene_name = str_replace(names(fa_sub), "^.*gene_symbol:", "") %>%
    str_replace("\\s+.*$", ""),
  stringsAsFactors = FALSE) %>%
  distinct()
# No version number
fa_tr2g_no_version <- fa_tr2g_expected %>%
  mutate(transcript = str_replace(transcript, "\\.\\d+$", ""),
    gene = str_replace(gene, "\\.\\d+$", ""))

write.csv(fa_tr2g_expected, file = "./inst/testdata/fa_tr2g_expected.csv",
  quote = FALSE, row.names = FALSE)
write.csv(fa_tr2g_no_version, file = "./inst/testdata/fa_tr2g_no_version.csv",
  quote = FALSE, row.names = FALSE)

# Non-ENS gene ID
fa_dm_sub <- fa_dm[1:5]
writeXStringSet(fa_dm_sub, "./inst/testdata/fasta_dm_test.fasta")
fa_tr2g_dm <- data.frame(transcript = str_extract(names(fa_dm_sub), "^[a-zA-Z\\d-\\.]+"),
  gene = str_replace(names(fa_dm_sub), "^.*gene:", "") %>%
    str_replace("\\s+.*$", ""),
  gene_name = str_replace(names(fa_dm_sub), "^.*gene_symbol:", "") %>%
    str_replace("\\s+.*$", ""),
  stringsAsFactors = FALSE) %>%
  distinct()
write.csv(fa_tr2g_dm, file = "./inst/testdata/fa_tr2g_dm.csv",
  quote = FALSE, row.names = FALSE)
