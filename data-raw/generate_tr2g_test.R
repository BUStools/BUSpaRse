# Generate toy gtf, gff3, and fasta files for testing
library(biomartr)
library(Biostrings)
library(plyranges)
library(dplyr)
library(tidyr)
getGTF(organism = "Homo sapiens", path = ".")
getGFF(organism = "Homo sapiens", path = ".")
# fasta file
getRNA(db = "ensembl", organism = "Homo sapiens", path = ".")

gtf <- read_gff("./Homo_sapiens.GRCh38.95_ensembl.gtf.gz")
gff3 <- plyranges::read_gff3("./Homo_sapiens.GRCh38.95_ensembl.gff3.gz")
fa <- readDNAStringSet("./Homo_sapiens.GRCh38.ncrna.fa.gz")

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
  select(gene = gene_id, transcript = transcript_id, gene_name) %>% 
  filter(complete.cases(.)) %>% 
  distinct() %>% 
  arrange(gene)
write.csv(tr2g_expected, "./inst/testdata/tr2g_expected.csv", row.names = FALSE,
          quote = FALSE)
# With version number
tr2g_expected_version <- mcols(gtf_sub) %>% 
  as.data.frame() %>% 
  select(gene = gene_id, transcript = transcript_id, gene_name,
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
fa_tr2g_expected <- data.frame(gene = c("ENSG00000222870.1", "ENSG00000252023.1",
                                        "ENSG00000252039.1", "ENSG00000223198.1",
                                        "ENSG00000207340.1"),
                               transcript = c("ENST00000410938.1", "ENST00000516214.1",
                                              "ENST00000516230.1", "ENST00000411266.1",
                                              "ENST00000384610.1"),
                               gene_name = c("RNU6-1328P", "RNU6-581P", "RNU6-287P",
                                             "RNU2-22P", "RNVU1-1"))
write.csv(fa_tr2g_expected, file = "./inst/testdata/fa_tr2g_expected.csv",
          quote = FALSE, row.names = FALSE)
