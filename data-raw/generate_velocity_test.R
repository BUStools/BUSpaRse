# Here I make toy datasets to test RNA velocity files.
# Cases to consider:
# 1. Multiple exons, all longer than L-1
# 2. One exon that is too short; it's not the first or the last exon
# 3. One exon that is too short; it is the last exon
# 4. Two adjacent exons that are too short
# 5. Only one exon
# 6. Two isoforms, due to an alternatively spliced exon
# 7. Two isoforms, due to an alternative splice site
# 8. Two isoforms, due to an alternative transcription start site
# 9. Some of the special cases, on the minus strand
# 10. More than one gene on the same chromosome
# Use L = 16 to see whether we get the correct results for exons between
# L-1 and 2(L-1)

library(GenomicRanges)
library(Biostrings)
library(plyranges)
library(tidyverse)
library(stringi)

# Make toy genomic ranges---------------------
gr_df <- tribble(
  ~seqnames, ~start, ~end, ~strand, ~gene_id, ~transcript_id, ~exon_number,
  "chr1",    5,      30,   "+",     "A",      "A1",           1,
  "chr1",    61,     90,   "+",     "A",      "A1",           2,
  "chr1",    141,    160,  "+",     "A",      "A1",           3,
  "chr1",    205,    230,  "+",     "Aa",     "A1a",          1,
  "chr1",    261,    290,  "+",     "Aa",     "A1a",          2,
  "chr1",    341,    360,  "+",     "Aa",     "A1a",          3,
  "chr2",    5,      30,   "+",     "B",      "B1",           1,
  "chr2",    61,     65,   "+",     "B",      "B1",           2,
  "chr2",    141,    160,  "+",     "B",      "B1",           3,
  "chr3",    5,      30,   "+",     "C",      "C1",           1,
  "chr3",    61,     90,   "+",     "C",      "C1",           2,
  "chr3",    141,    145,  "+",     "C",      "C1",           3,
  "chr4",    5,      90,   "+",     "D",      "D1",           1,
  "chr5",    5,      30,   "+",     "E",      "E1",           1,
  "chr5",    61,     65,   "+",     "E",      "E1",           2,
  "chr5",    85,     90,   "+",     "E",      "E1",           3,
  "chr5",    141,    160,  "+",     "E",      "E1",           4,
  "chr6",    5,      30,   "+",     "F",      "F1",           1,
  "chr6",    61,     90,   "+",     "F",      "F1",           2,
  "chr6",    141,    160,  "+",     "F",      "F1",           3,
  "chr6",    5,      30,   "+",     "F",      "F2",           1,
  "chr6",    141,    160,  "+",     "F",      "F2",           2,
  "chr7",    5,      30,   "+",     "G",      "G1",           1,
  "chr7",    61,     90,   "+",     "G",      "G1",           2,
  "chr7",    141,    160,  "+",     "G",      "G1",           3,
  "chr7",    5,      20,   "+",     "G",      "G2",           1,
  "chr7",    61,     90,   "+",     "G",      "G2",           2,
  "chr7",    141,    160,  "+",     "G",      "G2",           3,
  "chr8",    5,      30,   "+",     "H",      "H1",           1,
  "chr8",    61,     90,   "+",     "H",      "H1",           2,
  "chr8",    141,    160,  "+",     "H",      "H1",           3,
  "chr8",    61,     90,   "+",     "H",      "H2",           1,
  "chr8",    141,    160,  "+",     "H",      "H2",           2,
  "chr10",   5,      20,   "+",     "K",      "K1",           1
)

gr_df <- gr_df %>%
  mutate(exon_id = paste(transcript_id, exon_number, sep = "-"))

# Minus strand
gr_minus <- gr_df
gr_minus$strand <- "-"
gr_minus$gene_id <- paste0(gr_df$gene_id, "m")
gr_minus$transcript_id <- paste0(gr_df$transcript_id, "m")
gr_minus$exon_number <- c(rep(3:1, 4), 1, 4:1, 3:1, 2:1, rep(3:1, 2), 3:1, 2:1, 1)
gr_minus <- gr_minus %>% 
  arrange(seqnames, desc(start), desc(end))
gr_minus <- gr_minus[c(1:6, 8:34, 7),]
gr <- makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE)
gr_m <- makeGRangesFromDataFrame(gr_minus, keep.extra.columns = TRUE)
gr_full <- c(gr, gr_m)
gr_full$type <- "exon"
# Save GTF file
write_gff(gr_full, "./inst/testdata/velocity_annot.gtf")
tr2g_tx <- mcols(gr_full) %>%
  as.data.frame() %>%
  dplyr::select(transcript = transcript_id, gene = gene_id) %>%
  dplyr::filter(!str_detect(transcript, "K")) %>% 
  distinct()

# Make toy genome--------------------------
set.seed(123)
flank5p <- sample(c("A", "T", "C", "G"), 4, replace = TRUE)
set.seed(456)
flank3p <- sample(c("A", "T", "C", "G"), 5, replace = TRUE)
set.seed(789)
intergene <- sample(c("A", "T", "C", "G"), 44, replace = TRUE)
exon1 <- rep("A", 26)
intron1 <- rep("T", 30)
exon2 <- rep("C", 30)
intron2 <- rep("T", 50)
exon3 <- rep("G", 20)
chr1 <- paste(c(flank5p, exon1, intron1, exon2, intron2, exon3,
  intergene, exon1, intron1, exon2, intron2, exon3, flank3p),
collapse = "")
exon2b <- rep("C", 5)
intron2b <- rep("T", 75)
geneB <- paste(c(flank5p, exon1, intron1, exon2b, intron2b, exon3, flank3p),
  collapse = "")

exon3c <- rep("G", 5)
geneC <- paste(c(flank5p, exon1, intron1, exon2, intron2, exon3c, flank3p),
  collapse = "")

geneD <- paste(c(flank5p, rep("A", 86), flank3p), collapse = "")

exon2e <- rep("C", 5)
intron2e <- rep("T", 19)
exon2e2 <- rep("C", 6)
geneE <- paste(c(flank5p, exon1, intron1, exon2e, intron2e, exon2e2, intron2, exon3,
  flank3p), collapse = "")

gn <- setNames(c(chr1, geneB, geneC, geneD, geneE, rep(chr1, 4)),
  paste0("chr", 1:9))
gn <- DNAStringSet(gn, use.names = TRUE)
seqinfo(gn) <- Seqinfo(paste("chr", 1:9), isCircular = rep(FALSE, 9))
writeXStringSet(gn, "./inst/testdata/velocity_genome.fa")

# Make toy transcripts--------------------------
tx1 <- paste(c(exon1, exon2, exon3), collapse = "")
txB <- paste(c(exon1, exon2b, exon3), collapse = "")
txC <- paste(c(exon1, exon2, exon3c), collapse = "")
txD <- paste(rep("A", 86), collapse = "")
txE <- paste(c(exon1, exon2e, exon2e2, exon3), collapse = "")
txF2 <- paste(c(exon1, exon3), collapse = "")
txG2 <- paste(c(rep("A", 16), exon2, exon3), collapse = "")
txH2 <- paste(c(exon2, exon3), collapse = "")
plus_tx <- setNames(c(tx1, tx1, txB, txC, txD, txE, tx1, txF2, tx1, txG2, tx1, txH2),
  unique(gr$transcript_id)[1:12])
plus_tx <- DNAStringSet(plus_tx, use.names = TRUE)
minus_tx <- reverseComplement(plus_tx)
names(minus_tx) <- unique(gr_m$transcript_id)[1:12]
tx <- c(plus_tx, minus_tx)
writeXStringSet(tx, "./inst/testdata/velocity_tx.fa")
writeLines(names(tx), "./inst/testdata/velocity_cdna_tx.txt")

# Make toy intron fasta, with L=11, so flanking region is 15 or shorter---------
## Don't collapse isoforms--------------------
ifl1 <- paste(c(rep("A", 15), intron1, rep("C", 15)), collapse = "")
ifl2 <- paste(c(rep("C", 15), intron2, rep("G", 15)), collapse = "")
iflB1 <- paste(c(rep("A", 15), intron1, exon2b), collapse = "")
iflB2 <- paste(c(exon2b, intron2b, rep("G", 15)), collapse = "")
iflC2 <- paste(c(rep("C", 15), intron2, exon3c), collapse = "")
iflE2 <- paste(c(exon2e, intron2e, exon2e2), collapse = "")
iflE3 <- paste(c(exon2e2, intron2, rep("G", 15)), collapse = "")
iflF2 <- paste(c(rep("A", 15), intron1, exon2, intron2, rep("G", 15)), collapse = "")
iflG2 <- paste(c(rep("A", 25), intron1, rep("C", 15)), collapse = "")
ins_plus <- setNames(c(ifl1, ifl2,
  ifl1, ifl2,
  iflB1, iflB2,
  ifl1, iflC2,
  iflB1, iflE2, iflE3,
  ifl1, ifl2, iflF2,
  ifl1, ifl2, iflG2, ifl2,
  ifl1, ifl2, ifl2),
c("A1-I", "A1-I1",
  "A1a-I", "A1a-I1",
  "B1-I", "B1-I1",
  "C1-I", "C1-I1",
  "E1-I", "E1-I1", "E1-I2",
  "F1-I", "F1-I1", "F2-I",
  "G1-I", "G1-I1", "G2-I", "G2-I1",
  "H1-I", "H1-I1", "H2-I"))
ins_plus <- DNAStringSet(ins_plus, use.names = TRUE)
ins_minus <- setNames(reverseComplement(ins_plus),
  c("A1m-I", "A1m-I1",
    "A1am-I", "A1am-I1",
    "B1m-I", "B1m-I1",
    "C1m-I", "C1m-I1",
    "E1m-I", "E1m-I1", "E1m-I2",
    "F1m-I", "F1m-I1", "F2m-I",
    "G1m-I", "G1m-I1", "G2m-I", "G2m-I1",
    "H1m-I", "H1m-I1", "H2m-I"))
ins <- c(ins_plus, ins_minus)
ins <- ins[sort(names(ins))]
writeXStringSet(c(tx, ins), "./inst/testdata/velocity_introns.fa")
writeLines(names(ins), "./inst/testdata/velocity_introns_capture.txt")
tr2g_intron <- tibble(transcript = names(ins),
  gene = str_remove(transcript, "-I(\\d+)?") %>%
    str_remove("\\d"))
write_tsv(rbind(tr2g_tx, tr2g_intron), "./inst/testdata/velocity_tr2g.tsv",
  col_names = FALSE)

## Collapse isoforms--------------------
ins_coll_plus <- setNames(c(ins_plus[1:13],
  rep(c(ifl1, ifl2), 2)),
c("A-I", "A-I1",
  "Aa-I", "Aa-I1",
  "B-I", "B-I1",
  "C-I", "C-I1",
  "E-I", "E-I1", "E-I2",
  "F-I", "F-I1",
  "G-I", "G-I1",
  "H-I", "H-I1"))
ins_coll_plus <- DNAStringSet(ins_coll_plus, use.names = TRUE)
ins_coll_minus <- setNames(reverseComplement(ins_coll_plus),
  c("Am-I", "Am-I1",
    "Aam-I", "Aam-I1",
    "Bm-I", "Bm-I1",
    "Cm-I", "Cm-I1",
    "Em-I", "Em-I1", "Em-I2",
    "Fm-I", "Fm-I1",
    "Gm-I", "Gm-I1",
    "Hm-I", "Hm-I1"))
ins_coll <- c(ins_coll_plus, ins_coll_minus)
ins_coll <- ins_coll[sort(names(ins_coll))]
writeXStringSet(c(tx, ins_coll), "./inst/testdata/velocity_introns_coll.fa")
writeLines(names(ins_coll), "./inst/testdata/velocity_introns_coll_capture.txt")
tr2g_intron_coll <- tibble(transcript = names(ins_coll),
  gene = str_remove(transcript, "-.+"))
write_tsv(rbind(tr2g_tx, tr2g_intron_coll), "./inst/testdata/velocity_tr2g_coll.tsv",
  col_names = FALSE)

# Make toy exon-exon junctions-----------------
j1 <- paste(c(exon1[1:15], exon2[1:15]), collapse = "")
j1b <- paste(c(exon1[1:15], exon2b), collapse = "")
j2 <- paste(c(exon2[1:15], exon3[1:15]), collapse = "")
j2b <- paste(c(exon2b, exon3[1:15]), collapse = "")
j2c <- paste(c(exon2[1:15], exon3c), collapse = "")
j1e <- paste(c(exon1[1:15], exon2e), collapse = "")
j2e <- paste(c(exon2e, exon2e2), collapse = "")
j3e <- paste(c(exon2e2, exon3[1:15]), collapse = "")
jf2 <- paste(c(exon1[1:15], exon3[1:15]), collapse = "")

junc_plus <- setNames(c(j1, j2,
  j1, j2,
  j1b, j2b,
  j1, j2c,
  j1e, j2e, j3e,
  j1, j2, jf2,
  j1, j2, j1, j2,
  j1, j2, j2),
c("A1-J", "A1-J1",
  "A1a-J", "A1a-J1",
  "B1-J", "B1-J1",
  "C1-J", "C1-J1",
  "E1-J", "E1-J1", "E1-J2",
  "F1-J", "F1-J1", "F2-J",
  "G1-J", "G1-J1", "G2-J", "G2-J1",
  "H1-J", "H1-J1", "H2-J"))
junc_plus <- DNAStringSet(junc_plus, use.names = TRUE)
junc_minus <- setNames(reverseComplement(junc_plus),
  c("A1m-J", "A1m-J1",
    "A1am-J", "A1am-J1",
    "B1m-J", "B1m-J1",
    "C1m-J", "C1m-J1",
    "E1m-J", "E1m-J1", "E1m-J2",
    "F1m-J", "F1m-J1", "F2m-J",
    "G1m-J", "G1m-J1", "G2m-J", "G2m-J1",
    "H1m-J", "H1m-J1", "H2m-J"))
juncs <- c(junc_plus, junc_minus)
writeXStringSet(c(juncs, tx[c("D1", "D1m")], ins), "./inst/testdata/junctions.fa")
# Should capture all exon-exon junctions and transcripts without introns
writeLines(c(names(juncs), "D1", "D1m"), "./inst/testdata/junction_capture.txt")
tr2g_junc <- tibble(transcript = c(names(juncs), "D1", "D1m"),
  gene = str_remove_all(transcript, "-J(\\d+)?") %>%
    str_remove_all("\\d"))
write_tsv(rbind(tr2g_junc, tr2g_intron), "./inst/testdata/junction_tr2g.tsv",
  col_names = FALSE)

## Collapse isoforms------------
writeXStringSet(c(juncs, tx[c("D1", "D1m")], ins_coll), "./inst/testdata/junction_coll.fa")
write_tsv(rbind(tr2g_junc, tr2g_intron_coll), "./inst/testdata/junction_tr2g_coll.tsv",
  col_names = FALSE)
