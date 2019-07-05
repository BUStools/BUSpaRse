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

library(GenomicRanges)
library(plyranges)
library(tibble)
library(stringi)
# Make toy genomic ranges---------------------
gr_df <- tribble(
  ~seqnames, ~start, ~end, ~strand, ~gene_id, ~transcript_id, ~exon_number,
  "chr1",    5,      30,   "+",     "A",      "A1",           1,
  "chr1",    61,     90,   "+",     "A",      "A1",           2,
  "chr1",    141,    160,  "+",     "A",      "A1",           3,
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
  "chr8",    141,    160,  "+",     "H",      "H2",           2
)

gr <- makeGRangesFromDataFrame(gr_df, keep.extra.columns = TRUE)

# Minus strand
gr_minus <- gr
strand(gr_minus) <- "-"
gr_minus$gene_id <- paste0(gr$gene_id, "m")
gr_minus$transcript_id <- paste0(gr$transcript_id, "m")
gr_minus$exon_number <- c(rep(3:1, 3), 1, 4:1, 3:1, 2:1, rep(3:1, 2), 3:1, 2:1)
gr_full <- c(gr, gr_minus)
gr_full$gene_version <- 1
gr_full$transcript_version <- 1
gr_full$type <- "exon"
# Save GTF file
write_gff(gr_full, "./inst/testdata/test_velocity.gtf")

# Make toy genome--------------------------
set.seed(123)
flank5p <- sample(c("A", "T", "C", "G"), 4, replace = TRUE)
set.seed(456)
flank3p <- sample(c("A", "T", "C", "G"), 5, replace = TRUE)
exon1 <- rep("A", 26)
intron1 <- rep("T", 30)
exon2 <- rep("C", 30)
intron2 <- rep("T", 50)
exon3 <- rep("G", 20)
gene1 <- paste(c(flank5p, exon1, intron1, exon2, intron2, exon3, flank3p),
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

gn <- setNames(c(gene1, geneB, geneC, geneD, geneE, rep(gene1, 3)),
               paste0("chr", 1:8))
gn <- DNAStringSet(gn, use.names = TRUE)
writeXStringSet(gn, "./inst/testdata/test_velocity_genome.fa")

# Make toy transcripts--------------------------
tx1 <- paste(c(exon1, exon2, exon3), collapse = "")
txB <- paste(c(exon1, exon2b, exon3), collapse = "")
txC <- paste(c(exon1, exon2, exon3c), collapse = "")
txD <- paste(rep("A", 86), collapse = "")
txE <- paste(c(exon1, exon2e, exon2e2, exon3), collapse = "")
txF2 <- paste(c(exon1, exon3), collapse = "")
txG2 <- paste(c(rep("A", 16), exon2, exon3), collapse = "")
txH2 <- paste(c(exon2, exon3), collapse = "")
plus_tx <- setNames(c(tx1, txB, txC, txD, txE, tx1, txF2, tx1, txG2, tx1, txH2),
                    unique(gr$transcript_id))
plus_tx <- DNAStringSet(plus_tx, use.names = TRUE)
minus_tx <- reverseComplement(plus_tx)
names(minus_tx) <- unique(gr_minus$transcript_id)
writeXStringSet(c(plus_tx, minus_tx), "./inst/testdata/test_velocity_tx.fa")

# Make toy intron fasta, with L=11, so flanking region is 10 or shorter---------
## With truncation, with isoforms separate---------------------
ifl1 <- paste(c(rep("A", 10), intron1, rep("C", 10)), collapse = "")
ifl2 <- paste(c(rep("C", 10), intron2, rep("G", 10)), collapse = "")
iflB1 <- paste(c(rep("A", 10), intron1, exon2b), collapse = "")
iflB2 <- paste(c(exon2b, intron2b, rep("G", 10)), collapse = "")
iflC2 <- paste(c(rep("C", 10), intron2, exon3c), collapse = "")
iflE2 <- paste(c(exon2e, intron2e, exon2e2), collapse = "")
iflE3 <- paste(c(exon2e2, intron2, rep("G", 10)), collapse = "")
iflF2 <- paste(c(rep("A", 10), intron1, exon2, intron2, rep("G", 10)), collapse = "")
iflG2 <- paste(c(rep("A", 20), intron1, rep("C", 10)), collapse = "")
ins_trunc_plus <- setNames(c(ifl1, ifl2, 
                             iflB1, iflB2, 
                             ifl1, iflC2,
                             iflB1, iflE2, iflE3,
                             ifl1, ifl2, iflF2,
                             ifl1, ifl2, iflG2, ifl2,
                             ifl1, ifl2, ifl2),
                           c("A1-I1", "A1-I2",
                             "B1-I1", "B1-I2",
                             "C1-I1", "C1-I2",
                             "E1-I1", "E1-I2", "E1-I3",
                             "F1-I1", "F1-I2", "F2-I1",
                             "G1-I1", "G1-I2", "G2-I1", "G2-I2",
                             "H1-I1", "H1-I2", "H2-I1"))
ins_trunc_plus <- DNAStringSet(ins_trunc_plus, use.names = TRUE)
ins_trunc_minus <- setNames(reverseComplement(ins_trunc_plus),
                            c("A1m-I2", "A1m-I1",
                              "B1m-I2", "B1m-I1",
                              "C1m-I2", "C1m-I1",
                              "E1m-I3", "E1m-I2", "E1m-I1",
                              "F1m-I2", "F1m-I1", "F2m-I1",
                              "G1m-I2", "G1m-I1", "G2m-I2", "G2m-I1",
                              "H1m-I2", "H1m-I1", "H2m-I1"))
writeXStringSet(c(ins_trunc_plus, ins_trunc_minus), 
                "./inst/testdata/velocity_introns_trunc.fa")

## No truncation; pretend short exons don't exist-----------------
iflBf <- paste(c(rep("A", 10), intron1, exon2b, intron2, rep("G", 10)), collapse = "")
iflEf <- paste(c(rep("A", 10), intron1, exon2e, intron2e, exon2e2, intron2, rep("G", 10)),
               collapse = "")
ins_inc_plus <- setNames(c(ifl1, ifl2, 
                           iflBf, 
                           ifl1, iflC2,
                           iflEf,
                           ifl1, ifl2, iflF2,
                           ifl1, ifl2, iflG2, ifl2,
                           ifl1, ifl2, ifl2),
                         c("A1-I1", "A1-I2",
                           "B1-I1",
                           "C1-I1", "C1-I2",
                           "E1-I1",
                           "F1-I1", "F1-I2", "F2-I1",
                           "G1-I1", "G1-I2", "G2-I1", "G2-I2",
                           "H1-I1", "H1-I2", "H2-I1"))
ins_inc_plus <- DNAStringSet(ins_inc_plus, use.names = TRUE)
ins_inc_minus <- setNames(reverseComplement(ins_inc_plus),
                          c("A1m-I2", "A1m-I1",
                            "B1m-I1",
                            "C1m-I2", "C1m-I1",
                            "E1m-I1", 
                            "F1m-I2", "F1m-I1", "F2m-I1",
                            "G1m-I2", "G1m-I1", "G2m-I2", "G2m-I1",
                            "H1m-I2", "H1m-I1", "H2m-I1"))
writeXStringSet(c(ins_inc_plus, ins_inc_minus), 
                "./inst/testdata/velocity_introns_inc.fa")

## Collapse isoforms--------------------
ins_coll_plus <- setNames(c(ins_trunc_plus[1:11],
                            rep(c(ifl1, ifl2), 2)),
                          c("A-I1", "A-I2",
                            "B-I1", "B-I2",
                            "C-I1", "C-I2",
                            "E-I1", "E-I2", "E-I3",
                            "F-I1", "F-I2",
                            "G-I1", "G-I2",
                            "H-I1", "H-I2"))
ins_coll_plus <- DNAStringSet(ins_coll_plus, use.names = TRUE)
ins_coll_minus <- setNames(reverseComplement(ins_coll_plus),
                           c("A-I1", "A-I2",
                             "B-I1", "B-I2",
                             "C-I1", "C-I2",
                             "E-I1", "E-I2", "E-I3",
                             "F-I1", "F-I2",
                             "G-I1", "G-I2",
                             "H-I1", "H-I2"))
writeXStringSet(c(ins_coll_plus, ins_coll_minus),
                "./inst/testdata/velocity_introns_collapse.fa")
