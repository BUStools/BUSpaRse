# BUSpaRse 1.1.4 (2020-04-23)
* Now `tr2g_*` functions can also extract transcriptomes.
* By default, `tr2g_*` functions write the `tr2g.tsv` to disk.
* Write filtered GTF and GFF3 files.
* List of available gene and transcript biotypes from Ensembl can be queried by `data()`.
* RefSeq GFF3 annotations can now be used in `tr2g_gff3`.

# BUSpaRse 1.1.3 (2020-04-06)
* Allow filtering gene and transcript biotypes in tr2g and get_velocity_files
* Allow removing scaffolds and haplotypes from genome annotations for tr2g

# BUSpaRse 1.1.2 (2020-02-06)
* Removed RcppParallel and data.table dependencies, and hence all C++ multithreading
* Added knee_plot function
* Changed default in read_count_output to tcc = FALSE

# BUSpaRse 1.1.1 (2020-01-29)
* Sort GRanges for RNA velocity to make sure that exons are in the right order;
`split` no longer sorts within each element of GRangesList.
* Export subset_annot, which makes sure that all chromosomes in the annotation
are present in the genome prior to transcriptome extraction.

# BUSpaRse 0.99.25 (2019-09-11)
* Added message to indicate when get_velocity_files is extracting exon-exon junctions.
* Restored "separate" to be default isoform_action in get_velocity_files.

# BUSpaRse 0.99.24 (2019-09-06)
* When sorting tr2g from file, now the file must be formatted in ways required by bustools.

# BUSpaRse 0.99.23 (2019-09-06)
* Previous two version bumps did not accompany changes; those were used to trigger rebuilds on Bioconductor.
* Added the functionality to use L-1 (L is read length) bases around exon-exon junction to better distinguish between spliced and unspliced transcripts for RNA velocity.
* Fixed serious problem with get_velocity_file that counted reads mapping to exons of length between L-1 and 2(L-1) as from unspliced transcript. This was done by an reimplementation of the method to get flanked intronic ranges.
* Changed default isoform_action to "collapse".
* Make sure that all transcripts in tr2g.tsv from get_velocity_files are in the transcriptome.

# BUSpaRse 0.99.20 (2019-08-26)
* Addressed Bioconductor review.

# BUSpaRse 0.99.19 (2019-07-23)
* Finished get_velocity_files to generate files required for kallisto | bustools
RNA velocity.

# BUSpaRse 0.99.0 (2019-06-21)
* Submitted package to Bioconductor for review
