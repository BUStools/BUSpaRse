# BUSpaRse 0.99.25 (2019-09-11)
* Added message to indicate when get_velocity_files is extracting exon-exon junctions.

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
