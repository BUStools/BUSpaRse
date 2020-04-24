<!-- badges: start -->
[![Travis build status](https://travis-ci.com/BUStools/BUSpaRse.svg?branch=master)](https://travis-ci.com/BUStools/BUSpaRse)
[![BioC status](http://www.bioconductor.org/shields/build/devel/bioc/BUSpaRse.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/BUSpaRse)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

# BUSpaRse

This package processes `bus` files generated from single-cell RNA-seq FASTQ files, e.g. using [kallisto](http://pachterlab.github.io/kallisto/). The [`bus` format](https://github.com/BUStools/BUS-format) is a table with 4 columns: **B**arcode, **U**MI, **S**et, and counts, that represent key information in single-cell RNA-seq datasets. See [this paper](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz279/5487510?redirectedFrom=fulltext) for more information about the `bus` format. A gene count matrix for a single-cell RNA-seq experiment can be generated with the `kallisto bus` command and the [bustools](https://bustools.github.io/) suite of programs many times faster than with other programs. 

The most recent version of `bustools` can convert `bus` files to the gene count and transcript compatibility count (TCC) matrices very efficiently. This package has an alternative implementation of the algorithm that converts `bus` files to gene count and TCC matrices. This implementation is much less efficient (though still many times faster than, e.g., Cell Ranger). The purpose of this implementation is to facilitate experimentation with new algorithms or to adapt the methods for other applications. The implementation in this package is written in Rcpp, which is easier to work with than pure C++ code and requires less expertise of C++.

A file mapping transcripts to genes is required to convert the `bus` file to a gene count matrix, either with `bustools` or with this package. This package contains functions that produces this file or data frame, by directly querying Ensembl, by parsing GTF or GFF3 files, by extracting information from `TxDb` or `EnsDb` gene annotation resources from Bioconductor, or by parsing sequence names of fasta files of transcriptomes downloaded from Ensembl. This package can query Ensembl for not only vertebrates (i.e. www.ensembl.org), but also [plants](plants.ensembl.org), [fungi](fungi.ensembl.org), [invertebrates](metazoa.ensembl.org), and [protists](protists.ensembl.org). Now the functions used to map transcript to genes can also filter by biotypes and only keep standard chromosomes, and extract filtered transcriptomes.

This package can also generate the files required for running RNA velocity with `kallisto` and `bustools`, including a fasta file with not only the transcriptome but also appropriately flanked intronic sequences, lists of transcripts and introns to be captured, and a file mapping transcripts and introns to genes. For spliced transcripts, you may either use the cDNA sequences, or exon-exon junctions, for pseudoalignment. Using exon-exon junctions should more unambiguously distinguish between spliced and unspliced transcripts, since unspliced transcripts also have exonic sequences. 

## Example
See [the vignettes](https://bustools.github.io/BUS_notebooks_R/index.html) for examples of using `kallisto bus`, `bustools`, and `BUSpaRse` on real data. The vignettes contain a complete walk-through, starting with downloading the FASTQ files for an experiment and ending with an analysis. Google Colab version of those vignettes can be found [here](https://www.kallistobus.tools/tutorials). Also see `browseVignettes("BUSpaRse")` for vignettes for using `BUSpaRse` to get gene count matrix and for extracting filtered transcriptomes with `tr2g_*` functions.

## Installation

You can install development version of BUSpaRse with:

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("BUStools/BUSpaRse")
```

The release version can be installed from Bioconductor, or the development version with the `version = "devel"` argument:

```r
BiocManager::install("BUSpaRse")
```
