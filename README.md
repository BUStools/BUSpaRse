<!-- badges: start -->
[![Travis build status](https://travis-ci.com/BUStools/BUSpaRse.svg?branch=master)](https://travis-ci.com/BUStools/BUSpaRse)
[![BioC status](http://www.bioconductor.org/shields/build/devel/bioc/BUSpaRse.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/BUSpaRse)
[![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

# BUSpaRse

This package processes `bus` files generated from single-cell RNA-seq FASTQ files, e.g. using [kallisto](http://pachterlab.github.io/kallisto/). The [`bus` format](https://github.com/BUStools/BUS-format) is a table with 4 columns: **B**arcode, **U**MI, **S**et, and counts, that represent key information in single-cell RNA-seq datasets. See [this paper](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz279/5487510?redirectedFrom=fulltext) for more information about the `bus` format. A gene count matrix for a single-cell RNA-seq experiment can be generated with the `kallisto bus` command and the [bustools](https://bustools.github.io/) suite of programs many times faster than with other programs. 

The most recent version of `bustools` can convert `bus` files to the gene count and transcript compatibility count (TCC) matrices very efficiently. This package has an alternative implementation of the algorithm that converts `bus` files to gene count and TCC matrices. This implementation is much less efficient (though still many times faster than, e.g., Cell Ranger). The purpose of this implementation is to facilitate experimentation with new algorithms or to adapt the methods for other applications. The implementation in this package is written in Rcpp, which is easier to work with than pure C++ code and requires less expertise of C++.

A file mapping transcripts to genes is required to convert the `bus` file to a gene count matrix, either with `bustools` or with this package. This package contains functions that produces this file or data frame, by directly querying Ensembl, by parsing GTF or GFF3 files, by extracting information from `TxDb` or `EnsDb` gene annotation resources from Bioconductor, or by parsing sequence names of fasta files of transcriptomes downloaded from Ensembl. This package can query Ensembl for not only vertebrates (i.e. www.ensembl.org), but also [plants](plants.ensembl.org), [fungi](fungi.ensembl.org), [invertebrates](metazoa.ensembl.org), and [protists](protists.ensembl.org). 

This package can also generate the files required for running RNA velocity with `kallisto` and `bustools`, including a fasta file with not only the transcriptome but also appropriately flanked intronic sequences, lists of transcripts and introns to be captured, and a file mapping transcripts and introns to genes. For spliced transcripts, you may either use the cDNA sequences, or exon-exon junctions, for pseudoalignment. Using exon-exon junctions should more unambiguously distinguish between spliced and unspliced transcripts, since unspliced transcripts also have exonic sequences. 

## Example
See [the vignettes](https://bustools.github.io/BUS_notebooks_R/index.html) for examples of using `kallisto bus`, `bustools`, and `BUSpaRse` on real data. The vignettes contain a complete walk-through, starting with downloading the FASTQ files for an experiment and ending with an analysis. 

## Installation

You can install BUSpaRse with:

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("BUStools/BUSpaRse")
```

Or if you are using the development version of Bioconductor, you can install BUSpaRse with:

```r
BiocManager::install("BUSpaRse")
```

The first release of this package will be available when Bioconductor 3.10 is released.

### Installation note for MacOS
First of all, install Xcode command line tool. 

**Mojave**:

Go to https://developer.apple.com/download/more/ and search for "command line". Then Download "Command line tool for MacOS 10.14". Once dmg is downloaded install the package[^1].

**Older versions of MacOS**: 

Execute the command `xcode-select --install` on Terminal[^2].

_In case you encounter this error during installation_:

```
clang: error: unsupported option '-fopenmp'
```

This is what to do to resolve it:

This package contains compiled code, and a compiler that supports OpenMP is required to compile this package. However, the default clang that comes with MacOS does not support OpenMP. MacOS users using R 3.5 and above should download and install the appropriate version of the Clang compiler from [this webpage from CRAN](https://cran.r-project.org/bin/macosx/tools/), which has OpenMP enabled. R 3.5 no longer works with Clang 4, which was used for R 3.4. The Clang from CRAN are precompiled, so you can directly install them with the graphic installer rather than compile them from source. Any compiler with OpenMP support can be used instead of the Clang 6.0 from CRAN, but in general, you need to compile the compiler from source.

Then, if the file `~/.R/Makevars` does not exist, in the terminal, go to your home directory by `cd`, use `mkdir .R` to create the `.R` directory, and type `vim Makevars` to create and start editing the file. If it already exists, then type `vim Makevars` to edit it.

Alternatively, if you are uncomfortable with the command line, this can be done in RStudio. First use `file.exists("~/.R/Makevars")` to check if `~/.R/Makevars` exists. Then use `dir.exists("~/.R")` to check that if the `~/.R` directory exists. If it does not, then use `dir.create("~/.R")` to create the directory. Then use `file.create("~/.R/Makevars")` to create that file. Then navigate to that file in the Files pane in RStudio, open that file in RStudio, and edit it.

For instance, for R 3.5.x, add the following to the `~/.R/Makevars` file:

```
CC=/usr/local/clang6/bin/clang
SHLIB_CXXLD=/usr/local/clang6/bin/clang++
CXX= /usr/local/clang6/bin/clang++  -Wall
CXX1X= /usr/local/clang6/bin/clang++
CXX98= /usr/local/clang6/bin/clang++
CXX11= /usr/local/clang6/bin/clang++
CXX14= /usr/local/clang6/bin/clang++
CXX17= /usr/local/clang6/bin/clang++
LDFLAGS=-L/usr/local/clang6/lib

```

Above is the default path where this Clang 6.0 is installed (for R 3.5.x). Please change it if Clang 6.0 is installed in a custom path, or if a different version of Clang or another compiler with OpenMP support is used. This will tell R to use the compiler that has OpenMP enabled. Then restart the R session and reinstall this package.

[^1]: StackOverflow question https://stackoverflow.com/questions/52513927/installing-xcode-command-line-tools

[^2]: https://teuder.github.io/rcpp4everyone_en/020_install.html
