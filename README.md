[![Travis build status](https://travis-ci.com/BUStools/BUSpaRse.svg?branch=master)](https://travis-ci.com/BUStools/BUSpaRse)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/BUStools/BUSpaRse?branch=master&svg=true)](https://ci.appveyor.com/project/BUStools/BUSpaRse)

# BUSpaRse

This package processes `bus` files generated from raw single-cell RNA-seq data, e.g. using [kallisto](http://pachterlab.github.io/kallisto/). The [`bus` format](https://github.com/BUStools/BUS-format) is a table with 4 columns: **B**arbode, **U**MI, **S**et, and counts, that represent scRNA-seq datasets. See [this paper](https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/btz279/5487510?redirectedFrom=fulltext) for more information about the `bus` format. The gene count matrix can be generated with `kallisto bus` and `bustools` from FASTQ files many times faster than CellRanger while giving similar results. 

The most recent version of `bustools` can convert `bus` files to the gene count and transcript compatibility count (TCC) matrices very efficiently. This package has an alternative implementation of the algorithm that converts `bus` files to gene count and TCC matrices, though this implementation is much less efficient (though still many times faster than CellRanger). The purpose of this implementation is the ease for users to alter code to explore other algorithms for this purpose or to adapt the code for a somewhat different purpose. The implementation in this package is written in Rcpp, which is easier to work with than pure C++ code and requires less expertise of C++.

A file mapping transcripts to genes is required to convert the `bus` file to a gene count matrix, either with `bustools` or with this package. This package contains functions that produces this file or data frame, by directly querying Ensembl, by parsing GTF or GFF3 files, or by parsing sequence names of FASTA files of transcriptomes downloaded from Ensembl. This package can query Ensembl for not only vertebrates (i.e. www.ensembl.org), but also [plants](plants.ensembl.org), [fungi](fungi.ensembl.org), [invertebrates](metazoa.ensembl.org), and [protists](protists.ensembl.org). 

## Example
See [the vignettes](https://bustools.github.io/BUS_notebooks_R/index.html) for examples of using `kallisto bus`, `bustools`, and `BUSpaRse` on real data. The vignettes contain a complete walk-through, starting with downloading the FASTQ files for an experiment and ending with an analysis. There are currently vignettes for 10x v2 and 10x v3 technology datasets (Drop-seq, inDrops, and CEL-seq2 coming soon). 

## Installation

You can install BUSpaRse with:

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("BUStools/BUSpaRse")
```

The package will be available on Bioconductor shortly, but before its first Bioconductor release, **the user interface may change substantially**, and the changes will be reflected in the vignettes. If you reinstall this package prior to the first Bioconductor release, please check for changes in the user interface.


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
