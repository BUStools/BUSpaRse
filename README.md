[![Travis build status](https://travis-ci.com/BUStools/BUSpaRse.svg?branch=master)](https://travis-ci.com/BUStools/BUSpaRse)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/BUStools/BUSpaRse?branch=master&svg=true)](https://ci.appveyor.com/project/BUStools/BUSpaRse)

# BUSpaRse

This package processes `bus` files generated from raw single-cell RNA-seq data, e.g. using [kallisto](http://pachterlab.github.io/kallisto/). The [`bus` format](https://github.com/BUStools/BUS-format) is a table with 4 columns: **B**arbode, **U**MI, **S**et, and counts, that represent scRNA-seq datasets. The main goal of this package is to parse `bus` format in R into a sparse matrix that can be used for downstream analysis, such as with `Seurat` and `SingleCellExperiment`. See [this paper](https://www.biorxiv.org/content/early/2018/11/21/472571) for more information about the `bus` format.

## Example
See [the vignettes](https://bustools.github.io/BUS_notebooks_R/index.html) for examples of using `BUSpaRse`. The vignettes contain a complete walk-through, starting with downloading of the FASTQ files for an experiment and ending with an analysis. There are currently vignettes for 10x v2 and 10x v3 technology datasets (Drop-seq, inDrops, and CEL-seq2 coming soon). 

## Installation

You can install BUSpaRse with:

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("BUStools/BUSpaRse")
```

**Note**: This is work in progress. Though this package is currently incomplete, the functionality available in the master branch has been tested and can be used. More features are coming up in the `devel` branch, and when they become stable, they will be introduced to the master branch. Please do not install this package from the `devel` branch, which is unstable.

The package will be available on Bioconductor shortly, but before its first Bioconductor release, **the user interface may change substantially**, and the changes will be reflected in the vignettes. If you reinstall this package prior to the first Bioconductor release, please check for changes in the user interface.

Features currently available in the master branch:

* Create gene count matrix from BUS output (sorted and converted to text) in one step for organisms available on Ensembl. Only works if kallisto index was created from Ensembl transcriptomes.
* Create gene count matrix in 3 separate steps for organisms not available on Ensembl:
  - Map transcripts to genes (`tr2g_*` functions)
  - Determine the genes each equivalence class maps to (`EC2gene`)
  - Create the gene count matrix (`make_sparse_matrix`)
* Extract transcript and gene IDs from the following file formats with `tr2g_*` functions:
  - Ensembl (directly query biomart with species names)
  - FASTA files of transcriptome, only works for Ensembl transcriptomes
  - GTF files
  - GFF3 files

Upcoming features in `devel` branch:

* Create transcript compatibility count (TCC) matrices
* Multithreaded processing of BUS output when generating the sparse matrix

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

This package contains compiled code, and a compiler that supports OpenMP is required to compile this package. However, the default clang that comes with MacOS does not support OpenMP. MacOS users using R 3.5 should download and install Clang 6.0 and gfortran 6.1 compilers from [this webpage from CRAN](https://cran.r-project.org/bin/macosx/tools/), which has OpenMP enabled. R 3.5 no longer works with Clang 4, which was used for R 3.4. The Clang 6.0 and gfortran 6.1 from CRAN are precompiled, so you can directly install them with the graphic installer rather than compile them from source. Any compiler with OpenMP support can be used instead of the Clang 6.0 from CRAN, but in general, you need to compile the compiler from source.

Then, if the file `~/.R/Makevars` does not exist, in the terminal, go to your home directory by `cd`, use `mkdir .R` to create the `.R` directory, and type `vim Makevars` to create and start editing the file. If it already exists, then type `vim Makevars` to edit it.

Alternatively, if you are uncomfortable with the command line, this can be done in RStudio. First use `file.exists("~/.R/Makevars")` to check if `~/.R/Makevars` exists. Then use `dir.exists("~/.R")` to check that if the `~/.R` directory exists. If it does not, then use `dir.create("~/.R")` to create the directory. Then use `file.create("~/.R/Makevars")` to create that file. Then navigate to that file in the Files pane in RStudio, open that file in RStudio, and edit it.

Add the following to the `~/.R/Makevars` file:

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

Above is the default path where this Clang 6.0 is installed. Please change it if Clang 6.0 is installed in a custom path, or if a different compiler with OpenMP support is used. This will tell R to use the compiler that has OpenMP enabled. Then restart the R session and reinstall this package.

[^1]: StackOverflow question https://stackoverflow.com/questions/52513927/installing-xcode-command-line-tools
[^2]: https://teuder.github.io/rcpp4everyone_en/020_install.html
