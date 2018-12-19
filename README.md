[![Travis build status](https://travis-ci.com/BUStools/BUSpaRse.svg?branch=master)](https://travis-ci.com/BUStools/BUSpaRse)
[![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/BUStools/BUSpaRse?branch=master&svg=true)](https://ci.appveyor.com/project/BUStools/BUSpaRse)

# BUSpaRse

This package processes `bus` files generated from raw single-cell RNA-seq data, e.g. using [kallisto](http://pachterlab.github.io/kallisto/). The [`bus` format](https://github.com/BUStools/BUS-format) is a table with 4 columns: **B**arbode, **U**MI, **S**et, and counts, that represent scRNA-seq datasets. The main goal of this package is to parse `bus` format in R into a sparse matrix that can be used for downstream analysis, such as with `Seurat` and `SingleCellExperiment`. See [this paper](https://www.biorxiv.org/content/early/2018/11/21/472571) for more information about the `bus` format.

## Installation

You can install BUSpaRse with:

``` r
if (!require(devtools)) install.packages("devtools")
devtools::install_github("BUStools/BUSpaRse")
```

This is work in progress. The package will be available on Bioconductor shortly.

Installation note: This package contains compiled code. MacOS users using R 3.5 should download and install Clang 6.0 and gfortran 6.1 compilers from [this webpage](https://cran.r-project.org/bin/macosx/tools/). R 3.5 no longer works with Clang 4, which was used for R 3.4.

## Example
See [the vignettes](https://bustools.github.io/BUS_notebooks_R/index.html) for examples of using `BUSpaRse`. The vignettes contain a complete walk-through, starting with downloading of the FASTQ files for an experiment and ending with an analysis. There are currently vignettes for 10x v2 and 10x v3 technology datasets (Drop-seq, inDrops, and CEL-seq2 coming soon). 
