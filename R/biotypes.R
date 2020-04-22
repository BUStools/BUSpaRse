#' Gene biotypes from Ensembl
#' 
#' These vectors are here to make it easy to look up which biotypes are 
#' available for Ensembl without having to parse GTF and fasta files every time.
#' See \href{http://www.ensembl.org/info/genome/genebuild/biotypes.html}{this page from Ensembl}
#' for what the biotypes mean.
#' 
#' @format A character vector with all Ensembl gene biotypes.
#' @source This vector is all the unique gene biotypes in the Ensembl version 99
#' human GTF file.
#' @seealso ensembl_tx_biotypes cellranger_biotypes
"ensembl_gene_biotypes"

#' Transcript biotypes from Ensembl
#' 
#' These vectors are here to make it easy to look up which biotypes are 
#' available for Ensembl without having to parse GTF and fasta files every time.
#' See \href{http://www.ensembl.org/info/genome/genebuild/biotypes.html}{this page from Ensembl}
#' for what the biotypes mean.
#' 
#' @format A character vector with all Ensembl transcript biotypes.
#' @source This vector is all the unique transcript biotypes in the Ensembl 
#' version 99 human GTF file.
#' @seealso ensembl_gene_biotypes cellranger_biotypes
"ensembl_tx_biotypes"

#' Cell Ranger gene biotypes
#' 
#' In the GRCh38 Cell Ranger reference package, an Ensembl human GTF file is
#' filtered by gene biotypes. This vector includes all gene biotypes included by
#' this Cell Ranger reference package. 
#' 
#' @format A character vector with all Cell Ranger reference package gene
#' biotypes.
#' @source \url{https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references}
#' @seealso ensembl_gene_biotypes ensembl_tx_biotypes
"cellranger_biotypes"
