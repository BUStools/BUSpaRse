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

#' Tags in the attributes field of Ensembl GTF files
#' 
#' These are the column names of the `mcols` when the Ensembl GTF file is read
#' into R as a `GRanges`, including `gene_id`, `transcript_id`, `gene_biotype`,
#' `transcript_biotype`, `description`, and so on.
#' 
#' @format A character vector
#' @source Ensembl version 99 human GTF file
"ensembl_gtf_mcols"

# Tags in the attributes field of Ensembl GFF files
#' 
#' These are the column names of the `mcols` when the Ensembl GTF file is read
#' into R as a `GRanges`, including `gene_id`, `transcript_id`, `biotype`,
#' `description`, and so on, and the mandatory tags like `ID`, `Name`, and 
#' `Parent`.
#' 
#' @format A character vector
#' @source Ensembl version 99 human GFF3 file
"ensembl_gff_mcols"

#' Tags in the attributes field of RefSeq GFF files
#' 
#' These are the column names of the `mcols` when the Ensembl GTF file is read
#' into R as a `GRanges`, including `gene`, `transcript_id`, `gene_biotype`,
#' `description`, and so on, and the mandatory tags like `ID`, `Name`, and 
#' `Parent`.
#' 
#' @format A character vector
#' @source Ensembl version 99 human GTF file
"refseq_gff_mcols"
