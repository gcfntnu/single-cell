#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(Seurat))

if (!require(Seurat)) {
    stop("This script needs to have Seurat installed.")
}
if (!require(argparse)) {
    install.packages("argparse", repos="http://cran.rstudio.com")
    library("argparse")
}


if (!require(Matrix)) {
    install.packages("Matrix", repos="http://cran.rstudio.com")
    library("Matrix")
}


parser <- ArgumentParser()
parser$add_argument("-i", "--input", default="data/tmp/test_sample/outs/filtered_gene_bc_matrices", nargs="+",
                    help="Directory containing the matrix.mtx.gz, genes.tsv.gz, and barcodes.tsv.gz files provided by 10X. Multiple dirs allowed.")
parser$add_argument("-s", "--sample-info", dest="samples", default=NULL,
                    help="Sample info. Tab delimited file, needs column with `sample_id`")
parser$add_argument("-o", "--output", default="data/processed/test_sample/seurat_obj.rds",
                    help="Output filename (rds file)")
parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))
if (args$verbose == TRUE) options(echo=TRUE)

if (length(args$input)> 1){
    input.names <- NULL
    for (dir in args$input){
        parent.dir <- dirname(dir)
        if (basename(parent.dir) == "outs"){
            parent.dir <- dirname(parent.dir)
        }
        input.names <- c(input.names, basename(parent.dir))
    }
    names(args$input) <- input.names
}

data <- Read10X(args$input)
if (!is.null(args$samples)){
    meta.info <- read.delim(args$samples, sep="\t", row.names=1)
} else{
    meta.info <- NULL
}

seurat.obj  <- CreateSeuratObject(data, meta.data = meta.info, assay="RNA", names.delim="_")
saveRDS(seurat.obj, args$output)
