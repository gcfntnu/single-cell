#/usr/bin/env Rscript
suppressPackageStartupMessages(require(Seurat))

if (!require(argparse)) {
    install.packages("argparse", repos="http://cran.rstudio.com")
    library("argparse")
}
if (!require(Seurat)) {
    install.packages("Seurat", repos="http://cran.rstudio.com")
    library("Seurat")
}

if (!require(Matrix)) {
    install.packages("Matrix", repos="http://cran.rstudio.com")
    library("Matrix")
}


parser <- ArgumentParser()
parser$add_argument("-i", "--input", default="data/tmp/test_sample/outs/filtered_gene_bc_matrices",
                    help="Directory containing the matrix.mtx, genes.tsv, and barcodes.tsv files provided by 10X. A list can be given in order to load several data directories")
parser$add_argument("-s", "--sample-info",
                    help="Sample info. Tab delimited file, needs column with `sample_id`")
parser$add_argument("-o", "--output", default="data/processed/test_sample/seurat_obj.rds",
                    help="Output filename (rds file)")
parser$add_argument("-v", "--verbose", action="store_true",
                    default=FALSE, help="Print extra output")

args <- parser$parse_args(args=commandArgs(TRUE))
if (args$verbose == TRUE) options(echo=TRUE)

data <- Read10X(args$input)
seurat.obj  <- CreateSeuratObject(raw.data=data)
saveRDS(seurat.obj, out.fn)
