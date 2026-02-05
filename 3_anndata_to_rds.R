library(Seurat)
library(sceasy)
library(reticulate)
library(optparse)


use_python("/home6/jccarado/miniforge3/envs/scvi-env/bin/python", required = TRUE)
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input Seurat RDS file", metavar="character"), 
  make_option(c("-o", "--output"), type="character", default=NULL, help="Output h5ad file", metavar="character")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

adata_file <- opt$input
output <- opt$output

sceasy::convertFormat(adata_file, from = "anndata", to = "seurat", main_layer = "counts", outFile = paste0(output, ".rds"))
