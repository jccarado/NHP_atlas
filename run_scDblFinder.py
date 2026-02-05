import os
import seaborn as sns
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
import scib
import sys
import doubletdetection
import random
import anndata as ad
import pandas as pd
import scrublet as scr
from scipy.sparse import issparse
import celldex
import singlecellexperiment as sce
import singler
import re
import scipy.sparse
from scipy.sparse import csr_matrix
import celltypist
import subprocess
import rpy2.robjects as ro
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
import anndata2ri
import rpy2.robjects.packages as rpackages
import scvi
import rpy2


def scDblFinder(adata):
    seurat = rpackages.importr('Seurat')
    scdblfinder = rpackages.importr("scDblFinder")
    tidyr = rpackages.importr("tidyr")
    X = adata.X.toarray() if hasattr(adata.X, "toarray") else adata.X
    obs = adata.obs
    var = adata.var
    with localconverter(default_converter + numpy2ri.converter):
      ro.globalenv['X'] = X.T
    with localconverter(default_converter + pandas2ri.converter):
      ro.globalenv['obs'] = obs
    with localconverter(default_converter + pandas2ri.converter):
      ro.globalenv['var'] = var
    ro.r('''
    library(SummarizedExperiment)
    library(Seurat)
    library(scDblFinder)
    library(SingleCellExperiment)
    library(dplyr)
    
    # Create SummarizedExperiment
    adata <- SummarizedExperiment::SummarizedExperiment(
        assays = list(X = as.matrix(X)),
        colData = obs,
        rowData = var
    )
    
    counts <- SummarizedExperiment::assay(adata, "X")
    seurat_obj <- CreateSeuratObject(counts)
    options(future.globals.maxSize = 5e9)
    sce <- as.SingleCellExperiment(seurat_obj)
    sce <- scDblFinder(sce, dbr.sd = 1)
    meta_scdblfinder <- sce@colData@listData %>%
        as.data.frame() %>%
        dplyr::select(starts_with('scDblFinder'))
    rownames(meta_scdblfinder) <- rownames(sce@colData)
    seurat_obj <- AddMetaData(object = seurat_obj, metadata = meta_scdblfinder)
    ''')
    df_labels = ro.r('as.vector(seurat_obj@meta.data[[grep("scDblFinder.class", colnames(seurat_obj@meta.data), value=TRUE)]])')
    df_scores = ro.r('as.vector(seurat_obj@meta.data[[grep("scDblFinder.score", colnames(seurat_obj@meta.data), value=TRUE)]])')
    adata.obs['scDblFinder_status'] = list(df_labels)
    adata.obs['scDblFinder_score'] = list(df_scores)
    return adata

file = sys.argv[1]
adata = sc.read_h5ad(file)
sample_name = adata.obs["batch_id"].iloc[0]
adata = scDblFinder(adata)
doublets = adata[adata.obs["scDblFinder_status"] == "doublet", :].copy()
scores = adata.obs["scDblFinder_score"].tolist()
barcodes = doublets.obs.index.tolist()
os.chdir("/peng_2/peng_lab/ngs_data/private/rhesus_atlas/Doublets/scDblFinder")
with open(f"{sample_name}_doublets.txt", "w") as f:
     for barcode in barcodes:
         f.write(f"{barcode}\n")
         
with open(f"{sample_name}_scores.txt", "w") as f:
  for score in scores:
    f.write(f"{score}\n")

