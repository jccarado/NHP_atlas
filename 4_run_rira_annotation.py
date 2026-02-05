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
from rpy2.robjects import pandas2ri
import scvi
import rpy2.robjects as ro
from rpy2 import robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri, default_converter
from rpy2.robjects.conversion import localconverter
import anndata2ri
import rpy2.robjects.packages as rpackages
import scvi
import rpy2



batch_id = sys.argv[1]
directory = sys.argv[2]
directory = directory.strip()
batch_id = batch_id.strip()
rds = f"/peng_2/peng_lab/projects/rhesus_atlas/raw/{directory}/{batch_id}_downstream.rds"
ro.r.assign("file", rds)
print(batch_id, flush = True)
h5ad = f"/peng_2/peng_lab/projects/rhesus_atlas/raw/{directory}/{batch_id}_downstream.h5ad"
adata = sc.read_h5ad(h5ad)
counts = adata.layers["counts"].toarray()
data = adata.layers["data"].toarray()



with localconverter(default_converter + numpy2ri.converter):
    ro.globalenv['counts'] = counts
    ro.globalenv['data'] = data
    
ro.r(
'''
library(Seurat)
library(RIRA)
counts <- t(counts)
data <- t(data)
seurat_obj <- readRDS(file)
gene_names <- rownames(seurat_obj@assays$RNA)
cell_names <- colnames(seurat_obj@assays$RNA)
rownames(counts) <- gene_names
rownames(data) <- gene_names
colnames(counts) <- cell_names
colnames(data) <- cell_names
seurat_obj <- SetAssayData(object = seurat_obj,assay = "RNA",slot = "scale.data",new.data = as.matrix(seurat_obj@assays$RNA@counts))
seurat_obj <- SetAssayData(object = seurat_obj,assay = "RNA",slot = "counts",new.data = counts)
seurat_obj <- SetAssayData(object = seurat_obj,assay = "RNA",slot = "data",new.data = data)
seurat_obj@assays$RNA@counts <- as(seurat_obj@assays$RNA@counts, "dgCMatrix")
seurat_obj@assays$RNA@data <- as(seurat_obj@assays$RNA@data, "dgCMatrix")
png("/peng_2/peng_lab/projects/rhesus_atlas/raw/rira_predictions.png")
seurat_obj <- Classify_ImmuneCells(seurat_obj, retainProbabilityMatrix = TRUE)
seurat_obj <- Classify_TNK(seurat_obj, retainProbabilityMatrix = TRUE)
seurat_obj <- Classify_Myeloid(seurat_obj, retainProbabilityMatrix = TRUE)
dev.off()
'''
)
RIRA_predictions = ro.r('seurat_obj$RIRA_Immune_v2.predicted_labels')
adata.obs["RIRA_predictions"] = list(RIRA_predictions)
adata.obs["RIRA_TNK_predictions"] = ro.r('seurat_obj$RIRA_TNK_v2.predicted_labels')
adata.obs["RIRA_Myeloid_predictions"] = ro.r('seurat_obj$RIRA_Myeloid_v3.predicted_labels')
adata.write_h5ad(f"/peng_2/peng_lab/projects/rhesus_atlas/raw/{directory}/{batch_id}_downstream.h5ad")
os.chdir(f"/peng_2/peng_lab/projects/rhesus_atlas/raw/{directory}")
sc.pl.umap(adata, color = ["RIRA_predictions", "annotations", "dominant", "leiden"], ncols =2, save = f"_{batch_id}_rira_predictions.png")
