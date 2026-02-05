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
from rpy2.robjects import pandas2ri, conversion
import anndata2ri
import rpy2.robjects.packages as rpackages
import scvi
import rpy2


def doubletDetector(adata):
    clf = doubletdetection.BoostClassifier(n_iters=10, clustering_algorithm="leiden", standard_scaling=True,pseudocount=0.1,n_jobs=-1,)
    sample_name = adata.obs["batch_id"].iloc[0]
    doublets = clf.fit(adata.X).predict(p_thresh=1e-16, voter_thresh=0.5)
    doublet_score = clf.doublet_score()
    adata.obs["doublet_detector"] = doublets
    adata.obs["doublet_detector_score"] = doublet_score
    os.chdir("/peng_2/peng_lab/ngs_data/private/rhesus_atlas/Doublets/DoubletDetector")
    doubletdetection.plot.threshold(clf, save=f'{sample_name}_threshold_test.pdf', show=True, p_step=6)
    return adata

file = sys.argv[1]
adata = sc.read_h5ad(file)
sample_name = adata.obs["batch_id"].iloc[0]
adata = doubletDetector(adata)
doublets = adata[adata.obs["doublet_detector"] == 1.0, :].copy()
scores = adata.obs["doublet_detector_score"].tolist()
barcodes = doublets.obs.index.tolist()
os.chdir("/peng_2/peng_lab/ngs_data/private/rhesus_atlas/Doublets/DoubletDetector")
with open(f"{sample_name}_doublets.txt", "w") as f:
     for barcode in barcodes:
         f.write(f"{barcode}\n")
    
with open(f"{sample_name}_scores.txt", "w") as f:
  for score in scores:
      f.write(f"{score}\n")
