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
from rpy2.robjects import pandas2ri, numpy2ri
import anndata2ri
import rpy2.robjects.packages as rpackages
import scvi
import rpy2
from rpy2.robjects import r
import pyreadr


batch_id = sys.argv[1]
directory = sys.argv[2]
directory = directory.strip()
batch_id = batch_id.strip()

cc_genes_rhesus = pd.read_csv("/peng_2/peng_lab/tools/scVDJ/Rfunctions/CellCycle/regev_lab_cell_cycle_genes_orthologs.txt", sep="\t",header=0)
cc_genes_rhesus = cc_genes_rhesus[~cc_genes_rhesus["Rhesus"].isna()]
s_genes = cc_genes_rhesus.loc[cc_genes_rhesus["cycle"] == "S", "Rhesus"].tolist()
g2m_genes = cc_genes_rhesus.loc[cc_genes_rhesus["cycle"] == "G2M", "Rhesus"].tolist()
ref_data = celldex.fetch_reference("hpca", "2024-02-26", realize_assays=True)
VDJ_prefix_pattern = []
loc_metadata = pd.read_csv("/peng_2/peng_lab/ngs_data/private/rhesus_atlas/loc_genes_filter_out.txt", sep = "\t", header = 0)
loc_genes = loc_metadata["Gene"].tolist()
housekeeping_genes = pd.read_csv("/peng_2/peng_lab/ngs_data/private/rhesus_atlas/housekeeping_genes.txt", sep = "\t", header = 0)
housekeeping_genes = housekeeping_genes.iloc[:,0].tolist()

def read_barcodes(file_path):
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        try:
            df = pd.read_csv(file_path, header=None, sep="\t")
            if not df.empty:
                barcodes = df[0].astype(str).str.strip()
                return barcodes[barcodes != ''].tolist()
        except pd.errors.EmptyDataError:
            pass
    return []

blueprint_ref = celldex.fetch_reference("blueprint_encode", "2024-02-26", realize_assays=True)


with open("/peng_2/peng_lab/results/ongoing/Rhesus_TulaneCure/QC/rhesus_VDJ_prefixes.txt", 'r') as f:
        VDJ_prefixes = f.read().splitlines()
VDJ_prefix_pattern = VDJ_prefixes
os.chdir(f"/peng_2/peng_lab/ngs_data/private/rhesus_atlas/{directory}")
file = batch_id + ".h5ad"
adata = sc.read_h5ad(file)
dd_barcodes = read_barcodes(f"/peng_2/peng_lab/ngs_data/private/rhesus_atlas/Doublets/DoubletDetector/{batch_id}_doublets.txt")
scdbl_barcodes = read_barcodes(f"/peng_2/peng_lab/ngs_data/private/rhesus_atlas/Doublets/scDblFinder/{batch_id}_doublets.txt")
if len(dd_barcodes) == 0:
    doublet_barcodes = set(scdbl_barcodes)
elif len(scdbl_barcodes) == 0:
    doublet_barcodes = set(dd_barcodes)
else: 
    doublet_barcodes = set(dd_barcodes) & set(scdbl_barcodes)
print(f"Number of doublets: {len(doublet_barcodes)}")
adata.obs["dd_scores"] = pd.read_csv(f"/peng_2/peng_lab/ngs_data/private/rhesus_atlas/Doublets/DoubletDetector/{batch_id}_scores.txt", sep = "\t", header = None).iloc[:, 0].tolist()
adata.obs["scdbl_scores"] = pd.read_csv(f"/peng_2/peng_lab/ngs_data/private/rhesus_atlas/Doublets/scDblFinder/{batch_id}_scores.txt", sep = "\t", header = None).iloc[:, 0].tolist()
adata.obs["dd_doublet"] = "Singlet"
adata.obs["scdbl_doublet"] = "Singlet"
adata.obs.loc[dd_barcodes, "dd_doublet"] = "Doublet"
adata.obs.loc[scdbl_barcodes, "scdbl_doublet"] = "Doublet"
adata = adata[~adata.obs_names.isin(doublet_barcodes)].copy()
adata = adata[adata.obs["scdbl_scores"] < 0.75, :].copy()
adata = adata[adata.obs["dd_scores"] < 100, :].copy()
adata = adata[adata.obs["pct_counts_hb"] < 20, :].copy()
sce_adata = sce.SingleCellExperiment.from_anndata(adata)
mat = sce_adata.assay("X")


features = list(sce_adata.row_data.row_names)


results, integrated = singler.annotate_integrated(test_data = mat, ref_data = [blueprint_ref, ref_data], ref_labels = [blueprint_ref.get_column_data().column("label.main"), ref_data.get_column_data().column("label.main")], test_features = features)


predicted_celltypes = integrated["best_label"]
adata.obs['annotations'] = predicted_celltypes
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="batch_id", min_disp = 0.5)
var_features = adata.var_names[adata.var['highly_variable']]
ribosomal_genes = [g for g in var_features if re.search(r"^r[PR]-|^RP|rRNA|114674534|109910387", g, re.IGNORECASE)]
mitochondrial_genes = [g for g in var_features if re.search(r"^MT-", g, re.IGNORECASE)]
VDJ_genes = [g for g in var_features for pattern in VDJ_prefix_pattern if re.search(pattern, g, re.IGNORECASE)] 
microRNA_genes = [g for g in var_features if re.search(r"^MIR", g, re.IGNORECASE)]
genes_to_remove = list(set(ribosomal_genes + mitochondrial_genes + VDJ_genes + microRNA_genes + loc_genes + housekeeping_genes))
var_features_filtered = [g for g in var_features if g not in genes_to_remove]
adata.var['highly_variable'] = adata.var_names.isin(var_features_filtered)



adata.raw = adata
adata.layers["data"] = adata.X.copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, use_highly_variable=True, n_comps=50, svd_solver="arpack")
sc.tl.score_genes_cell_cycle(adata,s_genes=s_genes, g2m_genes=g2m_genes, copy=False)
adata.layers["scale.data"] = adata.X
sc.pp.neighbors(adata)
sc.tl.leiden(adata, resolution = 1.2)
sc.tl.umap(adata)
clusters = adata.obs["leiden"]
cell_types = adata.obs["annotations"]
for cluster in clusters.unique():
    predominant_cell_types =[]
    cells_in_cluster = cell_types[clusters == cluster]
    predominant_cell_types.append(cells_in_cluster.value_counts().idxmax())
    if predominant_cell_types[0] in ["B_cell", "Pre-B_cell_CD34-", "Pro-B_cell_CD34+", "B-cells"]:
        adata.obs.loc[adata.obs["leiden"] == cluster, "dominant"] = "B_cell"
    elif predominant_cell_types[0] in ["T_cells", "NK_cell", "CD4+ T-cells", "CD8+ T-cells", "NK cells"]:
        adata.obs.loc[adata.obs["leiden"] == cluster, "dominant"] = "T_cell/NK_cell"
    elif predominant_cell_types[0] in ["Monocyte", "DC", "Macrophage", "GMP", "CMP", "Pro-Myelocyte", "Monocytes", "Macrophages"]:
        adata.obs.loc[adata.obs["leiden"] == cluster, "dominant"] = "Monocytes"
os.chdir(f"/peng_2/peng_lab/projects/rhesus_atlas/raw/{directory}")
sc.pl.umap(adata, color = ["annotations", "dominant", "leiden"], ncols = 3, legend_loc = "on data", legend_fontsize = 4, save = f"_{batch_id}_celltypist_annotations.png")
adata.X = adata.layers["counts"].toarray()
adata.write_h5ad(f"/peng_2/peng_lab/ngs_data/private/rhesus_atlas/{directory}/{batch_id}_downstream.h5ad")


