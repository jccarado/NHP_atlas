import scanpy as sc
import os
import celltypist
from celltypist import models
import seaborn as sns
import matplotlib.pyplot as plt

sc.settings.set_figure_params(figsize=(3,3))
os.chdir("/peng_2/peng_lab/projects/rhesus_atlas/ongoing/atlas/scVI")
adata = sc.read_h5ad("b_cells_atlas_integrated.h5ad")
adata.obs_names_make_unique()

adata.X = adata.layers["counts"].toarray()
sc.pp.normalize_total(adata, target_sum = 1e4)
sc.pp.log1p(adata)

predictions_immune_low = celltypist.annotate(adata, model = "Immune_All_Low.pkl", majority_voting = True)
adata = predictions_immune_low.to_adata()
import pandas as pd

cluster_key = 'leiden'
label_key = 'predicted_labels'

sc.tl.leiden(adata, resolution = 4.0)

# dominant CellTypist label per cluster
cluster_to_label = (
    adata.obs
    .groupby(cluster_key)[label_key]
    .agg(lambda x: x.value_counts().idxmax())
)

# map back to cells
adata.obs['cell_type_cluster'] = adata.obs[cluster_key].map(cluster_to_label)

sc.pl.umap(adata, color = "cell_type_cluster", legend_loc = "on data", legend_fontsize = 4, save = "_celltypist.png")
sc.pl.umap(adata, color = "conf_score", save = "_conf_scores.png")

adata.write_h5ad("b_cells_atlas_integrated.h5ad")



