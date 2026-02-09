import os
import seaborn as sns
import scanpy as sc
import anndata
import scib
import scvi
import sys
import seaborn as sns
import ray
import torch
from ray import tune
import pandas as pd
import matplotlib.pyplot as plt
from scvi import autotune
import re

cell_type = sys.argv[1]
adata = sc.read_h5ad(f"/peng_2/peng_lab/results/ongoing/rhesus_atlas/atlas/unintegrated_processed_{cell_type}.h5ad")
adata_subset = adata[:, adata.var['highly_variable']].copy()
scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="batch_id")
model_scvi = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model_scvi.view_anndata_setup()
model_scvi.train(check_val_every_n_epoch=1)
history = model_scvi.history

ig, axs = plt.subplots(5, 1, figsize=(8, 18))  # 5 rows, 1 column

metrics_to_plot = list(history.keys())  # Replace with desired keys if needed

fig, axs = plt.subplots(len(metrics_to_plot), 1, figsize=(8, 4 * len(metrics_to_plot)))

for i, key in enumerate(metrics_to_plot):
    axs[i].plot(history[key])
    axs[i].set_title(key)
    axs[i].set_xlabel("Epoch")
    axs[i].set_ylabel(key)

plt.tight_layout()
plt.savefig(f"{cell_type}_integration.pdf")


adata.obsm["X_scVI"] = model_scvi.get_latent_representation()
sc.pp.neighbors(adata, use_rep = "X_scVI")
sc.tl.leiden(adata, flavor = "igraph", resolution = 0.8)
sc.tl.umap(adata)
adata.write_h5ad(f"/peng_2/peng_lab/results/ongoing/rhesus_atlas/atlas/scVI/{cell_type}_atlas_integrated.h5ad")




scvi.model.SCVI.setup_anndata(adata_subset, layer="counts", batch_key="batch_id", continuous_covariate_keys=["pct_counts_mt", "pct_counts_ribo", "pct_counts_housekeeping"])
model_scvi = scvi.model.SCVI(adata_subset, n_layers=2, n_latent=30, gene_likelihood="nb")
model_scvi.view_anndata_setup()
model_scvi.train(max_epochs = 100, early_stopping = True)

adata_subset.obsm["X_scVI"] = model_scvi.get_latent_representation()

cell_type_barcodes = pd.read_csv("/peng_2/peng_lab/ngs_data/private/rhesus_atlas/barcodes_metadata.txt", sep = "\t")
temp_obs = adata.obs.reset_index()
temp_obs.rename(columns={'index': 'barcode'}, inplace=True)
merged = temp_obs.merge(
    cell_type_barcodes[['barcode', 'batch_id', 'cell_type']], 
    on=['barcode', 'batch_id'], 
    how='left'
)


merged.set_index('barcode', inplace=True)

adata.obs = merged
scanvi_model = scvi.model.SCANVI.from_scvi_model(
    model_scvi,
    adata=adata,
    labels_key="cell_type",
    unlabeled_category="Unknown",
)

scanvi_model.train(max_epochs=100)
SCANVI_LATENT_KEY = "X_scANVI"
adata.obsm[SCANVI_LATENT_KEY] = scanvi_model.get_latent_representation(adata)
sc.pp.neighbors(adata_subset, use_rep = "X_scANVI")
sc.tl.leiden(adata_subset, flavor = "igraph", resolution = 0.8)
sc.tl.umap(adata_subset, min_dist=0.3)

adata_subset.write_h5ad(f"/peng_2/peng_lab/results/ongoing/rhesus_atlas/atlas/scVI/{cell_type}_atlas_integrated_highly_variable.h5ad")
