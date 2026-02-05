import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import os
import pandas as pd
import numpy as np

os.chdir("/peng_2/peng_lab/projects/rhesus_atlas/ongoing/atlas/scVI")

adata = sc.read_h5ad("nk_t_cells_atlas_integrated.h5ad")
# pick one representative naive cell
root_cell = adata.obs_names[
    adata.obs['cell_type'] == 'Quiescent'
][0]

adata.uns['iroot'] = adata.obs_names.get_loc(root_cell)

sc.tl.diffmap(adata)
sc.tl.dpt(adata, n_dcs=10)

plt.figure(figsize=(6,4))
sns.scatterplot(
    x=adata.obs['dpt_pseudotime'],
    y=adata.obs['conf_score'],
    s=10,
    alpha=0.4
)
plt.xlabel('T cell pseudotime')
plt.ylabel('CellTypist confidence')
plt.savefig("diffusion_dpt.png")
plt.close()


plt.figure(figsize=(6,4))

sns.regplot(
    x=adata.obs['dpt_pseudotime'],
    y=adata.obs['conf_score'],
    scatter=True,
    lowess=True
)

plt.xlabel('T cell pseudotime')
plt.ylabel('CellTypist confidence')

plt.tight_layout()
plt.savefig('tcell_pseudotime_vs_celltypist_confidence.png', dpi=300)
plt.close()

sc.pl.umap(adata, color = "dpt_pseudotime", save = "_dpt.png")