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
import anndata2ri
import rpy2.robjects.packages as rpackages
import scvi
import rpy2
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri, pandas2ri, conversion
import scipy.sparse as sp
import seaborn as sns
from scipy.stats import median_abs_deviation


# Load R libraries
SoupX = importr('SoupX')
Seurat = importr('Seurat')

def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

batch_id = sys.argv[1]
directory = sys.argv[2]
directory = directory.strip()
batch_id = batch_id.strip()
ig_directory = "10X_GEX_CR9_0_1_IMGT_IgTCR"
if directory == "Rhesus_DukeU19Allo_Luo1":
  ig_directory = "10X_GEX_pooled_CR9_0_1_IMGT_IgTCR"
elif directory == "Rhesus_TulaneCure":
  ig_directory = "10X_GEX5P_CR9_0_1_IMGT_IgTCR"
filePath_raw = f"/peng_2/peng_lab/results/ongoing/{directory}/cellranger/{ig_directory}/{batch_id}/outs/raw_feature_bc_matrix"
filePath_filtered = f"/peng_2/peng_lab/results/ongoing/{directory}/cellranger/{ig_directory}/{batch_id}/outs/filtered_feature_bc_matrix"

ro.globalenv['raw_data_dir'] = filePath_raw
ro.globalenv['filtered_data_dir'] = filePath_filtered
ro.r('raw_counts <- Seurat::Read10X(raw_data_dir)')
ro.r('filtered_counts <- Seurat::Read10X(filtered_data_dir)')

#Runs DropletUtils and corrects for ambient mRNA contamination#

#1. emptyDrops() classifies barcodes as empty or not based on raw data
#2. filtered_counts contains barcodes that are confirmed not empty by EmptyDroplets
#3. SoupX is ran with the droplet object being the raw counts (before empty droplet removal) and cell object being the filtered (after)
#4. Must run clustering and PCA on pre-Seurat Object before adjusting counts
ro.r('''
library(Seurat)
library(DropletUtils)
library(SoupX)
soup <- SoupX::SoupChannel(tod = raw_counts, toc = filtered_counts)
seurat_obj <- CreateSeuratObject(counts = filtered_counts)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
clusters <- seurat_obj$seurat_clusters
names(clusters) <- colnames(seurat_obj)
soup <- SoupX::setClusters(soup, clusters)
''')


ro.r('''png("soup.png")''')
ro.r('soup <- SoupX::autoEstCont(soup, forceAccept=TRUE)')
ro.r('corrected_counts <- SoupX::adjustCounts(soup, roundToInt=TRUE)')
ro.r('''dev.off()''')

corrected_counts = ro.r['corrected_counts']
x = np.array(corrected_counts.do_slot('x'))        
i = np.array(corrected_counts.do_slot('i'))         
p = np.array(corrected_counts.do_slot('p'))         
dims = tuple(corrected_counts.do_slot('Dim'))    

corrected_counts_csr = sp.csc_matrix((x, i, p), shape=dims)
corrected_counts_csr = corrected_counts_csr.tocsr() 
corrected_counts_csr = corrected_counts_csr.T
cell_names = list(ro.r('colnames(corrected_counts)'))
gene_names = list(ro.r('rownames(corrected_counts)'))

# Convert to AnnData
adata = sc.AnnData(
    X=corrected_counts_csr,
    obs=pd.DataFrame(index=cell_names),
    var=pd.DataFrame(index=gene_names)
)

root_1 = "/peng_2/peng_lab/projects/TulaneCure/ongoing"
mt_gene_patterns = pd.read_csv("/peng_2/peng_lab/results/ongoing/Rhesus_TulaneCure/QC/rhesusMitoGeneNames.txt", header=None)[0].apply(lambda x: f"^{x}$").str.cat(sep="|")
mt_gene_patterns = mt_gene_patterns.split("|")
mito_genes = []
loc_metadata = pd.read_csv("/peng_2/peng_lab/ngs_data/private/rhesus_atlas/loc_genes_filter_out.txt", sep = "\t", header = 0)
loc_genes = loc_metadata["Gene"].tolist()
loc_mito_genes = loc_metadata.loc[loc_metadata["Reason"] == "Mitochondrial", "Gene"].tolist()
loc_ribo_genes = loc_metadata.loc[loc_metadata["Reason"] == "Ribosomal", "Gene"].tolist()
housekeeping_genes = pd.read_csv("/peng_2/peng_lab/ngs_data/private/rhesus_atlas/housekeeping_genes.txt", sep = "\t", header = 0)
housekeeping_genes = housekeeping_genes.iloc[:,0].tolist()
ribo_genes = [gene for gene in loc_ribo_genes if gene in adata.var_names]



#Remove cells dynamically based on high mitochondrial content, high hemoglobin expression, low counts, zero housekeeping gene expression #
for gene in adata.var_names:
    for pattern in mt_gene_patterns:
        if re.search(pattern, gene):
            mito_genes.append(gene)
    for pattern in loc_mito_genes:
        if re.search(pattern, gene):
            mito_genes.append(gene)
adata.var_names = [f"MT-{gene}" if gene in mito_genes else gene for gene in adata.var_names]
adata.var['mt'] = adata.var_names.str.startswith('MT-') 
adata.var["ribo"] = adata.var_names.isin(ribo_genes)
hb_genes = ["HBA", "HBB", "HBD", "HBG1", "HBG2", 'HBQ1', 'HBP1']
hb_genes_present = [gene for gene in hb_genes if gene in adata.var_names]
adata.var["hb"] = [gene in hb_genes_present for gene in adata.var_names]
adata.var["housekeeping"] = adata.var_names.isin(housekeeping_genes)
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb", "housekeeping"], percent_top = [20,40,60,80,100],  inplace=True, log1p = True)

#Cell is classified as an outlier if any of the conditions are true


adata.obs["outlier"] = (
      is_outlier(adata, "log1p_total_counts", 5)
      | is_outlier(adata, "log1p_n_genes_by_counts", 5)
      | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
)

print("Number of cells pre-QC:", adata.obs.shape[0], flush=True)
print("Number of outliers:", adata.obs["outlier"].sum(), flush = True)


#MAD calculation for mitochondrial content
mt_percent_sorted = np.sort(adata.obs['pct_counts_mt'])
median = np.median(mt_percent_sorted)
mad = np.median(np.abs(mt_percent_sorted - median))
threshold = median + 3*mad
#If threshold is too small, make it 5% mitochondrial expression
threshold = max(5, threshold)
#If threshold is too large, make it 10% mitochondrial expression
threshold = min(10, threshold)

#==Plot log-log barcode plot to show outliers ===#


df = adata.obs.copy()
df['barcode'] = df.index
df['total_counts'] = adata.obs['total_counts']
df = df.sort_values('total_counts', ascending=False).reset_index(drop=True)
df['rank'] = np.arange(1, len(df) + 1)

# Prepare colors (optional: set color map)
palette = {False: 'blue', True: 'red'}

plt.figure(figsize=(10, 5))
# Plot using scatter, with log-log scale
sns.scatterplot(
    x='rank',
    y='total_counts',
    hue='outlier',
    data=df,
    palette=palette,
    s=8,
    linewidth=0,
    alpha=0.7,
    legend='full'
)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("Barcode rank (log scale)")
plt.ylabel("Total counts (log scale)")
plt.title(f"{batch_id} Log Barcode-Barcode Plot")
plt.legend(title="Outlier")
plt.tight_layout()
plt.savefig(f"/peng_2/peng_lab/projects/rhesus_atlas/raw/{directory}/figures/{batch_id}_log_plot.png")

#=============Removes outliers================
adata = adata[~adata.obs.outlier, :].copy()
adata = adata[adata.obs['pct_counts_mt'] < threshold, :].copy()

#==============Hard Treshold Filtering=============#
adata = adata[adata.obs["total_counts"] > 500, :].copy()
adata = adata[adata.obs["total_counts"] < 1e4, :].copy()
sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells=3)
mirna_mask = adata.var_names.str.upper().str.contains("MIR")
adata = adata[:, ~mirna_mask].copy()
adata = adata[adata.obs['pct_counts_hb'] < 20, :].copy()
adata = adata[adata.obs.n_genes_by_counts <= 2500, :]
adata = adata[adata.obs["pct_counts_housekeeping"] > 0, :].copy()
adata.obs["batch_id"] = batch_id
adata.obs["source"] = directory
adata.write_h5ad(f"/peng_2/peng_lab/projects/rhesus_atlas/raw/{directory}/{batch_id}.h5ad")
print("Number of cells post QC:", adata.obs.shape[0], flush = True)
os.chdir('..')
#==============Remove Prevalent Doublets=============#
subprocess.run(f'python run_scDblFinder.py /peng_2/peng_lab/ngs_data/private/rhesus_atlas/{directory}/{batch_id}.h5ad', shell = True)
subprocess.run(f'python run_doublet_detector.py /peng_2/peng_lab/ngs_data/private/rhesus_atlas/{directory}/{batch_id}.h5ad', shell = True)
