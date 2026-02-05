import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
from matplotlib import rcParams
import os
import re
import matplotlib.backends.backend_pdf
import sys


def viewVDJOnUMAP(adata, label):
    rcParams.update({'figure.figsize': (8,10), 'axes.labelsize': 5, 'axes.titlesize': 5, 'legend.title_fontsize':5, 'legend.fontsize':5, 'legend.markerscale': 0.25})
    fig, axs = plt.subplots(3,3,figsize=(8,10))
    fig.suptitle(f"{label}", fontsize=10)
    sc.pl.umap(adata, color = "leiden", ax=axs[0,0], show = False, legend_loc = None)
    axs[0,0].set_title("leiden")
    sc.pl.umap(adata, color='pair.TRAB', ax=axs[1,0], show=False, palette={'Both': 'purple', 'VDJ': 'red', 'VJ': 'blue', 'ND': 'grey'}, title='TRAB')
    sc.pl.umap(adata, color='pair.TRGD', ax=axs[1,1], show=False, palette={'Both': 'purple', 'VDJ': 'red', 'VJ': 'blue', 'ND': 'grey'}, title='TRGD')
    dual_colors = {'Y': 'red','N': 'grey'}
    sc.pl.umap(adata, color='pair.Dual', ax=axs[2,0], show=False, palette=dual_colors, title='Dual Expression')
    sc.pl.umap(adata, color = "pair.Ig", ax=axs[0,1], show = False, palette={'Both': 'purple', 'VDJ': 'red', 'VJ': 'blue', 'ND': 'grey'}, title='Ig')
    sc.pl.umap(adata, color = "C_region", ax=axs[0,2], show = False, title='C_region')
    fig.delaxes(axs[2,1])
    plt.tight_layout()
    plt.show()

def vdjHelper(adata, vDjPairs_sampleID):
    if vDjPairs_sampleID["paired"] == True:
        if vDjPairs_sampleID["pairTyp"] == "TRAB":
            adata.obs.loc[vDjPairs_sampleID["barcode"], 'pair.TRAB'] = "Both"
            adata.obs.loc[vDjPairs_sampleID["barcode"], 'C_region'] = vDjPairs_sampleID["Cregion.VDJ"]
        elif vDjPairs_sampleID["pairTyp"] == "TRGD":
            adata.obs.loc[vDjPairs_sampleID["barcode"], 'pair.TRGD'] = "Both"
            adata.obs.loc[vDjPairs_sampleID["barcode"], 'C_region'] = vDjPairs_sampleID["Cregion.VDJ"]
        elif vDjPairs_sampleID["pairTyp"] == "Ig":
            adata.obs.loc[vDjPairs_sampleID["barcode"], 'pair.Ig'] = "Both"
            adata.obs.loc[vDjPairs_sampleID["barcode"], 'C_region'] = vDjPairs_sampleID["Cregion.VDJ"]
    elif vDjPairs_sampleID["paired"] == False:
        if vDjPairs_sampleID["pairTyp"] == "TRAB":
            if vDjPairs_sampleID["is_cell.VDJ"]:
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'C_region'] = vDjPairs_sampleID["Cregion.VDJ"]
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'pair.TRAB'] = "VDJ"
            else:
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'C_region'] = vDjPairs_sampleID["Cregion.VJ"]
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'pair.TRAB'] = "VJ"
        elif vDjPairs_sampleID["pairTyp"] == "TRGD":
            if vDjPairs_sampleID["is_cell.VDJ"]:
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'C_region'] = vDjPairs_sampleID["Cregion.VDJ"]
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'pair.TRGD'] = "VDJ"
            else:
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'C_region'] = vDjPairs_sampleID["Cregion.VJ"]
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'pair.TRGD'] = "VJ"
        elif vDjPairs_sampleID["pairTyp"] == "Ig":
            if vDjPairs_sampleID["is_cell.VDJ"]:
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'C_region'] = vDjPairs_sampleID["Cregion.VDJ"]
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'pair.Ig'] = "VDJ"
            else:
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'C_region'] = vDjPairs_sampleID["Cregion.VJ"]
                adata.obs.loc[vDjPairs_sampleID["barcode"], 'pair.Ig'] = "VJ"
    return adata


batch_id = sys.argv[1]
directory = sys.argv[2]
adata = sc.read_h5ad(f"/peng_2/peng_lab/ngs_data/private/rhesus_atlas/{directory}/{batch_id}_downstream.h5ad")
outdir = "/peng_2/peng_lab/projects/rhesus_atlas/raw/VDJ_overlay"

if not os.path.exists(outdir):
    os.makedirs(outdir)
    
adata.obs['pair.TRAB'] = "ND"
adata.obs['pair.TRGD'] = "ND"
adata.obs['pair.Ig'] = "ND"
adata.obs['pair.Dual'] = "N"

adata.obs['pair.TRAB'] = adata.obs['pair.TRAB'].astype('category')
adata.obs['pair.TRGD'] = adata.obs['pair.TRGD'].astype('category')
adata.obs['pair.Ig'] = adata.obs['pair.Ig'].astype('category')
adata.obs['pair.Dual'] = adata.obs['pair.Dual'].astype('category')

adata.obs["pair.Dual"] = adata.obs["pair.Dual"].cat.add_categories(["Y"])
adata.obs["pair.TRAB"] = adata.obs["pair.TRAB"].cat.add_categories(["Both", "VDJ", "VJ"])
adata.obs["pair.TRGD"] = adata.obs["pair.TRGD"].cat.add_categories(["Both", "VDJ", "VJ"])
adata.obs["pair.Ig"] = adata.obs["pair.Ig"].cat.add_categories(["Both", "VDJ", "VJ"])
adata.obs["C_region"] = "ND"

adata.obs["barcode_batch"] = adata.obs.index + "_" + adata.obs["batch_id"].astype(str)

if directory != "Rhesus_DukeKwun":
    indir_VDJ = f"/peng_2/peng_lab/projects/{directory}/ongoing/contigFilter/combined_CR9_0_1"
    if directory == "Rhesus_TulaneCure":
        indir_VDJ = f"/peng_2/peng_lab/projects/{directory}/ongoing/contigFilter"
        full_metadata = pd.read_csv("/peng_2/peng_lab/ngs_data/private/Rhesus_TulaneCure/Tulane_meta_data.txt", sep="\t")
    elif directory == "Rhesus_UPR":
        full_metadata = pd.read_csv("/peng_2/peng_lab/ngs_data/private/Rhesus_UPR/10x_UPR_Library_Info.tsv", sep = "\t")
    elif directory == "Rhesus_DukeU19Allo_Luo1":
        full_metadata = pd.read_csv("/peng_2/peng_lab/ngs_data/private/Rhesus_DukeU19Allo_Luo1/DukeU19Allo_Luo1_Library_Info.txt", sep="\t")
    elif directory == "Rhesus_DukeU19Allo_Kwun1":
        full_metadata = pd.read_csv("/peng_2/peng_lab/ngs_data/private/Rhesus_DukeU19Allo_Kwun1/DukeU19Allo_Kwun1_Library_Info_J.txt", sep = "\t")
    vdj_metadata = full_metadata[(full_metadata['libTyp'] != "GEX")]
    prefix = "Single"
    vDjPairs = pd.read_csv(os.path.join(indir_VDJ, f"{prefix}_vDjPairs.txt"), sep='\t', header=0, quoting=3)
    vDjPairs.rename(columns={"sampleID": "sampleID_VDJ"}, inplace=True)

vdj_to_gex_file = pd.read_csv("/peng_2/peng_lab/ngs_data/private/rhesus_atlas/rhesus_vdj_metadata.txt", sep = "\t", header = 0)

if directory != "Rhesus_DukeKwun":
    for i in vdj_to_gex_file["VDJ"]:
        if vdj_to_gex_file.loc[vdj_to_gex_file["VDJ"] == i, "GEX"].values[0] == batch_id:
            print(i)
            vdj_id = i
            vDjPairs_subset = vDjPairs[vDjPairs["libraryID.VDJ"] == vdj_id].copy()
            vDjPairs_subset["sampleID"] = batch_id
            for j in range(len(adata.obs_names)):
                sampleID = adata.obs["barcode_batch"].iloc[j].split("_")[0]
                vDjPairs_sampleID = vDjPairs_subset.loc[(vDjPairs_subset["barcode"] == sampleID)]
                if vDjPairs_sampleID.empty:
                    continue
                if len(vDjPairs_sampleID.index) > 1:
                    print("More than one")
                    adata.obs.loc[adata.obs.index[j], 'pair.Dual'] = "Y"
                    for k in range(len(vDjPairs_sampleID.index)):
                        adata = vdjHelper(adata, vDjPairs_sampleID.iloc[k])
                else:
                    print("Just one")
                    adata.obs.loc[adata.obs.index[j], 'pair.Dual'] = "N"
                    adata = vdjHelper(adata, vDjPairs_sampleID.iloc[0])

pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join(outdir, f"{batch_id}_VDJ_overlay_CR9_0_1.pdf"))
fig = viewVDJOnUMAP(adata, batch_id)
pdf.savefig(fig)
pdf.close()
adata.write_h5ad(f"/peng_2/peng_lab/ngs_data/private/rhesus_atlas/{directory}/{batch_id}_vdj.h5ad")
