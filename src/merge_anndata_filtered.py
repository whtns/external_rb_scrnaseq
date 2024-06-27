#!/usr/bin/env python 

myproject = "full"

# imports------------------------------

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc
sc.settings.figdir = myproject
import scvi
from pathlib import Path
import anndata
import pandas as pd
from matplotlib.pyplot import rc_context
import os 
import glob
import numpy as np
from scvi.model.utils import mde
import pymde
import scrb
import seaborn as sns

adata_paths = {
  "SRR13884242" : "output/scanpy/SRR13884242_filtered_scanpy.h5ad",
  "SRR13884243" : "output/scanpy/SRR13884243_filtered_scanpy.h5ad",
  "SRR13884247" : "output/scanpy/SRR13884247_filtered_scanpy.h5ad",
  "SRR13884249" : "output/scanpy/SRR13884249_filtered_scanpy.h5ad",
  "SRR14800534" : "output/scanpy/SRR14800534_filtered_scanpy.h5ad",
  "SRR14800535" : "output/scanpy/SRR14800535_filtered_scanpy.h5ad",
  "SRR14800536" : "output/scanpy/SRR14800536_filtered_scanpy.h5ad",
  "SRR14800540" : "output/scanpy/SRR14800540_filtered_scanpy.h5ad",
  "SRR14800541" : "output/scanpy/SRR14800541_filtered_scanpy.h5ad",
  "SRR14800543" : "output/scanpy/SRR14800543_filtered_scanpy.h5ad",
  "SRR17960481" : "output/scanpy/SRR17960481_filtered_scanpy.h5ad",
  "SRR17960484" : "output/scanpy/SRR17960484_filtered_scanpy.h5ad"
}


# load adatas------------------------------

adatas = dict()
for k,v in adata_paths.items():
  adatas[k] = sc.read_h5ad(adata_paths[k])
  adatas[k].obs["sample_id"] = k
  adatas[k].obs.abbreviation = adatas[k].obs.abbreviation.astype('category')

combined_adata = anndata.concat(adatas, label="sample_id")
combined_adata.raw = combined_adata  # keep full dimension safe

combined_adata = combined_adata[combined_adata.obs['percent.mt'] < 5,:]

# sc.pp.highly_variable_genes(
#     combined_adata,
#     flavor="seurat_v3",
#     n_top_genes=2000,
#     layer="counts",
#     batch_key="sample_id",
#     subset=True,
# )

# sc.pp.highly_variable_genes(
#     combined_adata
# )

# scvi ------------------------------
combined_adata.obs = combined_adata.obs.drop("orig.ident", 1)
# combined_adata.write_h5ad("output/scanpy/combined_adata_filtered0.h5ad")

scvi.model.SCVI.setup_anndata(combined_adata, batch_key="sample_id")

vae = scvi.model.SCVI(combined_adata, n_layers=2, n_latent=30, gene_likelihood="nb")

vae.train()

vae.save("my_model/")

vae = scvi.model.SCVI.load("my_model/", adata=combined_adata)

combined_adata.obsm["X_scVI"] = vae.get_latent_representation()

# regress cell cycle genes -------------------------------
cell_cycle_genes = [x.strip() for x in open("data/regev_lab_cell_cycle_genes.txt")]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in combined_adata.var_names]

sc.tl.score_genes_cell_cycle(combined_adata, s_genes=s_genes, g2m_genes=g2m_genes)

# sc.pp.regress_out(combined_adata, ['S_score', 'G2M_score'])
# sc.pp.scale(combined_adata)

sc.pp.log1p(combined_adata)

sc.pp.neighbors(combined_adata, use_rep="X_scVI")
sc.tl.leiden(combined_adata, resolution = 0.3)
combined_adata.obsm["X_mde"] = mde(combined_adata.obsm["X_scVI"])

combined_adata.write_h5ad("output/scanpy/combined_adata_filtered.h5ad")

# proceeed ------------------------------

combined_adata = sc.read_h5ad("output/scanpy/combined_adata_filtered.h5ad")
sc.pp.log1p(combined_adata)

combined_adata.obs.scna = combined_adata.obs.scna.replace("", "none")

neg_cell_indices = np.where((combined_adata.obsm['X_mde'] < -8))
neg_cell_indices = list(set(neg_cell_indices[0].tolist()))

pos_cell_indices = np.where((combined_adata.obsm['X_mde'] > 6))
pos_cell_indices = list(set(pos_cell_indices[0].tolist()))

bad_cell_names = combined_adata.obs.index[neg_cell_indices + pos_cell_indices]
filtered_combined_adata = combined_adata[~combined_adata.obs.index.isin(bad_cell_names)].copy()

# plot genes ------------------------------
scrb.plot_my_embedding(filtered_combined_adata, "TFF1")
scrb.plot_my_embedding(filtered_combined_adata, "DST")
scrb.plot_my_embedding(filtered_combined_adata, "GAP43")
scrb.plot_my_embedding(filtered_combined_adata, "PRAME")
scrb.plot_my_embedding(filtered_combined_adata, "GSTP1")
scrb.plot_my_embedding(filtered_combined_adata, "ARR3")
scrb.plot_my_embedding(filtered_combined_adata, "NRXN1")

# percent.mt------------------------------
scrb.plot_my_embedding(filtered_combined_adata, "percent.mt")

# phase ------------------------------
scrb.plot_my_embedding(filtered_combined_adata, "phase")
scrb.plot_my_embedding(filtered_combined_adata, "S.Score")
scrb.plot_my_embedding(filtered_combined_adata, "G2M.Score")

# samples ------------------------------
scrb.plot_my_embedding(filtered_combined_adata, "sample_id")

# clusters------------------------------
sc.tl.leiden(filtered_combined_adata, resolution = 0.5)
scrb.plot_my_embedding(filtered_combined_adata, "leiden")
sc.tl.rank_genes_groups(filtered_combined_adata, 'leiden', method='t-test')
sc.tl.dendrogram(filtered_combined_adata, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(filtered_combined_adata, n_genes=3, save = "marker_genes_dotplot.pdf", swap_axes = True)
sc.pl.rank_genes_groups(filtered_combined_adata, n_genes=8, sharey=False, fontsize = 16, save = "marker_genes.pdf", ncols = 2)

g = sns.FacetGrid(filtered_combined_adata.obs, col="leiden", hue="leiden", col_wrap = 3)
g.map(sns.scatterplot, "S_score", "G2M_score")
g.savefig("results/combined_phase_facet.png") 

scrb.plot_my_embedding(filtered_combined_adata, "sample_id")

filtered_combined_adata.obs['clone_opt'] = filtered_combined_adata.obs['clone_opt'].astype('category')
scrb.plot_my_embedding(filtered_combined_adata, "clone_opt")

scrb.plot_my_embedding(filtered_combined_adata, "scna")

# preprocess ------------------------------

filtered_combined_adata.obs_names_make_unique()

sc.pp.filter_cells(filtered_combined_adata, min_genes=200)
sc.pp.filter_genes(filtered_combined_adata, min_cells=3)

filtered_combined_adata.var['mt'] = filtered_combined_adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(filtered_combined_adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

filtered_combined_adata = filtered_combined_adata[filtered_combined_adata.obs.n_genes_by_counts < 2500, :]
filtered_combined_adata = filtered_combined_adata[filtered_combined_adata.obs.pct_counts_mt < 30, :]

sc.pp.normalize_total(filtered_combined_adata, target_sum=1e4)

sc.pp.log1p(filtered_combined_adata)

sc.pp.highly_variable_genes(filtered_combined_adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

filtered_combined_adata = filtered_combined_adata[:, filtered_combined_adata.var.highly_variable]

sc.pp.regress_out(filtered_combined_adata, ['total_counts', 'pct_counts_mt'])

scrb.check_resolution(filtered_combined_adata, resolution = 0.3)



sc.pl.heatmap(filtered_combined_adata, markers, groupby='bulk_labels', swap_axes=True)


# sc.pp.scale(combined_adata, max_value=10)
# sc.tl.pca(combined_adata, svd_solver='arpack')
# sc.pl.pca_variance_ratio(combined_adata, log=True)
# sc.pp.neighbors(combined_adata, n_neighbors=10, n_pcs=40)
# sc.tl.leiden(combined_adata)
# combined_adata.write_h5ad("output/scanpy/combined_adata_filtered_processed.h5ad")
