#!/usr/bin/env python 

# imports------------------------------

myproject = "1q_or_6p"

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
 
adata_paths = {
  "SRR13884249" : "output/scanpy/SRR13884249_filtered_scanpy.h5ad", 
  "SRR14800534" : "output/scanpy/SRR14800534_filtered_scanpy.h5ad", 
  "SRR14800536" : "output/scanpy/SRR14800536_filtered_scanpy.h5ad", 
  "SRR14800535" : "output/scanpy/SRR14800535_filtered_scanpy.h5ad", 
  "SRR17960484" : "output/scanpy/SRR17960484_filtered_scanpy.h5ad", 
  "SRR13884247" : "output/scanpy/SRR13884247_filtered_scanpy.h5ad", 
  "SRR17960481" : "output/scanpy/SRR17960481_filtered_scanpy.h5ad"
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

vae.save("my_model_1q_or_6p/", overwrite = True)

vae = scvi.model.SCVI.load("my_model_1q_or_6p/", adata=combined_adata)

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

combined_adata.write_h5ad("output/scanpy/combined_adata_filtered_1q_or_6p.h5ad")

combined_adata = sc.read_h5ad("output/scanpy/combined_adata_filtered_1q_or_6p.h5ad")
sc.pp.log1p(combined_adata)

combined_adata.obs.scna = combined_adata.obs.scna.replace("", "none")

neg_cell_indices = np.where((combined_adata.obsm['X_mde'] < -8))
neg_cell_indices = list(set(neg_cell_indices[0].tolist()))

pos_cell_indices = np.where((combined_adata.obsm['X_mde'] > 6))
pos_cell_indices = list(set(pos_cell_indices[0].tolist()))

bad_cell_names = combined_adata.obs.index[neg_cell_indices + pos_cell_indices]
filtered_combined_adata = combined_adata[~combined_adata.obs.index.isin(bad_cell_names)].copy()

scrb.plot_my_embedding(filtered_combined_adata, "TFF1")
scrb.plot_my_embedding(filtered_combined_adata, "DST")
scrb.plot_my_embedding(filtered_combined_adata, "GAP43")
scrb.plot_my_embedding(filtered_combined_adata, "PRAME")
scrb.plot_my_embedding(filtered_combined_adata, "GSTP1")
scrb.plot_my_embedding(filtered_combined_adata, "ARR3")
scrb.plot_my_embedding(filtered_combined_adata, "NRXN1")

scrb.plot_my_embedding(filtered_combined_adata, "scna")
scrb.plot_my_embedding(filtered_combined_adata, "percent.mt")
scrb.plot_my_embedding(filtered_combined_adata, "phase")
scrb.plot_my_embedding(filtered_combined_adata, "sample_id")
sc.tl.leiden(filtered_combined_adata, resolution = 0.3)
scrb.plot_my_embedding(filtered_combined_adata, "leiden")
sc.tl.rank_genes_groups(filtered_combined_adata, 'leiden', method='t-test')
sc.tl.dendrogram(filtered_combined_adata, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(filtered_combined_adata, n_genes=10, save = "marker_genes_dotplot.pdf")
sc.pl.rank_genes_groups(filtered_combined_adata, n_genes=25, sharey=False, save = "marker_genes.pdf")


