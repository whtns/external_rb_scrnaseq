#!/usr/bin/env python 

# imports------------------------------

myproject = "1q_and_16q_regressed"

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
  "SRR14800534" : "output/scanpy/SRR14800534_filtered_scanpy.h5ad", 
  "SRR14800536" : "output/scanpy/SRR14800536_filtered_scanpy.h5ad", 
  "SRR14800535" : "output/scanpy/SRR14800535_filtered_scanpy.h5ad", 
}

# load adatas------------------------------

adatas = dict()
for k,v in adata_paths.items():
  adatas[k] = sc.read_h5ad(adata_paths[k])
  adatas[k].obs["sample_id"] = k
  adatas[k].obs.abbreviation = adatas[k].obs.abbreviation.astype('category')

combined_adata = anndata.concat(adatas, label="sample_id")
combined_adata.raw = combined_adata  # keep full dimension safe

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

# combined_adata = sc.read_h5ad("output/scanpy/combined_adata_filtered0.h5ad")

adata_cond = combined_adata.copy()

# condition integration on cell cycle effects  https://yoseflab.github.io/scvi-tools-reproducibility/scvi_covariates/ ------------------------------

cell_cycle_genes = [x.strip() for x in open("data/regev_lab_cell_cycle_genes.txt")]

nuisance_genes = [x for x in cell_cycle_genes if x in adata_cond.var.index.tolist()]

# then copy the expression of each nuisance gene into adata.obs where the key
# is the gene name
for g in nuisance_genes:
    exp = adata_cond[:,g].X
    adata_cond.obs[g] = exp.toarray()
    
# finally, remove the nuisance genes from the anndata
gene_subset = [g for g in combined_adata.var_names if g not in nuisance_genes]
adata_cond = adata_cond[:,gene_subset].copy()

# run setup_anndata with our list of nuisance genes as our continuous covariates
scvi.model.SCVI.setup_anndata(adata_cond,
                        batch_key = 'sample_id',
                        continuous_covariate_keys=nuisance_genes)

vae = scvi.model.SCVI(adata_cond, n_layers=2, n_latent=30, gene_likelihood="nb")

vae.train()

vae.save("my_model_1q_and_16q_regressed/", overwrite = True)

vae = scvi.model.SCVI.load("my_model_1q_and_16q_regressed/", adata=adata_cond)

adata_cond.obsm["X_scVI"] = vae.get_latent_representation()

# regress cell cycle genes -------------------------------
cell_cycle_genes = [x.strip() for x in open("data/regev_lab_cell_cycle_genes.txt")]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata_cond.var_names]

sc.tl.score_genes_cell_cycle(adata_cond, s_genes=s_genes, g2m_genes=g2m_genes)

# sc.pp.regress_out(adata_cond, ['S_score', 'G2M_score'])
# sc.pp.scale(adata_cond)

sc.pp.log1p(adata_cond)

sc.pp.neighbors(adata_cond, use_rep="X_scVI")
sc.tl.leiden(adata_cond, resolution = 0.3)
adata_cond.obsm["X_mde"] = mde(adata_cond.obsm["X_scVI"])

adata_cond.write_h5ad("output/scanpy/combined_adata_1q_and_16q_regressed.h5ad")

adata_cond = sc.read_h5ad("output/scanpy/combined_adata_1q_and_16q_regressed.h5ad")

sc.pp.log1p(adata_cond)

adata_cond.obs.scna = adata_cond.obs.scna.replace("", "none")

neg_cell_indices = np.where((adata_cond.obsm['X_mde'] < -8))
neg_cell_indices = list(set(neg_cell_indices[0].tolist()))

pos_cell_indices = np.where((adata_cond.obsm['X_mde'] > 6))
pos_cell_indices = list(set(pos_cell_indices[0].tolist()))


bad_cell_names = adata_cond.obs.index[neg_cell_indices + pos_cell_indices]
filtered_adata_cond = adata_cond[~adata_cond.obs.index.isin(bad_cell_names)].copy()


scrb.plot_my_embedding(filtered_adata_cond, "TFF1")

scrb.plot_my_embedding(filtered_adata_cond, "DST")
scrb.plot_my_embedding(filtered_adata_cond, "GAP43")
scrb.plot_my_embedding(filtered_adata_cond, "PRAME")
scrb.plot_my_embedding(filtered_adata_cond, "GSTP1")
scrb.plot_my_embedding(filtered_adata_cond, "MARCKSL1")

filtered_adata_cond = filtered_adata_cond[filtered_adata_cond.obs['percent.mt'] < 3,:]

scrb.plot_my_embedding(filtered_adata_cond, "percent.mt")
scrb.plot_my_embedding(filtered_adata_cond, "phase")
scrb.plot_my_embedding(filtered_adata_cond, "sample_id")
sc.tl.leiden(filtered_adata_cond, resolution = 0.3)
cells_w_enough_leiden = filtered_adata_cond.obs.groupby("leiden").filter(lambda x: len(x) > 10)
filtered_adata_cond = filtered_adata_cond[filtered_adata_cond.obs.index.isin(cells_w_enough_leiden.index), :]
scrb.plot_my_embedding(filtered_adata_cond, "leiden")
sc.tl.rank_genes_groups(filtered_adata_cond, 'leiden', method='wilcoxon')
sc.tl.dendrogram(filtered_adata_cond, groupby='leiden')
sc.pl.rank_genes_groups_dotplot(filtered_adata_cond, n_genes=10, save = "marker_genes_dotplot.pdf")
sc.pl.rank_genes_groups(filtered_adata_cond, n_genes=25, sharey=False, save = "marker_genes.pdf")

