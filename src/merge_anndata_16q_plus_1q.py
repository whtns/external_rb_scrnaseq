import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scanpy as sc
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


adata_16q_plus_1q = sc.read_h5ad("output/scanpy/16q_plus_1q_adata.h5ad")

sc.pp.log1p(adata_16q_plus_1q)

sc.tl.leiden(adata_16q_plus_1q, resolution = 0.4)

adata_16q_plus_1q.obsm["X_mde"] = mde(adata_16q_plus_1q.obsm["X_scVI"])

sc.tl.rank_genes_groups(adata_16q_plus_1q, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata_16q_plus_1q, n_genes=25, sharey=False, save = "marker_genes.pdf")

adata_16q_plus_1q.obs['cluster_short'] = adata_16q_plus_1q.obs.cluster.str.replace(".*_", "")

adata_16q_plus_1q.obs['clone_opt'] = adata_16q_plus_1q.obs.clone_opt.astype('category')

# write obs data to file for plotting R ggplot2 ------------------------------
adata_16q_plus_1q.obs.to_csv("results/adata_16q_plus_1q_metadata.csv")


top_genes = list()
for i in range(0,14):
  top_genes = top_genes + adata_16q_plus_1q.uns['rank_genes_groups']['names'].field(i)[0:1].tolist()


interesting_genes  = ["RXRG", "DEK", "KIF14", "SOX4", "NEK2"]

top_genes = top_genes  + interesting_genes

sc.pl.matrixplot(
    adata_16q_plus_1q,
    var_names = top_genes,
    groupby ='leiden',
    dendrogram = True,
    log = True,
    save = "matrixplot.pdf"
)


def plot_my_embedding(color):
  sc.pl.embedding(
    adata_16q_plus_1q,
    basis="X_mde",
    color=color,
    frameon=False,
    ncols=1,
    save = f'_{color}.pdf'
)


plot_my_embedding("sample_id")

plot_my_embedding("leiden")

plot_my_embedding("cluster_short")

plot_my_embedding("clone_opt")


# cell cycle phase ------------------------------
cl = "https://raw.githubusercontent.com/theislab/scanpy_usage/master/"
cl += "180209_cell_cycle/data/regev_lab_cell_cycle_genes.txt"
cc_prots = pd.read_csv(cl, header=None).squeeze()
s_genes = cc_prots[:43].tolist()
g2m_genes = cc_prots[43:].tolist()


sc.tl.score_genes_cell_cycle(adata_16q_plus_1q, s_genes = s_genes, g2m_genes = g2m_genes)

plot_my_embedding("phase")

sc.pl.embedding(
    adata_16q_plus_1q,
    basis="X_mde",
    color=top_genes,
    frameon=False,
    ncols=1,
    save = "_top_markers.png", 
    use_raw = False
)


# end ------------------------------

leiden_clusters = adata_16q_plus_1q.obs.leiden.unique().tolist()
leiden_clusters.sort()

os.makedirs("figures/violin_adata_16q_plus_1q")
os.makedirs("figures/rank_genes_groups_clone_opt_adata_16q_plus_1q")

for i in leiden_clusters:
  minor_adata = adata_16q_plus_1q[adata_16q_plus_1q.obs.leiden == i,:]
  sc.tl.rank_genes_groups(minor_adata, 'clone_opt', method='t-test')
  sc.pl.rank_genes_groups(minor_adata, n_genes=10, sharey=False, save = f'_adata_16q_plus_1q/{i}_marker_genes.pdf')
  top_genes = list()
  
  n_clones = len(minor_adata.obs.clone_opt.unique())
  
  print(n_clones)
  
  for j in range(0,n_clones):
    top_genes = top_genes + minor_adata.uns['rank_genes_groups']['names'].field(j)[0:1].tolist()

  print(top_genes)
  sc.pl.violin(
    minor_adata,
    keys=top_genes,
    groupby='clone_opt',
    rotation=90,
    log = False,
    use_raw = True,
    save = f'_adata_16q_plus_1q/{i}.pdf')