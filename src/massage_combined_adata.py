#!/usr/bin/env python

# imports ------------------------------

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

# functions ------------------------------
def plot_dend(sample_id, adata):
  adata.obs.cluster = adata.obs.cluster.astype('category')
  sc.tl.dendrogram(adata, groupby = "cluster")
  dend = sc.pl.dendrogram(adata, groupby = "cluster", orientation = "left")
  dend.figure.set_size_inches(18, 13)
  dend.figure.savefig(f'{sample_id}_dend.pdf', dpi=200)

def plot_my_embedding(myadata, color):
  sc.pl.embedding(
    myadata,
    basis="X_mde",
    color=color,
    frameon=True,
    ncols=1,
    save = f'_{color}.pdf'
)

adata = sc.read_h5ad("output/scanpy/combined_adata_filtered.h5ad")

plot_my_embedding(adata, color = "phase")

# view embedding

bad_cell_indices = np.where((adata.obsm['X_mde'] < -6))
bad_cell_indices = list(set(bad_cell_indices[0].tolist()))
bad_cell_names = adata.obs.index[bad_cell_indices]
filtered_adata = adata[~adata.obs.index.isin(bad_cell_names)].copy()

plot_my_embedding(filtered_adata, color = "phase")

plot_my_embedding(filtered_adata, color = "leiden")




