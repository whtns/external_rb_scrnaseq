#!/usr/bin/env python 

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import patchworklib as pw
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
import re
from plotnine import ggplot, geom_bar, aes, coord_flip
from matplotlib import cm, colors
from pathlib import Path
import scipy.stats as stats
import janitor


def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(l, key=alphanum_key)
  
def plot_my_embedding(adata_input, color, basis="X_mde", label="", **kwargs):
  myadata = adata_input.copy()
  
  mypal = cm.tab20.colors
  
  sc.pl.embedding(
    myadata,
    basis=basis,
    color=color,
    frameon=False,
    ncols=1,
    palette=mypal,
    save = f'_{color}{label}.pdf',
    show=False,
    **kwargs)

def scanpy_plot_meta_by_sample(input_adata, color_by="leiden", split_by="sample_id", label="", **kwargs):
  split_by_vals = input_adata.obs[split_by].unique().tolist()
  split_by_vals = [x for x in split_by_vals if str(x) != 'nan']
  split_by_vals = natural_sort(split_by_vals)
  color_by_vals = input_adata.obs[color_by].cat.categories
  base_pal = dict(zip(color_by_vals, cm.tab20.colors[0:len(color_by_vals)]))
  
  myadata = input_adata.copy()
  for i,v in enumerate(split_by_vals):
    myadata = input_adata.copy()
    myadata.obs.loc[myadata.obs['sample_id'] != v,color_by] = np.NaN
    myadata.obs[color_by] = myadata.obs[color_by].cat.remove_unused_categories()
    used_color_by_values = myadata.obs[color_by].cat.categories
    mini_pal = {k:base_pal[k] for k in used_color_by_values if k in base_pal}
    sc.pl.embedding(myadata, basis="X_mde", color = color_by, frameon=False, ncols=1, \
    use_raw = False, na_in_legend = False, na_color = 'white', \
    size=20, title=f'{v}_{color_by}', \
    save=f'{v}_{color_by}_{label}.png', palette=mini_pal, **kwargs)
  
  return base_pal

def add_GT_opt(adata):
  myadata = adata.copy()
  
  seu_meta_paths = Path('output/seurat').glob('*embeddings*')
  seu_metas = dict()
  for i in seu_meta_paths:
    k = i.name.replace("_embeddings.csv", "")
    seu_metas[k] = pd.read_csv(i)
    seu_metas[k]["sample_id"] = k
    
  seu_meta_df = pd.concat([v for k,v in seu_metas.items()])
  seu_meta_df['cell_short'] = seu_meta_df.cell.str.replace("\\-1", "")
  seu_meta_df = seu_meta_df.merge(myadata.obs, how = "right", on = ["cell_short", "sample_id"])
  seu_meta_df.index = myadata.obs.index
  
  myadata.obs['GT_opt'] = seu_meta_df.GT_opt.astype('category')
  
  myadata.obs.GT_opt = myadata.obs.GT_opt.cat.add_categories("original").fillna("original")
  
  myadata.obs['GT_opt'] = myadata.obs['GT_opt'].astype('category')
  
  return myadata


# load full dataset -------------------------------

# adata_all = sc.read_h5ad("output/scanpy/scvi_all.h5ad")
# adata_all.uns['log1p']["base"] = None
# sc.pp.log1p(adata_all)

# load filtered dataset -------------------------------

adata_filtered = sc.read_h5ad("output/scanpy/adata_all_filtered.h5ad")

cluster_dict = {"glia 1":"7", "glia 2":"9", "ganglion":"12"}

adata_filtered = adata_filtered[~adata_filtered.obs.leiden.isin(cluster_dict.values()),:]


ct_markers = ["PDE6A", "RHO", "NR2E3", "OPN1LW", "OPN1MW", "OPN1SW", "RLBP1", 
"APOE", "CLU", "GFAP", "HLA-DPA1", "HLA-DRA", "C1QA", "VSX2", 
"TMEM215", "VSX1", "SNCG", "SLC17A6", "RBPMS", "CALB1", "CHAT", 
"GAD2"]

sc.pl.dotplot(adata_filtered, ct_markers, groupby='leiden', save = "test.pdf")

plot_my_embedding(adata_filtered, color = "leiden")

# plot_my_embedding(new_adata_filtered, color = "CD24")


adata_filtered.uns['log1p']["base"] = None
sc.pp.log1p(adata_filtered)

# adata_filtered_low_mito = adata_filtered[adata_filtered.obs.pct_counts_mt < 1, :]
# plot_my_embedding(adata_filtered, color="pct_counts_mt", vmax=5, title = "Percent mito reads")



# plot by sample ------------------------------

new_adata_filtered = add_GT_opt(adata_filtered)

# add clone simplifications 
def add_clone_simplifications(new_adata_filtered):
  
  clone_simplifications = pd.read_csv("results/new_adata_filtered_meta.csv", index_col=0)[["scna"]]
  # clone_simplifications.rename({'...1' : "cells"}, axis = 1, inplace = True)
  # clone_simplifications.index  = clone_simplifications.cells
  
  clone_simplifications = clone_simplifications.merge(new_adata_filtered.obs, left_index = True, right_index=True)
  
  clone_simplifications.scna = clone_simplifications.scna.fillna('none')
  
  new_adata_filtered.obs['scna'] = clone_simplifications.scna.astype('category')
  
  return new_adata_filtered

new_adata_filtered = add_clone_simplifications(new_adata_filtered)

# regress cell cycle genes -------------------------------
cell_cycle_genes = [x.strip() for x in open("data/regev_lab_cell_cycle_genes.txt")]

s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in new_adata_filtered.var_names]

sc.tl.score_genes_cell_cycle(new_adata_filtered, s_genes=s_genes, g2m_genes=g2m_genes)

sc.pp.regress_out(new_adata_filtered, ['S_score', 'G2M_score'])
sc.pp.scale(new_adata_filtered)

# new_adata_filtered.write_h5ad("output/scanpy/adata_all_filtered_regressed.h5ad")

regressed_adata = sc.read_h5ad("output/scanpy/adata_all_filtered_regressed.h5ad")

sc.tl.leiden(regressed_adata, resolution = 0.2)

plot_my_embedding(regressed_adata, color = "leiden")



def make_dotplot(adata_filtered, vmax = 40, label = 0.2):
  
  n_clusters = len(adata_filtered.uns['rank_genes_groups']['names'].dtype.names)
  
  top_genes = dict()
  
  cluster_names = adata_filtered.uns['rank_genes_groups']['names'].dtype.names
  
  for i in cluster_names:
    top_genes[i] = list(adata_filtered.uns['rank_genes_groups']['names'][i])[0:5]
  
  markers = list()
  for k,v in top_genes.items():
    markers.extend(v)
  
  # regex = re.compile(r'^MT-.*')
  # filtered_markers = [ele for ele in markers if not regex.match(ele)]
  # 
  # regex = re.compile(r'^HES6')
  # filtered_markers = [ele for ele in filtered_markers if not regex.match(ele)]
  
  # return top_genes
  # interesting_genes  = ["RXRG", "DEK", "KIF14", "SOX4", "NEK2"]
  # top_genes = top_genes  + interesting_genes
  
  sc.pl.dotplot(
    adata_filtered,
    var_names = markers,
    groupby='leiden',
    dendrogram=False,
    swap_axes = True,
    vmax = vmax,
    save = f'{label}_matrixplot.pdf')

def check_resolution(adata_filtered, resolution = 0.2):
  
  sc.tl.leiden(adata_filtered, resolution = resolution)
  cells_w_enough_leiden = adata_filtered.obs.groupby("leiden").filter(lambda x: len(x) > 10)
  adata_filtered = adata_filtered[adata_filtered.obs.index.isin(cells_w_enough_leiden.index), :]
  plot_my_embedding(adata_filtered, "leiden", title = "inter-sample clusters", label = f'_{resolution}')
  plot_my_embedding(adata_filtered, "clone_opt", title = "SCNA clones", label = f'_{resolution}')
  plot_my_embedding(adata_filtered, "scna", title = "SCNAs", label = f'_{resolution}')
  plot_my_embedding(adata_filtered, "phase", title = "cell cycle phase", label = f'_{resolution}')
  
  meta="leiden"
  
  scanpy_plot_meta_by_sample(adata_filtered, color_by="leiden", split_by = "sample_id", label=resolution)
  scanpy_plot_meta_by_sample(adata_filtered, "clone_opt", "sample_id", label=resolution)
  scanpy_plot_meta_by_sample(adata_filtered, "scna", "sample_id", label=resolution)
  scanpy_plot_meta_by_sample(adata_filtered, "phase", "sample_id", label=resolution)
  
  font = {'size'   : 12}
  
  plt.rc('font', **font)
  plt.rcParams["figure.figsize"] = (3,3)
  sc.tl.rank_genes_groups(adata_filtered, groupby='leiden', method='wilcoxon')
  sc.pl.rank_genes_groups(adata_filtered, n_genes=10, fontsize = 12, sharey=False, ncols = 2, save = f'{resolution}_leiden.pdf')
  
  make_dotplot(adata_filtered, label = resolution)
  
  metadata_path = f'results/adata_filtered_metadata_{resolution}.csv'
  
  adata_filtered.obs.to_csv(metadata_path)
  
  return(adata_filtered)

# dropped_samples = ["SRR14800543"]
# dropped_clones = [1.0]
# # 
# new_adata_filtered = new_adata_filtered[~(new_adata_filtered.obs.sample_id.isin(dropped_samples) & new_adata_filtered.obs.clone_opt.isin(dropped_clones)),:]

# check_resolution(new_adata_filtered, resolution = 0.25)

# scanpy_plot_meta_by_sample(new_adata_filtered_025, "cluster_short", "sample_id")
# scanpy_plot_meta_by_sample(new_adata_filtered, "clone_opt", "sample_id")
# scanpy_plot_meta_by_sample(new_adata_filtered, "scna", "sample_id")
# scanpy_plot_meta_by_sample(new_adata_filtered, "phase", "sample_id")

check_resolution(regressed_adata, resolution = 0.2)
