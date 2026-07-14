#!/usr/bin/env python
"""Diagnostic plots for a BigSur run on one CellRanger sample.

Produces individual PNGs plus a combined PDF in results/bigsur/<sample>_diagnostics/.

Usage:
    python scripts/bigsur_diagnostics.py SRX10264523
"""

import argparse
import os
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib.backends.backend_pdf import PdfPages
from scipy.cluster.hierarchy import leaves_list, linkage

import bigsur
from bigsur_contaminants import GROUPS, drop_contaminants

warnings.filterwarnings("ignore")

CELLRANGER = "/project2/cobrinik_1090/external_rb_scrnaseq_proj/output/cellranger"
OUTROOT = "/project2/cobrinik_1090/external_rb_scrnaseq_proj/results/bigsur"


def load_counts(sample, min_cells=3, exclude=None):
    mtx = os.path.join(CELLRANGER, sample, "outs", "filtered_feature_bc_matrix")
    adata = sc.read_10x_mtx(mtx, var_names="gene_symbols", make_unique=True)
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    if exclude is not None:
        # Must match the runner's filter, or the panels won't line up with the
        # correlation h5ad they are plotted against.
        adata = drop_contaminants(adata, tuple(exclude) or GROUPS, verbose=False)
    return adata


def panel_mcfano_vs_mean(ax, var, mean_expr):
    hv = var.highly_variable.values
    ax.scatter(mean_expr[~hv], var.mc_Fano.values[~hv], s=3, c="0.75",
               label=f"not HVG (n={(~hv).sum():,})", rasterized=True)
    ax.scatter(mean_expr[hv], var.mc_Fano.values[hv], s=4, c="#c2410c",
               label=f"HVG (n={hv.sum():,})", rasterized=True)
    ax.axhline(1.0, ls="--", lw=1, c="#0369a1")
    ax.text(ax.get_xlim()[0], 1.05, " mcFano = 1 (pure Poisson)", fontsize=7,
            c="#0369a1", va="bottom")
    top = var.sort_values("mc_Fano", ascending=False).head(12)
    for g in top.index:
        i = var.index.get_loc(g)
        ax.annotate(g, (mean_expr[i], var.mc_Fano.values[i]), fontsize=6,
                    xytext=(3, 0), textcoords="offset points")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("mean UMI per cell"); ax.set_ylabel("mcFano")
    ax.set_title("A. mcFano vs expression\n(the feature-selection fit)", fontsize=9, loc="left")
    ax.legend(fontsize=6, frameon=False, loc="upper left")


def panel_fano_vs_mcfano(ax, var, raw_fano):
    hv = var.highly_variable.values
    ax.scatter(raw_fano[~hv], var.mc_Fano.values[~hv], s=3, c="0.75", rasterized=True)
    ax.scatter(raw_fano[hv], var.mc_Fano.values[hv], s=4, c="#c2410c", rasterized=True)
    lim = [min(raw_fano.min(), 1e-2), max(raw_fano.max(), var.mc_Fano.max())]
    ax.plot(lim, lim, ls="--", lw=1, c="0.4")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("uncorrected Fano (var/mean)"); ax.set_ylabel("mcFano")
    ax.set_title("B. Effect of BigSur's correction\n(points below the line are deflated)",
                 fontsize=9, loc="left")


def panel_pval_hist(ax, var):
    ax.hist(var.p_value.values, bins=50, color="#0369a1", edgecolor="white", lw=0.3)
    ax.set_xlabel("raw p-value"); ax.set_ylabel("genes")
    n_sig = int((var.FDR_adj_pvalue < 0.05).sum())
    ax.set_title(f"C. p-value calibration\n{n_sig:,} genes at FDR < 0.05", fontsize=9, loc="left")
    ax.set_yscale("log")


def panel_hvg_overlap(ax, adata, var, n_top):
    a2 = adata.copy()
    sc.pp.highly_variable_genes(a2, n_top_genes=n_top, flavor="seurat_v3")
    scanpy_hvg = set(a2.var_names[a2.var.highly_variable])
    bigsur_hvg = set(var.index[var.highly_variable])
    both = bigsur_hvg & scanpy_hvg
    ax.bar(["BigSur\nonly", "both", "scanpy\nonly"],
           [len(bigsur_hvg - scanpy_hvg), len(both), len(scanpy_hvg - bigsur_hvg)],
           color=["#c2410c", "#7c3aed", "#0369a1"])
    for i, v in enumerate([len(bigsur_hvg - scanpy_hvg), len(both),
                           len(scanpy_hvg - bigsur_hvg)]):
        ax.text(i, v, f"{v:,}", ha="center", va="bottom", fontsize=8)
    ax.set_ylabel("genes")
    jac = len(both) / len(bigsur_hvg | scanpy_hvg)
    ax.set_title(f"D. BigSur vs scanpy seurat_v3\n(top {n_top:,} each; Jaccard {jac:.2f})",
                 fontsize=9, loc="left")
    return bigsur_hvg, scanpy_hvg


def panel_mcpcc_dist(ax, mc, pv):
    tril = np.tril_indices_from(mc, k=-1)
    vals, pvals = mc[tril], pv[tril]
    sig = pvals < 0.05
    ax.hist(vals, bins=100, color="0.8", label=f"all pairs (n={vals.size:,})")
    ax.hist(vals[sig], bins=100, color="#c2410c",
            label=f"BH < 0.05 (n={sig.sum():,})")
    ax.set_yscale("log")
    ax.set_xlabel("mcPCC"); ax.set_ylabel("gene pairs")
    ax.set_title("E. Correlation strength\n(significant pairs are the tails)",
                 fontsize=9, loc="left")
    ax.legend(fontsize=6, frameon=False)


def panel_corr_heatmap(ax, mc, pv, names, n=60):
    sig = (pv < 0.05) & (pv > 0)
    full_sig = sig + sig.T
    degree = full_sig.sum(0)
    top = np.argsort(-degree)[:n]
    full_mc = mc + mc.T
    np.fill_diagonal(full_mc, 1.0)
    sub = full_mc[np.ix_(top, top)]
    order = leaves_list(linkage(sub, method="average"))
    sub = sub[np.ix_(order, order)]
    lab = names[top][order]
    im = ax.imshow(sub, cmap="RdBu_r", vmin=-1, vmax=1, interpolation="nearest")
    ax.set_xticks(range(n)); ax.set_xticklabels(lab, rotation=90, fontsize=3.5)
    ax.set_yticks(range(n)); ax.set_yticklabels(lab, fontsize=3.5)
    ax.set_title(f"F. mcPCC modules\n(top {n} genes by significant partners)",
                 fontsize=9, loc="left")
    plt.colorbar(im, ax=ax, fraction=0.046, label="mcPCC")
    return degree


def panel_degree(ax, degree):
    ax.hist(degree, bins=60, color="#0369a1", edgecolor="white", lw=0.3)
    ax.set_xlabel("significantly correlated partners per gene")
    ax.set_ylabel("genes")
    ax.set_title(f"G. Correlation network degree\nmedian {np.median(degree):.0f} partners",
                 fontsize=9, loc="left")


def umap_from_genes(adata, genes, key):
    a = adata[:, list(genes)].copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.scale(a, max_value=10)
    sc.tl.pca(a, n_comps=30)
    sc.pp.neighbors(a, n_neighbors=15)
    sc.tl.umap(a)
    sc.tl.leiden(a, key_added=key, resolution=1.0, flavor="igraph", n_iterations=2,
                 directed=False)
    return a


def main():
    p = argparse.ArgumentParser()
    p.add_argument("sample")
    p.add_argument("--min-cells", type=int, default=3)
    p.add_argument("--outroot", default=OUTROOT)
    p.add_argument("--exclude-contaminants", nargs="*", choices=GROUPS, default=None,
                   metavar="GROUP",
                   help="must match whatever bigsur_cellranger.py was run with")
    args = p.parse_args()

    outdir = os.path.join(args.outroot, f"{args.sample}_diagnostics")
    os.makedirs(outdir, exist_ok=True)

    print("loading counts", flush=True)
    adata = load_counts(args.sample, args.min_cells, args.exclude_contaminants)
    X = adata.X.toarray().astype(np.float64)
    mean_expr = X.mean(0)
    raw_fano = X.var(0) / np.maximum(mean_expr, 1e-12)

    print("feature selection", flush=True)
    bigsur.mcfano_feature_selection(adata, layer="X", min_mcfano_cutoff=0.9,
                                    p_val_cutoff=0.05, verbose=0)
    var = adata.var.copy()
    n_hvg = int(var.highly_variable.sum())
    print(f"{n_hvg} HVGs, CV={adata.uns['CV_for_mc_Fano_fit']:.3f}", flush=True)

    corr_h5 = os.path.join(args.outroot, f"{args.sample}_correlations.h5ad")
    print(f"loading {corr_h5}", flush=True)
    sub = sc.read_h5ad(corr_h5)
    mc = sub.varm["mcPCCs"].toarray()
    pv = sub.varm["BH-corrected p-values of mcPCCs"].toarray()
    names = sub.var_names.values

    figs = {}

    fig, ax = plt.subplots(figsize=(5.5, 4.5)); panel_mcfano_vs_mean(ax, var, mean_expr)
    fig.tight_layout(); figs["A_mcfano_vs_mean"] = fig

    fig, ax = plt.subplots(figsize=(5.5, 4.5)); panel_fano_vs_mcfano(ax, var, raw_fano)
    fig.tight_layout(); figs["B_fano_correction"] = fig

    fig, ax = plt.subplots(figsize=(5.5, 4.5)); panel_pval_hist(ax, var)
    fig.tight_layout(); figs["C_pvalue_calibration"] = fig

    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    bigsur_hvg, scanpy_hvg = panel_hvg_overlap(ax, adata, var, n_hvg)
    fig.tight_layout(); figs["D_hvg_overlap"] = fig

    fig, ax = plt.subplots(figsize=(5.5, 4.5)); panel_mcpcc_dist(ax, mc, pv)
    fig.tight_layout(); figs["E_mcpcc_distribution"] = fig

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    degree = panel_corr_heatmap(ax, mc, pv, names)
    fig.tight_layout(); figs["F_mcpcc_modules"] = fig

    fig, ax = plt.subplots(figsize=(5.5, 4.5)); panel_degree(ax, degree)
    fig.tight_layout(); figs["G_network_degree"] = fig

    print("UMAP on BigSur HVGs", flush=True)
    ab = umap_from_genes(adata, bigsur_hvg, "leiden_bigsur")
    print("UMAP on scanpy HVGs", flush=True)
    asc = umap_from_genes(adata, scanpy_hvg, "leiden_scanpy")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    sc.pl.umap(ab, color="leiden_bigsur", ax=axes[0], show=False, frameon=False,
               legend_loc="on data", legend_fontsize=6, title="")
    axes[0].set_title(f"H. UMAP / leiden on BigSur HVGs (n={len(bigsur_hvg):,})",
                      fontsize=9, loc="left")
    sc.pl.umap(asc, color="leiden_scanpy", ax=axes[1], show=False, frameon=False,
               legend_loc="on data", legend_fontsize=6, title="")
    axes[1].set_title(f"I. UMAP / leiden on scanpy HVGs (n={len(scanpy_hvg):,})",
                      fontsize=9, loc="left")
    fig.tight_layout(); figs["H_umap_comparison"] = fig

    # Marker expression on the BigSur UMAP: the modules the correlations flagged.
    markers = [g for g in ["C1QB", "CD74", "SPP1", "RCVRN", "CRX", "MKI67"]
               if g in adata.var_names]
    # Log-normalize for display: raw UMIs are dominated by a few high-count cells,
    # which flattens the colour scale.
    norm = adata.copy()
    sc.pp.normalize_total(norm, target_sum=1e4)
    sc.pp.log1p(norm)
    ab_full = ab.copy()
    # Suffix the obs columns: a bare gene name would be ambiguous with .var_names.
    cols = {g: f"{g}_lognorm" for g in markers}
    ab_full.obs = ab_full.obs.join(
        pd.DataFrame(np.asarray(norm[:, markers].X.todense()), index=norm.obs_names,
                     columns=[cols[g] for g in markers]))
    fig, axes = plt.subplots(2, 3, figsize=(11, 7))
    for axm, g in zip(axes.ravel(), markers):
        sc.pl.umap(ab_full, color=cols[g], ax=axm, show=False, frameon=False,
                   cmap="viridis", title=g, vmax="p99", colorbar_loc="right")
    fig.suptitle("J. Marker genes on the BigSur UMAP (log-normalized)", fontsize=9,
                 x=0.02, ha="left")
    fig.tight_layout(); figs["J_markers"] = fig

    pdf_path = os.path.join(outdir, f"{args.sample}_bigsur_diagnostics.pdf")
    with PdfPages(pdf_path) as pdf:
        for name, fig in figs.items():
            fig.savefig(os.path.join(outdir, f"{name}.png"), dpi=150)
            pdf.savefig(fig)
    print(f"wrote {pdf_path} and {len(figs)} PNGs to {outdir}", flush=True)


if __name__ == "__main__":
    main()
