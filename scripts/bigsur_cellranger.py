#!/usr/bin/env python
"""Run BigSur mcFano feature selection (and optionally correlations) on a
CellRanger filtered_feature_bc_matrix.

BigSur needs raw integer UMIs with no all-zero genes, which is exactly what
CellRanger's filtered matrix provides once empty genes are dropped.

Usage:
    python scripts/bigsur_cellranger.py SRX10264523 [--correlations] [--n-corr-genes 2000]
"""

import argparse
import os
import time

import numpy as np
import scanpy as sc

import bigsur
from bigsur_contaminants import GROUPS, drop_contaminants

CELLRANGER = "/project2/cobrinik_1090/external_rb_scrnaseq_proj/output/cellranger"
OUTDIR = "/project2/cobrinik_1090/external_rb_scrnaseq_proj/results/bigsur"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("sample")
    p.add_argument("--min-cells", type=int, default=3)
    p.add_argument("--correlations", action="store_true")
    p.add_argument("--n-corr-genes", type=int, default=2000,
                   help="correlations are genes x genes, so restrict to the top-N HVGs")
    p.add_argument("--max-trials", type=int, default=20000)
    p.add_argument("--outdir", default=OUTDIR)
    p.add_argument("--exclude-contaminants", nargs="*", choices=GROUPS, default=None,
                   metavar="GROUP",
                   help="drop contaminant families before feature selection; "
                        f"bare flag drops all of {', '.join(GROUPS)}")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    mtx_dir = os.path.join(CELLRANGER, args.sample, "outs", "filtered_feature_bc_matrix")

    adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True)
    adata.var_names_make_unique()
    print(f"{args.sample}: {adata.n_obs} cells x {adata.n_vars} genes (raw)", flush=True)

    # QC rejects all-zero genes; min_cells also keeps the dense matrix manageable.
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"after min_cells={args.min_cells}: {adata.n_vars} genes; "
          f"dense size {adata.n_obs * adata.n_vars * 8 / 1e9:.1f} GB", flush=True)

    if args.exclude_contaminants is not None:
        groups = tuple(args.exclude_contaminants) or GROUPS
        print(f"excluding contaminant families: {', '.join(groups)}", flush=True)
        adata = drop_contaminants(adata, groups)
        print(f"after contaminant filter: {adata.n_vars} genes", flush=True)

    t0 = time.time()
    bigsur.mcfano_feature_selection(adata, layer="X", min_mcfano_cutoff=0.9,
                                    p_val_cutoff=0.05, verbose=1)
    print(f"feature selection: {time.time() - t0:.1f}s | fitted CV "
          f"{adata.uns['CV_for_mc_Fano_fit']:.3f} | "
          f"{int(adata.var.highly_variable.sum())} HVGs", flush=True)

    fs = adata.var[["mc_Fano", "p_value", "FDR_adj_pvalue", "highly_variable"]]
    fs = fs.sort_values("mc_Fano", ascending=False)
    fs_path = os.path.join(args.outdir, f"{args.sample}_mcfano.csv")
    fs.to_csv(fs_path)
    print("top 15 genes by mcFano:\n", fs.head(15), flush=True)
    print(f"wrote {fs_path}", flush=True)

    if args.correlations:
        hvg = fs.index[fs.highly_variable][: args.n_corr_genes]
        sub = adata[:, hvg].copy()
        print(f"correlations on {sub.n_vars} genes "
              f"({sub.n_vars ** 2 * 8 / 1e9:.2f} GB output matrix)", flush=True)
        t0 = time.time()
        bigsur.calculate_correlations(sub, layer="X", max_trials=args.max_trials, verbose=1)
        print(f"correlations: {time.time() - t0:.1f}s", flush=True)

        h5_path = os.path.join(args.outdir, f"{args.sample}_correlations.h5ad")
        sub.write_h5ad(h5_path)
        pv = sub.varm["BH-corrected p-values of mcPCCs"].toarray()
        mc = sub.varm["mcPCCs"].toarray()
        sig = (pv > 0) & (pv < 0.05)
        print(f"{int(sig.sum())} significant gene pairs (BH < 0.05)", flush=True)
        if sig.sum():
            i, j = np.unravel_index(np.argsort(-np.abs(mc * sig), axis=None)[:10], mc.shape)
            names = sub.var_names
            print("top correlated pairs:", flush=True)
            for a, b in zip(i, j):
                print(f"  {names[a]:>12} ~ {names[b]:<12} mcPCC={mc[a, b]:+.3f}", flush=True)
        print(f"wrote {h5_path}", flush=True)


if __name__ == "__main__":
    main()
