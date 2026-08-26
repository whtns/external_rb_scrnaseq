#!/usr/bin/env python
"""Issue #40, Stage B. Select overdispersed genes + initialize cNMF via the mosaicMPI
API, bypassing the `mosaicmpi set-parameters` CLI.

The CLI's `-m default_topn` path is broken in mosaicmpi 2.5.0: it splits the method
name on "_" and compares to "top_n", but the choice is spelled "topn", so it forwards
an unsupported `topn=` kwarg and raises TypeError. It also defaults to selecting ~10k
overdispersed genes (min_score=1.0), far more than cNMF wants. Here we call
`Dataset.select_overdispersed_genes(top_n=...)` directly (the correct kwarg) and then
`initialize_cnmf(...)`, mirroring what the CLI does after gene selection.

Usage: python src/mosaic/set_parameters.py <output_dir> <name> <top_n> <k_first> <k_last> <n_iter> <seed> [beta_loss]
The working h5ad is read/written at <output_dir>/<name>/<name>.h5ad (model-odg's output).
"""
import sys
import mosaicmpi


def main(output_dir, name, top_n, k_first, k_last, n_iter, seed, beta_loss):
    h5ad = f"{output_dir}/{name}/{name}.h5ad"
    ds = mosaicmpi.Dataset.from_h5ad(h5ad)
    ds.select_overdispersed_genes(overdispersion_metric="odscore", top_n=int(top_n))
    n_sel = int(ds.adata.var["selected"].sum()) if "selected" in ds.adata.var else -1
    print(f"selected {n_sel} overdispersed genes (top_n={top_n})")
    ds.initialize_cnmf(
        cnmf_output_dir=output_dir, cnmf_name=name,
        kvals=range(int(k_first), int(k_last) + 1),
        n_iter=int(n_iter), beta_loss=beta_loss, seed=int(seed),
    )
    ds.write_h5ad(h5ad)
    print(f"initialized cNMF: k={k_first}..{k_last}, n_iter={n_iter}, seed={seed}, "
          f"beta_loss={beta_loss}")


if __name__ == "__main__":
    if len(sys.argv) not in (8, 9):
        raise SystemExit(__doc__)
    beta = sys.argv[8] if len(sys.argv) == 9 else "kullback-leibler"
    main(*sys.argv[1:8], beta)
