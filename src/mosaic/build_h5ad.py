#!/usr/bin/env python
"""Issue #40, Stage B. Build a mosaicMPI Dataset .h5ad from the sparse counts export.

Uses mosaicMPI's own ``Dataset.from_df`` constructor (NOT a hand-rolled AnnData): it
sets up the .X / normalization bookkeeping that ``model-odg`` expects (X=normalized,
raw.X=counts). We hand it RAW COUNTS with ``is_normalized=False`` so mosaicMPI does the
TPM normalization for overdispersed-gene selection itself.

Input matrix is features x cells (MatrixMarket); Dataset.from_df wants observations x
features, so we transpose to cells x genes. ``.X`` is kept DENSE -- mosaicMPI densifies
internally at every step and raises on sparse .X, so sparsify is NOT used (the sparse
export is only to avoid a giant dense-TSV round-trip).

Usage: python src/mosaic/build_h5ad.py <in_dir> <out_h5ad>
  <in_dir> holds counts.mtx, barcodes.txt, genes.txt, metadata.tsv (from export_counts.R)
"""
import sys
import numpy as np
import pandas as pd
import scipy.io
import mosaicmpi


def main(in_dir: str, out_h5ad: str) -> None:
    mtx = scipy.io.mmread(f"{in_dir}/counts.mtx")          # features x cells, sparse
    barcodes = [l.strip() for l in open(f"{in_dir}/barcodes.txt")]
    genes = [l.strip() for l in open(f"{in_dir}/genes.txt")]

    # -> observations (cells) x features (genes), dense DataFrame.
    counts = pd.DataFrame(
        np.asarray(mtx.todense()).T, index=barcodes, columns=genes, dtype="float32"
    )
    if not np.allclose(counts.values, np.round(counts.values)):
        raise SystemExit("counts.mtx is not integer-valued -- refusing (need raw counts)")

    meta = pd.read_csv(f"{in_dir}/metadata.tsv", sep="\t").set_index("barcode")
    meta = meta.reindex(counts.index)                      # align to cells, exact order

    ds = mosaicmpi.Dataset.from_df(data=counts, is_normalized=False, obs=meta)
    ds.write_h5ad(out_h5ad)

    print(f"wrote {out_h5ad}: {counts.shape[0]} cells x {counts.shape[1]} genes "
          f"(dense .X, raw counts, {meta.shape[1]} metadata cols)")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
