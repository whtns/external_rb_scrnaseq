#!/usr/bin/env python
"""Issue #40, Stage F. Export per-cell community usage from a mosaicMPI integration.

Loads the integration network and, for each requested dataset (sample), writes the
cells x communities usage matrix (barcode-indexed) plus the program->community map. This
is the load-bearing join for the "cross-tumor recurrent community" claim: in R we join
each community's per-cell usage onto the Seurat metadata by barcode and correlate against
hypoxia_score / S.Score / G2M.Score.

normalize=False -> raw community activity (what we correlate against the scores).

Usage: python src/mosaic/export_community_usage.py <network_pkl> <out_dir> <SRX> [<SRX> ...]
"""
import sys
import os
import mosaicmpi


def main(pkl: str, out_dir: str, datasets: list[str]) -> None:
    os.makedirs(out_dir, exist_ok=True)
    net = mosaicmpi.Network.from_pkl(pkl)

    # program -> community membership (program id = "dataset|k|program").
    with open(os.path.join(out_dir, "program_communities.tsv"), "w") as fh:
        fh.write("program\tcommunity\n")
        for prog, comm in net.program_communities.items():
            fh.write(f"{prog}\t{comm}\n")
    print(f"communities: {len(net.communities)}; programs placed: {len(net.program_communities)}")

    for srx in datasets:
        try:
            usage = net.get_community_usage(subset_datasets=srx, normalize=False)
        except Exception as e:  # noqa: BLE001 - dataset may have no surviving programs
            print(f"!! {srx}: get_community_usage failed ({e}); skipping")
            continue
        # usage index may be a (dataset, barcode) MultiIndex or plain barcodes.
        if usage.index.nlevels > 1:
            usage = usage.droplevel(list(range(usage.index.nlevels - 1)))
        out = os.path.join(out_dir, f"community_usage_{srx}.csv")
        usage.to_csv(out, index_label="barcode")
        print(f"wrote {out}: {usage.shape[0]} cells x {usage.shape[1]} communities")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
