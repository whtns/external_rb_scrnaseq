#!/usr/bin/env python
"""Summarize BigSur runs across all samples into one table.

Reads the per-sample mcFano CSVs plus the array-job logs and reports, per sample:
cells, HVGs, fitted CV, whether the CV fit fell back to the 0.5 default, the
number of significant gene pairs, the top mcFano genes, and flags for the
contamination programs (immunoglobulin, haemoglobin, crystallin, dissociation
stress) that dominate feature selection in some samples.

Usage:  python scripts/bigsur_cohort_summary.py [--logs 'logs/bigsur_all_10217585_*.log']
"""

import argparse
import glob
import os
import re

import pandas as pd

OUTROOT = "/project2/cobrinik_1090/external_rb_scrnaseq_proj/results/bigsur"

# Gene programs that, when they top the mcFano ranking, mean feature selection is
# being driven by contamination or a dissociation artifact rather than tumour biology.
PROGRAMS = {
    "immunoglobulin": re.compile(r"^(IGKC|IGHG\d|IGHA\d|IGHM|IGLC\d|JCHAIN)$"),
    "haemoglobin": re.compile(r"^HB[ABGDQZE]\d?$"),
    "crystallin": re.compile(r"^CRY[ABG][ABS]?\d?$"),
    "stress_IEG": re.compile(r"^(FOS|FOSB|JUN|JUNB|JUND|EGR1|DUSP1|ZFP36|GADD45B|HSPA1[AB]|SOCS3|IER2|ATF3)$"),
    "mito": re.compile(r"^MT-"),
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--logs", default="logs/bigsur_all_*.log")
    p.add_argument("--outroot", default=OUTROOT)
    p.add_argument("--top-n", type=int, default=50,
                   help="how many top mcFano genes to scan for the programs")
    args = p.parse_args()
    outroot = args.outroot

    # Pull cells / CV-fallback / significant-pair counts out of the run logs.
    log_info = {}
    for path in glob.glob(args.logs):
        txt = open(path).read()
        m = re.search(r"^(\S+): (\d+) cells x", txt, re.M)
        if not m:
            continue
        sample, cells = m.group(1), int(m.group(2))
        pairs = re.search(r"(\d+) significant gene pairs", txt)
        log_info[sample] = {
            "cells": cells,
            "cv_fit_failed": "CV cannot be fit in biological range" in txt,
            "sig_pairs": int(pairs.group(1)) if pairs else None,
        }

    rows = []
    for csv in sorted(glob.glob(os.path.join(outroot, "*_mcfano.csv"))):
        sample = os.path.basename(csv).replace("_mcfano.csv", "")
        var = pd.read_csv(csv, index_col=0)
        top = var.head(args.top_n).index

        row = {"sample": sample}
        row.update(log_info.get(sample, {}))
        row["genes_tested"] = len(var)
        row["hvgs"] = int(var.highly_variable.sum())
        row["fdr_sig_genes"] = int((var.FDR_adj_pvalue < 0.05).sum())
        row["max_mcfano"] = round(float(var.mc_Fano.max()), 1)
        row["top5"] = ",".join(var.head(5).index)
        for name, pat in PROGRAMS.items():
            row[f"n_{name}_in_top{args.top_n}"] = sum(bool(pat.match(g)) for g in top)
        rows.append(row)

    df = pd.DataFrame(rows)
    cols = ["sample", "cells", "genes_tested", "hvgs", "fdr_sig_genes", "cv_fit_failed",
            "sig_pairs", "max_mcfano"] + \
           [f"n_{k}_in_top{args.top_n}" for k in PROGRAMS] + ["top5"]
    df = df[[c for c in cols if c in df.columns]]

    out = os.path.join(outroot, "cohort_summary.csv")
    df.to_csv(out, index=False)

    pd.set_option("display.width", 250, "display.max_columns", 50)
    print(df.drop(columns=["top5"]).to_string(index=False))
    print()
    print(f"samples:                     {len(df)}")
    print(f"CV fit fell back to 0.5:     {int(df.cv_fit_failed.sum())}")
    print(f"median HVGs:                 {df.hvgs.median():.0f}")
    print(f"median significant pairs:    {df.sig_pairs.median():.0f}")
    for name in PROGRAMS:
        col = f"n_{name}_in_top{args.top_n}"
        n = int((df[col] > 0).sum())
        print(f"samples with {name:15s} in top {args.top_n}: {n}/{len(df)}")
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
