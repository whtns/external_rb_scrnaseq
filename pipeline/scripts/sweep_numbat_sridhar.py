#!/usr/bin/env python3
"""
Run numbat_sridhar over a tabular parameter grid and summarize outcomes.

Example:
  python scripts/sweep_numbat_sridhar.py \
    --grid config/numbat_sridhar_sweep_grid.tsv \
    --snakefile Snakefile \
    --configfile config.yaml \
    --cores 6
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import os
import re
import shlex
import subprocess
from pathlib import Path
from typing import Dict, List, Tuple

GRID_PARAMS = [
    "numbat_t",
    "gamma",
    "max_entropy",
    "min_LLR",
    "max_iter",
    "read_prop",
    "tau",
    "cell_ceiling",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Sweep numbat_sridhar parameters.")
    parser.add_argument("--grid", required=True, help="TSV with run_id, sample, and parameter columns.")
    parser.add_argument("--snakefile", default="Snakefile", help="Path to Snakefile.")
    parser.add_argument("--configfile", default="config.yaml", help="Path to Snakemake config file.")
    parser.add_argument("--cores", type=int, default=4, help="Snakemake cores per run.")
    parser.add_argument("--sweep-dir", default="numbat_sridhar_sweeps", help="Output subdirectory under pipeline output root.")
    parser.add_argument("--dry-run", action="store_true", help="Pass -n to snakemake.")
    parser.add_argument("--rerun-incomplete", action="store_true", help="Pass --rerun-incomplete to snakemake.")
    parser.add_argument("--keep-going", action="store_true", help="Pass --keep-going to snakemake.")
    return parser.parse_args()


def parse_config_value(config_path: Path, key: str, default: str) -> str:
    pattern = re.compile(rf"^\s*{re.escape(key)}\s*:\s*(.+?)\s*$")
    for raw in config_path.read_text().splitlines():
        line = raw.split("#", 1)[0].rstrip()
        if not line:
            continue
        m = pattern.match(line)
        if m:
            return m.group(1).strip().strip('"').strip("'")
    return default


def read_grid(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("Grid file has no header.")
        required = {"run_id", "sample", *GRID_PARAMS}
        missing = required - set(reader.fieldnames)
        if missing:
            raise ValueError(f"Grid missing required columns: {sorted(missing)}")
        rows = []
        for row in reader:
            if not any(v.strip() for v in row.values() if v is not None):
                continue
            rows.append({k: (v.strip() if v is not None else "") for k, v in row.items()})
    if not rows:
        raise ValueError("Grid has no data rows.")
    return rows


def parse_benchmark_seconds(path: Path) -> float | None:
    if not path.exists():
        return None
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    if len(lines) < 2:
        return None
    cols = re.split(r"\s+", lines[1].strip())
    if not cols:
        return None
    try:
        return float(cols[0])
    except ValueError:
        return None


def parse_numbat_log_metrics(log_path: Path) -> Dict[str, object]:
    metrics: Dict[str, object] = {
        "input_cells": None,
        "max_observed_iteration": None,
        "all_cells_succeeded": False,
        "no_cnv_warning": False,
        "elapsed_sec_from_log": None,
    }
    if not log_path.exists():
        return metrics

    txt = log_path.read_text(errors="replace")

    m_cells = re.search(r"Input metrics:\s*\n\s*(\d+)\s+cells", txt)
    if m_cells:
        metrics["input_cells"] = int(m_cells.group(1))

    iters = [int(x) for x in re.findall(r"Iteration\s+(\d+)", txt)]
    if iters:
        metrics["max_observed_iteration"] = max(iters)

    metrics["all_cells_succeeded"] = "All cells succeeded" in txt
    metrics["no_cnv_warning"] = "No CNV remains after filtering by entropy" in txt

    timestamps = re.findall(r"INFO \[([0-9\-: ]+)\]", txt)
    if len(timestamps) >= 2:
        try:
            t0 = dt.datetime.strptime(timestamps[0], "%Y-%m-%d %H:%M:%S")
            t1 = dt.datetime.strptime(timestamps[-1], "%Y-%m-%d %H:%M:%S")
            metrics["elapsed_sec_from_log"] = int((t1 - t0).total_seconds())
        except ValueError:
            pass

    return metrics


def detect_max_iter_from_files(run_sample_dir: Path) -> int | None:
    if not run_sample_dir.exists():
        return None
    max_iter = None
    for name in os.listdir(run_sample_dir):
        m = re.search(r"_(\d+)\.(tsv|tsv\.gz)$", name)
        if m:
            cur = int(m.group(1))
            max_iter = cur if max_iter is None else max(max_iter, cur)
    return max_iter


def read_rds_class(rds_path: Path) -> Tuple[str | None, bool | None]:
    if not rds_path.exists():
        return None, None
    code = (
        "x <- readRDS(commandArgs(TRUE)[1]); "
        "cat(class(x)[1], '\\n'); "
        "is_recovered <- isTRUE(is.list(x) && !is.null(x[['recovered']]) && x[['recovered']]); "
        "cat(ifelse(is_recovered, 'TRUE', 'FALSE'), '\\n')"
    )
    cmd = ["Rscript", "-e", code, str(rds_path)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        return None, None
    out = [ln.strip() for ln in proc.stdout.splitlines() if ln.strip()]
    if len(out) < 2:
        return None, None
    return out[0], out[1] == "TRUE"


def main() -> int:
    args = parse_args()

    snakefile = Path(args.snakefile).resolve()
    configfile = Path(args.configfile).resolve()
    grid_path = Path(args.grid).resolve()

    snake_dir = snakefile.parent
    output_rel = parse_config_value(configfile, "output", "../output/")
    output_root = (snake_dir / output_rel).resolve()

    rows = read_grid(grid_path)

    summary_rows: List[Dict[str, object]] = []

    for i, row in enumerate(rows, start=1):
        run_id = row["run_id"]
        sample = row["sample"]
        run_subdir = f"{args.sweep_dir}/{sample}/{run_id}"
        run_sample_dir = output_root / run_subdir / sample
        target_rds = output_root / run_subdir / f"{sample}_numbat.rds"

        overrides = [
            f"sample={sample}",
            f"numbat_sridhar_dirname={run_subdir}",
        ]
        overrides.extend(f"{k}={row[k]}" for k in GRID_PARAMS)

        cmd = [
            "snakemake",
            "-s",
            str(snakefile),
            "--configfile",
            str(configfile),
            "--cores",
            str(args.cores),
            str(target_rds),
            "--config",
            *overrides,
        ]

        if args.dry_run:
            cmd.insert(1, "-n")
        if args.rerun_incomplete:
            cmd.append("--rerun-incomplete")
        if args.keep_going:
            cmd.append("--keep-going")

        run_sample_dir.mkdir(parents=True, exist_ok=True)
        print(f"[{i}/{len(rows)}] Running {run_id} for {sample}")
        print("  " + " ".join(shlex.quote(x) for x in cmd))

        proc = subprocess.run(cmd, cwd=snake_dir, capture_output=True, text=True)

        (run_sample_dir / "sweep_stdout.log").write_text(proc.stdout)
        (run_sample_dir / "sweep_stderr.log").write_text(proc.stderr)

        numbat_log = run_sample_dir / "log.txt"
        bench_main = run_sample_dir / "benchmark_numbat.txt"

        log_metrics = parse_numbat_log_metrics(numbat_log)
        benchmark_sec = parse_benchmark_seconds(bench_main)
        max_iter_files = detect_max_iter_from_files(run_sample_dir)
        rds_class, recovered = read_rds_class(target_rds)

        summary = {
            "run_id": run_id,
            "sample": sample,
            "status": "success" if proc.returncode == 0 else "failed",
            "return_code": proc.returncode,
            "run_subdir": run_subdir,
            "target_rds": str(target_rds),
            "rds_exists": target_rds.exists(),
            "rds_class": rds_class,
            "rds_recovered": recovered,
            "benchmark_sec": benchmark_sec,
            "max_observed_iteration": log_metrics["max_observed_iteration"],
            "max_iteration_from_files": max_iter_files,
            "input_cells": log_metrics["input_cells"],
            "all_cells_succeeded": log_metrics["all_cells_succeeded"],
            "no_cnv_warning": log_metrics["no_cnv_warning"],
            "elapsed_sec_from_log": log_metrics["elapsed_sec_from_log"],
        }
        summary.update({k: row[k] for k in GRID_PARAMS})
        summary_rows.append(summary)

    # Sort: successful first, then non-recovered first, then no-cnv-warning last.
    summary_rows.sort(
        key=lambda r: (
            0 if r["status"] == "success" else 1,
            1 if r.get("rds_recovered") else 0,
            1 if r.get("no_cnv_warning") else 0,
            -(r.get("max_observed_iteration") or 0),
        )
    )

    sweep_root = output_root / args.sweep_dir
    sweep_root.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "run_id",
        "sample",
        "status",
        "return_code",
        *GRID_PARAMS,
        "run_subdir",
        "target_rds",
        "rds_exists",
        "rds_class",
        "rds_recovered",
        "benchmark_sec",
        "elapsed_sec_from_log",
        "input_cells",
        "max_observed_iteration",
        "max_iteration_from_files",
        "all_cells_succeeded",
        "no_cnv_warning",
    ]

    out_tsv = sweep_root / "summary.tsv"
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    out_csv = sweep_root / "summary.csv"
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summary_rows)

    print(f"Wrote summary: {out_tsv}")
    print(f"Wrote summary: {out_csv}")

    failures = sum(1 for r in summary_rows if r["status"] != "success")
    if failures:
        print(f"Completed with {failures} failed runs.")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
