#!/usr/bin/env python3
"""
Download FASTQ files for PRJNA728722 (SRR14483745–SRR14483975).

Strategy:
  1. Batch-query ffq --ftp for EBI/ENA FTP URLs (fast, no SRA quota).
  2. Download in parallel with wget (-c for resume support).
  3. Fall back to fasterq-dump (module load sratoolkit) if ffq finds no URL.

Usage:
  python3 download_prjna728722.py [options]

  --outdir / -o   Output directory for FASTQ files  (default: fastq)
  --threads / -t  Parallel wget downloads            (default: 4)
  --batch-size    SRR accessions per ffq call        (default: 50)
  --dry-run       Print URLs; do not download
  --fasterq-only  Skip ffq; use fasterq-dump for all accessions
"""

import argparse
import json
import logging
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

logging.basicConfig(
    format="%(asctime)s %(levelname)s %(message)s",
    datefmt="%H:%M:%S",
    level=logging.INFO,
)
log = logging.getLogger(__name__)

ACCESSIONS = (
    [f"SRR{i}" for i in range(14483745, 14483777)] +   # SRR14483745–14483776
    [f"SRR{i}" for i in range(14483781, 14483785)] +   # SRR14483781–14483784
    [f"SRR{i}" for i in range(14483937, 14483940)] +   # SRR14483937–14483939
    [f"SRR{i}" for i in range(14483958, 14483976)]     # SRR14483958–14483975
)  # 57 samples total


# ---------------------------------------------------------------------------
# ffq helpers
# ---------------------------------------------------------------------------

def run_ffq_batch(accessions: list[str]) -> list[dict]:
    """Call `ffq --ftp` for a batch of accessions; return list of file dicts."""
    cmd = ["ffq", "--ftp"] + accessions
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log.warning("ffq error: %s", result.stderr.strip()[:200])
        return []
    try:
        data = json.loads(result.stdout)
    except json.JSONDecodeError as exc:
        log.warning("ffq JSON parse error: %s", exc)
        return []
    return data if isinstance(data, list) else [data]


def collect_ftp_urls(accessions: list[str], batch_size: int) -> dict[str, list[str]]:
    """Return {accession: [ftp_url, ...]} for all accessions."""
    url_map: dict[str, list[str]] = {a: [] for a in accessions}
    n_batches = (len(accessions) + batch_size - 1) // batch_size

    for idx in range(0, len(accessions), batch_size):
        batch = accessions[idx : idx + batch_size]
        batch_num = idx // batch_size + 1
        log.info("ffq batch %d/%d  (%s … %s)", batch_num, n_batches, batch[0], batch[-1])
        entries = run_ffq_batch(batch)
        for entry in entries:
            url = entry.get("url", "")
            acc = entry.get("accession", "")
            if acc in url_map and url.startswith("ftp://"):
                url_map[acc].append(url)

    found = sum(1 for v in url_map.values() if v)
    log.info("ffq: FTP URLs found for %d / %d accessions", found, len(accessions))
    return url_map


# ---------------------------------------------------------------------------
# Download helpers
# ---------------------------------------------------------------------------

def wget_download(url: str, outdir: Path, retries: int = 3) -> bool:
    filename = url.rsplit("/", 1)[-1]
    outpath = outdir / filename
    if outpath.exists() and outpath.stat().st_size > 0:
        log.info("Skip (exists): %s", filename)
        return True
    for attempt in range(1, retries + 1):
        r = subprocess.run(
            ["wget", "-q", "-c", "--show-progress", "-P", str(outdir), url],
        )
        if r.returncode == 0:
            log.info("Downloaded: %s", filename)
            return True
        log.warning("wget attempt %d/%d failed: %s", attempt, retries, url)
    log.error("FAILED: %s", url)
    return False


def fasterq_download(accession: str, outdir: Path, threads: int = 6) -> bool:
    """Download via fasterq-dump (requires sratoolkit module to be loaded)."""
    existing = list(outdir.glob(f"{accession}*.fastq.gz"))
    if existing:
        log.info("Skip (exists): %s", accession)
        return True
    tmpdir = outdir / "tmp"
    tmpdir.mkdir(exist_ok=True)
    r = subprocess.run(
        [
            "fasterq-dump", accession,
            "--outdir", str(outdir),
            "--temp", str(tmpdir),
            "--threads", str(threads),
            "--split-files",
            "--gzip",
            "--progress",
        ],
    )
    if r.returncode == 0:
        log.info("fasterq-dump done: %s", accession)
        return True
    log.error("fasterq-dump FAILED: %s", accession)
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", "-o", default="fastq",
                        help="Output directory (default: fastq)")
    parser.add_argument("--threads", "-t", type=int, default=4,
                        help="Parallel wget downloads (default: 4)")
    parser.add_argument("--batch-size", "-b", type=int, default=50,
                        help="Accessions per ffq call (default: 50)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print URLs without downloading")
    parser.add_argument("--fasterq-only", action="store_true",
                        help="Skip ffq; use fasterq-dump for all accessions")
    parser.add_argument("--fasterq-threads", type=int, default=6,
                        help="Threads per fasterq-dump call (default: 6)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    log.info("Project: PRJNA728722  |  %d accessions  |  outdir: %s",
             len(ACCESSIONS), outdir)

    # ---- fasterq-dump-only mode ----------------------------------------
    if args.fasterq_only:
        log.info("Mode: fasterq-dump only")
        if args.dry_run:
            for acc in ACCESSIONS:
                print(acc)
            return
        failed = []
        for acc in ACCESSIONS:
            if not fasterq_download(acc, outdir, threads=args.fasterq_threads):
                failed.append(acc)
        _report_failures(failed)
        return

    # ---- ffq → wget mode -----------------------------------------------
    url_map = collect_ftp_urls(ACCESSIONS, batch_size=args.batch_size)

    # Accessions with no FTP URL fall back to fasterq-dump
    missing = [acc for acc, urls in url_map.items() if not urls]
    all_urls = [(acc, url) for acc, urls in url_map.items() for url in urls]

    # Save URL list for reference
    url_file = outdir / "ftp_urls.txt"
    with url_file.open("w") as fh:
        for _, url in all_urls:
            fh.write(url + "\n")
    log.info("Saved %d FTP URLs to %s", len(all_urls), url_file)

    if args.dry_run:
        for _, url in all_urls:
            print(url)
        if missing:
            log.info("No FTP URL (would use fasterq-dump): %s", missing)
        return

    # Parallel wget
    failed_urls: list[str] = []
    if all_urls:
        log.info("Downloading %d files with %d threads…", len(all_urls), args.threads)
        with ThreadPoolExecutor(max_workers=args.threads) as pool:
            futures = {
                pool.submit(wget_download, url, outdir): url
                for _, url in all_urls
            }
            for fut in as_completed(futures):
                if not fut.result():
                    failed_urls.append(futures[fut])

    # fasterq-dump fallback for accessions ffq couldn't resolve
    failed_fasterq: list[str] = []
    if missing:
        log.info("Falling back to fasterq-dump for %d accessions: %s…",
                 len(missing), missing[:3])
        for acc in missing:
            if not fasterq_download(acc, outdir, threads=args.fasterq_threads):
                failed_fasterq.append(acc)

    _report_failures(failed_urls + failed_fasterq)


def _report_failures(failed: list[str]) -> None:
    if failed:
        log.error("%d downloads failed:", len(failed))
        for item in failed:
            log.error("  %s", item)
        sys.exit(1)
    log.info("All downloads complete.")


if __name__ == "__main__":
    main()
