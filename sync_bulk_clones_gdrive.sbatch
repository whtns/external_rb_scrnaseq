#!/bin/bash
#SBATCH --job-name=sync_bulk_clones
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/sync_bulk_clones_gdrive_%j.log

# Push the 33 round-correct numbat bulk-clone panels (#42) to Drive.
# Page 1 is expression+allele, page 2 is the allele-only haplotype track (#41).
#
# New folder numbat_bulk_clones/{paper_retained,other}/ under the existing
# hypoxia_rebuilt destination, split the same way sample_summary is.
#
# ~321 MB across 33 files. `rclone copy` (never sync) so nothing already on
# Drive can be deleted by this job. Reads tar_read(paper_retained_samples)
# only -- no tar_make, so it takes no lock on _targets_r431.

set -uo pipefail

module load gcc/13.3.0
module load r/4.4.1
module load rclone/1.74.3

export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

command -v rclone >/dev/null || { echo "FATAL: rclone not on PATH"; exit 1; }
rclone listremotes | grep -q '^gdrive:' || { echo "FATAL: no gdrive: remote"; exit 1; }

Rscript -e '
suppressPackageStartupMessages(library(targets))
tar_config_set(store = "_targets_r431")
source("R/pipeline_constants.R")   # hypoxia_gdrive_destination

paths <- Sys.glob("results/numbat_bulk_clones/*.pdf")
paths <- paths[file.exists(paths)]
stopifnot(length(paths) > 0)
retained <- tar_read(paper_retained_samples)

cat("staging", length(paths), "bulk-clone panels\n")
cat("newest mtime:", format(max(file.mtime(paths))), "\n")

stage <- file.path(tempdir(), "bulk_clone_sync_stage")
unlink(stage, recursive = TRUE); dir.create(stage)
on.exit(unlink(stage, recursive = TRUE))

for (p in paths) {
  sid <- stringr::str_extract(basename(p), "SR[RX][0-9]+")
  grp <- if (!is.na(sid) && sid %in% retained) "paper_retained" else "other"
  d   <- file.path(stage, "numbat_bulk_clones", grp)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  file.symlink(normalizePath(p), file.path(d, basename(p)))
}

ret <- system2("rclone",
  c("copy", "--copy-links", "--stats-one-line", "--stats", "30s",
    shQuote(stage), hypoxia_gdrive_destination, "--progress"),
  stdout = TRUE, stderr = TRUE)
cat(ret, sep = "\n")
if (!is.null(attr(ret, "status")) && attr(ret, "status") != 0)
  stop("rclone failed")

cat("\nsynced", length(paths), "files to", hypoxia_gdrive_destination, "\n")
'
RC=$?
echo "### rclone stage exit: $RC"
exit $RC
