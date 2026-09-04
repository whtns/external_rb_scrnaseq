#!/bin/bash
#SBATCH --job-name=sync_summaries_gdrive
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/sync_summaries_gdrive_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/sync_summaries_gdrive_%j.err

# Push the 2026-09-03 sample_summaries rebuild (#41 allele page, #42 round-correct
# bulk panel, #43 depth row) to Drive. The copies there are all dated 2026-08-27,
# i.e. pre-fix.
#
# Deliberately NOT `sync_hypoxia_gdrive.sbatch` / tar_make(hypoxia_rebuilt_gdrive).
# That script's comment ("nothing upstream is touched") no longer holds: after the
# _targets.R edits for #41/#42/#44, tar_outdated() reports 72 targets upstream of
# hypoxia_rebuilt_gdrive as outdated, so a tar_make() would rebuild most of the
# pipeline to perform a file copy.
#
# This reads the already-built sample_summaries paths out of the store with
# tar_read() -- no tar_make(), so nothing is rebuilt, nothing is written to the
# store, and the one-job-per-store rule is not engaged.
#
# Staging mirrors hypoxia_rebuilt_gdrive exactly (R/pipeline_targets_figures.R:919-932):
# symlink each PDF under <plot type>/<paper_retained|other>/, then one
# `rclone copy --copy-links`. Destination paths therefore match the 33 files
# already on Drive and update them in place. rclone copy re-uploads only what
# changed and never deletes -- `sync` is not used.

module load r/4.4.1
module load rclone/1.74.3
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

command -v rclone >/dev/null || { echo "FATAL: rclone not on PATH"; exit 1; }
rclone listremotes | grep -q '^gdrive:' || { echo "FATAL: no gdrive: remote"; exit 1; }

export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

Rscript -e '
suppressPackageStartupMessages(library(targets))
tar_config_set(store = "_targets_r431")
source("R/pipeline_constants.R")   # hypoxia_gdrive_destination

paths    <- unique(unlist(tar_read(sample_summaries)))
paths    <- paths[nchar(paths) > 0 & !is.na(paths) & file.exists(paths)]
retained <- tar_read(paper_retained_samples)
stopifnot(length(paths) > 0)

cat("staging", length(paths), "sample summaries\n")
cat("newest mtime:", format(max(file.mtime(paths))), "\n")

stage <- file.path(tempdir(), "summary_sync_stage")
unlink(stage, recursive = TRUE); dir.create(stage)
on.exit(unlink(stage, recursive = TRUE))

for (p in paths) {
  sid <- stringr::str_extract(basename(p), "SR[RX][0-9]+")
  grp <- if (!is.na(sid) && sid %in% retained) "paper_retained" else "other"
  d   <- file.path(stage, "sample_summary", grp)
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

echo
echo "=== verification: Drive-side mtimes after sync ==="
rclone lsl "gdrive:rb_scrnaseq/hypoxia_rebuilt" | grep 'sample_summary/' | sort -k4

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"sync_summaries_gdrive ($SLURM_JOB_ID) finished: exit $RC\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi

exit $RC
