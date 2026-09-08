#!/bin/bash
#SBATCH --job-name=summaries_haplo
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_summaries_haplotype_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_summaries_haplotype_%j.err

# Rebuild sample_summaries so #41's haplotype track is actually visible.
#
# collate_sample_summary()'s annotate_panel() read page 1 only, so the
# allele-only haplotype page that numbat_bulk_panel.R has been emitting since
# #41 never reached the summaries. annotate_panel() now takes a page argument
# and a new row5b renders page 2. Verified: all 33 bulk-clone PDFs have exactly
# 2 pages.
#
# install_local, NOT load_all: sample_summaries runs on crew workers, which are
# separate processes and would otherwise pick up the installed (old) package.
#
# tar_invalidate is mandatory here -- targets does not track package function
# bodies, so the change above is invisible to the pipeline's own cue and
# sample_summaries would be reported up to date (AGENTS.md).
#
# ONE tar_make() against _targets_r431. Checked free before submitting.

set -uo pipefail

module load gcc/13.3.0
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

echo "=== installing numbatHelpers ==="
Rscript -e 'devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force = TRUE, upgrade = "never")' || exit 1

echo "=== confirming the installed copy has the page argument ==="
Rscript -e '
  src <- deparse(body(numbatHelpers::collate_sample_summary))
  stopifnot(any(grepl("page = 2L", src)), any(grepl("row5b", src)))
  cat("installed package carries the haplotype row\n")' || exit 1

echo "=== invalidating sample_summaries ==="
Rscript -e '
  suppressPackageStartupMessages(library(targets))
  tar_config_set(store = "_targets_r431")
  before <- tar_meta(names = "sample_summaries", fields = "time")$time
  tar_invalidate(sample_summaries)
  cat("invalidated; previous build time was", format(before), "\n")' || exit 1

echo "=== building ==="
Rscript run_targets_bg.R sample_summaries
RC=$?

echo "=== verifying the build time actually moved ==="
Rscript -e '
  suppressPackageStartupMessages(library(targets))
  tar_config_set(store = "_targets_r431")
  print(tar_meta(names = "sample_summaries", fields = c("time", "error")))'

echo "### exit: $RC"
exit $RC
