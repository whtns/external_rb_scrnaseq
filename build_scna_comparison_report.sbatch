#!/bin/bash
#SBATCH --job-name=scna_comparison
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/scna_comparison_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/scna_comparison_%j.err

# Two-clone collages + diffex for the tier-A RB SCNA contrasts, after extending
# scna_collage_samples and large_clone_comparisons.
# See docs/rb_scna_clone_comparison_report.md.
#
# ONE tar_make() against _targets_r431 -- store verified free before submitting.
#
# 64G/4cpu/12h deliberately: epyc-64 had 178 pending vs 3 running, and a
# 96G/8cpu/24h ask sat unscheduled behind them. The heavy work runs in
# crew-worker jobs with their own allocations, so the driver does not need more.

set -uo pipefail
module purge
module load gcc/13.3.0   # r/4.4.1 is nested under gcc in the Lmod hierarchy
module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers",
                                         force = TRUE, upgrade = "never", quiet = TRUE))
cat("numbatHelpers installed\n")' || { echo "FATAL: install failed"; exit 1; }

# The config targets are tar_file/plain targets keyed on file content, so the
# YAML edit invalidates them on its own; scna_collage_samples is a literal in
# the pipeline source and targets DOES track that. Invalidate the consumers
# explicitly anyway -- the collage branches carry cue(command = FALSE) in
# places, which suppresses re-running on a changed command.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
tar_invalidate(matches("^two_clone_res_collages_filtered"))
tar_invalidate(matches("^all_clone_res_collages_filtered"))
tar_invalidate(any_of(c("scna_collage_samples", "large_clone_comparisons",
                        "large_clone_comparisons_file", "clone_comparison_table")))
tar_invalidate(matches("^all_diffex_clones"))
cat("invalidated collage + diffex targets\n")' || { echo "FATAL: invalidation failed"; exit 1; }

for TGT in two_clone_res_collages_filtered_merged_1q \
           two_clone_res_collages_filtered_merged_2p \
           two_clone_res_collages_filtered_merged_6p \
           two_clone_res_collages_filtered_merged_16q \
           all_diffex_clones_for_each_cluster; do
  echo "=== building $TGT ==="
  Rscript run_targets_bg.R "$TGT"
  echo "--- $TGT returned $? (not authoritative; checked below) ---"
done

# run_targets_bg.R exits 0 even when tar_make() errors, so ask the store.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
n <- c(paste0("two_clone_res_collages_filtered_merged_", c("1q","2p","6p","16q")),
       "all_diffex_clones_for_each_cluster")
m <- tar_meta(any_of(n), fields = c("time","error"))
print(as.data.frame(m))
bad <- m$name[!is.na(m$error)]
if (length(bad)) { cat("FAILED:", paste(bad, collapse=", "), "\n"); quit(status=1) }
cat("all requested targets built clean\n")' || echo "### SOME TARGETS FAILED"

echo "=== collage PDFs on disk ==="
ls -la results/two_clone_scna_collage_filtered_merged_*.pdf 2>/dev/null

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"scna_comparison ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
