#!/bin/bash
#SBATCH --job-name=build_41_44
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/build_41_44_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/build_41_44_%j.err

# Production build for github #41, #42, #43, #44.
#
# TWO SEQUENTIAL tar_make() calls, never concurrent (AGENTS.md). Check nothing
# else is on the store first:
#   squeue -u $USER -o "%.10i %.22j %.8T %.10M %R"
#
# Why two stages rather than one target list: diploid_seu takes the default
# error = "stop", so when it failed in job 11742715 it halted the pipeline
# before sample_summaries finished -- #42 and #43 were held hostage by #44.
# Stage 1 now runs to completion on its own; stage 2 can fail without taking
# stage 1 with it.
#
# WARNING: diploid_seu OVERWRITES output/seurat/diploid_subsets/diploid_seu.rds,
# the hand-made 1.38 GB object from 2026-04-26. Backed up first to
# diploid_seu_handmade_20260426.rds.bak (md5 8c09dbb1685a1af08b6546fc8e415de8);
# the guard below refuses to run without it. Every downstream diploid figure
# will move.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

BAK=output/seurat/diploid_subsets/diploid_seu_handmade_20260426.rds.bak
[ -s "$BAK" ] || { echo "FATAL: backup $BAK missing; refusing to overwrite the hand-made diploid object"; exit 1; }
echo "backup present: $(ls -la $BAK)"

# install_local, not load_all: these targets run on crew workers, which are
# fresh R sessions that see only the installed package.
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers",
                                         force = TRUE, upgrade = "never", quiet = TRUE))
cat("numbatHelpers installed\n")' || { echo "FATAL: install failed"; exit 1; }

# targets does not track package function bodies. sample_summaries' COMMAND is
# unchanged but collate_sample_summary() is not, so it must be invalidated by
# hand. diploid_seu is invalidated because plot_celltype_predictions() changed.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
tar_invalidate(any_of(c("sample_summaries", "diploid_seu", "diploid_seu_composition")))
cat("invalidated the changed targets\n")' || { echo "FATAL: invalidation failed"; exit 1; }

echo
echo "############ stage 1: #42 / #41 / #43 ############"
Rscript run_targets_bg.R numbat_bulk_clone_pdfs sample_summaries
S1=$?
echo "### stage 1 exit: $S1"

echo
echo "############ stage 2: #44 (diploid rebuild) ############"
# Runs regardless of stage 1's exit: they are independent deliverables, and a
# stage-1 branch failure should not block the diploid work either.
Rscript run_targets_bg.R diploid_seu_composition
S2=$?
echo "### stage 2 exit: $S2"

echo
echo "=== stage 1: $S1   stage 2: $S2 ==="

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"build_41_44 ($SLURM_JOB_ID) done - stage1 exit $S1, stage2 exit $S2\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
