#!/bin/bash
#SBATCH --job-name=build_hypoxia_mito_box
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/build_hypoxia_mito_box_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/build_hypoxia_mito_box_%j.err

# Build hypoxia_mito_score_boxplots (the hypoxia-vs-mitochondrial decomposition
# of the composite score the split thresholds on) and re-sync the hypoxia
# deliverables to Drive so the new figure lands beside the composite boxplot.
#
# hypoxia_rebuilt_gdrive re-runs regardless: it gained a mito_score_boxplot
# entry, which changes its command. rclone skips byte-identical files, so the
# other ~705 already-synced files cost a listing, not a re-upload.
#
# One tar_make() at a time against _targets_r431 (AGENTS.md) -- check first:
#   squeue -u $USER -o "%.10i %.22j %.8T %.10M %R"

module load r/4.4.1
module load rclone/1.74.3   # needed by hypoxia_rebuilt_gdrive
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

command -v rclone >/dev/null || { echo "FATAL: rclone not on PATH"; exit 1; }
rclone listremotes | grep -q '^gdrive:' || { echo "FATAL: no gdrive: remote"; exit 1; }

# install_local, not load_all: hypoxia_mito_score_boxplots runs on a crew
# worker, which is a fresh R session that only sees the installed package.
Rscript -e '
suppressMessages(devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE))
cat("numbatHelpers installed\n")
' || { echo "FATAL: install failed"; exit 1; }

# targets does not track package function bodies, so a change inside
# plot_hypoxia_mito_score_boxplots() does not invalidate the target on its own.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
tar_invalidate(any_of("hypoxia_mito_score_boxplots"))
cat("invalidated hypoxia_mito_score_boxplots\n")
' || { echo "FATAL: invalidation failed"; exit 1; }

Rscript run_targets_bg.R hypoxia_mito_score_boxplots hypoxia_rebuilt_gdrive

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"build_hypoxia_mito_box ($SLURM_JOB_ID) finished\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi
