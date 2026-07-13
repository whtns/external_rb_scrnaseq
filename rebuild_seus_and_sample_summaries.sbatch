#!/bin/bash
#SBATCH --job-name=rebuild-seus-summaries
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_seus_and_sample_summaries_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_seus_and_sample_summaries_%j.err

# Rebuild unfiltered_seus, filtered_seus, seus_low_hypoxia, and sample_summaries.
# Steps:
#   1. Install updated numbatHelpers so crew workers pick up the latest code.
#   2. Invalidate all branches of the three cue-frozen Seurat targets.
#   3. tar_make(sample_summaries) — rebuilds the full chain through crew workers.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript - <<'EOF'
devtools::install_local(
  "/project2/cobrinik_1090/rpkgs/numbat_helpers",
  force = TRUE,
  upgrade = "never",
  quiet = TRUE
)
cat("numbatHelpers installed.\n"); flush(stdout())

targets::tar_config_set(store = "_targets_r431")

targets::tar_invalidate(
  names = tidyselect::starts_with("unfiltered_seus") |
          tidyselect::starts_with("filtered_seus") |
          tidyselect::starts_with("seus_low_hypoxia")
)
cat("Invalidated unfiltered_seus, filtered_seus, seus_low_hypoxia branches.\n"); flush(stdout())

targets::tar_make(
  names = tidyselect::all_of(c("sample_summaries", "large_numbat_pdfs")),
  reporter = "timestamp"
)
EOF
