#!/bin/bash
#SBATCH --job-name=rebuild-filt-seus-dict
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_filtered_seus_dict_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_filtered_seus_dict_%j.err

# Propagate the cone/rod >=3 dictionary fix into filtered_seus.
# Steps:
#   1. Install updated numbatHelpers so crew workers pick up the latest code.
#   2. Invalidate cluster_dictionary (targets cannot detect the in-package
#      function-body change on its own) AND the cue-frozen seu targets that
#      consume it (depend=FALSE, so they never auto-invalidate).
#   3. tar_make(filtered_seus) — rebuilds the chain through crew workers.

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
  names = tidyselect::starts_with("cluster_dictionary") |
          tidyselect::starts_with("unfiltered_seus") |
          tidyselect::starts_with("filtered_seus")
)
cat("Invalidated cluster_dictionary(+per_sample), unfiltered_seus, filtered_seus.\n")
flush(stdout())

targets::tar_make(
  names = tidyselect::all_of("filtered_seus"),
  reporter = "timestamp"
)
EOF
