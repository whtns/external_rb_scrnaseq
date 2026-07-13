#!/bin/bash
#SBATCH --job-name=rebuild-sample-summaries-vector-exp
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_sample_summaries_vector_exp_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_sample_summaries_vector_exp_%j.err

# Rebuild all sample_summaries with updated convert_numbat_pngs() that renders
# exp_roll_clust as a true vector PDF (base R pdf() device, no ragg/Cairo).
# Steps:
#   1. Install updated numbatHelpers so crew workers pick up the new code.
#   2. Invalidate all large/filtered/low_hypoxia_numbat_pdfs branches.
#   3. tar_make(sample_summaries) — rebuilds the full chain.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript - <<'EOF'
# Install so all crew workers get the updated convert_numbat_pngs()
devtools::install_local(
  "/project2/cobrinik_1090/rpkgs/numbat_helpers",
  force = TRUE,
  upgrade = "never",
  quiet = TRUE
)
cat("numbatHelpers installed.\n"); flush(stdout())

targets::tar_config_set(store = "_targets_r431")

# Invalidate all branches of the three numbat_pdfs targets so that
# convert_numbat_pngs() re-runs for every sample.
targets::tar_invalidate(
  names = tidyselect::starts_with("large_numbat_pdfs") |
          tidyselect::starts_with("filtered_numbat_pdfs") |
          tidyselect::starts_with("low_hypoxia_numbat_pdfs")
)
cat("Invalidated numbat_pdfs branches.\n"); flush(stdout())

targets::tar_make(
  names = tidyselect::all_of("sample_summaries"),
  reporter = "timestamp"
)
EOF
