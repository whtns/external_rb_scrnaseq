#!/bin/bash
#SBATCH --job-name=rebuild_hypoxia_summaries
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=48G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_summaries_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_summaries_%j.err

# Collate the freshly-rebuilt low/high/whole hypoxia collages into the persistent
# per-sample results/{sample}_hypoxia_summary.pdf. The heavy upstream (collages,
# score plots, hypoxia_seus) is already cached from job 9929429; this only runs
# collate_hypoxia_summary. numbatHelpers already has all fixes installed, so no
# document()/install here.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj
mkdir -p tmp   # the low/high collages live here as tempfiles that collate reads

# hypoxia_summaries has never been materialized, so tar_make builds it fresh;
# invalidate any pre-existing branches (0 is fine) so a re-run re-collates using
# the current collage tempfile paths. Its collage inputs stay cached (up to date).
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
meta <- tar_meta(fields="name", complete_only=FALSE)
to_inv <- meta$name[grepl("^hypoxia_summaries(_[0-9a-f]+)?$", meta$name)]
if (length(to_inv)) tar_invalidate(any_of(to_inv))
cat("invalidated:", length(to_inv), "hypoxia_summaries targets/branches\n")
' \
  || { echo "FATAL: invalidation step errored -- aborting before tar_make"; exit 1; }

Rscript run_targets_bg.R hypoxia_summaries
