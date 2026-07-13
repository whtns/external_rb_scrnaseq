#!/bin/bash
#SBATCH --job-name=rebuild_hmcoll_low
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hmcoll_low_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hmcoll_low_%j.err

# Rebuild heatmap_collages_low_hypoxia with the auto-generated cell-cycle phase
# cluster names (auto_phase_level / assign_auto_phase_clusters). plot_seu_marker_heatmap
# lives in numbatHelpers, so a function-body change is not tracked by targets;
# tar_invalidate the target then tar_make. Crew heavy workers need the installed
# package, so document + install_local first.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript -e 'devtools::document("/project2/cobrinik_1090/rpkgs/numbat_helpers")' \
  && Rscript -e 'devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force = TRUE, upgrade = "never")' \
  || { echo "FATAL: package build/install failed"; exit 1; }

# Fail fast unless the new helpers are in the installed package.
Rscript -e 'stopifnot(exists("assign_auto_phase_clusters", where = asNamespace("numbatHelpers")), exists("auto_phase_level", where = asNamespace("numbatHelpers")))' \
  || { echo "FATAL: auto-phase helpers not in installed package"; exit 1; }

# Invalidate all branches of heatmap_collages_low_hypoxia (function-body change,
# not command-tracked), then rebuild. Inputs (seus_low_hypoxia) stay cached.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store = "_targets_r431")
meta <- tar_meta(fields = "name", complete_only = FALSE)
to_inv <- meta$name[grepl("^heatmap_collages_low_hypoxia(_[0-9a-f]+)?$", meta$name)]
if (length(to_inv)) tar_invalidate(any_of(to_inv))
cat("invalidated:", length(to_inv), "heatmap_collages_low_hypoxia targets/branches\n")
' || { echo "FATAL: invalidation errored"; exit 1; }

Rscript run_targets_bg.R heatmap_collages_low_hypoxia
