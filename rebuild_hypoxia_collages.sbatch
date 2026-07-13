#!/bin/bash
#SBATCH --job-name=rebuild_hypoxia_collages
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_collages_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_collages_%j.err

# Rebuild ONLY the collage/heatmap plotting targets that failed in job 9926394 on
# degenerate hypoxia subsets, now that plot_seu_marker_heatmap /
# plot_hypoxia_gene_heatmap / find_diffex_clones_in_phase guard those cases.
# The expensive upstream (clone trees, clone-based diffex, sample_summaries) built
# successfully last run and is NOT touched here -- the fix is downstream of them.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# document() so NAMESPACE reflects any new @export; upgrade="never" so a failing
# dependency upgrade can't leave a stale build; verify a known symbol is present.
Rscript -e 'suppressMessages({devtools::document("/project2/cobrinik_1090/rpkgs/numbat_helpers", quiet=TRUE); devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, upgrade="never", quiet=TRUE)}); stopifnot(exists("myreadxl", envir=asNamespace("numbatHelpers"), inherits=FALSE)); cat("numbatHelpers installed; myreadxl present\n")' \
  || { echo "FATAL: numbatHelpers install/verify failed -- aborting before tar_make"; exit 1; }

# These plotting targets typically carry cue(command=FALSE), so a function-body
# change alone won't mark them outdated -- invalidate them explicitly (base target
# + all dynamic branches) so tar_make actually re-runs them with the fixed code.
Rscript -e '
suppressPackageStartupMessages(library(targets)); tar_config_set(store="_targets_r431")
meta <- tar_meta(fields="name", complete_only=FALSE)
bases <- c("heatmap_collages_low_hypoxia","heatmap_collages_high_hypoxia",
           "heatmap_collages_6p","annotated_heatmap_collages",
           "hypoxia_gene_heatmap_low_hypoxia","diffex_6p_g1")
pat <- paste0("^(", paste(bases, collapse="|"), ")(_[0-9a-f]+)?$")
to_inv <- meta$name[grepl(pat, meta$name)]
tar_invalidate(any_of(to_inv))
cat("invalidated:", length(to_inv), "collage/heatmap targets/branches\n")
stopifnot(length(to_inv) > 0)
' \
  || { echo "FATAL: invalidation step failed -- aborting before tar_make"; exit 1; }

Rscript run_targets_bg.R \
  heatmap_collages_low_hypoxia heatmap_collages_high_hypoxia heatmap_collages_6p \
  annotated_heatmap_collages hypoxia_gene_heatmap_low_hypoxia diffex_6p_g1
