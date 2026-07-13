#!/bin/bash
#SBATCH --job-name=rebuild_eof_plot
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_eof_plot_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_eof_plot_%j.err

# Regenerate ONLY the effect_of_filtering plot after changing plot_effect_of_filtering
# (plot_functions_40.R): the low/high hypoxia scna DimPlots are now hypoxia_score
# gradient FeaturePlots. The hypoxia seus themselves are unchanged (Jul 2 rebuild),
# so we invalidate only the plot targets -- NOT the hypoxia chain.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# install_local so crew workers (effect_of_filtering runs on workers) get the new
# plot code; load_all alone won't reach them.
Rscript -e 'devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, quiet=TRUE); cat("numbatHelpers installed\n")'

# Invalidate only effect_of_filtering (+ its branches) and the collated PDF.
Rscript -e '
suppressPackageStartupMessages(library(targets))
tar_config_set(store="_targets_r431")
meta <- tar_meta(fields = "name", complete_only = FALSE)
pat <- "^(effect_of_filtering|effect_of_filtering_collated)(_[0-9a-f]+)?$"
to_inv <- meta$name[grepl(pat, meta$name)]
tar_invalidate(any_of(to_inv))
cat("Invalidated", length(to_inv), "effect_of_filtering targets/branches\n")
'

Rscript run_targets_bg.R effect_of_filtering_collated
