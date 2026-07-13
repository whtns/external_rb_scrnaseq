#!/bin/bash
#SBATCH --job-name=effect_of_filtering
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_effect_of_filtering_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_effect_of_filtering_%j.err

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Install latest numbatHelpers so crew workers pick up the tryCatch fix in
# split_hypoxia_by_clusters (subset_seu_by_expression "No cells found" guard).
Rscript -e 'devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, quiet=TRUE); cat("numbatHelpers installed\n")'

# Invalidate only the 2 errored hypoxia_partition_paths branches and the
# downstream targets that depend on them. Completed branches stay cached.
Rscript -e '
suppressPackageStartupMessages(library(targets))
tar_config_set(store="_targets_r431")
meta <- tar_meta(fields=c(name,error), complete_only=FALSE)
err_branches <- meta$name[!is.na(meta$error) & nzchar(meta$error) &
                           grepl("^hypoxia_partition_paths_", meta$name)]
if (length(err_branches) > 0) {
  tar_invalidate(any_of(err_branches))
  cat("Invalidated", length(err_branches), "errored hypoxia_partition_paths branches\n")
}
tar_invalidate(any_of(c("seus_low_hypoxia","seus_high_hypoxia",
                         "hypoxia_split_log_collated","effect_of_filtering",
                         "effect_of_filtering_collated")))
cat("Invalidated downstream targets\n")
'

Rscript run_targets_bg.R effect_of_filtering_collated hypoxia_split_log_collated
