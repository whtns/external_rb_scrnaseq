#!/bin/bash
#SBATCH --job-name=rebuild_hypoxia_chain
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_chain_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_hypoxia_chain_%j.err

# Atomic rebuild of the ENTIRE hypoxia chain from the current filtered_seus.
#
# WHY THIS EXISTS (see AGENTS.md "path-string staleness" / memory
# seu_targets_depend_false_freeze): filtered_seus -> hypoxia_seus ->
# hypoxia_partition_paths -> seus_low/high_hypoxia -> effect_of_filtering all
# return *path strings*, so targets hashes the path text, not the file content.
# A regenerated filtered_seu.rds keeps the same path, so nothing downstream
# rebuilds on its own. The old rebuild_effect_of_filtering.sbatch invalidated
# only seus_low/high + the errored partition branches, leaving hypoxia_seus and
# the completed partition branches stale -> low/high counts drifted away from
# the new filtered counts (e.g. SRX11133593 filtered=4829 vs low+high=1744).
#
# This script invalidates the WHOLE chain as a unit so a partial refresh (the
# actual failure mode) is impossible. Always use THIS to refresh hypoxia after
# filtered_seus changes -- never invalidate a subset by hand.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# Install latest numbatHelpers so crew workers pick up the split_hypoxia_by_clusters
# tryCatch guard and the cone/rod dictionary fix.
Rscript -e 'devtools::install_local("/project2/cobrinik_1090/rpkgs/numbat_helpers", force=TRUE, quiet=TRUE); cat("numbatHelpers installed\n")'

# Invalidate the full hypoxia chain (parents + all dynamic branches) atomically.
Rscript -e '
suppressPackageStartupMessages(library(targets))
tar_config_set(store="_targets_r431")
meta <- tar_meta(fields = "name", complete_only = FALSE)
# Match each chain target and its "<name>_<branchhash>" children, but nothing else.
chain <- c("hypoxia_seus", "hypoxia_threshold_per_sample", "hypoxia_score_plots",
           "hypoxia_partition_paths", "seus_low_hypoxia", "seus_high_hypoxia",
           "hypoxia_split_log_collated", "effect_of_filtering",
           "effect_of_filtering_collated")
pat <- paste0("^(", paste(chain, collapse = "|"), ")(_[0-9a-f]+)?$")
to_inv <- meta$name[grepl(pat, meta$name)]
tar_invalidate(any_of(to_inv))
cat("Invalidated", length(to_inv), "hypoxia-chain targets/branches\n")
print(sort(intersect(chain, to_inv)))
'

# Rebuild the collated effect_of_filtering PDF (pulls the whole chain through) and
# the split decision-log collation.
Rscript run_targets_bg.R effect_of_filtering_collated hypoxia_split_log_collated
