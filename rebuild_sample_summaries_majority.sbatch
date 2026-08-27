#!/bin/bash
#SBATCH --job-name=samp_summ_maj
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=36:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_sample_summaries_maj_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/rebuild_sample_summaries_maj_%j.err

# Rebuild sample_summaries after the majority-of-rounds numbat RDS migration.
#
# Pass 1 (job 11381877) did NOT work: tar_outdated reported unfiltered_seus and
# filtered_seus outdated, but their depend = FALSE cue made tar_make dispatch
# zero branches, and sample_summaries skipped because its inputs are path
# strings. Clone trees/heatmaps/PDFs did rebuild off numbat_rds_files, leaving
# the store INCONSISTENT: new-round trees against old-round Seurat metadata.
#
# So invalidate the gated subgraph explicitly first, then build.
# See src/invalidate_numbat_downstream.R for the two gates.
#
# Heavy work runs on crew SLURM workers (epyc-64); this job only orchestrates.
# Per AGENTS.md: only one tar_make() may run against _targets_r431 at a time.

set -euo pipefail

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

module purge
module load gcc/13.3.0 r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

Rscript src/invalidate_numbat_downstream.R
Rscript run_targets_bg.R sample_summaries
