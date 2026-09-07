#!/bin/bash
#SBATCH --job-name=depth_clone_assoc
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/depth_clone_assoc_%j.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/depth_clone_assoc_%j.err

# Per-cell test of whether numbat clones track read depth (#43 follow-up).
#
# Read-only: opens output/numbat_sridhar/SRX*_numbat.rds and cell_qc_values in
# batch_hashes.sqlite, writes only results/ + docs/. No tar_make(), nothing
# written to _targets_r431, so the one-job-per-store rule is not engaged and this
# can run alongside anything.
#
# Stage 1 smokes two samples (one with a normal clone + 3 tumor clones, one
# 2-clone) and aborts the job if that fails, so a typo does not burn 33 RDS reads.

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"
cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

# load_all, not install_local: nothing here runs on a crew worker.
Rscript -e 'suppressMessages(devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers"))' \
  >/dev/null 2>&1 || echo "note: load_all check failed; script uses the installed package"

echo "############ stage 1: smoke (2 samples) ############"
Rscript src/diag_depth_clone_association.R SRX10264519 SRX22868104
S1=$?
echo "### smoke exit: $S1"
[ $S1 -eq 0 ] || { echo "FATAL: smoke failed; not running the cohort"; exit 1; }

echo
echo "############ stage 2: full cohort ############"
Rscript src/diag_depth_clone_association.R
S2=$?
echo "### cohort exit: $S2"

if [ -f "$HOME/.slack_webhook" ]; then
  curl -s -X POST -H 'Content-type: application/json' \
    --data "{\"text\":\"depth_clone_assoc ($SLURM_JOB_ID) finished: smoke $S1, cohort $S2\"}" \
    "$(cat "$HOME/.slack_webhook")" >/dev/null || true
fi

exit $S2
