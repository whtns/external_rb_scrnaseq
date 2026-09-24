#!/bin/bash
#SBATCH --job-name=nb_rds_t1e5
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --array=0-31%8
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_rds_t1e5_%A_%a.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_rds_t1e5_%A_%a.err

# Build *_numbat.rds from the t=1e-5 arm so the targets pipeline can consume it.
# Writes output/numbat_t1e5/<sample>_numbat.rds. Does NOT touch
# output/numbat_sridhar/ -- the production cohort stays exactly as it is, and can
# be returned to by reverting the one-line path change in R/pipeline_targets_inputs.R.
#
# ROUND POLICY: round_mode='final' -- the last round whose nine artefacts are all
# present. NOT the selected-round manifest. At t=1e-5 selection beats the final
# round on 1 of 34 samples and misses zero union events (vs +34 canonical events at
# t=1e-2), so there is nothing left to buy with the selecting-on-the-outcome bias.
# See docs/srx_vs_srr_numbat_settings.md section 8.
#
# SCOPE: SRX only. The 36 dirs under output/numbat_t1e5/ are the rerunnable SRX
# cohort; SRX11133592/93/94 were never in it (EXCLUDE in rerun_numbat_arm.sbatch)
# and no SRR sample is touched by this script or reachable from it.
#
# EXPECTED SKIPS -- 4 of 36 dirs have no complete round and are not built:
#   SRX10031191  numbat never ran (6 cells vs hardcoded min_cells = 10)
#   SRX11133590  "No CNV remains after filtering by LLR in pseudobulks"
#   SRX10031192  "No CNV remains after filtering by entropy in single cells"
#   SRX11133586  same
# That leaves a 32-sample cohort. The skip is deliberate (user decision
# 2026-09-21): these are honest "no CNV at this t" results, not crashes to patch.
#
# PREREQUISITE: SRX22868103's 240 GB rerun (job 12218137) must have landed, or it
# builds from its 3-round OOM remnant instead of the full 4 rounds.

set -uo pipefail

PROJ=/project2/cobrinik_1090/external_rb_scrnaseq_proj
SRC=$PROJ/output/numbat_t1e5
cd "$PROJ"

module load r/4.4.1
# Deliberately NOT setting R_LIBS: numbat 1.2.2 lives in ~/R/numbat-1.2.2 for the
# version replay and must never be the one that builds production objects.
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

# Buildable = has done.txt AND a round with all nine artefacts present.
mapfile -t SAMPLES < <(
  for d in "$SRC"/SRX*/; do
    s=$(basename "$d")
    [ -f "$d/done.txt" ] || continue
    for k in $(ls "$d"/segs_consensus_*.tsv 2>/dev/null \
               | sed 's/.*segs_consensus_\([0-9]*\)\.tsv/\1/' | sort -rn); do
      ok=1
      for f in segs_consensus_$k.tsv clone_post_$k.tsv joint_post_$k.tsv \
               exp_post_$k.tsv allele_post_$k.tsv geno_$k.tsv \
               treeML_$k.rds mut_graph_$k.rds tree_final_$k.rds; do
        [ -f "$d/$f" ] || { ok=0; break; }
      done
      [ $ok -eq 1 ] && { echo "$s"; break; }
    done
  done | sort
)

echo "buildable samples: ${#SAMPLES[@]}"
if [ "${SLURM_ARRAY_TASK_ID:-0}" -ge "${#SAMPLES[@]}" ]; then
  echo "array index ${SLURM_ARRAY_TASK_ID} beyond list; nothing to do"; exit 0
fi
S=${SAMPLES[${SLURM_ARRAY_TASK_ID:-0}]}

# Never let the 1.2.2 replay library build a production object.
NB=$(Rscript -e 'cat(as.character(packageVersion("numbat")))' 2>/dev/null)
[ "$NB" = "1.5.2" ] || { echo "FATAL: numbat $NB, expected 1.5.2" >&2; exit 1; }

echo "=== $S (numbat $NB) ==="
Rscript pipeline/scripts/process_numbat_rds.R \
  "done_file='$SRC/$S/done.txt'" \
  "round_mode='final'"
RC=$?

OUT=$SRC/${S}_numbat.rds
if [ $RC -eq 0 ] && [ -f "$OUT" ]; then
  echo "OK  $(du -h "$OUT" | cut -f1)  $OUT"
else
  echo "FAILED rc=$RC, rds present: $([ -f "$OUT" ] && echo yes || echo no)" >&2
fi
exit $RC
