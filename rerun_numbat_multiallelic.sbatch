#!/bin/bash
#SBATCH --job-name=nb_multiallelic
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --array=0-1
#SBATCH --output=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_multiallelic_%A_%a.log
#SBATCH --error=/project2/cobrinik_1090/external_rb_scrnaseq_proj/logs/nb_multiallelic_%A_%a.err

# Rerun numbat with multi_allelic = TRUE (numbat's own default; run_numbat.R
# previously hardcoded FALSE) to test whether 6p gains resolve differently.
#
# SRX10264524: chr6p amp segs by consensus round = 2, 2, 2, 0 -- the 6p gain was
#              called in rounds 1-3 and LOST in the final round, while
#              config/large_clone_comparisons.yaml annotates it 2_v_1_6p+.
# SRX14116944: 2, 1, 1, 1 -- one surviving 6p amp (0.3-49.5 Mb) that spans the
#              hardcoded MHC blacklist, so its breakpoints are undetermined.
#
# Params match the original snakemake run for these two samples (num_clones 5,
# t 1e-2, alpha 1e-3, gamma 20, max_entropy 0.7, min_LLR 5, tau 0.3, read_prop 1,
# ref sridhar_ref.rds) EXCEPT for three deliberate changes:
#
#   multi_allelic     FALSE -> TRUE   numbat's own default; enables bamp/bdel
#   check_convergence FALSE -> TRUE   the originals never tested convergence at
#                                     all, so max_iter was a blind budget
#
# max_iter stays at 4, matching the original runs, so segs_consensus_4.tsv here is
# directly comparable to the original round 4 with multi_allelic as the only
# difference. Iteration k does not depend on max_iter, so this is the same round 4
# an 8-iteration run would write.
#
# It was 8 for the first pass (job 10496424), to give convergence room to be
# reached. It never was: check_convergence never fired in 8 iterations, and rounds
# 5-8 added nothing but wander (non-neutral segs 23, 22, 25, 24, 26 for
# SRX10264524; 28, 28, 28, 30 for SRX14116944, with chr6p flickering 2,2,1,2).
# Freezing the diploid baseline did not help either -- see
# doc/diploid_baseline_pin_6p.md -- so the instability is in the segmentation, not
# in the iteration budget, and the extra four rounds cost ~3.5 h/sample to buy an
# equally unsettled answer. Round 4 is not convergence; it is a stated choice.
#
# Writes to output/numbat_multiallelic/<sample>/ -- the existing
# output/numbat_sridhar/ results are inputs to the targets pipeline and are NOT
# touched.

set -euo pipefail

SAMPLES=(SRX10264524 SRX14116944)
S=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

PROJ=/project2/cobrinik_1090/external_rb_scrnaseq_proj
OUTDIR=$PROJ/output/numbat_multiallelic/$S
mkdir -p "$OUTDIR"

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export TMPDIR=/scratch1/$USER/tmp
mkdir -p "$TMPDIR"

cd $PROJ

echo "=== $S : numbat rerun with multi_allelic=TRUE -> $OUTDIR ==="
echo "start: $(date)"

Rscript pipeline/scripts/run_numbat.R \
  "matrix_file='$PROJ/output/cellranger/$S/outs/filtered_feature_bc_matrix/matrix.mtx.gz'" \
  "seu_path='$PROJ/output/seurat/${S}_seu.rds'" \
  "allele_df='$PROJ/output/numbat/${S}_allele_counts.tsv.gz'" \
  "ref_path='/project2/cobrinik_1090/Homo_sapiens/numbat/sridhar_ref.rds'" \
  "out_dir='$OUTDIR'" \
  "init_k=5" \
  "t='1e-2'" \
  "alpha='1e-3'" \
  "gamma=20" \
  "max_entropy=0.7" \
  "min_LLR=5" \
  "max_iter=4" \
  "tau='0.3'" \
  "read_prop=1" \
  "multi_allelic=TRUE" \
  "check_convergence=TRUE" \
  "ncores=8"

echo "exit: $? ; end: $(date)"

# Immediate readout: chr6p calls per consensus round, old vs new.
echo "=== $S chr6p amp segments by consensus round ==="
for D in "$PROJ/output/numbat_sridhar/$S" "$OUTDIR"; do
  echo "--- $D ---"
  for f in "$D"/segs_consensus_*.tsv; do
    [ -f "$f" ] || continue
    n=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["CHROM"]==6 && $h["seg_start"]+0<58500000 && $h["cnv_state"]=="amp"{c++} END{print c+0}' "$f")
    b=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} $h["CHROM"]==6 && $h["seg_start"]+0<58500000 && ($h["cnv_state"]=="bamp"||$h["cnv_state"]=="bdel"){c++} END{print c+0}' "$f")
    echo "  $(basename "$f"): amp=$n  bamp/bdel=$b"
  done
done
