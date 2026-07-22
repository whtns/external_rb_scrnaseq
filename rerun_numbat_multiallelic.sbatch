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
#   max_iter          4     -> 8      room for convergence to be reached
#
# The iteration changes do not cost a separate arm: numbat writes
# segs_consensus_N.tsv for EVERY iteration, so segs_consensus_4.tsv from this run
# is directly comparable to the original round 4 (multi_allelic being the only
# difference at that round), while later rounds show where the run actually
# settles. One job answers both questions.
#
# Justification for touching iterations at all: SRX10264524's original run went
# 19 -> 12 -> 10 -> 11 non-neutral segments genome-wide -- still moving at the
# last round -- and its two round-3 chr6p amps sat at LLR 8.2 / 8.4 against
# min_LLR 5, i.e. marginal calls flickering across the threshold. Reproducing an
# unstable answer would not tell us much.
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
  "max_iter=8" \
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
