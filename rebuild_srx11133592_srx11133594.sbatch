#!/bin/bash
#SBATCH --job-name=rebuild-srx11133592-srx11133594
#SBATCH --account=cobrinik_1090
#SBATCH --partition=epyc-64
#SBATCH --qos=normal
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --output=logs/rebuild_srx11133592_srx11133594_%j.log
#SBATCH --error=logs/rebuild_srx11133592_srx11133594_%j.err

# Force-invalidate large/filtered/low_hypoxia numbat PDFs (built before re-copy from SRR),
# bulk_clones aggregates, and sample_summaries for SRX11133592/SRX11133594.
# large_numbat_pdfs branch suffix 1a5916df10f64d45=SRX11133592, 7b00063c93e5a274=SRX11133594

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript run_targets_bg.R \
  large_numbat_pdfs_1a5916df10f64d45 \
  large_numbat_pdfs_7b00063c93e5a274 \
  filtered_numbat_pdfs_1a5916df10f64d45 \
  filtered_numbat_pdfs_7b00063c93e5a274 \
  low_hypoxia_numbat_pdfs_1a5916df10f64d45 \
  low_hypoxia_numbat_pdfs_7b00063c93e5a274 \
  unfiltered_numbat_bulk_clones \
  filtered_numbat_bulk_clones \
  low_hypoxia_numbat_bulk_clones \
  sample_summaries_a4c3ddb220bdc62c \
  sample_summaries_f1ddced65b883049 \
  _targets_r431
