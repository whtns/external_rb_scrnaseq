#!/bin/bash
#SBATCH --job-name=hypoxia_heatmaps
#SBATCH --output=logs/hypoxia_gene_heatmaps_%j.out
#SBATCH --error=logs/hypoxia_gene_heatmaps_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=epyc-64

module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

cd /project2/cobrinik_1090/external_rb_scrnaseq_proj

Rscript - <<'EOF'
devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers")

paths <- c(
  "output/seurat/SRX10031194_hypoxia_low_seu.rds",
  "output/seurat/SRX10264517_hypoxia_low_seu.rds",
  "output/seurat/SRX10264518_hypoxia_low_seu.rds",
  "output/seurat/SRX10264519_hypoxia_low_seu.rds",
  "output/seurat/SRX10264520_hypoxia_low_seu.rds",
  "output/seurat/SRX10264521_hypoxia_low_seu.rds",
  "output/seurat/SRX10264522_hypoxia_low_seu.rds",
  "output/seurat/SRX10264523_hypoxia_low_seu.rds",
  "output/seurat/SRX10264524_hypoxia_low_seu.rds",
  "output/seurat/SRX10264525_hypoxia_low_seu.rds",
  "output/seurat/SRX10264526_hypoxia_low_seu.rds",
  "output/seurat/SRX10831280_hypoxia_low_seu.rds",
  "output/seurat/SRX10831281_hypoxia_low_seu.rds",
  "output/seurat/SRX10831282_hypoxia_low_seu.rds",
  "output/seurat/SRX10831283_hypoxia_low_seu.rds",
  "output/seurat/SRX10831286_hypoxia_low_seu.rds",
  "output/seurat/SRX10831287_hypoxia_low_seu.rds",
  "output/seurat/SRX11133585_hypoxia_low_seu.rds",
  "output/seurat/SRX11133587_hypoxia_low_seu.rds",
  "output/seurat/SRX11133588_hypoxia_low_seu.rds",
  "output/seurat/SRX11133589_hypoxia_low_seu.rds",
  "output/seurat/SRX11133592_hypoxia_low_seu.rds",
  "output/seurat/SRX11133593_hypoxia_low_seu.rds",
  "output/seurat/SRX11133594_hypoxia_low_seu.rds",
  "output/seurat/SRX14116944_hypoxia_low_seu.rds",
  "output/seurat/SRX14116945_hypoxia_low_seu.rds",
  "output/seurat/SRX14116946_hypoxia_low_seu.rds",
  "output/seurat/SRX14116947_hypoxia_low_seu.rds",
  "output/seurat/SRX14116948_hypoxia_low_seu.rds",
  "output/seurat/SRX22868102_hypoxia_low_seu.rds",
  "output/seurat/SRX22868103_hypoxia_low_seu.rds",
  "output/seurat/SRX22868104_hypoxia_low_seu.rds",
  "output/seurat/SRX22868105_hypoxia_low_seu.rds"
)

results <- purrr::map(paths, function(p) {
  message("Processing: ", p)
  tryCatch(
    plot_hypoxia_gene_heatmap(p),
    error = function(e) { message("ERROR for ", p, ": ", e$message); NA_character_ }
  )
})

message("Done. Output paths:")
purrr::walk(results, ~ message(.x))
EOF
