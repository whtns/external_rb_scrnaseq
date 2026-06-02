.libPaths(c("/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/.R_libs", .libPaths()))

pkgs <- c(
  # pinned: only data.table which has a dev version installed locally
  "data.table@1.18.2.1",

  # GitHub packages with SHAs
  "whtns/numbat@b9354cd",
  "cobriniklab/seuratTools",
  "mojaveazure/seurat-disk@9b89970eac2a3bd770e744f63c7763419486b14c",
  "CCICB/TCGAgistic@adf0782554c2d7916169cf81f69d19247047b7eb",

  # packages.R
  "targets",
  "tarchetypes",
  "crew",
  "future",
  "dotenv",
  "tidyverse",
  "fs",
  "readxl",
  "Seurat",
  "patchwork",
  "msigdbr",
  "glue",
  "ggplotify",
  "rlang",
  "clusterProfiler",
  "DOSE",
  "writexl",
  "magick",
  "EnhancedVolcano",
  "ggpubr",
  "plotly",
  "htmlwidgets",
  "igraph",
  "scales",
  "ggupset",
  "ggridges",
  "DESeq2",
  "dendextend",
  "ggbeeswarm",
  "pals",
  "clustifyr",
  "speckle",
  "ggraph",
  "org.Hs.eg.db",
  "plyranges",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "AnnotationHub",
  "maftools",
  "janitor",
  "plotgardener",
  "conflicted"
)

pak::lockfile_create(
  pkgs,
  lockfile = "/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/pak.lock"
)
