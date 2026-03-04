# Retinoblastoma scRNA-seq Analysis Pipeline

Single-cell RNA-seq analysis of retinoblastoma tumours with somatic copy-number
alteration (SCNA) stratification across four clonal subtypes: **1q+**, **2p+**,
**6p+**, and **16q−**.

---

## Quick start (current pipeline)

The main analysis runs with the R `{targets}` framework from the **project root**,
not from this directory.

```r
# inside R, from project root
library(targets)

tar_make()                                        # run all outdated targets
tar_make(names = fig_02)                          # run one specific target
tar_visnetwork()                                  # interactive dependency graph
tar_meta(fields = c("error", "warnings"))         # inspect failures
```

To run in a persistent session:

```bash
# from project root — run inside tmux for long jobs
tmux new -s targets
Rscript -e "targets::tar_make()"
```

Key pipeline files at the project root:

| File | Purpose |
|------|---------|
| `_targets.R` | Pipeline definition (all targets) |
| `packages.R` | `library()` calls loaded by `_targets.R` |
| `R/` | ~80 modular function files sourced automatically |
| `config/` | YAML configuration (clusters, markers, clone rules) |
| `batch_hashes.sqlite` | SQLite database tracking Seurat object metadata |
| `pipeline_overview.md` | Mermaid dependency diagram of the full pipeline |

---

## Pipeline stages

```
Configuration YAMLs
    └─► Seurat processing (unfiltered → filtered → regressed → hypoxia-scored)
            └─► Debranching by SCNA clone (1q / 2p / 6p / 16q)
                    ├─► Numbat CNV heatmaps & probability tracks
                    ├─► Differential expression (clone × cluster)
                    ├─► Integration across samples per SCNA type
                    ├─► Corresponding-state analysis (2p + 6p cross-SCNA)
                    ├─► Enrichment & oncoprint
                    └─► Manuscript figures (fig_02–10, fig_s04–s25)
```

Static branching via `tarchetypes::tar_map()` maps identical analyses over all
four SCNA types (see `scna_map_values` tibble at the top of `_targets.R`).

---

## Configuration files

All live in `config/` and are loaded as named lists by `_targets.R`:

| File | Contents |
|------|---------|
| `cluster_dictionary.yaml` | Cluster → cell-type annotation mappings |
| `celltype_markers.yaml` | Marker genes per cell type |
| `large_clone_comparisons.yaml` | Per-sample clone-pair comparisons |
| `large_clone_simplifications.yaml` | Simplification of clone arm labels |
| `large_filter_expression.yaml` | CNV score filtering thresholds |

---

## Metadata tracking (SQLite)

`batch_hashes.sqlite` at the project root stores Seurat object metadata across
pipeline runs. Schema is initialised automatically by `R/seu_metadata_db.R`:

| Table | Key | Contents |
|-------|-----|---------|
| `seurat_objects` | `filepath` | n_cells, sample_id, scna_type, processing stage, assay, timestamp |
| `cell_metadata` | `(filepath, column)` | dtype, n_unique, n_NA, JSON quantile/count summary |
| `cluster_composition` | `(filepath, cluster)` | cell counts and percentages per cluster |
| `qc_metrics` | `(filepath, metric)` | min/q25/median/mean/q75/max for nCount/nFeature/percent.mt |
| `hashes` | `filepath` | content hash (digest of cell barcodes) — legacy |
| `cluster_orders` | `file_id` | JSON cluster ordering — legacy |

The `seu_metadata_db` target re-populates these tables whenever any upstream
Seurat file changes (detected via `format = "file"` on `debranched_seus_*` targets).

---

## Key outputs

| Location | Contents |
|----------|---------|
| `output/seurat/` | Processed Seurat RDS files (filtered, debranched, integrated) |
| `output/numbat/` | Numbat CNV inference outputs per sample |
| `results/` | Final PDFs, CSVs, and RDS DE results (~22 GB, 2 600+ files) |
| `results/fig_*.pdf` | Manuscript figure panels |
| `results/fig_s*.pdf` | Supplementary figure panels |
| `results/table_*.xlsx` | Manuscript data tables |

---

## Requirements

R packages are declared in `packages.R`. Key dependencies:

- **Pipeline:** `targets`, `tarchetypes`, `crew`
- **Single-cell:** `Seurat`, `SeuratDisk`, `seuratTools`, `numbat`
- **DE & enrichment:** `DESeq2`, `clusterProfiler`, `DOSE`, `msigdbr`
- **Genomics:** `maftools`, `TCGAgistic`, `plyranges`, `org.Hs.eg.db`
- **Visualisation:** `patchwork`, `ggpubr`, `EnhancedVolcano`, `pals`, `plotgardener`
- **Data:** `DBI`, `RSQLite`, `jsonlite`, `digest`, `glue`, `fs`, `readxl`, `writexl`

---

## Legacy: ARMOR upstream pre-processing (this directory)

This `pipeline/` subdirectory contains the **ARMOR** Snakemake workflow used for
the original upstream pre-processing (STAR alignment → salmon quantification →
read QC). It is not part of the active analysis — its outputs feed into
`output/seurat/` and are treated as pre-computed inputs by `_targets.R`.

Retained for reproducibility. Re-running is only needed if raw FASTQs change.

```bash
# if you need to re-run the upstream pipeline:
cd pipeline/
conda activate ARMOR
snakemake -n --use-conda          # dry run
snakemake -j 6 --use-conda        # full run (run inside tmux)

# dependency graph
snakemake --rulegraph | dot -Tpdf > rulegraph.pdf
```

> If you get a `target cannot be locked` error, check that Snakemake is not
> already running in another session.
