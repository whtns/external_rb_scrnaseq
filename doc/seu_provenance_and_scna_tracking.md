# SEU Provenance and "scna" Metadata Tracking

## Overview

This document traces the provenance of Seurat (SEU) objects through the targets pipeline and tracks where and how the "scna" (somatic copy number alteration) metadata column is created and stored.

## Target Flow 📊

### Level 1: Input Loading
**File**: `pipeline_targets_inputs.R`

**Targets**:
- **numbat_rds_files** & **seus**: Loaded from disk via `tarchetypes::tar_files()`
  - Created by Snakemake pipeline (numbat_sridhar rules)
  - Raw Seurat objects without scna metadata at this stage
  - Tracked from `output/numbat_sridhar/` and `output/seurat/` directories

### Level 2: Metadata Enrichment
**File**: `pipeline_targets_seurat.R`

**Key Targets**:
- **unfiltered_seus** → calls `prep_unfiltered_seu(numbat_rds_files, ...)`
  - ✅ **FIRST adds "scna" metadata column**
- **filtered_seus** → calls `filter_cluster_save_seu(numbat_rds_files, seus, ...)`
  - ✅ Also adds "scna" metadata column
- **final_seus** → `set_final_seus(interesting_samples)`
  - Uses refined Seurat objects with scna metadata

### Level 3: SCNA-stratified Collections
**File**: `pipeline_targets_seurat.R`

These are file-tracking targets that organize Seurat objects by SCNA type:
- **debranched_seus_1q**: Seurat objects with 1q alterations (5 samples)
- **debranched_seus_2p**: Seurat objects with 2p alterations (6 samples)
- **debranched_seus_6p**: Seurat objects with 6p alterations (4 samples)
- **debranched_seus_16q**: Seurat objects with 16q alterations (3 samples)

Note: The SCNA suffix in the filename (e.g., `_1q.rds`) pre-filters objects, but metadata is formally added in Level 2.

### Level 4: Database Population
**File**: `pipeline_targets_seurat.R`

**Target**: **seu_metadata_db**
- Calls `bulk_extract_seu_metadata(scna_seus, ...)`
- Reads Seurat objects and populates SQLite database
- Deployment mode: "main" (single-threaded for database writes)
- References: `scna_seus` (combined SCNA-stratified collections)

## Where "scna" is First Added 🎯

### Source File
[R/plot_functions_3.R](../R/plot_functions_3.R) — Lines 50-68

### Functions
The "scna" metadata column is added by two functions:
1. **`prep_unfiltered_seu()`** — Creates unfiltered Seurat objects with scna metadata
2. **`filter_cluster_save_seu()`** — Creates filtered Seurat objects with scna metadata

### Implementation Details

**Input Data Source**:
- Numbat clone genotypes: `GT_opt` column from numbat post-processing RDS files
- Clone simplification config: `large_clone_simplifications` (YAML file)
  - Maps complex genotypes to categorical SCNA labels: "1q", "2p", "6p", "16q"

**Processing Steps** (from plot_functions_3.R):
```r
# Extract clone post-processing data from numbat output
nb_clone_post <- # numbat clone genotype table

# Build metadata table from clone genotypes
scna_metadata <- nb_clone_post %>%
  tibble::rownames_to_column("cell") %>%
  dplyr::rowwise() %>%
  # Map genotype to SCNA category
  dplyr::mutate(scna = simplify_gt_col(GT_opt, large_clone_simplifications)) %>%
  # Clean up SCNA labels
  dplyr::mutate(scna = wrap_scna_labels(scna)) %>%
  tidyr::unnest(scna) %>%
  dplyr::distinct(cell, .keep_all = TRUE) %>%
  tibble::column_to_rownames("cell") %>%
  dplyr::mutate(scna = as.character(scna))

# Match cells in metadata to actual Seurat object cells
scna_metadata <- scna_metadata[colnames(seu), ]

# Add to Seurat metadata
seu <- Seurat::AddMetaData(seu, scna_metadata)
```

**Result**:
- Each cell in the Seurat object receives one or more SCNA labels
- SCNA column is categorical: values in {"1q", "2p", "6p", "16q", or combinations}
- Missing values possible if clones lack clear SCNA assignments

## SQLite Database Tracking 📁

**Yes, "scna" is tracked in the sqlite database — but via two mechanisms:**

### Mechanism 1: seurat_objects Table
**Table**: `seurat_objects`  
**Relevant Column**: `scna_type`

**Important Note**: This is **inferred from filename**, not from metadata:
```r
scna_type <- case_when(
  grepl("_1q",  filepath) ~ "1q",
  grepl("_2p",  filepath) ~ "2p",
  grepl("_6p",  filepath) ~ "6p",
  grepl("_16q", filepath) ~ "16q",
  TRUE ~ NA_character_
)
```

**When populated**:
- Run 1 during database initialization via `init_seu_metadata_db()`
- Updated for each Seurat file by `extract_seu_metadata()`

### Mechanism 2: cell_metadata Table
**Table**: `cell_metadata`  
**Key Row**: One row per Seurat file × column name

**For the "scna" column**:
- **filepath**: Path to RDS file
- **column_name**: "scna" (the metadata column name)
- **dtype**: "character" (or "factor" depending on Seurat version)
- **n_cells**: Total cells in object
- **n_unique**: Count of unique scna values (e.g., 4 for {"1q", "2p", "6p", "16q"})
- **n_na**: Count of NA values (cells with no clear SCNA)
- **summary_json**: Frequency counts of top 50 values
  - Example: `{"1q": 2341, "2p": 1203, "16q": 891, "6p": 45, …}`

**When populated**:
- Called once per Seurat file by `extract_seu_metadata()` 
- Loops through all metadata columns and inserts/upserts into cell_metadata table

## Dependency Chain Summary

```
Snakemake RDS files 
("output/seurat/*filtered_seu.rds")
    ↓
numbat_rds_files (targets tracked files)
    ↓
unfiltered_seus / filtered_seus 
    ↓ (call prep_unfiltered_seu / filter_cluster_save_seu)
    ↓ 
    ← ← ← [ADD "scna" metadata column HERE] ← ← ←
    ↓
debranched_seus / final_seus / scna_seus
    (SCNA-stratified collections)
    ↓
seu_metadata_db 
    (calls bulk_extract_seu_metadata)
    ↓ (reads "scna" column from metadata)
    ↓
sqlite: 
  - seurat_objects.scna_type (inferred from filename)
  - cell_metadata[filepath × "scna"] (summary stats + frequencies)
```

## Key Files

| File | Purpose |
|------|---------|
| `pipeline_targets_inputs.R` | Load SEU files from disk |
| `pipeline_targets_seurat.R` | Enrich metadata, track by SCNA, populate database |
| `R/plot_functions_3.R` | **ADD "scna" column** (prep_unfiltered_seu, filter_cluster_save_seu) |
| `R/seu_metadata_db.R` | SQLite database schema and population logic |
| `_targets.R` | Main pipeline definition |

## Configuration Files

| File | Purpose |
|------|---------|
| `config/large_clone_simplifications.yaml` | Maps genotypes → SCNA labels |
| `batch_hashes.sqlite` | SQLite database storing SEU metadata |

## Sample-Level Metadata

**Samples and their SCNA content** (from _targets.R):

```r
rb_scna_samples <- list(
  "1q"  = c("SRR13884246", "SRR13884249", "SRR14800534", "SRR14800535", "SRR14800536"),
  "2p"  = c("SRR13884246", "SRR13884247", "SRR13884248", "SRR13884249", "SRR17960481", "SRR17960484"),
  "6p"  = c("SRR13884247", "SRR13884248", "SRR17960484"),
  "16q" = c("SRR14800534", "SRR14800535", "SRR14800536")
)
```

Note: Some samples carry multiple SCNA types (e.g., SRR13884246 has both 1q and 2p alterations).

## Data Quality Considerations

1. **Clones without clear SCNA**: Cells may have NA in the "scna" column if their clone lacks a clear SCNA assignment
2. **Dual alterations**: Cells can belong to clones with multiple concurrent SCNAs (e.g., "1q+2p")
3. **Filename vs. metadata mismatch**: The `scna_type` in seurat_objects table (inferred from filename) may differ from actual cell-level "scna" values if filtering has occurred
4. **Database staleness**: The SQLite database is only recomputed when debranched_seus_*q files change; manual edits to Seurat metadata may not be reflected until targets are invalidated

## Related Targets

- **regressed_seus**: Regression-adjusted versions of final_seus; inherits scna metadata
- **chosen_resolution_seus**: Phase cluster assignments per scna stratum
- **scna_seus**: Union of all SCNA-stratified objects, used for database population
- **filtering_recurrent_heatmaps**, **regressed_recurrent_heatmaps**: Group by scna for visualization

## Annotated Guide To pipeline_targets_seurat.R

This section documents the structure of [R/pipeline_targets_seurat.R](R/pipeline_targets_seurat.R) in execution order, with the intent of making maintenance and dependency tracing easier.

### 1. File-Tracked SCNA Seurat Inputs

These targets use tar_files so file changes invalidate downstream targets:

- debranched_seus_1q
- debranched_seus_2p
- debranched_seus_6p
- debranched_seus_16q

Purpose:
- Define per-SCNA debranched Seurat object paths.
- Serve as the canonical SCNA-stratified Seurat inputs for analysis and database population.

### 2. File-Tracked Integrated Seurat Inputs

Also tracked with tar_files:

- integrated_seus_1q
- integrated_seus_16q
- integrated_seus_2p
- integrated_seus_6p

Purpose:
- Register integrated Seurat objects used in integration-level figures and downstream comparison targets.

### 3. Per-Sample Seurat Processing Core

Core processing and preparation targets:

- filtered_large_plot_files
- unfiltered_seus
- filtered_seus
- final_seus
- cc_plots_wo_arms

Purpose:
- Build sample-level Seurat derivatives.
- Apply filtering and metadata augmentation.
- Materialize final sample set used by most downstream targets.

SCNA metadata note:
- scna is first added during unfiltered_seus and filtered_seus generation through functions documented in [R/plot_functions_3.R](R/plot_functions_3.R).

### 4. Debranched Aggregation And SCNA Collections

Aggregation and mapping targets:

- debranched_seu_files
- debranched_seus
- overall_seus
- scna_seus

Purpose:
- Establish path-level and named-vector abstractions over debranched Seurat objects.
- Build the combined SCNA collection consumed by sqlite extraction and many plotting targets.

### 5. SQLite Metadata Layer

Database and resolution metadata targets:

- seu_metadata_db
- resolution_dictionary
- chosen_resolution_seus

Purpose:
- Persist Seurat metadata summaries to sqlite.
- Retrieve and apply resolution mappings for phase cluster assignment.

Database write behavior:
- seu_metadata_db and resolution_dictionary run with deployment = main to avoid worker-level sqlite write contention.

### 6. Recurrent Marker Heatmap Track

Marker recurrence and heatmap targets:

- filtered_recurrent_genes
- filtered_recurrent_heatmaps
- filtered_recurrent_heatmap_file
- combined_recurrent_filtered_heatmap
- regressed_recurrent_genes
- regressed_recurrent_heatmaps
- regressed_recurrent_heatmap_file

Purpose:
- Produce recurrent marker outputs under filtered and regressed workflows.
- Emit merged PDF artifacts for figure assembly.

### 7. Regression Diagnostics Track

Regression-focused targets:

- regressed_seus
- effect_of_filtering
- regression_effect_plots
- fig_03_05
- regression_ora_plots
- cin_score_plots

Purpose:
- Quantify effects of filtering and regression.
- Generate diagnostics and enrichment views.

### 8. QC And Numbat Heatmap Track

QC/heatmap targets:

- fig_s03a_low_hypoxia_plots
- fig_s03a
- fig_s13
- filtered_numbat_heatmaps_file
- filtered_large_scna_prob_file

Purpose:
- Generate unfiltered and filtered Numbat heatmap outputs.
- Build merged artifacts used in supplemental reporting.

### 9. Heatmap Collage Track

Collage targets:

- heatmap_collages
- heatmap_collages_6p
- annotated_heatmap_collages
- collage_compilation
- collage_compilation_all_resolutions

Purpose:
- Produce per-sample and merged heatmap collage outputs across resolutions and SCNA strata.

### 10. Hypoxia Track

Hypoxia-oriented targets:

- hypoxia_seus
- hypoxia_score_plots
- seus_low_hypoxia
- seus_high_hypoxia
- heatmap_collages_hypoxia
- hypoxia_effect_plots
- silhouette_plots

Purpose:
- Score hypoxia, subset by hypoxia state, and visualize hypoxia-linked structure.

### 11. Integrated Seurat Objects And Integrated Collages

Integrated object/collage targets:

- integrated_seu_16q
- integrated_seus
- annotated_integrated_heatmap_collages

Purpose:
- Provide integrated objects used by integration-level visual outputs.

### 12. Static Branching For Clustrees

Static branching block:

- clustree_map via tar_map over debranched_map_values
- clustrees via tar_combine
- clustree_compilation

Purpose:
- Build per-sample clustree outputs as static branches.
- Merge branch outputs into a single clustree compilation artifact.

### 13. SCNA-Mapped Hypoxia Branching

Branching targets via tar_map over SCNA labels:

- hypoxia_seus_1q
- hypoxia_seus_2p
- hypoxia_seus_6p
- hypoxia_seus_16q

Purpose:
- Derive SCNA-specific low-hypoxia Seurat subsets from the global low-hypoxia set.

### 14. SCNA-Mapped Low-Hypoxia Integration

Branching targets via tar_map over SCNA labels:

- integrated_seu_low_hypoxia_1q
- integrated_seu_low_hypoxia_2p
- integrated_seu_low_hypoxia_6p
- integrated_seu_low_hypoxia_16q

Plus summary targets:

- hypoxia_low_integrated_heatmap_collages
- hypoxia_low_integrated_heatmap_collages0

Purpose:
- Perform and visualize SCNA-specific integration constrained to low-hypoxia cells.

### 15. Maintenance Notes

- The file mixes tar_files, tar_target, tar_map, and tar_combine patterns; this is expected and intentional.
- scna-sensitive outputs depend on the correctness of metadata augmentation in the preprocessing stage.
- sqlite outputs reflect tracked file dependencies; manual edits outside tracked paths will not invalidate targets unless upstream files/commands change.
