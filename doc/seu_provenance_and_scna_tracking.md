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
