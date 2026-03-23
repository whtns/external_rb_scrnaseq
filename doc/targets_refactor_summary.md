# `_targets.R` Refactor: Splitting into Thematic Pipeline Files

**Date:** 2026-03-19

## Overview

The monolithic `_targets.R` (≈1940 lines) was split into five thematic pipeline files plus a constants file, leaving `_targets.R` as a slim 23-line orchestrator. The goal was to make target definitions easy to browse and maintain by grouping them by concern.

## File Structure After Refactor

```
_targets.R                          # slim orchestrator (23 lines + preserved comments)
R/
  pipeline_constants.R              # shared constants (scna_map_values, output_plot_extensions)
  pipeline_targets_inputs.R         # defines pipeline_targets_inputs
  pipeline_targets_seurat.R         # defines pipeline_targets_seurat
  pipeline_targets_diffex.R         # defines pipeline_targets_diffex
  pipeline_targets_integration.R    # defines pipeline_targets_integration
  pipeline_targets_figures.R        # defines pipeline_targets_figures
```

## The Orchestrator `_targets.R`

```r
suppressPackageStartupMessages(source("./packages.R"))
lapply(list.files("./R", full.names = TRUE), source)

tar_option_set(
  memory = "transient",
  garbage_collection = TRUE,
  error = "continue",
  controller = crew_controller_local(workers = 4)
)

tar_plan(
  !!!pipeline_targets_inputs,
  !!!pipeline_targets_seurat,
  !!!pipeline_targets_diffex,
  !!!pipeline_targets_integration,
  !!!pipeline_targets_figures
)
```

`lapply(list.files("./R", full.names = TRUE), source)` auto-sources all R files including the new pipeline files. The `!!!` (rlang splice operator) inside `tar_plan()` unpacks each list into individual target arguments.

## Pipeline File Contents

### `R/pipeline_targets_inputs.R`
File-tracking targets, sample definitions (`interesting_samples`, `debranched_ids`, `rb_scna_samples`), config loading, cluster dictionaries, oncoprint settings files, QC/metadata, gene lists (`celltype_markers`, `interesting_genes`, `subtype_markers`, `liu_lu_supp_data`, `mps`, `stemness_markers`).

### `R/pipeline_targets_seurat.R`
Seurat object file targets (`debranched_seus_1q/2p/6p/16q`, `integrated_seus_1q/2p/6p/16q`), Seurat processing pipeline, clustrees, hypoxia targets via `tarchetypes::tar_map()` over `scna_map_values`, integration targets, `integrated_seus`, heatmap collages, `resolution_dictionary`, `chosen_resolution_seus`.

### `R/pipeline_targets_diffex.R`
All/cis/trans differential expression + volcanos + tables, clustree diffex (`clustree_tables`, `clustree_diffexes`, `clustree_cis/trans/all_changes`), `num_diffex_clone_scna_tally`, fig_07a/b (1q+ integration cluster diffex), fig_08a/b (16q- integration cluster diffex).

### `R/pipeline_targets_integration.R`
SCNA collages (`collages_2p/6p/16q`, `fig_s04_09`, `collected_scna_collages`), corresponding-state seu paths and dictionaries, corresponding-state 2p and 6p diffex/volcanos/enrichments, clone CC plots (`clone_cc_plots_by_scna_1q/16q`), phase diffex (2p/6p G1), oncoprint targets (`unfiltered_oncoprint_input_by_scna`, `oncoprint_input_by_scna`, per-cluster variants, by-region), enrichment over oncoprints (GO-BP + Hallmark), rod/celltype scoring, pseudobulk subtype scores, fig_09, fig_10.

### `R/pipeline_targets_figures.R`
`figures_and_tables` aggregate target (large named list referencing all main outputs), study metadata and QC figures (`fig_01`, `fig_s02/s03/s04`, `table_s01/s02/s03/s07/s09/s10`), marker plots (violins, heatmaps, featureplots), `tarchetypes::tar_map()` for per-SCNA cluster violin plots and factor heatmaps, clone distribution plots and pearls, large Numbat plots, clone trees, integration figures (`fig_03/04/05`, `fig_04_v2/07/09`, `fig_s07/s08/s10/s12/s20/s23/s25`), `pipeline_notification`.

## `R/pipeline_constants.R`

Defines shared constants used across pipeline files:

```r
output_plot_extensions <- c("pdf", "png")

scna_map_values <- tibble::tibble(
  scna     = c("1q", "2p", "6p", "16q"),
  seus_sym = rlang::syms(c("debranched_seus_1q", "debranched_seus_2p",
                            "debranched_seus_6p", "debranched_seus_16q")),
  var_y    = c("phase_level", "clusters", "clusters", "phase_level")
)
```

`scna_map_values` drives `tarchetypes::tar_map()` calls in `pipeline_targets_seurat.R` and `pipeline_targets_figures.R`. The `scna` column becomes the suffix appended to template target names (e.g., template name `hypoxia_seus` inside `tar_map(names="scna")` produces `hypoxia_seus_1q`, `hypoxia_seus_2p`, etc.).

## Key Technical Notes

### `tar_map()` template names are not real target names
Inside `tarchetypes::tar_map(names = "scna", tar_target(foo, ...))`, the name `foo` is a template — only the suffixed variants `foo_1q`, `foo_2p`, `foo_6p`, `foo_16q` are registered as targets. Apparent grep hits for `tar_target(foo` do not indicate a duplicate.

### Combining `tar_files()` blocks with `list()` and `tar_map()`
The pattern used in `pipeline_targets_seurat.R`:
```r
pipeline_targets_seurat <- c(
  tar_files_block_1,   # tarchetypes::tar_files() returns a list
  list(tar_target(a, ...), tar_target(b, ...)),  # wrap singles in list()
  tar_map(...)         # also returns a list
)
```
`tar_plan()` uses `unlist(recursive = FALSE)` internally to flatten one level of nesting.

### Duplicate target handling in the original `_targets.R`
The original file had several duplicate target definitions (`rb_scna_samples`, `oncoprint_settings*`, `subtype_markers`, `divergent_cluster_file`, `clustree_tables`). In the refactored version, each target appears in exactly one pipeline file.

### Dependency name fix
The original `_targets.R` had `tally_num_diffex(oncoprint_input_by_scna_unfiltered)` — a nonexistent target name. The correct name used everywhere else is `unfiltered_oncoprint_input_by_scna`. Fixed in `pipeline_targets_diffex.R`.

## Testing the Refactored Pipeline

### Step 1 — Parse check (fast, catches syntax and name errors)
```r
targets::tar_manifest(store = "_targets_r431")
```
Returns a data frame of all targets if the pipeline loads without errors. If it errors, the message identifies which target has a problem.

### Step 2 — Dependency graph (catches missing upstream dependencies)
```r
targets::tar_visnetwork(store = "_targets_r431")
```
Opens an interactive dependency graph in the viewer. Targets with missing dependencies will be visually disconnected or highlighted.

### Step 3 — Check for outdated targets
```r
targets::tar_outdated(store = "_targets_r431")
```
Lists targets that need to be rebuilt. Since the refactor only reorganizes definitions without changing code, all previously-built targets should remain up to date.

### Step 4 — Smoke test a lightweight target
```r
targets::tar_make(
  names = "interesting_samples",
  store = "_targets_r431"
)
```
Runs a single cheap target end-to-end to confirm the pipeline executes correctly.

### Duplicate target name check
```bash
grep -h "tar_target(" R/pipeline_targets_*.R \
  | grep -oP "tar_target\(\K[a-zA-Z_0-9]+" \
  | sort | uniq -d
```
Any output indicates a target name defined in more than one pipeline file.
