# Unmade Targets — Priority List

Generated from `_targets.R` analysis (2026-03-05).

---

## Priority 1 — Fix Broken Function Reference

### `fig_s11` / `fig_s23` — Unknown function placeholder
- **Issue:** Both targets call `not_sure_what_this_does(integrated_seus_2p/6p, cluster_orders, "results/fig_s11.pdf")`
- **Action:** Identify the correct function name and replace the placeholder
- **Blocks:** `fig_s11`, `fig_s23`, and any downstream figure compilation

---

## Priority 2 — Hypoxia Stratification Targets

These targets are commented out and need input target corrections before enabling.

### `heatmap_collages_hypoxia_low`
- **Issue:** Currently an identical copy of `heatmap_collages_hypoxia`; needs to use `seus_low_hypoxia`
- **Fix:** Update `pattern = map(seus_low_hypoxia)` and `iteration = "list"`
- **Location:** ~lines 830–834

### `heatmap_collages_hypoxia_high`
- **Issue:** Currently an identical copy of `heatmap_collages_hypoxia`; needs to use `seus_high_hypoxia`
- **Fix:** Update `pattern = map(seus_high_hypoxia)` and `iteration = "list"`
- **Location:** ~lines 835–839

---

## Priority 3 — 1q Collage / Integrated Hypoxia

These targets depend on intermediate targets that are not yet defined.

### `collages_1q`
- **Issue:** References undefined targets `hypoxia_1q`, `hypoxia_1q_low`, `hypoxia_1q_high`
- **Fix:** Define those intermediates (likely derived from `debranched_seus_1q` or via `tar_map`); note that low/high labels were previously inverted
- **Location:** ~lines 925–933

### `integrated_seu_1q_low_hypoxia_heatmaps`
- **Issue:** References undefined target `integrated_seu_1q_low_hypoxia`
- **Fix:** Define `integrated_seu_1q_low_hypoxia` analogous to `integrated_seu_low_hypoxia`; then call `plot_seu_marker_heatmap_by_scna()`
- **Location:** ~lines 1060–1062

---

## Priority 4 — Enrichment Plot Targets

These targets are referenced in the `figures_and_tables` aggregation list but commented out.

### `oncoprint_enrich_clones_plots_gobp`
- **Issue:** Listed in `figures_and_tables` but commented out; active hallmark equivalent exists
- **Action:** Implement GO-BP enrichment clone plots or remove from `figures_and_tables`

### `oncoprint_enrich_clusters_plots_gobp`
- **Issue:** Same as above for cluster-level GO-BP enrichment
- **Action:** Implement or remove from `figures_and_tables`

---

## ~~Priority 5 — Malformed / Incomplete Target~~ (resolved)

### ~~`large_cluster_markers`~~ — removed
- **Resolution:** Was a duplicate of the active `cluster_markers` target (line 291) — same function, same arguments. Commented-out block removed from `_targets.R`.

---

## Background: Commented-out Data Entries

The following are not targets themselves but affect pipeline completeness:

- **`corresponding_seus`** (~lines 1209–1210): Several `.rds` file paths are commented out (e.g., `SRR13884246_branch_5_filtered_seu_2p.rds`). Verify if the omitted samples are intentionally excluded.
- **`corresponding_states_dictionary`** (~lines 1224–1238): Multiple state comparison rows are commented out. Review whether these comparisons are needed for downstream figures.
