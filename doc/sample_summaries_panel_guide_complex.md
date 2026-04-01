# Sample Summary Panels

Each `results/{SRR}_summary.pdf` is produced by `collate_sample_summary` in
[R/plot_functions_52.R](../R/plot_functions_52.R), called as a branched target
mapped over `filtered_clone_tree_files` in
[R/pipeline_targets_figures.R](../R/pipeline_targets_figures.R).

---

## Layout

```
┌─────────────────────────────────────────────────────┐
│  Clone tree (1..n)  |  Segment tree  (Row 1)        │
├──────────────┬──────────────────────────────────────┤
│              │  Unfiltered heatmap  (Row 2)          │
│  Karyogram   ├──────────────────────────────────────┤
│              │  Filtered heatmap    (Row 3)          │
├──────────────┴──────────────────────────────────────┤
│  Unfiltered expression  |  Filtered expression       │
├─────────────────────────────────────────────────────┤
│  Unfiltered bulk clones |  Filtered bulk clones      │
└─────────────────────────────────────────────────────┘
```

---

## Row 1 — Clone trees and segment tree

**Targets:** `filtered_clone_tree_files`, `clone_trees_segments_files`

One panel is shown per tree file associated with the sample. Both target types
call `save_clone_tree_from_path` on the filtered Seurat object and the
corresponding Numbat RDS.

- **Clone tree** — phylogenetic tree inferred by Numbat with SCNA clone labels
  simplified by `large_clone_simplifications`. Multiple panels appear when the
  sample has branched subclones. Source: `filtered_clone_tree_files`.
- **Segment tree** — same tree topology but annotated with raw Numbat segment
  labels (no SCNA simplification, `large_clone_simplifications = NULL`).
  Source: `clone_trees_segments_files`.

---

## Rows 2–3 — Karyogram (left column) and Numbat heatmaps (right column)

The karyogram spans the full height of both heatmap rows. The heatmaps are
stacked in the right column.

### Karyogram

**File:** `results/{SRR}_karyogram.pdf`
**Target:** `ideogram_res_s06a` → combined into `fig_s06a`

Chromosomal ideogram showing copy-number state per cell across all chromosomes,
produced by `make_rb_scna_ideograms`. Colours indicate gain (red), loss (blue),
balanced gain, and LOH.

### Unfiltered heatmap (Row 2)

**Target:** `fig_s03a_unfiltered_plots`

Produced by `make_numbat_heatmaps(unfiltered_seus, ...)` with
`extension = "_unfiltered"`. Shows Numbat SCNA probability heatmap across cells
for the **unfiltered** Seurat object (`p_min = 0.9`). X-axis labels show
segment names.

### Filtered heatmap (Row 3)

**Target:** `fig_s03a_plots`

Produced by `make_numbat_heatmaps(filtered_seus, ...)` with
`extension = "_filtered"`. Same heatmap type as above but for the **filtered**
Seurat object — cells failing QC or non-tumour cells have been removed.

> If neither target has output for the sample, the function falls back to raw
> Numbat PDFs on disk at
> `results/numbat_sridhar/{SRR}/{SRR}_unfiltered.pdf` /
> `{SRR}_filtered.pdf`.

---

## Row 4 — Smoothed expression heatmap (unfiltered | filtered)

**Targets:** `large_numbat_expression` (left), `filtered_numbat_expression` (right)

Both retrieve `exp_roll_clust.pdf` from Numbat output via
`retrieve_numbat_plot_type`. This is the rolling-window smoothed gene expression
heatmap with cells ordered by Numbat clone assignment.

- **Unfiltered** — from `results/numbat_sridhar/{SRR}/exp_roll_clust.pdf`
  (large unfiltered Numbat run).
- **Filtered** — from `results/{SRR}/exp_roll_clust.pdf` (filtered Numbat run).
  Only present for samples with a filtered RDS in
  `output/numbat_sridhar_filtered/`.

---

## Row 5 — Bulk clone plots (unfiltered | filtered)

**Targets:** `large_numbat_bulk_clones` (left), `filtered_numbat_bulk_clones` (right)

Both retrieve `bulk_clones_final.pdf` from Numbat output via
`retrieve_numbat_plot_type`. Shows allele frequency and read depth across
chromosomes for each inferred clone — the diagnostic plot Numbat uses to
evaluate clone assignments.

- **Unfiltered** — from `results/numbat_sridhar/{SRR}/bulk_clones_final.pdf`.
- **Filtered** — from `results/{SRR}/bulk_clones_final.pdf`. Only present for
  samples with a filtered Numbat run.
