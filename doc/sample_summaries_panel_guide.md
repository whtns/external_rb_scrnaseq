# Sample Summary Panels

Each per-sample summary PDF is assembled from multiple upstream analysis outputs
and displays the following panels.

---

## Layout

```
┌──────────────────────────────────────────┬─────────────┐
│  Clone tree (1..n)  |  Segment tree      │  Karyogram  │  Row 1
├──────────────────────────────┬───────────┴─────────────┤
│  Unfiltered heatmap          │  Filtered heatmap        │
├──────────────────────────────┼──────────────────────────┤
│  Unfiltered expression       │  Filtered expression     │
├──────────────────────────────┼──────────────────────────┤
│  Unfiltered bulk clones      │  Filtered bulk clones    │
└──────────────────────────────┴──────────────────────────┘
         Left column (unfiltered)   Right column (filtered)
```

Row 1 spans full width. All remaining content is arranged in two persistent
columns: **unfiltered** on the left and **filtered** on the right. Within each
column, panels are stacked vertically in the order: heatmap → expression →
bulk clones.

---

## Row 1 — Clone trees, segment tree, and karyogram

One panel per tree file associated with the sample.

- **Clone tree** — Phylogenetic tree inferred by Numbat showing the
  evolutionary relationships between tumour subclones. SCNA clone labels are
  simplified to major arm-level events. Multiple panels appear when the sample
  has branched subclones.

- **Segment tree** — Same tree topology as the clone tree but annotated with
  raw Numbat segment labels without any simplification applied.

- **Karyogram** — Chromosomal ideogram rotated 90° showing copy-number state
  per cell across all chromosomes. Colours indicate gain (red), loss (blue),
  balanced gain, and LOH. Placed to the right of the trees. Produced by the
  ideogram pipeline step combining Numbat posterior estimates across all
  chromosome arms.

---

## Left column — Unfiltered panels

Derived from the unfiltered Seurat object and the corresponding unfiltered
Numbat run. All cells passing initial QC are included; tumour/normal separation
has not yet been applied.

### Unfiltered heatmap

Numbat SCNA posterior probability heatmap across cells. Each column is a cell;
rows show SCNA segments. Segment names are labelled on the x-axis.

### Unfiltered expression

Rolling-window smoothed gene expression heatmap with cells ordered by Numbat
clone assignment. Rows are genes; columns are cells grouped by inferred clone.

### Unfiltered bulk clones

Numbat diagnostic plot showing allele frequency and read depth across all
chromosomes, broken down by inferred clone. Used to evaluate the quality of
clone assignments.

---

## Right column — Filtered panels

Derived from the filtered Seurat object and the corresponding filtered Numbat
run. Non-tumour cells and low-quality cells have been removed. Differences
between the left and right columns reflect the effect of cell filtering on the
inferred copy-number landscape.

These panels are absent for samples that did not undergo a separate filtered
Numbat run.

### Filtered heatmap

Same as the unfiltered heatmap but restricted to the filtered cell set.

### Filtered expression

Same as the unfiltered expression heatmap but restricted to the filtered cell
set.

### Filtered bulk clones

Same as the unfiltered bulk clone plot but restricted to the filtered cell set.

---

> **Note:** if neither filtered nor unfiltered heatmap was produced for a
> sample by the pipeline, the function falls back to raw Numbat output files
> stored on disk.
