# Cell Filtering Sequence

Filtering is applied in `filter_cluster_save_seu()` (`R/plot_functions_3.R:20`) on the
`filtered_seus` target.  Each step operates on the cells surviving the previous step.

## Steps

### 1. Numbat clone confidence (NA removal)
Cells without a confident clone assignment (`clone_opt == NA`) from the numbat
`clone_post` table are dropped immediately after numbat metadata is joined.

- **Criterion**: `!is.na(seu$clone_opt)`
- **Source**: `R/plot_functions_3.R:38`

### 2. Bulk QC metrics
Hard thresholds on three standard QC fields computed from the raw count matrix.

| Metric | Threshold |
|---|---|
| `percent.mt` | < 10 % |
| `nCount_gene` | > 1 000 |
| `nFeature_gene` | > 1 000 |

- **Function**: `filter_sample_qc()`
- **Source**: `R/plot_functions_3.R:90–92`

### 3. Cluster-level annotation removal
Clusters flagged `remove == "1"` in the per-sample `cluster_dictionary` (curated
from marker-gene review) are dropped by `gene_snn_res.0.2` cluster ID.

- **Source**: `R/plot_functions_3.R:94–97`, `R/plot_functions_3.R:103`

### 4. MALAT1 contamination cluster removal
If any cluster is annotated `"MALAT1"` (a well-known ambient RNA / low-quality cell
marker), all cells in that cluster are additionally removed.

- **Source**: `R/plot_functions_3.R:104–106`

### 5. Manual cell-barcode exclusion (`cells_to_remove_final.xlsx`)
Individual cell barcodes listed in `data/cells_to_remove_final.xlsx` are removed.
The file is an Excel workbook with one sheet per sample, named by SRR ID; each sheet
contains at least a `cell` column of barcodes.

**Origin**: the file represents manual curation — cells that passed all automated
filters but were flagged during interactive inspection (e.g. marker expression
review, visual outlier detection in UMAP).  It is read at pipeline-definition time
into the `cells_to_remove` target (`R/pipeline_targets_inputs.R:200`) via
`read_cells_to_remove()` (`R/plot_functions_28.R:1`).

- **Source**: `R/plot_functions_3.R:98–109`

## Downstream partitioning (post-filtering)

After `filtered_seus` → `final_seus`, cells are further partitioned (not dropped) by
hypoxia score for parallel analysis tracks:

| Target | Criterion | Threshold |
|---|---|---|
| `seus_low_hypoxia` | `hypoxia_score <= hypoxia_threshold` | 0.5 (set in `R/pipeline_constants.R`) |
| `seus_high_hypoxia` | `hypoxia_score > hypoxia_threshold` | 0.5 |

These subsets are created by `subset_seu_by_expression()`
(`R/plot_functions_19.R:572`) and are used for SCNA-stratified analyses but do not
feed back into the main `final_seus` track.

## Summary of cell loss

A `plot_filtering_timeline()` call at the end of `filter_cluster_save_seu()` records
cell counts at each stage (pre-QC, post-SCNA, post-QC, post-cluster) and writes a
PDF to `results/{sample_id}_filtering_timeline_new.pdf`.
