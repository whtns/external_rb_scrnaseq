# The consensus diploid object

**Issues:** [#44](https://github.com/whtns/external_rb_scrnaseq/issues/44)
(build a target), [#15](https://github.com/whtns/external_rb_scrnaseq/issues/15)
(define candidates), [#16](https://github.com/whtns/external_rb_scrnaseq/issues/16)
(exclude 484)
**Targets:** `diploid_panel_samples_file` → `diploid_panel_sample_ids` →
`diploid_seu` → `diploid_seu_composition`
**Sample table:** `data/diploid_panel_samples.csv`

## What was wrong

Three separate problems, all of which had to be fixed before the object could be
rebuilt at all.

### 1. `diploid_seu` had no build command

```r
tar_target(diploid_seu,
  "output/seurat/diploid_subsets/diploid_seu.rds",
  format = "file")
```

A bare path. Nothing in the pipeline created that file. The object every diploid
figure depended on was **1.38 GB, dated 2026-04-26**, made by hand, and
predating the majority-round numbat rebuild — so its clone assignments came from
a superseded round selection.

### 2. `ks_diploid_seu` was dead, and pointed somewhere else

It called `assemble_diploid_seu()` **without `out_path`**, so it would have
written to the function's default `output/seurat/diploid_seu.rds` — a *different*
path from the one `diploid_seu` read. That file does not exist. The target was
also referenced nowhere else in the repo.

### 3. `assemble_diploid_seu()` could never have succeeded

At `numbat_helpers/R/plot_functions_3.R:664`:

```r
cone_mask <- tolower(seu$type) == "cones"    # plural
if (!any(cone_mask)) return(NULL)
```

`seu$type` is assigned by `plot_celltype_predictions()`
(`plot_functions_43.R:32-134`), whose reference is built with
`filter(type %in% str_remove(celltypes, "s$"))` — the **singular** spellings.
Confirmed against the reference data itself: the distinct `type` values in
`data/plae_pseudobulk_counts.csv` include `'Cone'` and `'Rod'`, and **no plural
form appears anywhere in the file**.

So `cone_mask` was always `FALSE`, every sample returned `NULL`, `cone_seus` was
empty after `compact()`, and the function then errored —
`integration_workflow(list())` on the integrate path, or `seus_list[[1]]`
subscript-out-of-bounds on the other.

This is why the target was dead and the object was made by hand. Fixed to
`tolower(seu$type) %in% c("cone", "cones")`, tolerant of both rather than
swapping one bug for its mirror image.

## What replaces it

Sample selection moves out of an inline `grepl` into
`data/diploid_panel_samples.csv` — `sample_id`, `include`, `percent_rod_cells`,
`reason` — which is what #15 asks for and makes #16 a one-row edit.

**22 of 31 samples included.** Exclusions:

| sample | reason |
|---|---|
| SRX11133585, SRX11133588, SRX14116948 | rod-rich (31%, 16%, 8% rod cells) |
| SRX14116944 | sample 484 (SRR17960484) — **#16** |
| SRX10031194, SRX10264517, SRX10264518, SRX10264523, SRX14116946 | carried over from the previous hard-coded list; **original reason not recorded — flagged REVIEW** |

Two things need your judgement:

- The five carried-over exclusions have **no recorded rationale** anywhere in the
  code or git history. They are preserved so behaviour does not silently change,
  but each is marked `REVIEW` in the table.
- `SRX11133592/93/94` are currently **included**. AGENTS.md holds them out of
  numbat rebuilds, but that restriction is about writes to
  `output/numbat_sridhar/`, and building a diploid Seurat object does not write
  there. Also marked `REVIEW` rather than silently decided either way.

### Resolving #16

484 = SRR17960484 = **SRX14116944** (`data/combined_metadata.tsv:10`). Note it
was *not* in the old hard-coded exclusion list, and because the on-disk object
was hand-made there was no way to tell whether it was in there — that is
precisely the reproducibility problem #44 describes.

`diploid_seu_composition` (`results/diploid_seu_composition.csv`) now reports
samples and cell counts, so the question is answerable by reading a CSV instead
of loading 1.4 GB.

## What is unchanged

The body of `assemble_diploid_seu()` beyond the mask fix: diploid cells are
`is.na(seu$scna) | seu$scna == ""`, then cell-type prediction, then cones only,
then `seuratTools::integration_workflow(resolution = c(0.2, 0.4))`. It writes via
`add_hash_metadata()`, which also upserts `(filepath, hash, n_cells)` into
`batch_hashes.sqlite`.

Consumers keep the same interface (`diploid_seu` is still a file path):
`low_hypoxia_diploid_merged_seu`, the high-hypoxia equivalent,
`plot_diploid_seu_umaps()`, `fig_single_sample_panels_with_diploid`.

## Expect the object to change

This will be the **first programmatic build**. The April object's criteria are
unknown and unreproducible, so the new object will differ and every downstream
diploid figure will move. Read `results/diploid_seu_composition.csv` before
trusting the new figures, and keep the old file until you have.
