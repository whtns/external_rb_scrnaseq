# Continuous hypoxia / mitochondrial column annotations on the filtered SCNA collages

Adds two continuous per-cell scores as heatmap **column annotations**, in their own
colours, to both filtered collage families:

- `two_clone_scna_collage_filtered` (`two_clone_res_collages_filtered_{1q,2p,16q}`)
- `all_clone_scna_collage_filtered` (`all_clone_res_collages_filtered_{1q,2p,16q}`)

The low-hypoxia twins (`two_clone_scna_collage`) are deliberately unchanged — they
keep the cell-cycle-only annotation.

## What is drawn

| annotation      | column          | colour ramp        | meaning |
|-----------------|-----------------|--------------------|---------|
| `G2M.Score`     | `G2M.Score`     | white → `#F8766D`  | unchanged |
| `S.Score`       | `S.Score`       | white → `#00BFC4`  | unchanged |
| `hypoxia_score` | `hypoxia_score` | white → `#762A83` (purple) | composite score from `add_hypoxia_score()`: `rescale(rowMeans(hypoxia, MT), 0..1)` — the score the hypoxia split thresholds on |
| `mito_score`    | `mito_score`    | white → `#8C510A` (brown)  | mitochondrial module score, **natural orientation** (dark = high mitochondrial content) |

`mito_score` is `-MT`. `add_hypoxia_score()` stores `MT` already multiplied by `-1`
so that averaging it with `hypoxia` yields the composite; annotating with `MT`
directly would read backwards — palest exactly where mitochondrial content is
highest.

`G2M.Score` / `S.Score` are pinned to the two colours `scales::hue_pal()` gave them
when they were the only continuous annotations. Without the pin, adding two more
numeric annotations shifts the positional palette and silently recolours the
cell-cycle rows, so an annotated collage would no longer match its un-annotated
twin panel for panel.

## Where the scores come from

The `*_filtered_seu.rds` objects these collages are built from carry **none** of
`hypoxia`, `MT`, `hypoxia_score`, `S.Score`, `G2M.Score` — only `percent.mt`
(verified with `src/diag_hypoxia_score_cols.R`). Only the hypoxia-split objects
have been through `add_hypoxia_score()`. So the scores are computed by
`.ensure_hypoxia_scores()`, the sibling of the existing `.ensure_cc_scores()`.

**Computed once per sample, on the whole object, before any subsetting.** Two
reasons this matters:

1. `hypoxia_score` is rescaled to `[0, 1]` over the cells it is computed on.
   Scoring a two-clone subset would place its 0 and 1 at that subset's extremes,
   so no two panels — and neither the two-clone nor the all-clone view of the same
   sample — would share a scale.
2. The same object feeds both builders once per SCNA of interest.
   `add_hypoxia_score()` pulls all of MSigDB and runs `AddModuleScore()` over every
   cell; uncached, that is four or more identical recomputations per sample
   (~1 min each).

## Cache: the `cell_scores` table

Cached in `batch_hashes.sqlite`, alongside the other per-cell numerics in
`cell_qc_values`:

```sql
CREATE TABLE IF NOT EXISTS cell_scores (
  filepath      TEXT,   -- source *_filtered_seu.rds; the cache key
  cell          TEXT,
  sample_id     TEXT,
  hypoxia       REAL,   -- HALLMARK_HYPOXIA module score
  mt            REAL,   -- MT as add_hypoxia_score() stores it (already negated)
  hypoxia_score REAL,   -- rescaled composite
  recorded_at   TEXT,
  PRIMARY KEY (filepath, cell)
)
```

Accessors: `save_cell_scores_to_db()` / `read_cell_scores_from_db()`
(`R/metadata_functions_1.R`), both going through `connect_hash_db()` + `db_retry()`
per AGENTS.md — never a raw `DBI::dbConnect()`.

Two correctness guards:

- **Writes are `DELETE` + append inside one transaction**, so a reader arriving
  mid-write sees either the old row set or the new one, never a half-replaced
  mixture of two rescalings.
- **A cache is used only if it covers every cell the caller asks about.** A partial
  hit returns `NULL` and the caller rescores. Rows written from one cell set carry
  a `hypoxia_score` rescaling that does not apply to a different one, so "close
  enough" is wrong here.

## Files changed

| file | change |
|------|--------|
| `numbat_helpers/R/plot_functions_19.R` | `plot_seu_marker_heatmap(score_annotations=, score_annotation_hues=)`; NA/constant-safe `numeric_col_fun()`; `cc_groupby` now drops absent columns instead of keeping them |
| `numbat_helpers/R/plot_functions_53.R` | `seu_complex_heatmap(numeric_col_hues=)` — pin a named numeric annotation's colour instead of taking its positional `hue_pal()` slot |
| `numbat_helpers/R/scna_two_clone_res_collages.R` | `.ensure_hypoxia_scores()`; `score_annotations` passthrough |
| `numbat_helpers/R/scna_all_clone_res_collages.R` | `score_annotations` passthrough |
| `numbat_helpers/R/metadata_functions_1.R` | `save_cell_scores_to_db()` / `read_cell_scores_from_db()` |
| `numbat_helpers/R/seu_metadata_db.R` | `cell_scores` in the schema + `init_seu_metadata_db()` |
| `R/pipeline_targets_integration.R` | both filtered collage targets pass `score_annotations = c("hypoxia_score", "mito_score")` |

## Verification

- `src/diag_hypoxia_score_cols.R` — which score columns each object stage carries.
- `src/smoke_two_clone_score_annotations.R` (`smoke_two_clone_score_annotations.sbatch`,
  `debug` partition) — cold scoring + cache write, warm read round-trips
  identically, and both builders render a PDF with the annotations.

## Rebuild

`rebuild_filtered_collage_score_annotations.sbatch` invalidates only the two
filtered collage families and rebuilds them. One `tar_make()` at a time against
`_targets_r431` — check `squeue -u $USER` first.
