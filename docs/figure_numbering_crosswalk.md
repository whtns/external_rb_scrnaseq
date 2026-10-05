# Figure and table numbering crosswalk

Draft: `doc/kstachelek_rb_scna_tx_effects_paper_2026-06-28.docx`. Pipeline labels:
`figure_order` in `R/pipeline_constants.R`, which only fills the `document_label`
column of the figures manifest (`R/pipeline_targets_figures.R:74`). Built 2026-10-04.

Three sources are compared for each figure:

- **Text:** the citations in the body text.
- **Caption:** the caption in the "Figures and Tables" and "Supplemental Data" sections.
- **Pipeline:** the label in `figure_order`.

Embedded images were identified by looking at them (`word/media/imageN.png`).

## Main figures

| Content | Text cites | Caption | Image | Pipeline | Problem |
|---|---|---|---|---|---|
| Concept + workflow schematic | 1A (l.43, 45), 1b (l.46, 50), **1a** (l.108, same workflow as l.46) | **none** | none | none | Fig. 1 has no caption or image; l.108 says 1a for the panel l.46 calls 1b |
| RB SCNA frequency in TCGA | "Supp Fig 'RB SCNAs in Other Cancers'" (l.45) | **"Figure S1"**, placed among the main figures | image2 | Fig. S1 | Supplementary figure sitting in the main-figure block |
| Single tumor: Numbat, states, subclone distribution | Fig. 2A–F | Figure 2 | image3 | Fig. 2 (`fig_single_sample_panels`) | OK. `fig_single_sample_panels_with_diploid` = "Fig. 2 (diploid)" |
| Single tumor: cluster markers | Fig. 3A, 3b, 3 | Figure 3 | image4 | none (comment says partly `fig_01`) | No pipeline target |
| 1q+ integration | Fig. 4, 4A–E | Figure 4 | image5 | Fig. 4 (`fig_1q_integrated`) | OK. `fig_1q_integrated_v2` = "Fig. 4v2", `fig_1q_16q_combined` = "Fig. 4-7" |
| 16q- integration | Fig. 5A, 5C–E | Figure 5 | image6 | Fig. 5 (`fig_16q_integrated`) | OK |
| Kooi candidate drivers | "Supplementary Fig. S1" (l.45, cites Kooi 2016) | "[Figure 6 omitted … removed]" | **image7 still embedded** | none | Caption says removed, image still there; l.45 calls it S1 |
| (none) | none | none | none | **Fig. 7a/7b** (1q cluster diffex), **8a/8b** (16q), **9** (2p corresponding clusters), **10** (6p corresponding clusters) | Pipeline claims main-figure numbers that don't exist in the draft |
| Copy of Fig. 2 in "Extra" | Fig. 5C, 5E (l.243–244) | **"Figure 5"** (para 386) | image29 (old Fig. 2) | none | A second "Figure 5" caption; leftover |

## Supplementary figures

| Caption no. | Content | Text cites | Pipeline | Problem |
|---|---|---|---|---|
| S1 (twice) | TCGA SCNA frequency | see above | S1 `fig_tcga_scna_frequency` | Duplicate S1 caption |
| S2 | TCGA GISTIC | none | S2 `fig_tcga_gistic` | OK |
| **S1 (again)** | QC of scRNA-seq samples | **"Figure S3"** (l.49) | **S4** `fig_study_cell_stats` | Three different numbers |
| **S2 (again)** | Expression-based SCNA profiles / hierarchical clustering (11 tumors) | none ("S5-7" in l.51?) | S3b `fig_numbat_expression_smoothed` | Duplicate S2 |
| S3 | Numbat SCNA heatmaps (+ phylogenies) | l.51 "(Fig. 2, Fig. S5-7)" | S3 `fig_numbat_heatmaps` | l.51 should be S2–S3, not S5–7 |
| S4 | Removed non-tumor clusters | S4 (l.50, 126) | none (`fig_study_cell_stats` mislabelled S4) | Text and caption agree; pipeline doesn't |
| S5 | 1q+ per sample, integrated | S5 | S5 `fig_1q_sample_specific_integrated` | OK |
| S6 | 1q+ per sample, not integrated | S6 | none; **`fig_karyograms` = "S6a"** | Karyograms have no caption anywhere |
| S7 | 1q+ alternative resolutions | S7 | none | OK |
| S8 | 1q+ clone diffex / enrichment | S8, S8A–B | S8 `fig_1q_clone_diffex_within_clusters` | OK |
| S9 | 1q+ G1 MT vs G1 DIFF/SUB2 enrichment | S9 | S9 `fig_1q_cluster_diffex_of_interest` | OK |
| S10 | 16q- candidate drivers | S10 | S10 `fig_16q_clone_diffex_within_clusters` | OK |
| S11 | 16q- alternative resolutions | S11 | none | OK |
| S12 | 16q- per sample, integrated | S12 | S12 `fig_16q_sample_specific_integrated` | OK |
| S13 | 16q- per sample, not integrated | S13 | none (`fig_numbat_heatmaps_permissive` used to be S13) | OK |
| S14 | 2p+ integration | S14 | S14 `fig_2p_integrated` | OK |
| S15 | 2p+ alternative resolutions | S15 (l.102 uses it for per-sample results) | none | l.102 also cites "**Fig 4.9E**" (thesis leftover) |
| S16 | 2p+ per sample, integrated | S16B | S16 `fig_2p_sample_specific_integrated` | OK |
| S17 | 2p+ g1 cluster diffex | S17 | none | OK |
| S18 | 2p+ candidate drivers | S18 | none (`fig_2p_clone_diffex_within_clusters` = "S-2p-clde") | Probably the same figure |
| S19 | 6p+ boundaries | none | none | |
| S20 (twice) | 6p+ per sample, not integrated | S20 | none; `fig_2p_sample_specific_unintegrated` = **"S4.9"** | Duplicate S20 caption; "S4.9" is a thesis number |
| S21 | 6p+ enriched terms | S21 | none | OK |
| S22 | 6p+ cis DE | none | none | |
| S23 | cNMF | S23 | none (mosaicMPI, outside the pipeline) | OK |
| S24 | Numbat clone phylogenies (image27) | S24 (l.87); **"Fig. 1b"** (l.56); "clone_tree_collage" (l.55, 56); "phylogenies of non-informative samples" (l.53) | **"Fig. 1b"** `clone_tree_collage` | Four names for one figure |
| S25 | Subtype violins | S25 | S25 `fig_subtype_markers` | OK |
| S26 | CENPF violins | none | none | |
| none | Merged replicate phylogenies | "Supp Fig. clone_tree_collage_of_merged_replicates" (l.51) | "Fig. 1b-rep" | Not made yet |
| none | Regression diagnostics | none | "Fig. 3-5" | Draft |

## Tables

| Caption | Text cites | Pipeline | Problem |
|---|---|---|---|
| S1 TCGA frequencies | "Supp Fig … + Table" (l.45) | none | |
| S2 datasets | S2 (l.49, 122, 133) | S2 `table_sample_metadata` | OK, but l.122 means Taylor et al.'s Table S2 |
| **X2** clone % per cluster (A 1q+, B 16q-, C 2p+) | "Table **S4b**" (l.83, 89), "**S4bB**" (l.85), "**S4B**" (l.102) | "X2", "X2 (16q-)", "X2 (2p+)" | Same table, five names; collides with S4 = QC |
| S3 rod proportions | "Table S3" (l.50, cited for "Extended Data 1" Numbat output) | S3 `table_rod_rich_samples` | l.50 citation looks wrong |
| S4 tumor QC | S4 (l.50, retina contamination) | S4 `table_qc_stats` | OK |
| S5 SRA metadata | none | none | |
| S6 Numbat segments | none | none | |
| S7 cis/trans DE | none | none | |
| S8 removed clusters | S8 (l.105, but for candidate drivers) | S8 `table_removed_clusters` | l.105 should be S9 |
| S9 candidate drivers | none | none | |
| none | "Table #" (l.50, % cells excluded) | none | Placeholder; probably S8 |

## Applied 2026-10-04

Decisions: renumber supplementary figures by first citation; keep the Kooi figure as a
supplementary figure; Fig. 1 is still to be made; target names stay unchanged.

**Review docx** (`doc/..._review_t1e5.docx`, author "Claude (review)"): 62 tracked edits
to citations and captions, plus 8 "[figure numbering]" comments. Rejecting this pass
gives back the previous text exactly. Validation against the untouched original passes,
including the tracked-change author check.

**Pipeline** (`R/pipeline_constants.R`): `figure_order` label values updated. Names are
unchanged, so only the figures-manifest target goes out of date.

| New | Old caption | Content | Pipeline target |
|---|---|---|---|
| S1 | S1 | TCGA SCNA frequency | `fig_tcga_scna_frequency` |
| S2 | S2 | TCGA GISTIC | `fig_tcga_gistic` |
| S3 | "Fig. 6 omitted" | Kooi candidate drivers | none |
| S4 | S1 (dup) | QC | `fig_study_cell_stats` |
| S5 | S4 | Removed non-tumor clusters | none |
| S6 | S3 | Numbat heatmaps (+ karyograms) | `fig_numbat_heatmaps`, `fig_karyograms` |
| S7 | S2 (dup) | Expression-based SCNA profiles | `fig_numbat_expression_smoothed` |
| S8 | none | Merged-replicate phylogenies (not made) | `clone_tree_collage_of_merged_replicates` |
| S9 | S24 | Clone phylogeny collage | `clone_tree_collage` |
| S10 | S7 | 1q+ alternative resolutions | none |
| S11 | S11 | 16q- alternative resolutions | none |
| S12 | S15 | 2p+ alternative resolutions | none |
| S13 | S10 | 16q- candidate drivers | `fig_16q_clone_diffex_within_clusters` |
| S14 | S12 | 16q- per sample, integrated | `fig_16q_sample_specific_integrated` |
| S15 | S13 | 16q- per sample, not integrated | none |
| S16 | S5 | 1q+ per sample, integrated | `fig_1q_sample_specific_integrated` |
| S17 | S6 | 1q+ per sample, not integrated | none |
| S18 | S25 | Subtype violins | `fig_subtype_markers` |
| S19 | S8 | 1q+ clone diffex | `fig_1q_clone_diffex_within_clusters` |
| S20 | S9 | 1q+ G1-state enrichment | `fig_1q_cluster_diffex_of_interest` |
| S21 | S14 | 2p+ integration | `fig_2p_integrated` |
| S22 | S16 | 2p+ per sample, integrated (also replaces thesis "Fig 4.9E") | `fig_2p_sample_specific_integrated` |
| S23 | S17 | 2p+ g1 cluster diffex | none |
| S24 | S18 | 2p+ candidate drivers | `fig_2p_clone_diffex_within_clusters` |
| S25 | S23 | cNMF | none (mosaicMPI) |
| S26 | S19 | 6p+ boundaries (uncited) | none |
| S27 | S20 (x2) | 6p+ per sample, not integrated | none |
| S28 | S21 | 6p+ enriched terms | none |
| S29 | S22 | 6p+ cis DE (uncited) | none |
| S30 | S26 | CENPF violins (uncited) | none |

Tables: only Table X2 was renumbered, to **S10** (A 1q+, B 16q-, C 2p+), replacing
"S4b", "S4bB" and "S4B". "Table S8" for candidate drivers became S9. "Table #" and the
"Table S3" citation for Numbat output have comments.

Still open (comments in the docx):

- Fig. 1 panels (1A / 1a / 1b).
- Captions are still in the old order; re-sort them.
- Move the TCGA and Kooi images into Supplemental Data.
- Duplicate S1 and S27 captions.
- The "Extra" section's old Fig. 2 copy captioned "Figure 5".
- Pipeline Fig. 7–10 relabelled "Draft (was …)".
