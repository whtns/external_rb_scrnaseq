# Pipeline-level constants used across pipeline_targets_*.R files.
# Defined outside tar_plan() so they are available at pipeline-definition time.

output_plot_extensions <- c(
  "dimplot.pdf",
  "merged_marker.pdf",
  "combined_marker.pdf",
  "phylo_probability.pdf",
  "clone_distribution.pdf"
)

# Values tibble for tarchetypes::tar_map() over the 4 SCNA types.
# Columns:
#   scna      - SCNA label used in target names and as list index key
#   seus_sym  - symbol of the corresponding low-hypoxia Seurat path target
#   var_y     - y-axis grouping for clone_pearls plots
scna_map_values <- tibble::tibble(
  scna     = c("1q",         "2p",       "6p",       "16q"),
  seus_sym = rlang::syms(c(
    "hypoxia_seus_1q", "hypoxia_seus_2p",
    "hypoxia_seus_6p", "hypoxia_seus_16q"
  )),
  var_y    = c("phase_level", "clusters", "clusters", "phase_level")
)

# rclone destination for figures_and_tables sync.
# Configure with: module load rclone && rclone config
gdrive_destination <- "gdrive:rb_scrnaseq/figures_and_tables"

# rclone destination for the rebuilt hypoxia downstream deliverables
# (per-sample hypoxia summaries, low-hypoxia clone trees, numbat/hypoxia-gene
# heatmaps, annotated collages). See hypoxia_rebuilt_gdrive target.
hypoxia_gdrive_destination <- "gdrive:rb_scrnaseq/hypoxia_rebuilt"

# Document-order mapping: semantic target name → display label.
# Edit this vector to renumber figures without touching target definitions.
figure_order <- c(
  # Main figures (per figure_and_table_captions.txt)
  fig_single_sample_panels               = "Fig. 2",
  # Fig. 3 = cluster marker gene analysis in single pilot tumor — covered by fig_01 (partial)
  fig_1q_integrated                      = "Fig. 4",
  fig_16q_integrated                     = "Fig. 5",
  fig_1q_cluster_diffex_integrated       = "Fig. 7a",
  fig_1q_cluster_diffex_unintegrated     = "Fig. 7b",
  fig_16q_cluster_diffex_integrated      = "Fig. 8a",
  fig_16q_cluster_diffex_unintegrated    = "Fig. 8b",
  fig_2p_corresponding_clusters          = "Fig. 9",
  fig_6p_corresponding_clusters          = "Fig. 10",
  # Supplemental figures (document order per figure_and_table_captions.txt)
  fig_tcga_scna_frequency                = "Fig. S1",
  fig_tcga_gistic                        = "Fig. S2",
  fig_numbat_heatmaps                    = "Fig. S3",
  fig_numbat_expression_smoothed         = "Fig. S3b",      # sub-panel of S3; no direct manuscript cite
  fig_study_cell_stats                   = "Fig. S4",
  fig_1q_sample_specific_integrated      = "Fig. S5",
  fig_karyograms                         = "Fig. S6a",
  # Fig. S6  = 1q+ sample-specific without integration — no pipeline target
  # Fig. S7  = alt Louvain resolutions for integrated 1q+ — no pipeline target
  fig_1q_clone_diffex_within_clusters    = "Fig. S8",
  fig_1q_cluster_diffex_of_interest      = "Fig. S9",
  fig_16q_clone_diffex_within_clusters   = "Fig. S10",
  # Fig. S11 = alt resolutions for integrated 16q- — no pipeline target
  fig_16q_sample_specific_integrated     = "Fig. S12",
  # Fig. S13 = 16q- sample-specific without integration — no pipeline target
  fig_2p_integrated                      = "Fig. S14",
  # Fig. S15 = alt resolutions for integrated 2p+ — no pipeline target
  fig_2p_sample_specific_integrated      = "Fig. S16",
  # Fig. S17 = 2p+ cluster DE — no pipeline target
  # Fig. S18 = 2p+ candidate drivers — no pipeline target
  # Fig. S19 = 6p+ SCNA boundaries — no pipeline target
  # Fig. S20 = 6p+ sample-specific without integration — no pipeline target
  # Fig. S21 = 6p+ enriched terms — no pipeline target
  # Fig. S22 = 6p+ DE in cis — no pipeline target
  fig_subtype_markers                    = "Fig. S25",
  # Draft / unassigned pipeline figures (document position not yet confirmed)
  fig_6p_integrated                      = "Fig. 6p-draft",   # was Fig. 5; captions Fig. 5 = 16q-
  fig_numbat_heatmaps_permissive         = "Fig. S-nbt-perm", # was Fig. S13; S13 = 16q- unintegrated
  fig_2p_clone_diffex_within_clusters    = "Fig. S-2p-clde",  # was Fig. S20; S20 = 6p+ unintegrated
  fig_6p_sample_specific_integrated      = "Fig. S-6p-si",    # no confirmed position in captions
  clone_tree_collage                     = "Fig. 1b",   # static PNG; also cited as Fig. S24
  clone_tree_collage_of_merged_replicates = "Fig. 1b-rep", # static PNG (not yet created)
  fig_2p_sample_specific_unintegrated    = "Fig. S4.9",
  fig_1q_integrated_v2                   = "Fig. 4v2",
  fig_1q_16q_combined                    = "Fig. 4-7",
  fig_regression_diagnostics             = "Fig. 3-5",
  fig_single_sample_panels_with_diploid  = "Fig. 2 (diploid)",
  # Tables (per figure_and_table_captions.txt)
  table_sample_metadata                  = "Table S2",
  table_rod_rich_samples                 = "Table S3",          # S3 = rod cell proportions
  table_qc_stats                         = "Table S4",          # S4 = Tumor QC
  table_removed_clusters                 = "Table S8",          # S8 = tally of removed clusters
  table_1q_clone_per_cluster             = "Table X2",          # X2 = clone per cluster (1q+ subtable A)
  table_16q_clone_per_cluster            = "Table X2 (16q-)",   # X2 subtable B
  table_2p_clone_per_cluster             = "Table X2 (2p+)"     # X2 subtable C
)

# Values tibble for tarchetypes::tar_map() over debranched sample/branch IDs.
# Columns:
#   id       - sample/branch label used in target names (e.g. clustree_SRX10264519)
#   seu_path - path to the filtered Seurat RDS for that sample/branch
debranched_map_values <- tibble::tibble(
  id = c(
    "SRX10264519",          "SRX10264520",
    "SRX10264523",
    "SRX10264524",
    "SRX10264525",          "SRX10264526",
    "SRX11133594",          "SRX11133593",          "SRX11133592",
    "SRX11133588",
    "SRX11133587",
    "SRX11133585",
    "SRX14116947",          "SRX14116944",
    "SRX22868105",
    "SRX22868102"
  )
) |>
  dplyr::mutate(seu_path = paste0("output/seurat/", id, "_filtered_seu.rds"))

## ------------------------------------------------------------------------------------ ##
## Which numbat run the pipeline consumes.
##
## Switched from "output/numbat_sridhar/" (t = 1e-2) to the t = 1e-5 arm on
## 2026-09-21. On 32 matched samples t=1e-5 raised final-round canonical RB events
## 82 -> 104 and collapsed the cross-round instability gap from 38% to 2%, with
## canonical segments 2.5x longer and 3x the LLR; better or equal on 31 of 32.
## See docs/srx_vs_srr_numbat_settings.md appendices C-D and pipeline/config.yaml.
##
## SCOPE: this directory holds SRX samples only. The 32 SRR samples cannot be
## rerun (their cellranger/seurat/allele inputs were never migrated from the old
## workstation) and AGENTS.md forbids rebuilding them, so pointing here
## deliberately narrows the cohort to SRX rather than mixing two values of t
## across the SRR/SRX split -- which is also a study/batch split, so the confound
## would land exactly on the comparison of interest.
##
## Four samples have no complete round at t=1e-5 and are absent here:
## SRX10031191 (numbat never ran, 6 cells), SRX11133590 (no CNV above min_LLR),
## SRX10031192 and SRX11133586 (no CNV surviving the max_entropy filter).
##
## THREE MIXED-PROVENANCE EXCEPTIONS (user decision 2026-09-21). SRX11133592/93/94
## are held out of every rerun by project policy, so they have no t=1e-5 object.
## Their existing objects were COPIED IN (sources in output/numbat_sridhar/ were
## read only, never modified) so they stay in the cohort. They do NOT match the
## other 31 samples, and not in the way you would guess:
##
##   sample        numbat   t      alpha   min_LLR  max_entropy
##   (the other 31) 1.5.2   1e-5   1e-4    5        0.7
##   SRX11133592    1.2.2   1e-5   1e-4    2        0.5   <- legacy workstation run
##   SRX11133593    1.5.2   1e-2   1e-3    5        0.7   <- the only t mismatch
##   SRX11133594    1.2.2   1e-5   1e-4    2        0.5   <- legacy workstation run
##
## So 92/94 match on t but are numbat 1.2.2 with looser thresholds; only 93 differs
## on t. Per-sample parameters for every member of the cohort are in
## results/numbat_cohort_provenance.csv -- consult it before any cross-sample
## claim, and exclude these three from anything comparing call rates or extents.
##
## REVERT: set this back to "output/numbat_sridhar/" to return to the production
## t = 1e-2, 71-sample cohort. Nothing else needs to change.
## ------------------------------------------------------------------------------------ ##
## NOT a variable, deliberately. It was one until 2026-09-22, and the rebuild
## (job 12258888) died with "object 'numbat_rds_dir' not found": the path is used
## inside tarchetypes::tar_files(), whose command is evaluated on a crew worker
## that could not resolve the global, taking 49 targets down with it. The literal
## string is inlined instead -- the pattern that worked before -- at:
##     R/pipeline_targets_inputs.R  (4 call sites)
##     R/pipeline_targets_qc.R      (1 call site)
## To switch cohorts, change the string at all five. Keep them in sync.
