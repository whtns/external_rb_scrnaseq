# Targets Summary

This document lists each target in the pipeline, its dependencies, and a brief description (inferred from code and comments).


## figures_and_tables

## Additional Targets (from line 1241 onward)

### fig_s12
- **Dependencies:** integrated_seus_16q, cluster_orders
- **Description:** Sample-specific analyses of tumors with 16q- subclones after integration.

### fig_s04_08
- **Dependencies:** integrated_seus_2p, cluster_orders
- **Description:** Sample-specific analyses of tumors with 2p+ subclones after integration.

### fig_s11
- **Dependencies:** integrated_seus_2p, cluster_orders
- **Description:** Sample-specific analyses of tumors with 2p+ subclones after integration (alternate figure).

### fig_s23
- **Dependencies:** integrated_seus_6p, cluster_orders
- **Description:** Sample-specific analyses of tumors with 6p+ subclones after integration.

### fig_s25
- **Dependencies:** subtype_markers
- **Description:** Plots for subtype marker analysis (figure S25).

### fig_s0x
- **Dependencies:** None
- **Description:** Karyogram plots (figure S0x).

### fig_03
- **Dependencies:** cluster_orders
- **Description:** Alternative resolutions for integrated 16q- analysis (figure 03).

### fig_04_07
- **Dependencies:** cluster_orders
- **Description:** Alternative resolutions for integrated 16q- analysis (figure 04_07).

### fig_04
- **Dependencies:** integrated_seus_2p, cluster_orders, large_clone_comparisons
- **Description:** Plots for integrated 2p analysis (figure 04).

### fig_05
- **Dependencies:** integrated_seus_6p, cluster_orders, large_clone_comparisons
- **Description:** Plots for integrated 6p analysis (figure 05).

### fig_s08
- **Dependencies:** Not specified
- **Description:** Additional figure S08 (details in code).

### fig_s10
- **Dependencies:** None
- **Description:** Differential expression comparisons between integrated 16q clones within clusters of interest (figure S10).

### fig_s20
- **Dependencies:** Not specified
- **Description:** Additional figure S20 (details in code).

### integrated_seus_1q
- **Dependencies:** None
- **Description:** Integrated Seurat objects for 1q+ subclones after integration.

### integrated_seus_16q
- **Dependencies:** None
- **Description:** Integrated Seurat objects for 16q- subclones after integration.

### integrated_seus_2p
- **Dependencies:** None
- **Description:** Integrated Seurat objects for 2p+ subclones after integration.

### integrated_seus_6p
- **Dependencies:** None
- **Description:** Integrated Seurat objects for 6p+ subclones after integration.

### fig_07a_input
- **Dependencies:** Not specified
- **Description:** Input for figure 07a (details in code).

### fig_07a
- **Dependencies:** Not specified
- **Description:** Figure 07a (details in code).

### fig_07b_input
- **Dependencies:** Not specified
- **Description:** Input for figure 07b (details in code).

### fig_07b
- **Dependencies:** Not specified
- **Description:** Figure 07b (details in code).

### fig_s09
- **Dependencies:** None
- **Description:** Differential expression comparisons between integrated 1q clusters of interest (figure S09).

### fig_08a_input
- **Dependencies:** Not specified
- **Description:** Input for figure 08a (details in code).

### fig_08a
- **Dependencies:** Not specified
- **Description:** Figure 08a (details in code).

### fig_08b_input
- **Dependencies:** Not specified
- **Description:** Input for figure 08b (details in code).

### fig_08b
- **Dependencies:** Not specified
- **Description:** Figure 08b (details in code).

### fig_09
- **Dependencies:** corresponding_seus_2p, corresponding_seus, corresponding_clusters_diffex, corresponding_clusters_enrichments
- **Description:** Differential expression and recurrence analysis for 2p (figure 09).

### fig_10
- **Dependencies:** corresponding_seus_6p, corresponding_seus, corresponding_clusters_diffex, corresponding_clusters_enrichments
- **Description:** Differential expression and recurrence analysis for 6p (figure 10).

### corresponding_seus
- **Dependencies:** None
- **Description:** Seurat objects corresponding to sample-specific analyses of tumors with 16q- subclones without integration.

### corresponding_states_dictionary
- **Dependencies:** Not specified
- **Description:** State dictionary for corresponding analyses (details in code).

### states_dictionary_2p
- **Dependencies:** Not specified
- **Description:** State dictionary for 2p analyses (details in code).

### corresponding_clusters_diffex_2p
- **Dependencies:** states_dictionary_2p, large_clone_comparisons, numbat_rds_files
- **Description:** Finds differential expression clusters between corresponding states for 2p.

### corresponding_clusters_volcanos_2p
- **Dependencies:** corresponding_clusters_diffex_2p
- **Description:** Plots volcano plots for corresponding clusters (2p).

### corresponding_clusters_enrichments_2p
- **Dependencies:** corresponding_clusters_diffex_2p
- **Description:** Plots enrichment for corresponding clusters (2p).

### states_dictionary_6p
- **Dependencies:** Not specified
- **Description:** State dictionary for 6p analyses (details in code).

### corresponding_state_6p_seus
- **Dependencies:** Not specified
- **Description:** Seurat objects for corresponding state 6p analyses.

### corresponding_clusters_diffex_6p
- **Dependencies:** Not specified
- **Description:** Finds differential expression clusters for 6p (details in code).

### corresponding_clusters_volcanos_6p
- **Dependencies:** Not specified
- **Description:** Plots volcano plots for corresponding clusters (6p).

### corresponding_clusters_enrichments_6p
- **Dependencies:** Not specified
- **Description:** Plots enrichment for corresponding clusters (6p).

### corresponding_clusters_volcanos
- **Dependencies:** Not specified
- **Description:** Plots volcano plots for corresponding clusters (general).

### corresponding_clusters_heatmaps
- **Dependencies:** Not specified
- **Description:** Plots heatmaps for corresponding clusters (general).

### corresponding_clusters_enrichments
- **Dependencies:** Not specified
- **Description:** Plots enrichment for corresponding clusters (general).

### clone_cc_plots_by_scna_1q
- **Dependencies:** debranched_seus_1q, large_clone_comparisons
- **Description:** Plots clone cluster comparisons by SCNA for 1q.

### clone_cc_plots_by_scna_16q
- **Dependencies:** debranched_seus_16q, large_clone_comparisons
- **Description:** Plots clone cluster comparisons by SCNA for 16q.

### enrichment_2p_g1
- **Dependencies:** diffex_2p_g1
- **Description:** Compares enrichment for 2p G1.

### enrichment_6p_g1
- **Dependencies:** diffex_6p_g1
- **Description:** Compares enrichment for 6p G1.

- **Dependencies:** clustree_compilation, collage_compilation, table_all_diffex_clones, oncoprint_plots, cluster_comparisons_by_phase_for_disctinct_clones
- **Description:** Collects all major figures and tables for reporting.

## plae_ref
- **Dependencies:** None
- **Description:** Reference data generated by `generate_plae_ref()`.

## numbat_rds_files
- **Dependencies:** interesting_samples
- **Description:** List of Numbat RDS files for selected samples.

## seus
- **Dependencies:** interesting_samples
- **Description:** List of Seurat objects for selected samples.

## celltype_markers
- **Dependencies:** None
- **Description:** Marker genes for major retinal cell types.

## cells_to_remove
- **Dependencies:** None
- **Description:** Loads cell IDs to remove from analysis.

## cluster_dictionary
- **Dependencies:** None
- **Description:** Loads cluster dictionary for annotation.

## hallmark_gene_sets
- **Dependencies:** None
- **Description:** Loads MSigDB hallmark gene sets.

## interesting_samples
- **Dependencies:** None
- **Description:** List of sample IDs considered interesting for analysis.

## debranched_ids
- **Dependencies:** None
- **Description:** List of sample IDs with debranched clones.

## large_clone_comparisons
- **Dependencies:** None
- **Description:** Mapping of sample IDs to clone comparison groups.

## large_clone_simplifications
- **Dependencies:** None
- **Description:** Simplified clone labels for each sample.

## large_filter_expressions
- **Dependencies:** None
- **Description:** Filtering expressions for clones in each sample.

## clone_comparison_table
- **Dependencies:** large_clone_comparisons
- **Description:** Tabulates clone comparisons for reporting.

## whole_pseudobulks
- **Dependencies:** numbat_rds_files, subtype_markers
- **Description:** Scores whole pseudobulks for each sample.

## unfiltered_derived_pseudobulk_subtype_scores
- **Dependencies:** unfiltered_seus
- **Description:** Derives subtype scores from unfiltered Seurat objects.

## filtered_derived_pseudobulk_subtype_scores
- **Dependencies:** final_seus
- **Description:** Derives subtype scores from filtered Seurat objects.

## study_cell_stats
- **Dependencies:** None
- **Description:** Collects study-level cell statistics.

## fig_s15
- **Dependencies:** study_cell_stats
- **Description:** Plots study metadata.

## table_s01
- **Dependencies:** study_cell_stats
- **Description:** Table of UMI, genes detected, mito %.

## table_s02
- **Dependencies:** None
- **Description:** Table of detailed sample metadata.

## table_s03
- **Dependencies:** cluster_dictionary
- **Description:** Table of clusters removed with marker genes.

## table_s04_fig_s16
- **Dependencies:** None
- **Description:** RB SCNA frequency in TCGA by cancer type.

## table_s07
- **Dependencies:** None
- **Description:** For each 1q+ sample, percentage of clone in cluster.

## table_s09
- **Dependencies:** None
- **Description:** For each 16q- sample, percentage of clone in cluster.

## table_s10
- **Dependencies:** None
- **Description:** For each 16q- sample, percentage of clone in cluster.

## total_metadata
- **Dependencies:** None
- **Description:** Loads total metadata from TSV file.

## excluded_samples
- **Dependencies:** None
- **Description:** List of samples excluded from analysis.

## interesting_genes
- **Dependencies:** None
- **Description:** List of genes of interest for analysis.

## subtype_violins
- **Dependencies:** final_seus, numbat_rds_files, large_clone_simplifications, subtype_markers
- **Description:** Generates violin plots of subtype scores for each sample.

## subtype_violin_files
- **Dependencies:** subtype_violins
- **Description:** Combines violin plots into a single PDF.

## patchwork_phase_plot
- **Dependencies:** final_seus
- **Description:** Plots phase distribution of all samples by SCNA.

## mp_heatmaps
- **Dependencies:** final_seus, mps["Cancer"]
- **Description:** Generates heatmaps for cancer-related markers.

## integrated_seus_16q
- **Dependencies:** None
- **Description:** Integrated Seurat objects for 16q- subclones after integration.

## fig_s12
- **Dependencies:** integrated_seus_16q, cluster_orders
- **Description:** Sample-specific analyses of tumors with 16q- subclones after integration.

## integrated_seus_1q
- **Dependencies:** None
- **Description:** Integrated Seurat objects for 1q+ subclones after integration.

## fig_s07
- **Dependencies:** integrated_seus_1q, cluster_orders
- **Description:** Sample-specific analyses of tumors with 1q+ subclones after integration.

## integrated_seus_2p
- **Dependencies:** None
- **Description:** Integrated Seurat objects for 2p+ subclones after integration.

## fig_s04_08
- **Dependencies:** integrated_seus_2p, cluster_orders
- **Description:** Sample-specific analyses of tumors with 2p+ subclones after integration.

## integrated_seus_6p
- **Dependencies:** None
- **Description:** Integrated Seurat objects for 6p+ subclones after integration.

## fig_s23
- **Dependencies:** integrated_seus_6p, cluster_orders
- **Description:** Sample-specific analyses of tumors with 6p+ subclones after integration.

---

_Note: This summary now includes all targets from the full pipeline, including integrated_seus_16q and related analysis targets._
