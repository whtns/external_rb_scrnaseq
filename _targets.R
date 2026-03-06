
# All R scripts in ./R/ are sourced below, including packages.R and functions.R if present.
## Load your packages, e.g. library(targets).
suppressPackageStartupMessages(source("./packages.R"))
## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

# Uncomment and use these for debugging as needed:

tar_option_set(
  memory = "transient",
  garbage_collection = TRUE,
  error = "continue"  # "null" silently propagates NULLs; "continue" marks failed targets
                      # and lets you inspect errors via tar_meta(fields = "error")
  # controller = crew_controller_local(workers = 4)
)

output_plot_extensions <- c(
  "dimplot.pdf",
  "merged_marker.pdf",
  "combined_marker.pdf",
  "phylo_probability.pdf",
  "clone_distribution.pdf"
)

# Values tibble for tarchetypes::tar_map() over the 4 SCNA types.
# Defined here (outside tar_plan) so it is available at pipeline definition time.
# Columns:
#   scna      - SCNA label used in target names and as list index key
#   seus_sym  - symbol of the corresponding debranched Seurat path target
#   var_y     - y-axis grouping for clone_pearls plots
scna_map_values <- tibble::tibble(
  scna     = c("1q",               "2p",               "6p",               "16q"),
  seus_sym = rlang::syms(c("debranched_seus_1q", "debranched_seus_2p",
                            "debranched_seus_6p", "debranched_seus_16q")),
  var_y    = c("phase_level",      "clusters",         "clusters",         "phase_level")
)

## tar_plan supports drake-style targets and also tar_target()
tar_plan(

  # Aggregate all main figures and tables into a single list for final output tracking.
  tar_target(
    figures_and_tables,
    list(
      fig_s02_table_s04, # tcga frequency
      fig_s03, # gistic plots
      fig_s04, # study cell stats
      fig_s05, # smoothed expression
      fig_s06a, # karyograms,
      fig_02, # integrated 1q fig for all 1q+ samples thesis fig 4.3
      fig_s07, # Sample-specific analyses of tumors with 1q+ subclones after integration.
      fig_s11,
      fig_s12, # Sample-specific analyses of tumors with 2p+ subclones after integration.
      fig_s04_08, # Sample-specific analyses of tumors with 16q- subclones after integration.
      clustree_compilation,
      collage_compilation,
      table_all_diffex_clones,
      # oncoprint_enrich_clones_plots_gobp,
      # oncoprint_enrich_clusters_plots_gobp,
      oncoprint_plots, 
      cluster_comparisons_by_phase_for_disctinct_clones
    )
  ),

  # Purpose: Track file path input/output for all numbat rds files.
  tarchetypes::tar_files(all_numbat_rds_files, retrieve_numbat_rds_files("output/numbat_sridhar/"), format = "file"),
  # Purpose: Prepare Seurat object path(s) for all seus.
  tarchetypes::tar_files(all_seus, retrieve_seus("output/seurat/"), format = "file"),
  # Purpose: Track file path input/output for numbat rds files.
  tarchetypes::tar_files(numbat_rds_files, retrieve_numbat_rds_files("output/numbat_sridhar/", interesting_samples), format = "file"),
  # Purpose: Track file path input/output for numbat rds filtered files.
  tarchetypes::tar_files(numbat_rds_filtered_files, retrieve_numbat_rds_files("output/numbat_sridhar_filtered/", interesting_samples), format = "file"),
  # Purpose: Build dependency target for seus.
  tarchetypes::tar_files(seus, retrieve_seus("output/seurat/", interesting_samples), format = "file"),

  # Read cluster dictionary mapping cluster IDs to cell type labels from TSV.
  tar_target(cluster_dictionary, read_cluster_dictionary("data/cluster_dictionary.tsv")),

 # Load clone comparison pairs (which clones to contrast) from YAML config.
 tar_target(large_clone_comparisons,
   yaml::read_yaml(here::here("config/large_clone_comparisons.yaml"))
 ),
 # Load clone label simplification mappings from YAML config.
 tar_target(large_clone_simplifications,
   yaml::read_yaml(here::here("config/large_clone_simplifications.yaml"))
 ),
  # Load cell filter expressions for removing unwanted clone populations from YAML.
  tar_target(large_filter_expressions,
    yaml::read_yaml(here::here("config/large_filter_expression.yaml"))
  ),
  
  # Define marker gene lists for major retinal cell types used in annotation.
  tar_target(
    celltype_markers,
    list(
      "cones" = c("PDE6H", "ARR3", "MYL4", "RCVRN", "GNAT2"),
      "rods" = c("GNGT1", "RHO", "CNGA1", "SAG", "PDE6G"),
      "Müller glia" = c("TF", "NPVF", "RLBP1", "RRH", "RGR"),
      "retinal astrocytes" = c("TRH", "MLC1", "NDP", "PAX5", "PAX2"),
      "microglia" = c("RGS1", "C1QB", "C1QA", "C1QC", "CCL3"),
      "bipolar cells" = c("CABP5", "TMEM215", "GABRA1", "CA10", "PCP2"),
      "retinal ganglion cells" = c("NEFL", "NEFM", "PVALB", "SNCG", "STMN2"),
      "amacrine cells" = c("C1QL2", "GAD2", "CARTPT", "GAD1", "SLC32A1"),
      "RPE" = c("TTR", "RPE65", "RGR", "RLBP1", "BEST1")
      # "neurogenic" = c("PTTG1", "SCG3", "SHD", "SOX11", "ONECUT2", "VSX1", "ASCL1", "TCF4", "SOX4", "RORB", "GADD45G", "MYBL1", "MFAP4", "HES6", "NEUROD1", "C8orf46", "OTX2", "OLIG2", "INSM1", "PCBP4", "NEUROD4", "FAM131C"),
      # "Neurogenic_supp" = c("HES2", "PKIB", "PCDH17", "ATOH7", "DLX1", "DLX2", "ELAVL4", "SOX11", "PENK", "GAL", "ONECUT2", "GADD45A", "HES6", "ISL1", "ONECUT1", "ONECUT3", "FOXP4", "NEUROD1", "ASCL1", "GADD45G", "DLL4", "NEUROD4", "SOX4", "DLL1", "C8orf46", "OLIG2", "NEUROG2", "OTX2", "HES4", "PTF1A", "SSTR2", "BTG2", "RGS16")
    )
  ),

  # setup ------------------------------

  # Read list of low-quality cells to exclude from analysis from Excel file.
  tar_target(cells_to_remove, read_cells_to_remove("data/cells_to_remove_final.xlsx")),
  # Fetch MSigDB Hallmark gene sets for human from msigdbr.
  tar_target(hallmark_gene_sets, msigdbr(species = "Homo sapiens", category = "H")),
  # Define the vector of SRR accession IDs selected for downstream analysis.
  tar_target(interesting_samples,
    c(
      "SRR13884242",
      "SRR13884243",
      "SRR13884244",
      "SRR13884245",
      "SRR13884246",
      "SRR13884247",
      "SRR13884248",
      "SRR13884249",
      "SRR14800534",
      "SRR14800535",
      "SRR14800536",
      "SRR14800540",
      "SRR14800541",
      "SRR14800543",
      "SRR17960481",
      "SRR17960484",
      "SRR27187899",
      "SRR27187900",
      "SRR27187901",
      "SRR27187902"
    )
  ),
  # Define paths to the unfiltered Seurat RDS files for all samples.
  tar_target(original_seus,
    c(
      "output/seurat/SRR13884242_unfiltered_seu.rds",
      "output/seurat/SRR13884243_unfiltered_seu.rds",
      "output/seurat/SRR13884246_unfiltered_seu.rds",
      "output/seurat/SRR13884247_unfiltered_seu.rds",
      "output/seurat/SRR13884248_unfiltered_seu.rds",
      "output/seurat/SRR13884249_unfiltered_seu.rds",
      "output/seurat/SRR14800534_unfiltered_seu.rds",
      "output/seurat/SRR14800535_unfiltered_seu.rds",
      "output/seurat/SRR14800536_unfiltered_seu.rds",
      "output/seurat/SRR14800540_unfiltered_seu.rds",
      "output/seurat/SRR14800541_unfiltered_seu.rds",
      "output/seurat/SRR14800543_unfiltered_seu.rds",
      "output/seurat/SRR17960481_unfiltered_seu.rds",
      "output/seurat/SRR17960484_unfiltered_seu.rds",
      "output/seurat/SRR27187899_unfiltered_seu.rds",
      "output/seurat/SRR27187902_unfiltered_seu.rds"
    )
  ),
  # Samples excluded from analysis (lack a diploid ancestral clone).
    tar_target(
    excluded_samples,
    c(
      "SRR14800541" # no diploid ancestral clone; only diploid are B2M and APOE
    )
  ),
  # Define IDs for per-branch debranched Seurat objects across all samples.
  tar_target(debranched_ids,
    c(
      "SRR13884242", "SRR13884243",
      "SRR13884246_branch_5", "SRR13884246_branch_6",
      "SRR13884247_branch_6", "SRR13884247_branch_4", "SRR13884247_branch_5",
      "SRR13884248", "SRR13884249",
      "SRR14800534", "SRR14800535", "SRR14800536",
      "SRR14800540_branch_2", "SRR14800540_branch_3",
      "SRR14800541_branch_4", "SRR14800541_branch_7",
      "SRR14800543_branch_3", "SRR14800543_branch_4",
      "SRR17960481", "SRR17960484",
      "SRR27187899",
      "SRR27187902_branch_3", "SRR27187902_branch_4"
    )
  ),

  # Tabulate clone comparison pairs into a summary table from the comparisons config.
  tar_target(
    clone_comparison_table,
    tabulate_clone_comparisons(large_clone_comparisons)
  ),

  # total metadata ------------------------------
  # Collect per-sample cell count and QC statistics across all studies.
  tar_target(study_cell_stats, collect_study_metadata()),

  # Plot study-level metadata summary for supplemental figure S04.
  tar_target(fig_s04, plot_study_metadata(study_cell_stats)),

  # Purpose: Generate output for table s01.
  tar_target(table_s01, make_table_s01(study_cell_stats)), # umi, genes detected, mito %

  # Purpose: Generate output for table s02.
  tar_target(table_s02, make_table_s02()), # detailed sample metadata

  # Purpose: Generate output for table s03.
  tar_target(table_s03, make_table_s03(cluster_dictionary)), # tally clusters removed (filtered out) with marker genes
  
  # Purpose: Generate output for fig s02 table s04.
  tar_target(fig_s02_table_s04, rb_scna_frequency_in_tcga_by_cancer_type()), # RB SCNA frequency in TCGA by cancer type (Taylor et al. 2018).
  
  # Purpose: Generate output for fig s03.
  tar_target(fig_s03, make_gistic_plot()), # RB SCNA frequency in TCGA by cancer type (Taylor et al. 2018).
  # Purpose: Generate output for table s07.
  tar_target(table_s07, make_table_s07()), # For each 1q+ sample, percentage of clone in cluster.
  # Purpose: Generate output for table s09.
  tar_target(table_s09, make_table_s09(seu_path = "output/seurat/integrated_16q/integrated_seu_16q_complete.rds", table_path = "results/table_s09.csv")), # For each 16q- sample, percentage of clone in cluster.
  # Purpose: Generate output for table s10.
  tar_target(table_s10, make_table_s10(seu_path = "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", table_path = "results/table_s10.csv")), # For each 16q- sample, percentage of clone in cluster.
  
  # Read sample-level metadata table from TSV.
  tar_target(total_metadata, read_tsv("data/metadata.tsv")),

  # Define per-SCNA and per-subtype lists of interesting candidate genes.
  tar_target(
    interesting_genes,
    list(
      "1q" = c("ANP32E", "ASH1L", "ASPM", "CENPF", "CKS1B", "CNIH4", "CRABP2", "ENAH", "ENO1", "IPO9", "KIF14", "NEK2", "NENF", "NUF2", "PSMD4", "RXRG", "UBE2T"),
      "2p" = c("CEP68", "MEIS1", "MDH1", "OST4", "PDIA6", "POLE4", "RAB1A", "SNRPG", "SOX11"),
      "6p" = c("ARID1B", "CASP8AP2", "CLIC1", "CUTA", "DDX39B", "DEK", "DST", "GLO1", "HDAC2", "HDDC2", "HMGA1", "LIN28B", "MCM3", "MNDA", "PRDM1", "SOX4"),
      "16q" = c("CHD9", "CNGB1", "CYBA", "MT1E", "MT1X", "MT2A"),
      "S1" = c("ARR3", "CRX", "PDC"),
      "S2" = c("EBF3", "DCX", "ROBO1", "SOX11", "GAP43", "PCDHB10", "STMN2", "NEFM", "POU4F2", "TFF1", "CD24"),
      "interesting" = c("SETD6", "RBL2", "DDX19A", "EZH2", "RB1", "MDM4", "ESRRG"),
      "cell_cycle" = c("CCNA2", "CCNB1")
    )
  ),

  # plot subtype scores
  # Score each sample for subtype signatures and generate violin plots per sample.
  tar_target(subtype_violins,
    score_and_vlnplot_seu(final_seus, numbat_rds_files, large_clone_simplifications, subtype_markers),
    pattern = map(final_seus),
    iteration = "list"
  ),

  # Combine all per-sample subtype violin PDFs into a single merged PDF.
  tar_target(
  	subtype_violin_files,
  	qpdf::pdf_combine(unlist(subtype_violins), "results/scna_vlns.pdf")
  ),
  # Plot cell cycle phase distribution across all samples grouped by SCNA type.
  tar_target(
    patchwork_phase_plot,
    plot_phase_distribution_of_all_samples_by_scna(final_seus)
  ),

  # Score cells for cancer meta-program signatures and plot per-sample heatmaps.
  tar_target(mp_heatmaps,
    score_and_heatmap_seu(final_seus, mps[["Cancer"]], leiden_cluster_file = "results/adata_filtered_metadata_0.2.csv"),
    pattern = map(final_seus),
    iteration = "list"
  ),

  # Plot violin plots of interesting genes across samples grouped by SCNA clone.
  tar_target(
    marker_gene_vlnplots_by_clone,
    plot_putative_marker_across_samples(unlist(interesting_genes), scna_seus, plot_type = VlnPlot, group_by = "scna")
  ),

  # Plot violin plots of interesting genes across all SCNA samples grouped by cluster.
  tar_target(
  	cluster_vlns_all,
    plot_putative_marker_across_samples(unlist(interesting_genes), scna_seus, plot_type = VlnPlot, group_by = "clusters")
  ),
  
  # cluster_vlns_1q / _2p / _6p / _16q — generated via tar_map over scna_map_values
  tarchetypes::tar_map(
    values = scna_map_values[, c("scna", "seus_sym")],
    names  = "scna",
    # Purpose: Build dependency target for cluster vlns.
    tar_target(cluster_vlns,
      plot_putative_marker_across_samples(
        interesting_genes[[scna]], seus_sym,
        plot_type = VlnPlot, group_by = "clusters", extension = scna
      )
    )
  ),
  
  # Plot violin plots of S2-subtype marker genes across samples by cluster.
  tar_target(
  	cluster_vlns_s2,
  	plot_putative_marker_across_samples(interesting_genes[["S2"]], scna_seus, plot_type = VlnPlot, group_by = "clusters", extension = "s2")
  ),

  # Plot violin plots of S1-subtype marker genes across samples by cluster.
  tar_target(
  	cluster_vlns_s1,
  	plot_putative_marker_across_samples(interesting_genes[["S1"]], scna_seus, plot_type = VlnPlot, group_by = "clusters", extension = "s1")
  ),
  
  # factor_heatmaps_{scna} + factor_heatmaps_combined_{scna} via tar_map
  # Creates targets: factor_heatmaps_1q/2p/6p/16q  (dynamic branches, one per sample)
  #              and factor_heatmaps_combined_1q/2p/6p/16q  (per-SCNA merged PDFs)
  tarchetypes::tar_map(
    values = scna_map_values[, c("scna", "seus_sym")],
    names  = "scna",
    # Purpose: Build dependency target for factor heatmaps.
    tar_target(factor_heatmaps,
      plot_seu_gene_heatmap(seus_sym, large_clone_comparisons, scna_of_interest = scna, w = 8, h = 12),
      pattern   = map(seus_sym),
      iteration = "list"
    ),
    # Purpose: Build dependency target for factor heatmaps combined.
    tar_target(factor_heatmaps_combined,
      qpdf::pdf_combine(factor_heatmaps, paste0("results/factor_heatmaps_", scna, ".pdf"))
    )
  ),
  # Combine all four SCNA factor heatmap PDFs into the fig_s16 supplemental PDF.
  tar_target(all_factor_heatmaps,
    qpdf::pdf_combine(
      c(factor_heatmaps_combined_1q, factor_heatmaps_combined_2p,
        factor_heatmaps_combined_6p, factor_heatmaps_combined_16q),
      "results/fig_s16.pdf"
    )
  ),
  
  # Plot feature plots of interesting genes across samples grouped by cluster.
  tar_target(
    marker_gene_featureplots_by_cluster,
    plot_putative_marker_across_samples(unlist(interesting_genes), scna_seus, plot_type = FeaturePlot, group_by = "clusters")
  ),
  
  # Purpose: Build dependency target for cluster markers.
  tar_target(
    cluster_markers,
    collect_all_markers(numbat_rds_files, "results/clusters.xlsx")
  ),

  # Purpose: Build dependency target for large clusters.
  tar_target(
    large_clusters,
    collect_clusters_from_seus(final_seus)
  ),

  # Purpose: Create PDF artifact(s) for large numbat pdfs.
  tar_target(large_numbat_pdfs,
    convert_numbat_pngs(numbat_rds_files),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),

  # Purpose: Track file path input/output for branch dictionary file.
  tarchetypes::tar_file(
    branch_dictionary_file,
    "data/branch_dictionary.csv"
  ),
  # Purpose: Build dependency target for branch dictionary.
  tar_target(
    branch_dictionary,
    pull_branches(branch_dictionary_file)
  ),

  # Purpose: Track file path input/output for cluster file.
  tar_target(cluster_file,
  					 "data/scna_cluster_order.csv", 
  					 format = "file_fast"),
  
  # Purpose: Build dependency target for cluster orders.
  tar_target(cluster_orders,
    pull_cluster_orders(cluster_file)
  ),
  
  # Purpose: Prepare Seurat object path(s) for overall seus.
  tar_target(overall_seus,
  					 merge_orders(debranched_seus,
  					 						 final_seus)
  ),
  
  # Purpose: Track file path input/output for cluster comparisons file.
  tarchetypes::tar_file(
  	cluster_comparisons_file,
  	"data/cluster_comparisons_by_phase_for_disctinct_clones.csv"
  ),
  
  # Purpose: Build dependency target for cluster comparisons.
  tar_target(cluster_comparisons,
  					 pull_cluster_comparisons(cluster_comparisons_file)),
  
  # Purpose: Compute clone-related output for cluster comparisons by phase for disctinct clones.
  tar_target(cluster_comparisons_by_phase_for_disctinct_clones,
  					 make_cluster_comparisons_by_phase_for_disctinct_clones(cluster_comparisons, overall_seus),
  					 pattern = map(cluster_comparisons),
  					 iteration = "list"
  ),
  
  # Purpose: Track file path input/output for filtered large plot files.
  tar_target(filtered_large_plot_files,
    make_numbat_plot_files(numbat_rds_files, seus, cluster_dictionary, large_filter_expressions, large_clone_simplifications, extension = "_filtered"),
    pattern = map(numbat_rds_files, seus),
    iteration = "list"
  ),
  
  # Purpose: Prepare Seurat object path(s) for unfiltered seus.
  tar_target(unfiltered_seus,
    prep_unfiltered_seu(numbat_rds_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions, extension = "_unfiltered"),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),
  

  
  # Purpose: Prepare Seurat object path(s) for filtered seus.
  tar_target(filtered_seus,
  filter_cluster_save_seu(numbat_rds_files, seus, cluster_dictionary, large_clone_simplifications, filter_expressions = NULL, cells_to_remove, extension = "", leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv"),
  pattern = map(numbat_rds_files),
  iteration = "list"
  ),

  # Purpose: Prepare Seurat object path(s) for final seus.
  tar_target(
    final_seus,
    set_final_seus(interesting_samples)
  ),

  
  # Purpose: Build dependency target for cc plots wo arms.
  tar_target(
  	cc_plots_wo_arms,
  	plot_phase_wo_arm(final_seus)
  ),
  
  # Purpose: Track file path input/output for debranched seu files.
  tar_file(debranched_seu_files,
    c(
      "SRR13884242"         = "output/seurat/SRR13884242_filtered_seu.rds",
      "SRR13884243"         = "output/seurat/SRR13884243_filtered_seu.rds",
      "SRR13884246_branch_5" = "output/seurat/SRR13884246_branch_5_filtered_seu.rds",
      "SRR13884246_branch_6" = "output/seurat/SRR13884246_branch_6_filtered_seu.rds",
      "SRR13884247_branch_6" = "output/seurat/SRR13884247_branch_6_filtered_seu.rds",
      "SRR13884247_branch_4" = "output/seurat/SRR13884247_branch_4_filtered_seu.rds",
      "SRR13884247_branch_5" = "output/seurat/SRR13884247_branch_5_filtered_seu.rds",
      "SRR13884248"         = "output/seurat/SRR13884248_filtered_seu.rds",
      "SRR13884249"         = "output/seurat/SRR13884249_filtered_seu.rds",
      "SRR14800534"         = "output/seurat/SRR14800534_filtered_seu.rds",
      "SRR14800535"         = "output/seurat/SRR14800535_filtered_seu.rds",
      "SRR14800536"         = "output/seurat/SRR14800536_filtered_seu.rds",
      "SRR14800540_branch_2" = "output/seurat/SRR14800540_branch_2_filtered_seu.rds",
      "SRR14800540_branch_3" = "output/seurat/SRR14800540_branch_3_filtered_seu.rds",
      "SRR14800541_branch_4" = "output/seurat/SRR14800541_branch_4_filtered_seu.rds",
      "SRR14800541_branch_7" = "output/seurat/SRR14800541_branch_7_filtered_seu.rds",
      "SRR14800543_branch_3" = "output/seurat/SRR14800543_branch_3_filtered_seu.rds",
      "SRR14800543_branch_4" = "output/seurat/SRR14800543_branch_4_filtered_seu.rds",
      "SRR17960481"         = "output/seurat/SRR17960481_filtered_seu.rds",
      "SRR17960484"         = "output/seurat/SRR17960484_filtered_seu.rds",
      "SRR27187899"         = "output/seurat/SRR27187899_filtered_seu.rds",
      "SRR27187902_branch_3" = "output/seurat/SRR27187902_branch_3_filtered_seu.rds",
      "SRR27187902_branch_4" = "output/seurat/SRR27187902_branch_4_filtered_seu.rds"
    )
  ),

  # Purpose: Prepare Seurat object path(s) for debranched seus.
  tar_target(debranched_seus,
    set_names(debranched_seu_files, str_remove(fs::path_file(debranched_seu_files), "_filtered_seu.rds"))
  ),
  
  
  # Purpose: Build dependency target for filtered recurrent genes.
  tar_target(
    filtered_recurrent_genes,
    pull_common_markers(final_seus, mps[["Cancer"]])
  ),
  # Purpose: Build dependency target for filtered recurrent heatmaps.
  tar_target(filtered_recurrent_heatmaps,
    heatmap_marker_genes(
      final_seus,
      filtered_recurrent_genes,
      subtype_markers,
      "filtered_",
      marker_col = "SCT_snn_res.0.2",
      group.by = c("SCT_snn_res.0.2", "Phase", "scna"),
      col_arrangement = c("SCT_snn_res.0.2", "Phase", "scna")
    ),
    pattern = map(final_seus),
    iteration = "list"
  ),
  # Purpose: Track file path input/output for filtered recurrent heatmap file.
  tar_target(
    filtered_recurrent_heatmap_file,
    qpdf::pdf_combine(filtered_recurrent_heatmaps, "results/filtered_recurrent_heatmaps.pdf")
  ),
  # Purpose: Build dependency target for combined recurrent filtered heatmap.
  tar_target(
    combined_recurrent_filtered_heatmap,
    heatmap_marker_genes(
      "output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds",
      filtered_recurrent_genes,
      subtype_markers,
      "filtered_",
      label = "SRR14800534_SRR14800535_SRR14800536",
      marker_col = "SCT_snn_res.0.2",
      group.by = c("integrated_snn_res.0.2", "Phase", "scna"),
      col_arrangement = c("integrated_snn_res.0.2", "Phase", "scna")
    )
  ),
  # Purpose: Build dependency target for regressed recurrent genes.
  tar_target(
    regressed_recurrent_genes,
    pull_common_markers(regressed_seus, mps[["Cancer"]])
  ),
  # Purpose: Build dependency target for regressed recurrent heatmaps.
  tar_target(regressed_recurrent_heatmaps,
    heatmap_marker_genes(
      regressed_seus,
      regressed_recurrent_genes,
      subtype_markers,
      "regressed_",
      marker_col = "SCT_snn_res.0.2",
      group.by = c("SCT_snn_res.0.2", "Phase", "scna"),
      col_arrangement = c("SCT_snn_res.0.2", "scna")
    ),
    pattern = map(regressed_seus),
    iteration = "list"
  ),
  
  # Purpose: Track file path input/output for regressed recurrent heatmap file.
  tar_target(regressed_recurrent_heatmap_file,
    qpdf::pdf_combine(regressed_recurrent_heatmaps, "results/regressed_recurrent_heatmaps.pdf")
  ),

  # Purpose: Prepare Seurat object path(s) for regressed seus.
  tar_target(regressed_seus,
    regress_filtered_seu(final_seus),
    pattern = map(final_seus),
    iteration = "list"
  ),

  # Purpose: Build dependency target for effect of filtering.
  tar_target(effect_of_filtering,
    plot_effect_of_filtering(unfiltered_seus, final_seus),
    pattern = map(unfiltered_seus),
    iteration = "list"
  ),

  # Purpose: Generate plot set for regression effect plots.
  tar_target(regression_effect_plots,
    plot_effect_of_regression(final_seus, regressed_seus, w = 18, h = 12),
    pattern = map(regressed_seus),
    iteration = "list"
  ),
  
  # regression diagnostics ------------------------------
  tar_target(fig_s02_04, qpdf::pdf_combine(unlist(regression_effect_plots), "results/fig_s02_04.pdf")),
  
  # Purpose: Generate plot set for regression ora plots.
  tar_target(regression_ora_plots,
    ora_effect_of_regression(final_seus, regressed_seus),
    pattern = map(regressed_seus),
    iteration = "list"
  ),
  # Purpose: Generate plot set for cin score plots.
  tar_target(cin_score_plots,
    score_chrom_instability(unfiltered_seus),
    pattern = map(unfiltered_seus),
    iteration = "list"
  ),
  
  # rb scna ideograms
  tar_target(fig_s06a,
  					 plot_fig_s06a()),

    # unfiltered_numbat_heatmaps
  tar_target(fig_s03a_plots,
    make_numbat_heatmaps(original_seus, numbat_rds_files, p_min = 0.9, line_width = 0.1, extension = "_unfiltered"),
    pattern = map(original_seus),
    iteration = "list"
  ),
  
  # Purpose: Generate output for fig s03a.
  tar_target(fig_s03a, qpdf::pdf_combine(unlist(fig_s03a_plots), "results/unfiltered_heatmaps.pdf")),
  
  # Purpose: Build dependency target for debranched samples.
  tar_target(
    debranched_samples,
    str_remove(fs::path_file(debranched_seus), "_filtered_seu.rds")
  ),
  
  # Purpose: Build dependency target for clustrees.
  tar_target(clustrees,
    make_clustrees_for_sample(debranched_seus, mylabel = debranched_samples, assay = "SCT"),
    pattern = map(debranched_seus, debranched_samples),
    iteration = "list"
  ),
  
  # Purpose: Build dependency target for clustree compilation.
  tarchetypes::tar_files(clustree_compilation, qpdf::pdf_combine(unlist(clustrees), "results/clustrees.pdf"), format = "file"),
  
  # filtered_numbat_heatmaps
  tar_target(fig_s13,
    make_numbat_heatmaps(final_seus, numbat_rds_files, p_min = 0.5, line_width = 0.1, extension = "_filtered"),
    pattern = map(numbat_rds_files),
    iteration = "list",
    cue = tar_cue("always")  # force re-render each run: heatmap layout depends on
                              # session graphics device state, not captured by hash
  ),
  
  # Purpose: Track file path input/output for filtered numbat heatmaps file.
  tar_target(filtered_numbat_heatmaps_file, qpdf::pdf_combine(map_chr(fig_s13, 1), "results/filtered_heatmaps.pdf")),
  
  # Purpose: Track file path input/output for filtered large scna prob file.
  tar_target(filtered_large_scna_prob_file, qpdf::pdf_combine(map_chr(fig_s13, 2), "results/filtered_scna_probabilities.pdf")),
  
  # Purpose: Create PDF artifact(s) for large montage pdfs.
  tar_target(
    large_montage_pdfs,
    make_pdf_montages(filtered_large_plot_files, large_heatmaps)
  ),
  # Purpose: Create PDF artifact(s) for expr heatmap pdfs.
  tar_target(
    expr_heatmap_pdfs,
    make_expression_heatmap_comparison(large_numbat_pdfs, filtered_numbat_heatmaps)
  ),

  # all ------------------------------
  tar_target(all_diffex_clones,
    find_diffex_clones(debranched_seus, numbat_rds_files, large_clone_comparisons, location = "all"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  # Purpose: Compute differential expression results for all diffex clones for each cluster.
  tar_target(all_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(debranched_seus, numbat_rds_files, large_clone_comparisons, cluster_orders, location = "all"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),

  # Purpose: Build dependency target for volcano thresholds all.
  tar_target(
    volcano_thresholds_all,
    list(
      "SRR13884242" = 4,
      "SRR13884243" = 4,
      "SRR13884247" = 4,
      "SRR13884249" = 4,
      "SRR14800534" = 2,
      "SRR14800535" = 3,
      "SRR14800536" = 3,
      "SRR14800540" = 4,
      "SRR14800541" = 3,
      "SRR14800543" = 3,
      "SRR17960481" = 3,
      "SRR17960484" = 3,
      "SRR27187899" = 3,
      # "SRR27187900" = 3,
      "SRR27187902" = 3
    )
  ),
  # Purpose: Build dependency target for volcano all.
  tar_target(
    volcano_all,
    make_volcano_diffex_clones(
      all_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_all.pdf",
      all_diffex_clones,
      "results/diffex_bw_clones_all.pdf"
    ),
  ),
  # Purpose: Generate output for table all diffex clones.
  tar_target(
    table_all_diffex_clones,
    tabulate_diffex_clones(
      all_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_all.xlsx",
      "results/diffex_bw_clones_per_cluster_all_by_chr.xlsx",
      all_diffex_clones,
      "results/diffex_bw_clones_all.xlsx",
      "results/diffex_bw_clones_all_by_chr.xlsx"
    ),
  ),

  # cis ------------------------------
  tar_target(cis_diffex_clones,
    find_diffex_clones(debranched_seus, numbat_rds_files, large_clone_comparisons, location = "cis"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  
  # Purpose: Compute differential expression results for cis diffex clones for each cluster.
  tar_target(cis_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(debranched_seus, numbat_rds_files, large_clone_comparisons, cluster_orders, location = "cis"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  
  
  # Purpose: Build dependency target for volcano thresholds in segment.
  tar_target(
    volcano_thresholds_in_segment,
    list(
      "SRR13884242" = 4,
      "SRR13884243" = 4,
      "SRR13884247" = 4,
      "SRR13884249" = 4,
      "SRR14800534" = 2,
      "SRR14800535" = 3,
      "SRR14800536" = 3,
      "SRR14800540" = 4,
      "SRR14800541" = 3,
      "SRR14800543" = 3,
      "SRR17960481" = 3,
      "SRR17960484" = 3
    )
  ),
  # Purpose: Build dependency target for volcano cis.
  tar_target(
    volcano_cis,
    make_volcano_diffex_clones(
      cis_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_in_segment.pdf",
      cis_diffex_clones,
      "results/diffex_bw_clones_large_in_segment.pdf"
    ),
  ),
  # Purpose: Generate output for table cis diffex clones.
  tar_target(
    table_cis_diffex_clones,
    tabulate_diffex_clones(
      cis_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_in_segment.xlsx",
      "results/diffex_bw_clones_per_cluster_large_in_segment_by_chr.xlsx",
      cis_diffex_clones,
      "results/diffex_bw_clones_large_in_segment.xlsx",
      "results/diffex_bw_clones_large_in_segment_by_chr.xlsx"
    ),
  ),

  # trans ------------------------------
  tar_target(trans_diffex_clones,
    find_diffex_clones(debranched_seus, numbat_rds_files, large_clone_comparisons, location = "out_of_segment"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  # Purpose: Compute differential expression results for trans diffex clones for each cluster.
  tar_target(trans_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(debranched_seus, numbat_rds_files, large_clone_comparisons, cluster_orders, location = "out_of_segment"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  # Purpose: Build dependency target for volcano trans.
  tar_target(
    volcano_trans,
    make_volcano_diffex_clones(
      trans_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_out_of_segment.pdf",
      trans_diffex_clones,
      "results/diffex_bw_clones_large_out_of_segment.pdf"
    ),
  ),

  # kooi candidates
  tar_target(
    found_kooi_candidates,
    tally_kooi_candidates(cis_diffex_clones = "results/diffex_bw_clones_large_in_segment_by_chr.xlsx", trans_diffex_clones = "results/diffex_bw_clones_large_out_of_segment_by_chr.xlsx")
  ),
  # Purpose: Generate output for table trans diffex clones.
  tar_target(
    table_trans_diffex_clones,
    tabulate_diffex_clones(
      trans_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_out_of_segment.xlsx",
      "results/diffex_bw_clones_per_cluster_large_out_of_segment_by_chr.xlsx",
      trans_diffex_clones,
      "results/diffex_bw_clones_large_out_of_segment.xlsx",
      "results/diffex_bw_clones_large_out_of_segment_by_chr.xlsx"
    ),
  ),

  # Purpose: Prepare integrated-analysis input/output for integrated seu 16q.
  tar_target(
    integrated_seu_16q,
    readRDS("output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds")
  ),

  # Purpose: Build dependency target for large numbat expression.
  tar_target(
    large_numbat_expression,
    retrieve_numbat_plot_type(large_numbat_pdfs, "exp_roll_clust.pdf")
  ),
  # Purpose: Generate output for fig s05.
  tar_target(fig_s05,
    qpdf::pdf_combine(large_numbat_expression, "results/numbat_expression.pdf")
  ),
  # Purpose: Compute clone-related output for clone trees old.
  tar_target(clone_trees_old,
    retrieve_numbat_plot_type(filtered_large_plot_files, "tree_filtered.pdf")
  ),
  
  # Purpose: Generate plot set for paired clone distribution plots.
  tar_target(paired_clone_distribution_plots,
    calculate_clone_distribution(scna_seus, cluster_orders, pairwise = TRUE),
    iteration = "list"
  ),
  
  # clone_pearls_1q/2p/6p/16q — generated via tar_map over scna_map_values
  # var_y: 1q and 16q use "phase_level"; 2p and 6p use "clusters"
  tarchetypes::tar_map(
    values = scna_map_values,
    names  = "scna",
    # Purpose: Compute clone-related output for clone pearls.
    tar_target(clone_pearls,
      plot_clone_pearls(seus_sym, var_y = var_y),
      pattern   = map(seus_sym),
      iteration = "list"
    )
  ),

  # Purpose: Build dependency target for cluster terms.
  tar_target(cluster_terms,
    enrich_by_cluster(scna_seus),
    pattern   = map(scna_seus),
    iteration = "list"
  ),
  
  # Purpose: Generate plot set for unpaired clone distribution plots.
  tar_target(unpaired_clone_distribution_plots,
    calculate_clone_distribution(final_seus, cluster_orders, pairwise = FALSE),
    pattern = map(final_seus),
    iteration = "list"
  ),
  
  # Purpose: Build dependency target for cluster scorecard.
  tar_target(
    cluster_scorecard,
    score_clusters_up_down(paired_clone_distribution_plots, interesting_samples)
  ),
  
  # Purpose: Build dependency target for clusters and markers.
  tar_target(clusters_and_markers,
    plot_seu_clusters_and_markers(debranched_seus, cluster_orders),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  
  # Purpose: Track file path input/output for debranched clone tree files.
  tar_target(debranched_clone_tree_files,
    save_clone_tree_from_path(debranched_seus, numbat_rds_files, large_clone_simplifications, label = "_debranched_clone_tree", legend = FALSE, horizontal = FALSE),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  
  # Purpose: Compute clone-related output for debranched clone trees.
  tar_target(debranched_clone_trees,
  					 plot_clone_tree_from_path(debranched_seus, numbat_rds_files, large_clone_simplifications, label = "_debranched_clone_tree", legend = FALSE, horizontal = FALSE),
  					 pattern = map(debranched_seus),
  					 iteration = "list"
  ),

  # Purpose: Compute clone-related output for clone trees.
  tar_target(clone_trees,
  					 plot_clone_tree_from_path(seus, numbat_rds_files, large_clone_simplifications, label = "_debranched_clone_tree", legend = FALSE, horizontal = FALSE),
  					 pattern = map(seus),
  					 iteration = "list"
  ),
  
  # Purpose: Track file path input/output for debranched clone trees file.
  tar_target(
  	debranched_clone_trees_file,
  	qpdf::pdf_combine(debranched_clone_tree_files, "results/debranched_clone_trees.pdf")
  ),
  
  # Purpose: Track file path input/output for final clone tree files.
  tar_target(final_clone_tree_files,
    save_clone_tree_from_path(final_seus, numbat_rds_files, large_clone_simplifications, label = "_clone_tree", legend = FALSE, horizontal = FALSE),
    pattern = map(final_seus),
    iteration = "list"
  ),
  
  # Purpose: Track file path input/output for final clone trees file.
  tar_target(
  	final_clone_trees_file,
  	qpdf::pdf_combine(final_clone_trees_files, "results/final_clone_trees.pdf")
  ),

  # Purpose: Build dependency target for heatmap collages.
  tar_target(heatmap_collages,
    plot_seu_marker_heatmap_all_resolutions(debranched_seus, numbat_rds_files, large_clone_simplifications),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  
  # Purpose: Build dependency target for heatmap collages 6p.
  tar_target(heatmap_collages_6p,
  					 plot_seu_marker_heatmap_all_resolutions(debranched_seus_6p, numbat_rds_files, large_clone_simplifications),
  					 pattern = map(debranched_seus_6p),
  					 iteration = "list"
  ),
  
  # Purpose: Build dependency target for annotated heatmap collages.
  tar_target(annotated_heatmap_collages,
    plot_seu_marker_heatmap(debranched_seus, cluster_orders, numbat_rds_files, large_clone_simplifications),
    pattern = map(debranched_seus),
    iteration = "list"
  ),

  # hypoxia ------------------------------
  tar_target(hypoxia_seus,
    load_and_save_hypoxia_score(filtered_seus),
    pattern = map(filtered_seus),
    iteration = "list"
  ),

  # Purpose: Generate plot set for hypoxia score plots.
  tar_target(
    hypoxia_score_plots,
    plot_hypoxia_score(hypoxia_seus, threshold = 0.5),
    pattern = map(hypoxia_seus),
    iteration = "list"
  ),

  # Purpose: Compute hypoxia-related outputs for seus low hypoxia.
  tar_target(seus_low_hypoxia,
    subset_seu_by_expression(hypoxia_seus, run_hypoxia_clustering = TRUE, hypoxia_expr = "hypoxia_score <= 0.5", slug = "hypoxia_low"),
    pattern = map(hypoxia_seus),
    iteration = "list"
  ),

  # Purpose: Compute hypoxia-related outputs for seus high hypoxia.
  tar_target(seus_high_hypoxia,
    subset_seu_by_expression(hypoxia_seus, run_hypoxia_clustering = TRUE, hypoxia_expr = "hypoxia_score > 0.5", slug = "hypoxia_high"),
    pattern = map(hypoxia_seus),
    iteration = "list"
  ),

  # Purpose: Compute hypoxia-related outputs for heatmap collages hypoxia.
  tar_target(heatmap_collages_hypoxia,
    plot_seu_marker_heatmap(hypoxia_seus, nb_paths = numbat_rds_files, clone_simplifications = large_clone_simplifications, tmp_plot_path = TRUE),
    pattern = map(hypoxia_seus),
    iteration = "list"
  ),

  # TODO: heatmap_collages_hypoxia_low, heatmap_collages_hypoxia_low_clusters,
  # and heatmap_collages_hypoxia_high were identical copies of heatmap_collages_hypoxia.
  # These need to be updated to use seus_low_hypoxia / seus_high_hypoxia respectively.
  # tar_target(heatmap_collages_hypoxia_low,
  #   plot_seu_marker_heatmap(seus_low_hypoxia, nb_paths = numbat_rds_files, clone_simplifications = large_clone_simplifications, tmp_plot_path = TRUE),
  #   pattern = map(seus_low_hypoxia),
  #   iteration = "list"
  # ),
  # tar_target(heatmap_collages_hypoxia_high,
  #   plot_seu_marker_heatmap(seus_high_hypoxia, nb_paths = numbat_rds_files, clone_simplifications = large_clone_simplifications, tmp_plot_path = TRUE),
  #   pattern = map(seus_high_hypoxia),
  #   iteration = "list"
  # ),

  tar_target(hypoxia_effect_plots,
    list(
      qpdf::pdf_combine(unlist(hypoxia_score_plots), "01_hypoxia_score_plots.pdf"),
      qpdf::pdf_combine(unlist(heatmap_collages_hypoxia), "02_hypoxia_heatmap.pdf")
      # TODO: add low/high hypoxia heatmaps once those targets are fixed
    ),
  ),
  
  # Purpose: Generate plot set for silhouette plots.
  tar_target(silhouette_plots,
  					 calc_silhouette(debranched_seus),
  					 pattern = map(debranched_seus),
  					 iteration = "list"
  					 ),
  
  # format = "file" so that changes to the RDS files on disk invalidate
  # downstream targets (e.g. factor heatmaps, diffex, clone pearls).
  tarchetypes::tar_files(debranched_seus_1q,
    c(
      "SRR13884246" = "output/seurat/SRR13884246_branch_5_filtered_seu_1q.rds",
      "SRR13884249" = "output/seurat/SRR13884249_filtered_seu_1q.rds",
      "SRR14800534" = "output/seurat/SRR14800534_filtered_seu_1q.rds",
      "SRR14800535" = "output/seurat/SRR14800535_filtered_seu_1q.rds",
      "SRR14800536" = "output/seurat/SRR14800536_filtered_seu_1q.rds"
    )
  ),

  # Purpose: Build dependency target for debranched seus 2p.
  tarchetypes::tar_files(debranched_seus_2p,
    c(
      "SRR13884246" = "output/seurat/SRR13884246_branch_5_filtered_seu_2p.rds",
      "SRR13884247" = "output/seurat/SRR13884247_branch_6_filtered_seu.rds",
      "SRR13884248" = "output/seurat/SRR13884248_filtered_seu_2p.rds",
      # "SRR13884249" = "output/seurat/SRR13884249_filtered_seu_2p.rds",
      # "SRR17960481" = "output/seurat/SRR17960481_filtered_seu.rds",
      "SRR17960484" = "output/seurat/SRR17960484_filtered_seu_2p.rds" # too many changes at once to conclude anything
    )
  ),

  # Purpose: Build dependency target for debranched seus 6p.
  tarchetypes::tar_files(debranched_seus_6p,
    c(
      "SRR13884247" = "output/seurat/SRR13884247_filtered_seu.rds",
      "SRR13884248" = "output/seurat/SRR13884248_filtered_seu_6p.rds",
      "SRR17960484" = "output/seurat/SRR17960484_filtered_seu_6p.rds",
      "SRR27187899" = "output/seurat/SRR27187899_filtered_seu.rds"
    )
  ),

  # Purpose: Build dependency target for debranched seus 16q.
  tarchetypes::tar_files(debranched_seus_16q,
    c(
      "SRR14800534" = "output/seurat/SRR14800534_filtered_seu_16q.rds",
      "SRR14800535" = "output/seurat/SRR14800535_filtered_seu_16q.rds",
      "SRR14800536" = "output/seurat/SRR14800536_filtered_seu_16q.rds"
    )
  ),
  
  # Purpose: Prepare Seurat object path(s) for scna seus.
  tar_target(scna_seus,
    unlist(list(
      "1q"  = debranched_seus_1q,
      "2p"  = debranched_seus_2p,
      "6p"  = debranched_seus_6p,
      "16q" = debranched_seus_16q
    ))
  ),

  # Populate SQLite metadata tables for all SCNA-stratified Seurat objects.
  # Upserts four tables: seurat_objects, cell_metadata, cluster_composition,
  # qc_metrics. Re-runs automatically when any debranched_seus_* file changes
  # (propagated via format = "file" on the upstream file-path targets).
  tar_target(seu_metadata_db,
    bulk_extract_seu_metadata(scna_seus, sqlite_path = "batch_hashes.sqlite"),
    deployment = "main"  # DB writes must run on the main process, not a worker
  ),

  # medium resolution (sample-specific SCNA collages) ------------------------------

  # TODO: collages_1q and hypoxia-stratified variants below reference
  # hypoxia_1q / hypoxia_1q_low / hypoxia_1q_high which are not defined in
  # this pipeline. Use hypoxia_seus_1q (defined below) or debranched_seus_1q.
  # Names were also inverted (hypoxia_low used hypoxia_1q_high and vice versa).
  # Uncomment and fix once the input targets are defined.

  # tar_target(collages_1q,
  #   plot_seu_marker_heatmap_by_scna(
  #     unlist(hypoxia_seus_1q), cluster_orders, numbat_rds_files,
  #     large_clone_simplifications, rb_scna_samples = rb_scna_samples,
  #     large_clone_comparisons = large_clone_comparisons, scna_of_interest = "1q"
  #   ),
  #   pattern = map(hypoxia_seus_1q),
  #   iteration = "list"
  # ),

  tar_target(collages_2p,
    plot_seu_marker_heatmap_by_scna(
      unlist(debranched_seus_2p), cluster_orders, numbat_rds_files, large_clone_simplifications,
      rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p"
    ),
    pattern = map(debranched_seus_2p),
    iteration = "list"
  ),
  # Sample-specific analyses of tumors with 2p+ subclones without integration (Fig. S4.9)
  tar_target(fig_s04_09,
    qpdf::pdf_combine(collages_2p, "results/fig_s04_09.pdf")
  ),

  # integrated 1q fig for all 1q+ samples (thesis fig 4.3)
  tar_target(fig_02,
    plot_integrated_1q_fig(cluster_orders)
  ),

  # Sample-specific analyses of tumors with 1q+ subclones after integration.
  tar_target(fig_s07,
    sample_specific_analyses_of_tumors_with_scna_subclones_after_integration(integrated_seus_1q, cluster_orders, "results/fig_s07.pdf")
  ),

  # Sample-specific analyses of tumors with 16q- subclones after integration.
  tar_target(fig_s12,
    sample_specific_analyses_of_tumors_with_scna_subclones_after_integration(integrated_seus_16q, cluster_orders, "results/fig_s12.pdf")
  ),

  # Sample-specific analyses of tumors with 2p+ subclones after integration.
  tar_target(fig_s04_08,
    sample_specific_analyses_of_tumors_with_scna_subclones_after_integration(integrated_seus_2p, cluster_orders, "results/fig_s04_08.pdf")
  ),

  # TODO: replace not_sure_what_this_does with the correct function name
  tar_target(fig_s11,
    not_sure_what_this_does(integrated_seus_2p, cluster_orders, "results/fig_s11.pdf")
  ),

  # Purpose: Generate output for fig s23.
  tar_target(fig_s23,
    not_sure_what_this_does(integrated_seus_6p, cluster_orders, "results/fig_s23.pdf")
  ),

  # Purpose: Generate output for fig s25.
  tar_target(fig_s25,
    plot_fig_s25(subtype_markers = subtype_markers)
  ),

  # Alternative resolutions for integrated 16q- analysis
  tar_target(fig_03,
    plot_fig_03(cluster_orders)
  ),

  # Purpose: Generate output for fig 04 07.
  tar_target(fig_04_07,
    plot_fig_04_07(cluster_orders)
  ),

  # Purpose: Generate output for fig 04.
  tar_target(fig_04,
    plot_fig_04_05(
      c("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds",
        "output/seurat/SRR13884247_branch_6_filtered_seu.rds"),
      corresponding_clusters_enrichments[[6]],
      plot_path = "results/fig_04.pdf",
      integrated_seu_paths = integrated_seus_2p,
      cluster_order = cluster_orders,
      large_clone_comparisons = large_clone_comparisons,
      scna_of_interest = "2p", width = 14, height = 10
    )
  ),

  # Purpose: Generate output for fig 05.
  tar_target(fig_05,
    plot_fig_04_05(
      "output/seurat/integrated_6p/integrated_seu_6p_duo.rds",
      corresponding_clusters_enrichments[[7]],
      integrated_seu_paths = integrated_seus_6p,
      plot_path = "results/fig_05.pdf",
      cluster_order = cluster_orders,
      large_clone_comparisons = large_clone_comparisons,
      scna_of_interest = "6p", width = 18, height = 10
    )
  ),

  # Purpose: Generate output for fig s08.
  tar_target(fig_s08,  # Differential expression between integrated 1q clones within clusters
    plot_fig_s08()
  ),

  # Purpose: Generate output for fig s10.
  tar_target(fig_s10,  # Differential expression between integrated 16q clones within clusters
    plot_fig_s10()
  ),

  # Purpose: Generate output for fig s20.
  tar_target(fig_s20,  # Differential expression between integrated 2p clones within clusters
    plot_fig_s20()
  ),

  # Purpose: Prepare integrated-analysis input/output for integrated seus 1q.
  tarchetypes::tar_files(integrated_seus_1q,
    c(
      "SRR13884249_filtered_seu.rds" = "output/seurat/integrated_1q/SRR13884249_integrated_1q_filtered_seu.rds",
      "SRR14800534_filtered_seu.rds" = "output/seurat/integrated_1q/SRR14800534_integrated_1q_filtered_seu.rds",
      "SRR14800535_filtered_seu.rds" = "output/seurat/integrated_1q/SRR14800535_integrated_1q_filtered_seu.rds",
      "SRR14800536_filtered_seu.rds" = "output/seurat/integrated_1q/SRR14800536_integrated_1q_filtered_seu.rds"
    )
  ),

  tarchetypes::tar_map(
    values = tibble::tibble(
      scna        = c("1q", "2p", "6p", "16q"),
      hypoxia_sym = rlang::syms(c("hypoxia_seus_1q", "hypoxia_seus_2p",
                                  "hypoxia_seus_6p", "hypoxia_seus_16q"))
    ),
    names = "scna",
    # Purpose: Compute hypoxia-related outputs for integrated seu low hypoxia.
    tar_target(integrated_seu_low_hypoxia,
      integration_by_scna_clones(
        hypoxia_sym, scna_of_interest = scna,
        large_clone_comparisons, filter_expr = "hypoxia_score <= 0.5"
      )
    )
  ),

  # Purpose: Compute hypoxia-related outputs for hypoxia low integrated heatmap collages.
  tar_target(hypoxia_low_integrated_heatmap_collages,
    plot_seu_marker_heatmap_integrated(integrated_seu_low_hypoxia_1q)
  ),

  # Purpose: Compute hypoxia-related outputs for hypoxia low integrated heatmap collages0.
  tar_target(hypoxia_low_integrated_heatmap_collages0,
    plot_integrated_1q_fig_low_hypoxia(integrated_seu_low_hypoxia_1q)
  ),

# tar_target(integrated_seu_1q_low_hypoxia_heatmaps,
# 					 plot_seu_marker_heatmap_by_scna(integrated_seu_1q_low_hypoxia, scna_of_interest = "1q")
# ),

  # tar_target(collages_6p,
  # 					 plot_seu_marker_heatmap_by_scna(unlist(debranched_seus_6p), cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "6p"),
  # 					 pattern = map(debranched_seus_6p),
  # 					 iteration = "list"
  # ),


  tarchetypes::tar_files(integrated_seus_16q,
    c(
      "SRR14800534_filtered_seu.rds" = "output/seurat/integrated_16q/SRR14800534_integrated_16q_filtered_seu.rds",
      "SRR14800535_filtered_seu.rds" = "output/seurat/integrated_16q/SRR14800535_integrated_16q_filtered_seu.rds",
      "SRR14800536_filtered_seu.rds" = "output/seurat/integrated_16q/SRR14800536_integrated_16q_filtered_seu.rds"
    )
  ),

  # Purpose: Prepare integrated-analysis input/output for integrated seus 2p.
  tarchetypes::tar_files(integrated_seus_2p,
    c(
      "SRR13884248_integrated_2p_filtered_seu.rds" = "output/seurat/integrated_2p/SRR13884248_integrated_2p_filtered_seu.rds",
      "SRR17960484_integrated_2p_filtered_seu.rds" = "output/seurat/integrated_2p/SRR17960484_integrated_2p_filtered_seu.rds"
    )
  ),

  # Purpose: Prepare integrated-analysis input/output for integrated seus 6p.
  tarchetypes::tar_files(integrated_seus_6p,
    c(
      "SRR13884248_integrated_6p_filtered_seu.rds" = "output/seurat/integrated_6p/SRR13884248_integrated_6p_filtered_seu.rds",
      "SRR17960484_integrated_6p_filtered_seu.rds" = "output/seurat/integrated_6p/SRR17960484_integrated_6p_filtered_seu.rds"
    )
  ),

  # Purpose: Generate output for fig 07a input.
  tar_target(fig_07a_input,
    find_diffex_bw_clones_for_each_cluster(integrated_seus_1q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "1q+"),
    pattern = map(integrated_seus_1q),
    iteration = "list"
  ),

  # Purpose: Generate output for fig 07a.
  tar_target(fig_07a,  # 1q+ cluster diffex after integration
    plot_fig_07_08(fig_07a_input, plot_path = "results/fig_07a.pdf", plot_title = "fig_07a: 1q+ cluster diffex after integration", height = 5, width = 4),
    iteration = "list"
  ),

  # Purpose: Generate output for fig 07b input.
  tar_target(fig_07b_input,
    find_diffex_bw_clones_for_each_cluster(debranched_seus_1q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "1q+"),
    pattern = map(debranched_seus_1q),
    iteration = "list"
  ),

  # Purpose: Generate output for fig 07b.
  tar_target(fig_07b,  # 1q+ cluster diffex without integration
    plot_fig_07_08(fig_07b_input, plot_path = "results/fig_07b.pdf", plot_title = "fig_07b: 1q+ cluster diffex without integration", height = 5, width = 4),
    iteration = "list"
  ),

  # Purpose: Generate output for fig s09.
  tar_target(fig_s09,  # Differential expression between integrated 1q clusters of interest
    plot_fig_s09()
  ),

  # Purpose: Generate output for fig 08a input.
  tar_target(fig_08a_input,
    find_diffex_bw_clones_for_each_cluster(integrated_seus_16q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "16q-"),
    pattern = map(integrated_seus_16q),
    iteration = "list"
  ),

  # Purpose: Generate output for fig 08a.
  tar_target(fig_08a,  # 16q- cluster diffex after integration
    plot_fig_07_08(fig_08a_input, plot_title = "fig_08a: 16q- cluster diffex after integration", plot_path = "results/fig_08a.pdf", height = 5, width = 6),
    iteration = "list"
  ),

  # Purpose: Generate output for fig 08b input.
  tar_target(fig_08b_input,
    find_diffex_bw_clones_for_each_cluster(debranched_seus_16q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "16q-"),
    pattern = map(debranched_seus_16q),
    iteration = "list"
  ),

  # Purpose: Generate output for fig 08b.
  tar_target(fig_08b,  # 16q- cluster diffex without integration
    plot_fig_07_08(fig_08b_input, plot_title = "fig_08b: 16q- cluster diffex without integration", plot_path = "results/fig_08b.pdf", p_adj_threshold = 0.05, height = 5, width = 6),
    iteration = "list"
  ),

# Purpose: Build dependency target for corresponding seus 2p.
tar_target(
	corresponding_seus_2p,
	c(
		# "output/seurat/SRR13884246_branch_5_filtered_seu_2p.rds", 
		# "SRR13884247_branch_6_filtered_seu" = "output/seurat/SRR13884247_branch_6_filtered_seu.rds", 
		"SRR13884248_filtered_seu_2p" = "output/seurat/SRR13884248_filtered_seu_2p.rds", 
		"SRR17960484_filtered_seu_2p.rds" = "output/seurat/SRR17960484_filtered_seu_2p.rds"
		# "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds"
	)
),

  # Purpose: Generate output for fig 09.
  tar_target(fig_09,
    plot_fig_09_10(
      corresponding_seus_2p, corresponding_seus,
      corresponding_clusters_diffex, corresponding_clusters_enrichments,
      recurrence_threshold = 3, plot_path = "results/fig_09.pdf",
      widths = rep(4, 3), heights = rep(8, 3),
      common_seus = c("SRR13884248_filtered_seu_2p.rds", "SRR17960484_filtered_seu_2p.rds")
    )
  ),

  # Purpose: Build dependency target for corresponding seus 6p.
  tar_target(corresponding_seus_6p,
    c(
      "SRR13884248_filtered_seu_6p.rds" = "output/seurat/SRR13884248_filtered_seu_6p.rds",
      "SRR17960484_filtered_seu_6p.rds" = "output/seurat/SRR17960484_filtered_seu_6p.rds"
    )
  ),

  # Purpose: Generate output for fig 10.
  tar_target(fig_10,
    plot_fig_09_10(
      corresponding_seus_6p, corresponding_seus,
      corresponding_clusters_diffex, corresponding_clusters_enrichments,
      recurrence_threshold = 2, plot_path = "results/fig_10.pdf",
      widths = rep(4, 3), heights = c(12, 4, 12),
      common_seus = c("SRR13884247_filtered_seu_6p.rds", "SRR17960484_filtered_seu_6p.rds")
    )
  ),

  # Purpose: Build dependency target for collages 6p.
  tar_target(collages_6p,
    plot_seu_marker_heatmap_by_scna(
      unlist(debranched_seus_6p), cluster_orders, numbat_rds_files, large_clone_simplifications,
      rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "6p"
    ),
    pattern = map(debranched_seus_6p),
    iteration = "list"
  ),

  # Sample-specific analyses of tumors with 16q- subclones without integration.
  tar_target(collages_16q,
    plot_seu_marker_heatmap_by_scna(
      unlist(debranched_seus_16q), cluster_orders, numbat_rds_files, large_clone_simplifications,
      rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "16q"
    ),
    pattern = map(debranched_seus_16q),
    iteration = "list"
  ),

  # Purpose: Build dependency target for collected scna collages.
  tar_target(collected_scna_collages,
  					 qpdf::pdf_combine(
  					 	unlist(list(collages_1q, collages_2p, collages_6p, collages_16q)), 
  					 	"results/fig_s17.pdf")
  					 ),

  # Purpose: Prepare Seurat object path(s) for corresponding seus.
  tar_target(corresponding_seus,
    list(
      # "output/seurat/SRR13884246_branch_5_filtered_seu_2p.rds",
      "output/seurat/SRR13884247_branch_6_filtered_seu.rds",
      "output/seurat/SRR13884248_filtered_seu_2p.rds",
      "output/seurat/SRR17960484_filtered_seu_2p.rds",
      "output/seurat/SRR13884247_filtered_seu.rds",
      "output/seurat/SRR17960484_filtered_seu_6p.rds",
      "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds",
      "output/seurat/integrated_6p/integrated_seu_6p_duo.rds",
      "output/seurat/SRR13884248_filtered_seu_6p.rds"
    )
  ),

# Purpose: Build dependency target for corresponding states dictionary.
tar_target(corresponding_states_dictionary,
					 make_corresponding_states_dictionary(tibble::tribble(
					 	~file_name, ~w_scna, ~wo_scna, ~scna_of_interest,
					 	# # "SRR13884246_branch_5_filtered_seu_2p.rds", "g1_1-g1_0-g1_3-g1_11", "g1_2-g1_5", "2p",
					 	# # "SRR13884246_branch_5_filtered_seu_2p.rds",     "s_g2_7-s_9",      "s_g2_8", "2p",
					 	
					 	"SRR13884247_branch_6_filtered_seu.rds", "g1_3-g1_6", "g1_5", "2p",
					 	# "SRR13884247_branch_6_filtered_seu.rds", "g1_3", "g1_5", "2p",
					 	# "SRR13884247_branch_6_filtered_seu.rds", "g1_6", "g1_5", "2p",
					 	
					 	"SRR13884248_filtered_seu_2p.rds",   "g1_1",    "g1_6-g1_7", "2p",
					 	# "SRR13884248_filtered_seu_2p.rds",   "g1_1",    "g1_6", "2p",
					 	# "SRR13884248_filtered_seu_2p.rds",   "g1_1",    "g1_7", "2p",
					 	
					 	"SRR17960484_filtered_seu_2p.rds",     "g1_3",    "g1_1-g1_5", "2p",
					 	# "SRR17960484_filtered_seu_2p.rds",     "g1_3",    "g1_1", "2p",
					 	# "SRR17960484_filtered_seu_2p.rds",     "g1_3",    "g1_5", "2p",
					 	# "SRR17960484_filtered_seu_2p.rds",     "s_6",      "s_7", "2p",
					 	
					 	"SRR13884247_filtered_seu.rds",   "g1_5-g1_7",    "g1_4-g1_0", "6p",
					 	
					 	"SRR17960484_filtered_seu_6p.rds",     "g1_0",      "g1_1", "6p",
					 	# "SRR17960484_filtered_seu_6p.rds",     "s_4",      "s_6", "6p",
					 	
					 	"seurat_2p_integrated_duo.rds",     "g1_1-g1_4",      "g1_9", "2p",
					 	"seurat_2p_integrated_duo.rds",     "g1_4",      "g1_9", "2p",
					 	"seurat_2p_integrated_duo.rds",     "g1_1",      "g1_9", "2p",
					 	"seurat_2p_integrated_duo.rds",     "g1_1",      "g1_4", "2p",
					 	
					 	"integrated_seu_6p_duo.rds",     "g1_4",      "g1_3", "6p",
					 	
					 	"SRR13884248_filtered_seu_6p.rds",   "g1_0",    "g1_2", "6p"
					 ))
),

  # corresponding states: 2p ------------------------------

  tar_target(states_dictionary_2p,
    make_corresponding_states_dictionary(tibble::tribble(
      ~file_name,                       ~w_scna,    ~wo_scna, ~scna_of_interest,
      "seurat_2p_integrated_duo.rds",   "g1_1-g1_4", "g1_9",  "2p",
      "seurat_2p_integrated_duo.rds",   "g1_4",      "g1_9",  "2p",
      "seurat_2p_integrated_duo.rds",   "g1_1",      "g1_9",  "2p",
      "seurat_2p_integrated_duo.rds",   "g1_1",      "g1_4",  "2p"
    ))
  ),

  # Purpose: Compute differential expression results for corresponding clusters diffex 2p.
  tar_target(corresponding_clusters_diffex_2p,
    find_diffex_clusters_between_corresponding_states(
      "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds",
      states_dictionary_2p, large_clone_comparisons,
      numbat_rds_files = numbat_rds_files, location = "all"
    )
  ),

  # Purpose: Build dependency target for corresponding clusters volcanos 2p.
  tar_target(corresponding_clusters_volcanos_2p,
    plot_corresponding_clusters_diffex_volcanos(
      corresponding_clusters_diffex_2p,
      "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds"
    )
  ),

  # Purpose: Run enrichment analysis for corresponding clusters enrichments 2p.
  tar_target(corresponding_clusters_enrichments_2p,
    plot_corresponding_enrichment(
      corresponding_clusters_diffex_2p,
      "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds"
    )
  ),

  # corresponding states: 6p ------------------------------

  tar_target(states_dictionary_6p,
    make_corresponding_states_dictionary(tibble::tribble(
      ~file_name,                        ~w_scna,     ~wo_scna,    ~scna_of_interest,
      "SRR13884247_filtered_seu.rds",    "g1_0-g1_1", "g1_4-g1_5", "6p",
      "SRR17960484_filtered_seu_6p.rds", "g1_0",      "g1_1",      "6p"
    ))
  ),

  # Purpose: Prepare Seurat object path(s) for corresponding state 6p seus.
  tar_target(corresponding_state_6p_seus,
    list(
      "output/seurat/SRR13884247_filtered_seu.rds",
      "output/seurat/SRR17960484_filtered_seu_6p.rds"
    )
  ),

  # Purpose: Compute differential expression results for corresponding clusters diffex 6p.
  tar_target(corresponding_clusters_diffex_6p,
    find_diffex_clusters_between_corresponding_states(
      unlist(corresponding_state_6p_seus), states_dictionary_6p,
      large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all"
    ),
    pattern = map(corresponding_state_6p_seus, states_dictionary_6p),
    iteration = "list"
  ),

  # Purpose: Build dependency target for corresponding clusters volcanos 6p.
  tar_target(corresponding_clusters_volcanos_6p,
    plot_corresponding_clusters_diffex_volcanos(
      corresponding_clusters_diffex_6p, unlist(corresponding_state_6p_seus)
    ),
    pattern = map(corresponding_clusters_diffex_6p, corresponding_state_6p_seus),
    iteration = "list"
  ),

  # Purpose: Run enrichment analysis for corresponding clusters enrichments 6p.
  tar_target(corresponding_clusters_enrichments_6p,
    plot_corresponding_enrichment(
      corresponding_clusters_diffex_6p, unlist(corresponding_state_6p_seus)
    ),
    pattern = map(corresponding_clusters_diffex_6p, corresponding_state_6p_seus),
    iteration = "list"
  ),

  # corresponding states: combined (2p + 6p) ------------------------------

  tar_target(corresponding_clusters_diffex,
    find_diffex_clusters_between_corresponding_states(
      unlist(corresponding_seus), corresponding_states_dictionary,
      large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all"
    ),
    pattern = map(corresponding_seus, corresponding_states_dictionary),
    iteration = "list"
  ),

  # Purpose: Build dependency target for corresponding clusters volcanos.
  tar_target(corresponding_clusters_volcanos,
    plot_corresponding_clusters_diffex_volcanos(corresponding_clusters_diffex, unlist(corresponding_seus)),
    pattern = map(corresponding_clusters_diffex, corresponding_seus),
    iteration = "list"
  ),

  # Purpose: Build dependency target for corresponding clusters heatmaps.
  tar_target(corresponding_clusters_heatmaps,
    plot_corresponding_clusters_diffex_heatmaps(
      corresponding_clusters_diffex, unlist(corresponding_seus),
      corresponding_states_dictionary, large_clone_comparisons,
      numbat_rds_files = numbat_rds_files, location = "all"
    ),
    pattern = map(corresponding_clusters_diffex, corresponding_seus, corresponding_states_dictionary),
    iteration = "list"
  ),

  # Purpose: Run enrichment analysis for corresponding clusters enrichments.
  tar_target(corresponding_clusters_enrichments,
    plot_corresponding_enrichment(corresponding_clusters_diffex, unlist(corresponding_seus)),
    pattern = map(corresponding_clusters_diffex, corresponding_seus),
    iteration = "list"
  ),

  # Purpose: Compute clone-related output for clone cc plots by scna 1q.
  tar_target(clone_cc_plots_by_scna_1q,
    clone_cc_plots_by_scna(debranched_seus_1q, scna_of_interest = "1q", large_clone_comparisons = large_clone_comparisons)
  ),

  # Purpose: Compute clone-related output for clone cc plots by scna 16q.
  tar_target(clone_cc_plots_by_scna_16q,
    clone_cc_plots_by_scna(debranched_seus_16q, scna_of_interest = "16q", large_clone_comparisons = large_clone_comparisons)
  ),

  
  # Purpose: Compute differential expression results for diffex 2p g1.
  tar_target(diffex_2p_g1, # 2p+ diffex; diffex 2p+
  					 find_diffex_clones_in_phase(corresponding_seus_2p, phase = "g1", scna_of_interest = "2p", numbat_rds_files, large_clone_comparisons, location = "all"),
  					 pattern = map(corresponding_seus_2p),
  					 iteration = "list"
  ),

# Purpose: Run enrichment analysis for enrichment 2p g1.
tar_target(enrichment_2p_g1,
					 compare_enrichment(diffex_2p_g1)
),
  
  # Purpose: Compute differential expression results for diffex 6p g1.
  tar_target(diffex_6p_g1,
  					 find_diffex_clones_in_phase(debranched_seus_6p, phase = "g1", scna_of_interest = "6p", numbat_rds_files, large_clone_comparisons, location = "all"),
  					 pattern = map(debranched_seus_6p),
  					 iteration = "list"
  ),

# Purpose: Run enrichment analysis for enrichment 6p g1.
tar_target(enrichment_6p_g1,
					 compare_enrichment(diffex_6p_g1)
),
  
  # Purpose: Prepare Seurat object path(s) for integrated seus.
  tar_target(integrated_seus,
  					 list(
  					 	"diploid_v_16q" = "output/seurat/SRR1_diploid_v_16q_filtered_seu.rds",
  					 	"16q_v_16q-1q" = "output/seurat/SRR1_16q_v_16q-1q_filtered_seu.rds",
  					 "diploid_v_16q-1q" = "output/seurat/SRR1_diploid_v_16q-1q_filtered_seu.rds"
  					 )
  					 ),
  
  # Purpose: Prepare integrated-analysis input/output for annotated integrated heatmap collages.
  tar_target(annotated_integrated_heatmap_collages,
  					 plot_seu_marker_heatmap_integrated(integrated_seus),
  					 pattern = map(integrated_seus),
  					 iteration = "list"
  ),

  # Purpose: Build dependency target for collage compilation.
  tar_target(
    collage_compilation,
    qpdf::pdf_combine(
    	annotated_heatmap_collages, 
    	"results/heatmap_collages.pdf")
  ),
  
  # Purpose: Build dependency target for collage compilation all resolutions.
  tar_target(
  	collage_compilation_all_resolutions,
  	qpdf::pdf_combine(heatmap_collages, "results/heatmap_collages_all_resolutions.pdf")
  ),
  
  # Purpose: Track file path input/output for divergent cluster file.
  tarchetypes::tar_file(divergent_cluster_file, "data/clustree_divergent_clusters.csv"),
  # Purpose: Build dependency target for clustree tables.
  tar_target(
    clustree_tables,
    pull_clustree_tables(clustrees, divergent_cluster_file)
  ),
  # Purpose: Compute differential expression results for clustree diffexes.
  tar_target(clustree_diffexes,
    find_all_diffex_from_clustree(clustree_tables, debranched_seus, clone_comparisons = large_clone_comparisons),
    pattern = map(clustree_tables),
  ),
  # Purpose: Build dependency target for clustree cis changes.
  tar_target(
    clustree_cis_changes,
    find_candidate_cis_in_clustree_diffexes(clustree_diffexes)
    # pattern = map(clustree_diffexes),
  ),
  # Purpose: Build dependency target for clustree trans changes.
  tar_target(
    clustree_trans_changes,
    find_candidate_trans_in_clustree_diffexes(clustree_diffexes)
    # pattern = map(clustree_diffexes),
  ),
  # Purpose: Build dependency target for clustree all changes.
  tar_target(
    clustree_all_changes,
    find_candidate_all_in_clustree_diffexes(clustree_diffexes)
    # pattern = map(clustree_diffexes),
  ),
  # Purpose: Track file path input/output for large numbat heatmap files.
  tar_target(large_numbat_heatmap_files, retrieve_numbat_plot_type(large_numbat_pdfs, "panel_2.pdf")),

  # Purpose: Create PDF artifact(s) for large numbat heatmap pdf.
  tar_target(large_numbat_heatmap_pdf, qpdf::pdf_combine(large_numbat_heatmap_files, "results/large_numbat_heatmaps.pdf")),

  # Purpose: Create PDF artifact(s) for large numbat sample pdfs.
  tar_target(large_numbat_sample_pdfs,
    reroute_done_to_results_pdf(numbat_rds_files, "_large"),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),
  # Purpose: Build dependency target for filtering timelines.
  tar_target(
    filtering_timelines,
    qpdf::pdf_combine(dir_ls("results", regexp = ".*SRR[0-9]*_filtering_timeline.pdf"), "results/filtering_timelines.pdf")
  ),

  # Purpose: Compute differential expression results for num diffex clone scna tally.
  tar_target(
    num_diffex_clone_scna_tally,
    tally_num_diffex(oncoprint_input_by_scna_unfiltered)
  ),
  
  # Purpose: Build dependency target for rb scna samples.
  tar_target(
  	rb_scna_samples,
  	list(
  	"1q" = c("SRR13884246", "SRR13884249", "SRR14800534", "SRR14800535", "SRR14800536"),
  	"2p"= c("SRR13884246", "SRR13884247", "SRR13884248", "SRR13884249", "SRR17960481", "SRR17960484"),
  	"6p"= c("SRR13884247", "SRR13884248", "SRR17960484"),
  	"16q" = c("SRR14800534", "SRR14800535", "SRR14800536")
  )),

  # hypoxia_seus_1q/2p/6p/16q — subset seus_low_hypoxia to each SCNA's samples
  tarchetypes::tar_map(
    values = scna_map_values[, "scna", drop = FALSE],
    names  = "scna",
    # Purpose: Prepare Seurat object path(s) for hypoxia seus.
    tar_target(hypoxia_seus,
      str_subset(unlist(seus_low_hypoxia), str_c(rb_scna_samples[[scna]], collapse = "|"))
    )
  ),

  # Purpose: Prepare oncoprint data/plots for unfiltered oncoprint input by scna.
  tar_target(
  	unfiltered_oncoprint_input_by_scna,
  	make_oncoprint_diffex(large_filter_expressions, cluster_dictionary, debranched_ids, cis_diffex_clones, trans_diffex_clones, all_diffex_clones, large_clone_comparisons, rb_scna_samples, n_slice = 20)
  ),

  # Purpose: Prepare oncoprint data/plots for oncoprint input by scna.
  tar_target(
  	oncoprint_input_by_scna,
    filter_oncoprint_diffex(unfiltered_oncoprint_input_by_scna, oncoprint_settings)
  ),
  
  # Purpose: Track file path input/output for oncoprint settings file.
  tar_file(oncoprint_settings_file, "data/oncoprint_settings.tsv"),
  
  # Purpose: Track file path input/output for oncoprint settings per cluster file.
  tar_file(oncoprint_settings_per_cluster_file, "data/oncoprint_settings_by_cluster.tsv"),
  
  # Purpose: Prepare oncoprint data/plots for oncoprint settings.
  tar_target(oncoprint_settings, read_tsv(oncoprint_settings_file)),
  
  # Purpose: Prepare oncoprint data/plots for oncoprint settings per cluster.
  tar_target(oncoprint_settings_per_cluster, read_tsv(oncoprint_settings_per_cluster_file)),

  # Purpose: Prepare oncoprint data/plots for unfiltered oncoprint input by scna for each cluster.
  tar_target(
    unfiltered_oncoprint_input_by_scna_for_each_cluster,
    make_oncoprint_diffex(large_filter_expressions, cluster_dictionary, debranched_ids, cis_diffex_clones_for_each_cluster, trans_diffex_clones_for_each_cluster, all_diffex_clones_for_each_cluster, large_clone_comparisons, rb_scna_samples, by_cluster = TRUE, n_slice = 20)
  ),
  
  # Purpose: Prepare oncoprint data/plots for oncoprint input by scna for each cluster.
  tar_target(
  	oncoprint_input_by_scna_for_each_cluster,
  	filter_oncoprint_diffex(unfiltered_oncoprint_input_by_scna_for_each_cluster, oncoprint_settings_by_cluster)
  ),
  
  # Purpose: Run enrichment analysis for oncoprint enrich clones gobp.
  tar_target(
    oncoprint_enrich_clones_gobp,
    enrich_oncoprints(large_filter_expressions,
      cluster_dictionary,
      debranched_ids,
      cis_diffex_clones,
      trans_diffex_clones,
      all_diffex_clones,
      large_clone_comparisons,
      gene_set = "gobp"
    )
  ),
  
  # hallmark ------------------------------
  
  tar_target(
    oncoprint_enrich_clones_hallmark,
    enrich_oncoprints(large_filter_expressions,
      cluster_dictionary,
      debranched_ids,
      cis_diffex_clones,
      trans_diffex_clones,
      all_diffex_clones,
      large_clone_comparisons,
      gene_set = "hallmark"
    )
  ),
  # Purpose: Build dependency target for rod rich samples.
  tar_target(rod_rich_samples,
    score_samples_for_rod_enrichment(numbat_rds_files),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),

  # Purpose: Generate output for table 06.
  tar_target(table_06,
    rod_rich_samples
  ),

  # Purpose: Build dependency target for celltype rich samples.
  tar_target(celltype_rich_samples,
    score_samples_for_celltype_enrichment(unfiltered_seus, final_seus, celltype_markers),
    pattern = map(unfiltered_seus),
    iteration = "list"
  ),
  # Purpose: Run enrichment analysis for oncoprint enrich clones plots gobp.
  tar_target(
    oncoprint_enrich_clones_plots_gobp,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_gobp,
      cis_plot_file = "results/enrichment/cis_enrichment_plots_by_clone_gobp.pdf",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_clone_gobp.pdf",
      cis_table_file = "results/enrichment/cis_enrichment_tables_by_clone_gobp.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_clone_gobp.xlsx"
    )
  ),
  # Purpose: Run enrichment analysis for oncoprint enrich clones plots hallmark.
  tar_target(
    oncoprint_enrich_clones_plots_hallmark,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_hallmark,
      cis_plot_file = "results/enrichment/cis_enrichment_plots_by_clone_hallmark.pdf",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_clone_hallmark.pdf",
      cis_table_file = "results/enrichment/cis_enrichment_tables_by_clone_hallmark.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_clone_hallmark.xlsx"
    )
  ),
  # Purpose: Build dependency target for subtype markers.
  tar_target(
    subtype_markers,
    pull_subtype_genes(supp_excel = "data/Liu_Radvanyi_2022_supp_data/41467_2021_25792_MOESM6_ESM.xlsx")
  ),
  # Purpose: Build dependency target for liu lu supp data.
  tar_target(
    liu_lu_supp_data,
    read_liu_lu_supp_tables()
  ),
  # Purpose: Build dependency target for mps.
  tar_target(
    mps,
    read_mps("/dataVolume/storage/Homo_sapiens/3ca/ITH_hallmarks/MPs_distribution/MP_list.RDS")
  ),
  # Purpose: Build dependency target for stemness markers.
  tar_target(
    stemness_markers,
    pull_stem_cell_markers()
  ),

  # enrichment by cluster msigdb hallmark
  tar_target(
    oncoprint_enrich_clusters_hallmark,
    enrich_oncoprints_clusters(large_filter_expressions,
      cluster_dictionary,
      debranched_ids,
      cis_diffex_clones_for_each_cluster,
      trans_diffex_clones_for_each_cluster,
      large_clone_comparisons,
      gene_set = "hallmark"
    )
  ),

  # enrichment plots and tables by cluster msigdb hallmark
  tar_target(
    oncoprint_enrich_clusters_plots_hallmark,
    compile_cis_trans_enrichment_recurrence_by_cluster(oncoprint_enrich_clusters_hallmark,
    																				cis_plot_file = "results/enrichment/cis_enrichment_plots_by_cluster_hallmark.pdf",
    																				trans_plot_file = "results/enrichment/trans_enrichment_plots_by_cluster_hallmark.pdf",
    																				cis_table_file = "results/enrichment/cis_enrichment_tables_by_cluster_hallmark.xlsx",
    																				trans_table_file = "results/enrichment/trans_enrichment_tables_by_cluster_hallmark.xlsx",
    																				by_cluster = TRUE
    )
  ),
  
  # Purpose: Prepare oncoprint data/plots for oncoprint input by region.
  tar_target(oncoprint_input_by_region,
    inspect_oncoprints(oncoprint_input_by_scna$cis, oncoprint_input_by_scna$trans, oncoprint_input_by_scna$all)
  ),

	# Purpose: Build dependency target for resolution dictionary.
	tar_target(resolution_dictionary,
						 list("SRR13884246_branch_5_filtered_seu_1q.rds"= "0",
						 		 "SRR13884249_filtered_seu_1q.rds"= "0",
						 		 "SRR14800534_filtered_seu_1q.rds"= "0",
						 		 "SRR14800535_filtered_seu_1q.rds"= "0",
						 		 "SRR14800536_filtered_seu_1q.rds"= "0",
						 		 "SRR13884246_branch_5_filtered_seu_2p.rds"= "0",
						 		 "SRR13884247_branch_6_filtered_seu.rds"= "-1",
						 		 "SRR13884248_filtered_seu_2p.rds"= "-1",
						 		 "SRR13884249_filtered_seu_2p.rds"= "-1",
						 		 "SRR17960481_filtered_seu.rds"= "0",
						 		 "SRR17960484_filtered_seu_2p.rds"= "-1",
						 		 "SRR13884247_filtered_seu.rds"= "-1",
						 		 "SRR13884248_filtered_seu_6p.rds"= "-1",
						 		 "SRR17960484_filtered_seu_6p.rds"= "0",
						 		 "SRR27187899_filtered_seu.rds"= "0",
						 		 "SRR14800534_filtered_seu_16q.rds"= "-1",
						 		 "SRR14800535_filtered_seu_16q.rds"= "0",
						 		 "SRR14800536_filtered_seu_16q.rds"= "0"
						 )
						 ),

	# Purpose: Prepare Seurat object path(s) for chosen resolution seus.
	tar_target(chosen_resolution_seus,
						 assign_designated_phase_clusters(scna_seus, cluster_orders, resolution_dictionary),
						 pattern = map(scna_seus),
						 iteration = "list"
	),
  
  # Purpose: Generate plot set for oncoprint plots.
  tar_target(oncoprint_plots, 
    make_oncoprint_plots(oncoprint_input_by_scna, debranched_clone_trees, oncoprint_settings, label = "_by_clone")
  ),
  
  # Purpose: Prepare oncoprint data/plots for oncoprint plots by cluster.
  tar_target(
  	oncoprint_plots_by_cluster,
  	make_oncoprint_plots(oncoprint_input_by_scna_for_each_cluster$cis, oncoprint_input_by_scna_for_each_cluster$trans, oncoprint_input_by_scna_for_each_cluster$all, debranched_clone_trees, oncoprint_settings_by_cluster, label = "_by_cluster", p_val_threshold = 1)
  ),
  
  # Purpose: Generate plot set for stachelek score plots.
  tar_target(stachelek_score_plots,
    score_stachelek(final_seus, oncoprint_input_by_scna),
    pattern = map(final_seus),
    iteration = "list"
  ),
  
  # pseudobulks -------------------------------
  
  tar_target(whole_pseudobulks,
    score_whole_pseudobulks(numbat_rds_files, subtype_markers),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),
  # Purpose: Build dependency target for unfiltered derived pseudobulk subtype scores.
  tar_target(
    unfiltered_derived_pseudobulk_subtype_scores,
    derive_pseudobulk_subtype_scores(unfiltered_seus),
  ),
  # Purpose: Build dependency target for filtered derived pseudobulk subtype scores.
  tar_target(
    filtered_derived_pseudobulk_subtype_scores,
    derive_pseudobulk_subtype_scores(final_seus),
  ),
    # Purpose: Build dependency target for pipeline notification.
    tar_target(
    pipeline_notification,
    {
      # depend on figures_and_tables so this runs after the main outputs are ready
      invisible(figures_and_tables)
      tryCatch({
        send_pipeline_notification(subject = "Targets pipeline completed", body = NULL)
        TRUE
      }, error = function(e) {
        warning("Failed to send pipeline notification: ", e$message)
        FALSE
      })
    },
  )

  # end of plan------------------------------
)

# =============================================================================
# Sample-level analysis notes (2023-03-07)
# =============================================================================
#
# INCLUDED SAMPLES (interesting_samples):
#   wu  SRR13884242: likely 16q GT; f31 sample; cluster 1 (prolif.) enriched GT 2;
#                    can attribute proliferation to 16q- acquisition in GT 2
#   wu  SRR13884243: likely 16q GT; biological replicate of SRR13884242;
#                    cluster 4 analogue of SRR13884242 c3
#   wu  SRR13884247: likely 16q GT; no GSEA diffex output; c5 marker DLK1 in
#                    miRNA cluster with MEG3
#   wu  SRR13884249: possible 1q/2p/16q GT; GT 1 (only 16q-) does not contribute
#                    to proliferating c2/c4 (TOP2A high); split by phase
#   yang SRR14800534: GT1 lacks 16q-; does not contribute to proliferating c2/c4;
#                     GT 2 lacks 1q—can't identify tx distinction b/w GT2/3
#   yang SRR14800535: 16q GT; GT 1 decreased in c1 (high TOP2A)
#   yang SRR14800536: possible 16q GT; GT1 decreased in c2/c3 (high G2M markers)
#   yang SRR14800540: clear 16q GT; three GTs with SCNAs; each progressively more
#                     proliferative; GT5/c4 markers: C1QA/B, CD74, HLA genes
#   yang SRR14800541: clear 1q/6p/16q GTs; GT1 not proliferating (no c2 G2M contribution)
#   yang SRR14800543: possible minor 16q/1q GT; dual/individual 1q+16q contribution
#                     with 13q CNLOH; MYCN marks clusters
#   field SRR17960481: 6p interesting; c2/c4 notable; c2 has mito genes; likely
#                      clonal 6p with PRs and stressed cells in GT 1
#   field SRR17960484: c1 enriched for wt GT 1; also has high Xist expression
#
# MAYBE (not in interesting_samples):
#   SRR13884240: possible 2p GT
#   SRR13884241: possible 1q GT
#   SRR13884244: possible 1q GT
#   SRR13884245: possible 1q GT
#   SRR13884246: possible 16q GT
#   SRR17960480: possible minor 16q- GT
#   SRR14800539: possible 16q GT
#
# EXCLUDED:
#   SRR14800537: possible 16q GT — excluded
#   SRR17960482: too complicated
#   SRR13884248: clear 6p (missing possible 2p in expression)
#   SRR17960483: cluster 6 maybe interesting
