## Load your packages, e.g. library(targets).
source("./packages.R")
source("./functions.R")
# plan(callr)
# debug(make_numbat_heatmaps)
# debug(make_numbat_plot_files)
# debug(diffex_cells)
# debug(filter_numbat_cells)
# debug(diffex_by_cluster)
# debug(enrich_diffex_by_cluster)
# debug(enrich_by_cluster)
# debug(find_diffex_bw_clusters_for_each_clone)

tar_option_set(memory = "transient", garbage_collection = TRUE)

output_plot_extensions <- c(
  "dimplot.pdf",
  "merged_marker.pdf",
  "combined_marker.pdf",
  "phylo_probability.pdf",
  "clone_distribution.pdf"
)

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

## tar_plan supports drake-style targets and also tar_target()
tar_plan(
  tar_target(
    paper_figures,
    copy_paper_data()
  ),

  # tarchetypes::tar_files(all_done_files, retrieve_done_files("output/numbat/"), format = "file"),

  tarchetypes::tar_files(large_done_files, retrieve_done_files("output/numbat_sridhar/", analyzed_samples), format = "file"),

  tar_target(
    celltype_markers,
    list(
      rods = c("PDE6A", "RHO", "NR2E3"),
      mature.cone = c("OPN1LW", "OPN1MW", "OPN1SW"),
      Müller.glia = c("RLBP1", "APOE", "CLU"),
      retinal.astrocytes = c("GFAP"),
      microglia = c("HLA-DPA1", "HLA-DRA", "C1QA"),
      bipolar.cells = c("VSX2", "TMEM215", "VSX1"),
      retinal.ganglion.cells = c("SNCG", "SLC17A6", "RBPMS"),
      amacrine.cells = c("CALB1", "CHAT", "GAD2")
    )
  ),

  # setup ------------------------------

  tar_target(cluster_dictionary, read_cluster_dictionary("data/cluster_dictionary.csv")),
  tar_target(hallmark_gene_sets, msigdbr(species = "Homo sapiens", category = "H")),
  tar_target(cluster_comparisons, list(
    "SRR13884240" = c(2, 3, 4, 5, 6, 7), # 2p
    "SRR13884241" = c(2, 3, 4, 5), # 2p
    "SRR13884242" = c(2, 3, 4), # 16q
    "SRR13884243" = c(2, 3), # 16q and 1q
    "SRR13884244" = c(2, 3), # nothing; should be 1q
    "SRR13884245" = c(2), # nothing; 19 i guess
    "SRR13884246" = c(2), # nothing; too complicated
    "SRR13884247" = c(2, 3, 4), # 6p; 1q gain only affects GT 4
    "SRR13884248" = c(2), # 6p
    "SRR13884249" = c(2, 3), # 1q; 2p only GT 3
    "SRR14800534" = c(3), # 1q only GT 3; maybe reconsider
    "SRR14800535" = c(2, 3), # 1q; maybe reconsider
    "SRR14800536" = c(2, 3, 4), # 16q
    "SRR14800537" = c(2, 3, 4, 5), # 16q
    "SRR14800539" = c(2), # nothing; 18 i guess
    "SRR14800540" = c(2), # 16q; just 2
    "SRR14800541" = c(2), # 16q is only in GT 2; 2,3,4,5 have 1q/2p/6p
    "SRR14800543" = c(3), # 1q/16q only in GT 3
    "SRR17960480" = c(2), # nothing; only 7p in GT 2
    "SRR17960481" = c(2, 3, 4), # 1q/2p/6p
    "SRR17960482" = c(2), # nothing; too complicated
    "SRR17960483" = c(2, 3), # nothing; clonal 1q/16q
    "SRR17960484" = c(4) # 16q only in GT 4; 2,3,4 have 1q; interesting sample
  )),
  tar_target(large_clone_comparisons, list(
    "SRR13884242" = list("2_v_1_1q+_16q-" = c("1b", "16b"), "3_v_2_8p+_11p+" = c("8a", "11a")), # 16q
    "SRR13884243" = list("2_v_1_1q+_16q-" = c("1b", "16b"), "3_v_2_8p+_11p+" = c("8a", "11a")), # 16q and 1q
    "SRR13884247" = list("2_v_1_6p+" = c("6c"), "3_v_2_17q+" = c("17d"), "4_v_3_10q+" = c("10b"), "5_v_3_1p-" = c("1a"), "6_v_2_2p+" = c("2b")), # 6p; 1q gain only affects GT 4
    "SRR13884249" = list("2_v_1_1q+" = c("1b"), "3_v_2_2p+" = c("2a")), # 1q; 2p only GT 3
    "SRR14800534" = list("2_v_1_16q-" = c("16c"), "3_v_2_1q+" = c("1b")), # 1q only GT 3; maybe reconsider
    "SRR14800535" = list("2_v_1_16q-" = c("16b"), "3_v_2_1q+" = c("1b")), # 1q only GT 3; maybe reconsider
    "SRR14800536" = list("2_v_1_16q-" = c("16b"), "3_v_2_1q+" = c("1b"), "4_v_3_15q+" = c("15b")), # 16q
    "SRR14800540" = list("2_v_1_6p+_16q-" = c("6b", "16e"), "3_v_2_1q+_7q+" = c("1d", "1e", "7b")), # 16q; just 2
    "SRR14800541" = list("2_v_1_1q+_6p+" = c("1f", "6a"), "3_v_2_16q-" = c("16d")), # 16q is only in GT 2; 2,3,4,5 have 1q/2p/6p
    "SRR14800543" = list("2_v_1_13loh" = c("13a"), "3_v_2_1q+_16q-" = c("1c", "16c")), # 1q/16q only in GT 3
    "SRR17960481" = list("2_v_1_1q+_6p+" = c("1b", "6a", "6b"), "3_v_2_2p+" = c("2a")), # 1q/2p/6p
    "SRR17960484" = list("2_v_1_1q+" = c("1c"), "3_v_2_6p+" = c("6a"), "4_v_3_2p+" = c("2a")) # 16q only in GT 4; 2,3,4 have 1q; interesting sample
  )),
  tar_target(large_clone_simplifications, list(
    "SRR13884242" = c("1q+" = "1b", "2p+" = "2b", "8p-" = "8a", "11p-" = "11a", "12p-" = "12a", "16q-" = "16b"), # 16q
    "SRR13884243" = c("1q+" = "1b", "2p+" = "2b", "8p-" = "8a", "11p-" = "11a", "12p-" = "12a", "16q-" = "16b", "8p-" = "8a", "11p-" = "11a"),
    "SRR13884247" = c("1p-" = "1a", "2p+" = "2b", "6p+" = "6c", "10q+" = "10b", "11q-" = "11b", "13q-" = "13b", "17q+" = "17d", "20p-" = "20a"),
    "SRR13884249" = c("1q+" = "1b", "2p+" = "2a", "13cnloh" = "13b", "16q-" = "16b"),
    "SRR14800534" = c("1q+" = "1b", "16q-" = "16c"), # 1q only GT 3; maybe reconsider
    "SRR14800535" = c("1q+" = "1b", "13cnloh" = "13a", "16q-" = "16b"), # 1q only GT 3; maybe reconsider
    "SRR14800536" = c("1q+" = "1b", "13-" = "13b", "15q+" = "15b", "16q-" = "16b", "19q-" = "19c"), # 16q
    "SRR14800540" = c("1q+" = "1d", "3q-" = "3b", "6q+" = "6b", "7q+" = "7b", "11q-" = "11b", "12p-" = "12f", "13-" = "13b", "16q-" = "16e", "20q+" = "20b"),
    "SRR14800541" = c("1q+" = "1f", "2p+" = "2a", "5p+" = "5b", "6p+" = "6a", "6q-" = "6f", "8p-" = "8a", "10p+" = "10a", "11q+" = "11c", "12p-" = "12b", "13cnloh" = "13e", "15cnloh" = "15e", "16q-" = "16d", "19p+" = "19a", "21+" = "21a", "22+" = "22b"),
    "SRR14800543" = c("1q+" = "1c", "12p+" = "12a", "13cnloh" = "13a", "16q-" = "16c", "17p-" = "17a", "18+" = "18b"), # 1q/16q only in GT 3
    "SRR17960481" = c("1q+" = "1b", "2p+" = "2a", "6p+" = "6a", "9q+" = "9b", "13cnloh" = "13b"), # 1q/2p/6p
    "SRR17960484" = c("1q+" = "1c", "2p+" = "2a", "5qcnloh" = "5e", "6p+" = "6a", "9qcnloh" = "9b", "10q+" = "10b", "11pcnloh" = "11a", "15+" = "15b", "16q-" = "16c") # 16q only in GT 4; 2,3,4 have 1q; interesting sample
  )),
  tar_target(
    clone_comparison_table,
    tabulate_clone_comparisons(large_clone_comparisons)
  ),
  tar_target(cluster_for_diffex, list(
    "SRR13884242" = 3, # 16q
    "SRR13884243" = 4, # 16q and 1q
    "SRR13884247" = 5, # 6p; 1q gain only affects GT 4
    "SRR13884249" = 2, # 1q; 2p only GT 3
    "SRR14800534" = 2, # 1q only GT 3; maybe reconsider
    "SRR14800535" = 1, # 1q; maybe reconsider
    "SRR14800536" = 2, # 16q
    "SRR14800540" = 4, # 16q; just 2
    "SRR14800541" = 1, # 16q is only in GT 2; 2,3,4,5 have 1q/2p/6p
    "SRR14800543" = 2, # 1q/16q only in GT 3
    "SRR17960481" = 2, # 1q/2p/6p
    "SRR17960484" = 1 # 16q only in GT 4; 2,3,4 have 1q; interesting sample
  )),
  tar_target(whole_pseudobulks,
    score_whole_pseudobulks(large_done_files, subtype_markers),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(
    derived_pseudobulk_subtype_scores,
    derive_pseudobulk_subtype_scores(large_done_files),
  ),
  tar_target(future_k_init, list(
    "SRR13884240" = 3,
    "SRR13884241" = 3,
    "SRR13884242" = 4,
    "SRR13884243" = 4,
    "SRR13884244" = 3,
    "SRR13884245" = 3,
    "SRR13884246" = 3,
    "SRR13884247" = 4,
    "SRR13884248" = 3,
    "SRR13884249" = 4,
    "SRR14800534" = 2,
    "SRR14800535" = 3,
    "SRR14800536" = 3,
    "SRR14800537" = 4,
    "SRR14800539" = 3,
    "SRR14800540" = 4,
    "SRR14800541" = 3,
    "SRR14800543" = 3,
    "SRR17960480" = 3,
    "SRR17960481" = 3,
    "SRR17960482" = 3,
    "SRR17960483" = 3,
    "SRR17960484" = 3
  )),
  tar_target(
    current_params,
    retrieve_snakemake_params(large_done_files),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(
    current_init_k,
    retrieve_current_param(current_params, "init_k")
  ),

  # wu metadata ------------------------------

  tar_target(wu_metadata, get_merged_metadata("output/scanpy/wu_merged/metadata.csv")),

  # yang metadata ------------------------------

  tar_target(yang_metadata, get_merged_metadata("output/scanpy/yang_merged/metadata.csv")),

  # field metadata ------------------------------

  tar_target(field_metadata, get_merged_metadata("output/scanpy/field_merged/metadata.csv")),

  # # collin metadata ------------------------------
  #
  # tar_target(collin_metadata, get_merged_metadata("output/scanpy/collin_merged/metadata.csv")),

  # total metadat ------------------------------
  tar_target(total_metadata, read_tsv("data/metadata.tsv")),
  tar_target(
    analyzed_samples, c(
      "SRR13884242",
      "SRR13884243",
      "SRR13884247",
      "SRR13884249",
      "SRR14800534",
      "SRR14800535",
      "SRR14800536",
      "SRR14800540",
      "SRR14800541",
      "SRR14800543",
      "SRR17960481",
      "SRR17960484"
    )
  ),

  tar_target(
    excluded_samples,
    c(
      "SRR14800541" # no diploid ancestral clone; only diploid are B2M and APOE
    )
  ),

  # large files ------------------------------

  tar_target(
    interesting_genes,
    unlist(list(
      "1q" = c("ANP32E", "ASH1L", "ASPM", "CENPF", "CKS1B", "CNIH4", "CRABP2", "ENAH", "ENO1", "IPO9", "KIF14", "NEK2", "NENF", "NUF2", "PSMD4", "RXRG", "UBE2T"),
      "2p" = c("CEP68", "MEIS1", "MDH1", "OST4", "PDIA6", "POLE4", "RAB1A", "SNRPG", "SOX11"),
      "6p" = c("ARID1B", "CASP8AP2", "CLIC1", "CUTA", "DDX39B", "DEK", "DST", "GLO1", "HDAC2", "HDDC2", "HMGA1", "LIN28B", "MCM3", "MNDA", "PRDM1", "SOX4"),
      "16q" = c("CHD9", "CNGB1", "CYBA", "MT1E", "MT1X", "MT2A"),
      "S1" = c("ARR3", "GUCA1C", "CRX", "PDC"),
      "S2" = c("EBF3", "DCX", "ROBO1", "SOX11", "GAP43", "PCDHB10", "STMN2", "NEFM", "POU4F2", "EBF1"),
      "interesting" = c("SETD6", "RBL2", "DDX19A", "EZH2", "RB1", "TFF1", "MDM4", "CD24"),
      "cell_cycle" = c("CCNA2", "CCNB1")
    ))
  ),

  # plot subtype scores
  tar_target(
    subtype_violins,
    score_and_vlnplot_seu(regressed_seus, subtype_markers, leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv"),
    pattern = map(regressed_seus),
    iteration = "list"
  ),

  # tar_target(
  #   mp_violins,
  #   score_and_vlnplot_seu(filtered_seus, mps[["Cancer"]], leiden_cluster_file = "results/adata_filtered_metadata_0.2.csv"),
  #   pattern = map(filtered_seus),
  #   iteration = "list"
  # ),

  tar_target(mp_heatmaps,
    score_and_heatmap_seu(filtered_seus, mps[["Cancer"]], leiden_cluster_file = "results/adata_filtered_metadata_0.2.csv"),
    pattern = map(filtered_seus),
    iteration = "list"
  ),

  # tar_target(
  #   stemness_violins,
  #   score_and_vlnplot_seu(filtered_seus, stemness_markers[["smith"]]),
  #   pattern = map(filtered_seus),
  #   iteration = "list"
  # ),

  tar_target(
    subtype_violin_files,
    compile_subtype_violins(analyzed_samples, subtype_violins)
  ),
  tar_target(
    marker_gene_featureplots_filtered,
    plot_putative_marker_across_samples(interesting_genes, unlist(filtered_seus), plot_type = FeaturePlot, group_by = "clone_opt", cluster_dictionary, extension = "filtered")
  ),

  # tar_target(
  #   marker_gene_featureplots_unfiltered,
  #   plot_putative_marker_across_samples(interesting_genes, unlist(unfiltered_seus), plot_type = FeaturePlot, group_by = "clone_opt", cluster_dictionary, extension = "unfiltered")
  # ),

  tar_target(
    marker_gene_vlnplots_by_clone,
    plot_putative_marker_across_samples(interesting_genes, unlist(filtered_seus), plot_type = VlnPlot, group_by = "clone_opt", cluster_dictionary, extension = "filtered")
  ),
  tar_target(
    marker_gene_vlnplots_by_cluster,
    plot_putative_marker_across_samples(interesting_genes, unlist(filtered_seus), plot_type = VlnPlot, group_by = "abbreviation", cluster_dictionary, extension = "filtered")
  ),
  tar_target(
    cluster_markers,
    collect_all_markers(large_done_files, "results/clusters.xlsx")
  ),
  tar_target(
    large_clusters,
    collect_clusters_from_seus(large_done_files)
  ),
  tar_target(
    large_filter_expressions,
    list(
      "SRR13884242" = c(
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "1c"',
        'clone_opt %in% c(2,3,4) & p_cnv < 0.7 & seg == "1c"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "16b"',
        'clone_opt %in% c(2,3,4) & p_cnv < 0.7 & seg == "16b"'
      ),
      "SRR13884243" = c(
        'clone_opt %in% c(2,3) & p_cnv < 0.7 & seg == "16b"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "16b"',
        'clone_opt %in% c(2,3) & p_cnv < 0.7 & seg == "1b"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "1b"'
      ),
      "SRR13884247" = c(
        'clone_opt %in% c(2,3,4,5,6) & p_cnv < 0.7 & seg == "6c"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "6c"',
        'clone_opt %in% c(1,2,3,4,5) & p_cnv > 0.3 & seg == "2b"',
        'clone_opt %in% c(6) & p_cnv < 0.7 & seg == "2b"'
      ),
      "SRR13884249" = c(
        'clone_opt %in% c(2,3) & p_cnv < 0.7 & seg == "1b"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "1b"',
        'clone_opt %in% c(1,2) & p_cnv > 0.3 & seg == "2a"',
        'clone_opt %in% c(3) & p_cnv < 0.7 & seg == "2a"'
      ),
      "SRR14800534" = c(
        'clone_opt %in% c(3) & p_cnv < 0.7 & seg == "1b"',
        'clone_opt %in% c(1,2) & p_cnv > 0.3 & seg == "1b"',
        'clone_opt %in% c(2,3) & p_cnv < 0.7 & seg == "16c"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "16c"'
      ),
      "SRR14800535" = c(
        'clone_opt %in% c(2,3) & p_cnv < 0.7 & seg == "16b"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "16b"',
        'clone_opt %in% c(3) & p_cnv < 0.7 & seg == "1b"',
        'clone_opt %in% c(1,2) & p_cnv > 0.3 & seg == "1b"'
      ),
      "SRR14800536" = c(
        'clone_opt %in% c(3,4) & p_cnv < 0.7 & seg == "1b"',
        'clone_opt %in% c(1,2) & p_cnv > 0.3 & seg == "1b"',
        'clone_opt %in% c(2,3,4) & p_cnv < 0.7 & seg == "16b"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "16b"'
      ),
      "SRR14800540" = c(
        'clone_opt %in% c(3) & p_cnv < 0.7 & seg == "1d"',
        'clone_opt %in% c(1,2) & p_cnv > 0.3 & seg == "1d"',
        'clone_opt %in% c(3) & p_cnv < 0.7 & seg == "1e"',
        'clone_opt %in% c(1,2) & p_cnv > 0.3 & seg == "1e"',
        'clone_opt %in% c(1,3) & p_cnv > 0.3 & seg == "6b"',
        'clone_opt %in% c(2) & p_cnv < 0.7 & seg == "6b"',
        'clone_opt %in% c(1,3) & p_cnv > 0.3 & seg == "16e"',
        'clone_opt %in% c(2) & p_cnv < 0.7 & seg == "16e"'
      ),
      "SRR14800541" = c(
        'clone_opt %in% c(2,3,4,5,6,7) & p_cnv < 0.7 & seg == "1f"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "1f"',
        'clone_opt %in% c(2,3,4,5,6,7) & p_cnv < 0.7 & seg == "2a"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "2a"',
        'clone_opt %in% c(2,3,4,5,6,7) & p_cnv < 0.7 & seg == "6a"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "6a"',
        'clone_opt %in% c(3,4) & p_cnv < 0.7 & seg == "16d"',
        'clone_opt %in% c(1,2,5,6,7) & p_cnv > 0.3 & seg == "16d"'
      ),
      "SRR14800543" = c(
        'clone_opt %in% c(3) & p_cnv < 0.7 & seg == "1c"',
        'clone_opt %in% c(1,2) & p_cnv > 0.3 & seg == "1c"',
        'clone_opt %in% c(1,2,4) & p_cnv > 0.3 & seg == "16c"',
        'clone_opt %in% c(3) & p_cnv < 0.7 & seg == "16c"'
      ),
      "SRR17960481" = c(
        'clone_opt %in% c(2,3) & p_cnv < 0.7 & seg == "1b"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "1b"'
      ),
      "SRR17960484" = c(
        'clone_opt %in% c(2,3,4,5) & p_cnv < 0.7 & seg == "1c"',
        'clone_opt %in% c(1) & p_cnv > 0.3 & seg == "1c"',
        'clone_opt %in% c(4,5) & p_cnv < 0.7 & seg == "2a"',
        'clone_opt %in% c(1,2,3) & p_cnv > 0.3 & seg == "2a"',
        'clone_opt %in% c(1,2) & p_cnv > 0.3 & seg == "6a"',
        'clone_opt %in% c(3,4,5) & p_cnv < 0.7 & seg == "6a"'
      )
    )
  ),
  tar_target(large_numbat_pdfs,
    convert_numbat_pngs(large_done_files),
    pattern = map(large_done_files),
    iteration = "list"
  ),

  # tar_target(
  #   large_plot_files,
  #   make_numbat_plot_files(large_done_files, cluster_dictionary, clone_simplifications = large_clone_simplifications),
  #   pattern = map(large_done_files),
  #   iteration = "list"
  # ),

  tar_target(filtered_large_plot_files,
    make_numbat_plot_files(large_done_files, cluster_dictionary, filter_expressions = large_filter_expressions, clone_simplifications = large_clone_simplifications, extension = "_filtered"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(filtered_seus,
    filter_cluster_save_seu(large_done_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions, extension = "_filtered", leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(regressed_seus,
    regress_filtered_seu(filtered_seus),
    pattern = map(filtered_seus),
    iteration = "list"
  ),
  tar_target(regression_effect_plots,
    plot_effects_of_cell_cycle_regression(large_done_files),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(unfiltered_seus,
    load_seu(large_done_files),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(large_heatmaps,
    make_numbat_heatmaps(large_done_files, p_min = 0.5, line_width = 0.1, large_filter_expressions, cluster_dictionary, extension = "_filtered"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(
    large_montage_pdfs,
    make_pdf_montages(filtered_large_plot_files, large_heatmaps)
  ),
  tar_target(
    expr_heatmap_pdfs,
    make_expression_heatmap_comparison(large_numbat_pdfs, large_heatmaps)
  ),

  # all ------------------------------
  tar_target(large_all_diffex_clones,
    find_diffex_clones(large_done_files, large_clone_comparisons, cluster_dictionary, location = "all"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(large_all_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(large_done_files, large_clone_comparisons, cluster_dictionary, location = "all"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(large_all_diffex_clones_for_each_phase,
    find_diffex_bw_clones_for_each_phase(large_done_files, large_clone_comparisons, cluster_dictionary, location = "all"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
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
      "SRR17960484" = 3
    )
  ),
  tar_target(
    volcano_large_all,
    make_volcano_diffex_clones(
      large_all_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_all.pdf",
      large_all_diffex_clones,
      "results/diffex_bw_clones_large_all.pdf",
      large_all_diffex_clones_for_each_phase,
      "results/diffex_bw_clones_per_phase_large_all.pdf"
    ),
  ),
  tar_target(
    table_large_all_diffex_clones,
    tabulate_diffex_clones(
      large_all_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_all.xlsx",
      "results/diffex_bw_clones_per_cluster_large_all_by_chr.xlsx",
      large_all_diffex_clones,
      "results/diffex_bw_clones_large_all.xlsx",
      "results/diffex_bw_clones_large_all_by_chr.xlsx",
      large_all_diffex_clones_for_each_phase,
      "results/diffex_bw_clones_per_phase_large_all.xlsx",
      "results/diffex_bw_clones_per_phase_large_all_by_chr.xlsx"
    ),
  ),

  # in segment ------------------------------
  tar_target(large_in_segment_diffex_clones,
    find_diffex_clones(large_done_files, large_clone_comparisons, cluster_dictionary, location = "in_segment"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(large_in_segment_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(large_done_files, large_clone_comparisons, cluster_dictionary, location = "in_segment"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(large_in_segment_diffex_clones_for_each_phase,
    find_diffex_bw_clones_for_each_phase(large_done_files, large_clone_comparisons, cluster_dictionary, location = "in_segment"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
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
  tar_target(
    volcano_large_in_segment,
    make_volcano_diffex_clones(
      large_in_segment_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_in_segment.pdf",
      large_in_segment_diffex_clones,
      "results/diffex_bw_clones_large_in_segment.pdf",
      large_in_segment_diffex_clones_for_each_phase,
      "results/diffex_bw_clones_per_phase_large_in_segment.pdf"
    ),
  ),
  tar_target(
    table_large_in_segment_diffex_clones,
    tabulate_diffex_clones(
      large_in_segment_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_in_segment.xlsx",
      "results/diffex_bw_clones_per_cluster_large_in_segment_by_chr.xlsx",
      large_in_segment_diffex_clones,
      "results/diffex_bw_clones_large_in_segment.xlsx",
      "results/diffex_bw_clones_large_in_segment_by_chr.xlsx",
      large_in_segment_diffex_clones_for_each_phase,
      "results/diffex_bw_clones_per_phase_large_in_segment.xlsx",
      "results/diffex_bw_clones_per_phase_large_in_segment_by_chr.xlsx"
    ),
  ),

  # out of segment ------------------------------
  tar_target(large_out_of_segment_diffex_clones,
    find_diffex_clones(large_done_files, large_clone_comparisons, cluster_dictionary, location = "out_of_segment"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(large_out_of_segment_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(large_done_files, large_clone_comparisons, cluster_dictionary, location = "out_of_segment"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(large_out_of_segment_diffex_clones_for_each_phase,
    find_diffex_bw_clones_for_each_phase(large_done_files, large_clone_comparisons, cluster_dictionary, location = "out_of_segment"),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(
    volcano_large_out_of_segment,
    make_volcano_diffex_clones(
      large_out_of_segment_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_out_of_segment.pdf",
      large_out_of_segment_diffex_clones,
      "results/diffex_bw_clones_large_out_of_segment.pdf",
      large_out_of_segment_diffex_clones_for_each_phase,
      "results/diffex_bw_clones_per_phase_large_out_of_segment.pdf"
    ),
  ),

  # kooi candidates
  tar_target(
    found_kooi_candidates,
    tally_kooi_candidates(cis_diffex_clones = "results/diffex_bw_clones_large_in_segment_by_chr.xlsx", trans_diffex_clones = "results/diffex_bw_clones_large_out_of_segment_by_chr.xlsx")
  ),
  tar_target(
    table_large_out_of_segment_diffex_clones,
    tabulate_diffex_clones(
      large_out_of_segment_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_out_of_segment.xlsx",
      "results/diffex_bw_clones_per_cluster_large_out_of_segment_by_chr.xlsx",
      large_out_of_segment_diffex_clones,
      "results/diffex_bw_clones_large_out_of_segment.xlsx",
      "results/diffex_bw_clones_large_out_of_segment_by_chr.xlsx",
      large_out_of_segment_diffex_clones_for_each_phase,
      "results/diffex_bw_clones_per_phase_large_out_of_segment.xlsx",
      "results/diffex_bw_clones_per_phase_large_out_of_segment_by_chr.xlsx"
    ),
  ),

  # # 16q
  # tar_target(integrated_seu_16q, seu_integrate_rbs(numbat_dir = "output/numbat_sridhar",
  #                                                  kept_samples = c('SRR14800534', 'SRR14800535', 'SRR14800536'),
  #                                                  cluster_dictionary)
  #            ),

  tar_target(
    integrated_seu_16q,
    readRDS("output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds")
  ),

  # tar_target(clone_comparisons_for_16q, list("2_v_1_16q-" = c("16b"), "3_v_2_1q+" = c("1b"))), #1q only GT 3; maybe reconsider
  #
  #
  #
  # # 16q/1q
  # tar_target(integrated_seu_1q_16q, seu_integrate_rbs(numbat_dir = "output/numbat_sridhar",
  #                                                     kept_samples = c('SRR13884242', 'SRR13884243', 'SRR13884249', 'SRR14800534', 'SRR14800535', 'SRR14800536', 'SRR14800543'),
  #                                                     cluster_dictionary)
  #            ),
  #
  # # 1q 6p
  # tar_target(integrated_seu_1q_6p, seu_integrate_rbs(numbat_dir = "output/numbat_sridhar",
  #                                                    kept_samples = c('SRR14800543', 'SRR17960481', 'SRR17960484'),
  #                                                    cluster_dictionary)
  #            ),


  # tar_target(
  #   large_in_segment_cluster_gse_plots,
  #   gse_plot_from_cluster_diffex(large_in_segment_diffex_clones_for_each_cluster),
  #   pattern = map(large_in_segment_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),
  #
  # tar_target(
  #   large_out_of_segment_cluster_gse_plots,
  #   gse_plot_from_cluster_diffex(large_out_of_segment_diffex_clones_for_each_cluster),
  #   pattern = map(large_out_of_segment_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),

  # tar_target(
  #   large_in_segment_cluster_gse_plots,
  #   gse_plot_from_clone_diffex(large_in_segment_diffex_clones_for_each_cluster),
  #   pattern = map(large_in_segment_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),
  #

  # tar_target(
  #   large_out_of_segment_cluster_gse_plots,
  #   gse_plot_from_clone_diffex(large_out_of_segment_diffex_clones_for_each_cluster),
  #   pattern = map(large_out_of_segment_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),

  tar_target(large_diffex_bw_clusters_for_each_clone,
    find_diffex_bw_clusters_for_each_clone(large_done_files, cluster_dictionary),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(
    large_numbat_expression,
    retrieve_numbat_plot_type(large_numbat_pdfs, "exp_roll_clust.pdf")
  ),
  tar_target(
    clone_trees,
    retrieve_numbat_plot_type(filtered_large_plot_files, "tree_filtered.pdf")
  ),
  tar_target(
    clone_distribution_plots,
    retrieve_numbat_plot_type(filtered_large_plot_files, "clone_distribution_filtered.pdf")
  ),

  # tar_target(large_clone_dist, retrieve_numbat_plot_type(large_plot_files, "exp_roll_clust.pdf")),

  # tar_target(large_numbat_expression_pdf, qpdf::pdf_combine(large_numbat_expression, "results/numbat_sridhar_large/large_numbat_expression.pdf")),
  #
  # tar_target(large_numbat_heatmap_files, retrieve_numbat_plot_type(large_numbat_pdfs, "panel_2.pdf")),
  #
  # tar_target(large_numbat_heatmap_pdf, qpdf::pdf_combine(large_numbat_heatmap_files, "results/numbat_sridhar_large/large_numbat_heatmaps.pdf")),
  #
  # tar_target(large_numbat_bulk_clones_final, retrieve_numbat_plot_type(large_numbat_pdfs, "bulk_clones_final.pdf")),
  #
  # tar_target(large_numbat_bulk_clones_final_pdf, qpdf::pdf_combine(large_numbat_bulk_clones_final, "results/numbat_sridhar_large/large_bulk_clones_final.pdf")),
  #
  # # large_numbat_sample_pdfs ------------------------------
  # tar_target(large_numbat_sample_pdfs,
  #   reroute_done_to_results_pdf(large_done_files, "_large"),
  #   pattern = map(large_done_files),
  #   iteration = "list",
  #   format = "file"
  # ),

  tar_target(large_numbat_sample_pdfs,
             reroute_done_to_results_pdf(large_done_files, "_large"),
             pattern = map(large_done_files),
             iteration = "list"
  ),

  tar_target(
    oncoprint_input_by_scna_unfiltered,
    make_oncoprint_diffex_unfiltered(large_filter_expressions, cluster_dictionary, analyzed_samples, large_in_segment_diffex_clones, large_out_of_segment_diffex_clones, large_clone_comparisons, n_slice = 20)
  ),

  # tar_target(
  #   oncoprint_input_by_scna_for_each_cluster_unfiltered,
  #   make_oncoprint_diffex_unfiltered(large_filter_expressions, cluster_dictionary, analyzed_samples, large_in_segment_diffex_clones_for_each_cluster, large_out_of_segment_diffex_clones_for_each_cluster, large_clone_comparisons, n_slice = 20)
  # ),

  tar_target(
    num_diffex_clone_scna_tally,
    tally_num_diffex(oncoprint_input_by_scna_unfiltered)
  ),
  #
  # tar_target(
  #   num_diffex_cluster_scna_tally,
  #   tally_num_diffex(oncoprint_input_by_scna_for_each_cluster_unfiltered)
  # ),

  tar_target(
    oncoprint_input_by_scna,
    make_oncoprint_diffex(large_filter_expressions, cluster_dictionary, analyzed_samples, large_in_segment_diffex_clones, large_out_of_segment_diffex_clones, large_clone_comparisons, n_slice = 20)
  ),

  # tar_target(
  #   oncoprint_input_by_scna_for_each_cluster,
  #   make_oncoprint_diffex(large_filter_expressions, cluster_dictionary, analyzed_samples, large_in_segment_diffex_clones_for_each_cluster, large_out_of_segment_diffex_clones_for_each_cluster, large_clone_comparisons, n_slice = 20)
  # ),

  # large_in_segment_diffex_clones_for_each_cluster

  tar_target(
    oncoprint_enrich_clones_gobp,
    enrich_oncoprints(large_filter_expressions,
      cluster_dictionary,
      analyzed_samples,
      large_in_segment_diffex_clones,
      large_out_of_segment_diffex_clones,
      large_clone_comparisons,
      gene_set = "gobp"
    )
  ),
  tar_target(
    oncoprint_enrich_clones_hallmark,
    enrich_oncoprints(large_filter_expressions,
      cluster_dictionary,
      analyzed_samples,
      large_in_segment_diffex_clones,
      large_out_of_segment_diffex_clones,
      large_clone_comparisons,
      gene_set = "hallmark"
    )
  ),
  tar_target(rod_rich_samples,
    score_samples_for_rod_enrichment(large_done_files),
    pattern = map(large_done_files),
    iteration = "list"
  ),
  tar_target(
    oncoprint_enrich_clones_plots_gobp,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_gobp,
      cis_plot_file = "results/enrichment/cis_enrichment_plots_by_clone_gobp",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_clone_gobp",
      cis_table_file = "results/enrichment/cis_enrichment_tables_by_clone_gobp.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_clone_gobp.xlsx"
    )
  ),
  tar_target(
    oncoprint_enrich_clones_plots_hallmark,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_hallmark,
      cis_plot_file = "results/enrichment/cis_enrichment_plots_by_clone_hallmark",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_clone_hallmark",
      cis_table_file = "results/enrichment/cis_enrichment_tables_by_clone_hallmark.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_clone_hallmark.xlsx"
    )
  ),
  tar_target(
    subtype_markers, pull_subtype_genes(supp_excel = "data/Liu_Radvanyi_2022_supp_data/41467_2021_25792_MOESM6_ESM.xlsx")
  ),
  tar_target(
    mps, read_mps("/dataVolume/storage/Homo_sapiens/3ca/ITH_hallmarks/MPs_distribution/MP_list.RDS")
  ),
  tar_target(
    stemness_markers,
    pull_stem_cell_markers()
  ),

  # phases enrichment ------------------------------
  # enrichment by phase msigdb gobp
  tar_target(
    oncoprint_enrich_phases_gobp,
    enrich_oncoprints_phases(large_filter_expressions,
      cluster_dictionary,
      analyzed_samples,
      large_in_segment_diffex_clones_for_each_phase,
      large_out_of_segment_diffex_clones_for_each_phase,
      large_clone_comparisons,
      gene_set = "gobp"
    )
  ),


  # enrichment by phase msigdb hallmark
  tar_target(
    oncoprint_enrich_phases_hallmark,
    enrich_oncoprints_phases(large_filter_expressions,
      cluster_dictionary,
      analyzed_samples,
      large_in_segment_diffex_clones_for_each_phase,
      large_out_of_segment_diffex_clones_for_each_phase,
      large_clone_comparisons,
      gene_set = "hallmark"
    )
  ),

  # enrichment plots and tables by phase gobp
  tar_target(
    oncoprint_enrich_phases_plots_gobp,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_phases_gobp,
      cis_plot_file = "results/enrichment/cis_enrichment_plots_by_phase_gobp",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_phase_gobp",
      cis_table_file = "results/enrichment/cis_enrichment_tables_by_phase_gobp.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_phase_gobp.xlsx",
      by_cluster = TRUE
    )
  ),

  # enrichment plots and tables by phase msigdb hallmark
  tar_target(
    oncoprint_enrich_phases_plots_hallmark,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_phases_hallmark,
      cis_plot_file = "results/enrichment/cis_enrichment_plots_by_phase_hallmark",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_phase_hallmark",
      cis_table_file = "results/enrichment/cis_enrichment_tables_by_phase_hallmark.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_phase_hallmark.xlsx",
      by_cluster = TRUE
    )
  ),

  # clusters enrichment ------------------------------
  # enrichment by cluster msigdb gobp
  tar_target(
    oncoprint_enrich_clusters_gobp,
    enrich_oncoprints_clusters(large_filter_expressions,
      cluster_dictionary,
      analyzed_samples,
      large_in_segment_diffex_clones_for_each_cluster,
      large_out_of_segment_diffex_clones_for_each_cluster,
      large_clone_comparisons,
      gene_set = "gobp"
    )
  ),


  # enrichment by cluster msigdb hallmark
  tar_target(
    oncoprint_enrich_clusters_hallmark,
    enrich_oncoprints_clusters(large_filter_expressions,
      cluster_dictionary,
      analyzed_samples,
      large_in_segment_diffex_clones_for_each_cluster,
      large_out_of_segment_diffex_clones_for_each_cluster,
      large_clone_comparisons,
      gene_set = "hallmark"
    )
  ),

  # enrichment plots and tables by cluster gobp
  tar_target(
    oncoprint_enrich_clusters_plots_gobp,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clusters_gobp,
      cis_plot_file = "results/enrichment/cis_enrichment_plots_by_cluster_gobp",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_cluster_gobp",
      cis_table_file = "results/enrichment/cis_enrichment_tables_by_cluster_gobp.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_cluster_gobp.xlsx",
      by_cluster = TRUE
    )
  ),

  # enrichment plots and tables by cluster msigdb hallmark
  tar_target(
    oncoprint_enrich_clusters_plots_hallmark,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clusters_hallmark,
      cis_plot_file = "results/enrichment/cis_enrichment_plots_by_cluster_hallmark",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_cluster_hallmark",
      cis_table_file = "results/enrichment/cis_enrichment_tables_by_cluster_hallmark.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_cluster_hallmark.xlsx",
      by_cluster = TRUE
    )
  ),
  tar_target(
    oncoprint_input_by_region,
    inspect_oncoprints(oncoprint_input_by_scna$cis, oncoprint_input_by_scna$trans)
  ),
  tar_target(
    oncoprint_plots,
    make_oncoprint_plots(oncoprint_input_by_scna$cis, oncoprint_input_by_scna$trans)
  ),

  #
  # tar_target(
  #   large_cluster_markers,
  #   collect_all_markers(large_done_files, "results/clusters.xlsx")
  # )
  #
  # # ),

  # end of plan------------------------------
)

# notes 2023-03-07 ------------------------------

# yes ------------------------------
# wu SRR13884242 (RB2 rep 1): likely 16q GT; f31 sample; cluster 1 (prolif.) has higher proportion of GT 2; can attribute this proliferation to acquisition of 16q- in GT 2
# wu SRR13884243 (RB2 rep 2): likely 16q GT; biological replication of SRR13884242; cluster 4 is an analogue of SRR13884242 c3 # nolint

# wu SRR13884247: likely 16q GT; no GSEA diffex output; cluster 5 is interesting; c5 marker DLK1 in amicroRNA cluster with MEG3
# wu SRR13884249: possible 1q/2p/16q GT; evidence that GT 1 (only 16q-) does not contribute to proliferating clusters c2 and c4 (TOP2A high); split by phase
# yang SRR14800534: GT1 lacks 16q-; does not contribute to proliferating clusters c2 and c4; GT 2 lacks 1q--can't identify tx distinction b/w gt2/3
# yang SRR14800535: 16q GT; GT 1 decreased in c1; c1 high TOP2A
# yang SRR14800536: possible 16q GT; GT1 decreased in c2/3 high G2M markers
# yang SRR14800540: clear 16q GT; very complicated tumor; three GTs with SCNAs; each is proliferating to a greater degree; confused about GT5 with high c4; markers are C1QA/B, CD74, HLA genes
# yang SRR14800541: clear 1q 6p and 16q GTs; GT1 not proliferating--no contribution to c2 with G2M markers
# yang SRR14800543: possible minor 16q/1q GT; can identify dual or individual contribution of 1q and 16q with 13q CNLOH; strange MYCN marker of clusters
# field SRR17960481: 6p interesting; cluster 2 and 4; cluster 2 has mito genes; not promising; likely clonal 6p with PRs and stressed cells composing GT 1
# field SRR17960484: cluster 1 enriched for wt GT 1; also has high Xist expression

# maybe ------------------------------
# SRR13884240: possible 2p GT
# SRR13884241: possible 1q GT
# SRR13884244: possible 1q GT
# SRR13884245: possible 1q GT
# SRR13884246: possible 16q GT
# SRR17960480: possible minor 16q- GT
# SRR14800539: possible 16q GT

# no ------------------------------
# SRR14800537: possible 16q GT
# SRR17960482: too complicated
# SRR13884248: clear 6p (missing possible 2p in expression)
# SRR17960483: cluster 6 maybe interesting
