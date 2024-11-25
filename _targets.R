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
# debug(calculate_clone_distribution)

tar_option_set(
  memory = "transient",
  garbage_collection = TRUE,
  error = "null"
)

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
    figures_and_tables,
    list(
      clustree_compilation,
      collage_compilation,
      table_all_diffex_clones,
      # oncoprint_enrich_clones_plots_gobp,
      # oncoprint_enrich_clusters_plots_gobp,
      oncoprint_plots, 
      cluster_comparisons_by_phase_for_disctinct_clones
      
    )
  ),
  tar_target(
    plae_ref,
    generate_plae_ref()
  ),
  tarchetypes::tar_files(numbat_rds_files, retrieve_numbat_rds_files("output/numbat_sridhar/", interesting_samples), format = "file"),
  tarchetypes::tar_files(seus, retrieve_seus("output/seurat/", interesting_samples), format = "file"),
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

  tar_target(cells_to_remove, read_cells_to_remove("data/cells_to_remove_final.xlsx")),
  tar_target(cluster_dictionary, read_cluster_dictionary("data/cluster_dictionary.tsv")),
  tar_target(hallmark_gene_sets, msigdbr(species = "Homo sapiens", category = "H")),
  tar_target(interesting_samples,
    c(
      "SRR13884242",
      "SRR13884243",
      # "SRR13884244",
      # "SRR13884245",
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
      # "SRR27187901",
      "SRR27187902"
    )
  ),
  tar_target(debranched_ids,
  					 c("SRR13884242", "SRR13884243", "SRR13884246_branch_5", "SRR13884246_branch_6", 
  					 	"SRR13884247_branch_6", "SRR13884247_branch_4", "SRR13884247_branch_5", 
  					 	"SRR13884248", "SRR13884249", "SRR14800534", "SRR14800535", "SRR14800536", 
  					 	"SRR14800540_branch_2", "SRR14800540_branch_3", "SRR14800541_branch_4", 
  					 	"SRR14800541_branch_7", "SRR14800543_branch_3", "SRR14800543_branch_4", 
  					 	"SRR17960481", "SRR17960484", "SRR27187899", "SRR27187902_branch_3", 
  					 	"SRR27187902_branch_4")
  ),
  
  tar_target(large_clone_comparisons,
    list(
      "SRR13884242" = list("2_v_1_1q+_16q-" = c("1b", "16b"), "3_v_2_8p-_11p-" = c("8a", "11a")), # 16q
      "SRR13884243" = list("2_v_1_1q+_16q-" = c("1b", "16b"), "3_v_2_8p-_11p-" = c("8a", "11a")), # 16q and 1q
      "SRR13884244" = list("2_v_1_1q+_16q-" = c("1b", "16a")), # 16q and 1q
      "SRR13884245" = list("2_v_1_1q+_16q-" = c("1b", "16b")), # 16q and 1q
      "SRR13884246" = list("2_v_1_1q+" = c("1b"), "5_v_4_2p+" = c("2b")), # 16q and 1q
      "SRR13884246_branch_5" = list("2_v_1_1q+" = c("1b"), "5_v_4_2p+" = c("2b")), # 16q and 1q
      "SRR13884246_branch_6" = list("2_v_1_1q+" = c("1b")), # 16q and 1q
      "SRR13884247" = list("2_v_1_6p+" = c("6c"), "3_v_2_17q+" = c("17d"), "4_v_3_10q+" = c("10b"), "5_v_3_1p-" = c("1a"), "6_v_2_2p+" = c("2b")), # 6p; 1q gain only affects GT 4
      "SRR13884247_branch_6" = list("2_v_1_6p+" = c("6c"), "6_v_2_2p+" = c("2b")), # 6p; 1q gain only affects GT 4
      "SRR13884247_branch_4" = list("2_v_1_6p+" = c("6c"), "3_v_2_17q+" = c("17d"), "4_v_3_10q+" = c("10b")), # 6p; 1q gain only affects GT 4
      "SRR13884247_branch_5" = list("2_v_1_6p+" = c("6c"), "3_v_2_17q+" = c("17d"), "5_v_3_1p-" = c("1a"), "6_v_2_2p+" = c("2b")), # 6p; 1q gain only affects GT 4
      "SRR13884248" = list("2_v_1_6p+" = c("6a"), "3_v_2_2p+" = c("2b")), # 16q and 1q
      "SRR13884249" = list("2_v_1_1q+" = c("1b"), "3_v_2_2p+" = c("2a")), # 1q; 2p only GT 3
      "SRR14800534" = list("2_v_1_16q-" = c("16c"), "3_v_2_1q+" = c("1b")), # 1q only GT 3; maybe reconsider
      "SRR14800535" = list("2_v_1_16q-" = c("16b"), "3_v_2_1q+" = c("1b")), # 1q only GT 3; maybe reconsider
      "SRR14800536" = list("2_v_1_16q-" = c("16b"), "3_v_2_1q+" = c("1b"), "4_v_3_15q+" = c("15b")), # 16q
      "SRR14800540_branch_2" = list("2_v_1_6p+_16q-" = c("6b", "16e")), # 16q; just 2
      "SRR14800540_branch_3" = list("3_v_2_1q+_7q+" = c("1d", "1e", "7b")), # 16q; just 2
      "SRR14800541_branch_4" = list("2_v_1_16q-" = c("1f", "6a"), "3_v_2_16q-" = c("16d")), # 16q is only in GT 2; 2,3,4,5 have 1q/2p/6p
      "SRR14800541_branch_7" = list("2_v_1_16q-" = c("1f", "6a")), # 16q is only in GT 2; 2,3,4,5 have 1q/2p/6p
      "SRR14800543" = list("2_v_1_13loh" = c("13a"), "3_v_2_1q+_16q-" = c("1c", "16c")), # 1q/16q only in GT 3
      "SRR14800543_branch_3" = list("2_v_1_13loh" = c("13a"), "3_v_2_1q+_16q-" = c("1c", "16c")), # 1q/16q only in GT 3
      "SRR14800543_branch_4" = list("2_v_1_13loh" = c("13a")), # 1q/16q only in GT 3
      "SRR17960481" = list("2_v_1_1q+_6p+" = c("1b", "6a", "6b"), "3_v_2_2p+" = c("2a")), # 1q/2p/6p
      "SRR17960484" = list("2_v_1_1q+" = c("1c"), "3_v_2_6p+" = c("6a"), "4_v_3_2p+" = c("2a")),
      "SRR27187899" = list("2_v_1_6p+" = c("6b"), "3_v_2_16q-" = c("16e")),
      # "SRR27187900" = list("3_v_2_1q+_6p+_16q-" = c("1d", "6b", "16b")), # not enough cells "2_v_1_1q+_6p+_16q-" = c("1d", "6b", "16b")
      # "SRR27187901" = list("2_v_1_1q+" = c("1c"), "3_v_2_6p+" = c("6a"), "4_v_3_2p+" = c("2a")),
      "SRR27187902_branch_3" = list("2_v_1_1q+_2p+" = c("1k", "2b")),
      "SRR27187902_branch_4" = list("2_v_1_1q+_2p+" = c("1k", "2b"), "4_v_2_16q-" = c("16a"))
    )
  ),
  tar_target(
    large_clone_simplifications,
    list(
      "SRR13884242" = c("1q+" = "1b", "2p+" = "2b", "8p-" = "8a", "11p-" = "11a", "12p-" = "12a", "16q-" = "16b"), # 16q
      "SRR13884243" = c("1q+" = "1b", "2p+" = "2b", "8p-" = "8a", "11p-" = "11a", "12p-" = "12a", "16q-" = "16b"),
      "SRR13884244" = c("1q+" = "1b", "5q+" = "5a", "5q-" = "5f", "12p-" = "12b", "12q+" = "12f", "16q-" = "16a", "19p+" = "19a"),
      "SRR13884245" = c("1q+" = "1b", "5p+" = "5a", "5q-" = "5c", "12p-" = "12b", "12q-" = "12g", "16q-" = "16b", "19p+" = "19a"),
      "SRR13884246" = c("1q+" = "1b", "2p+" = "2b", "8+" = "8a", "12p-" = "12a", "14q-" = "14d", "16p-" = "16b", "17p-" = "17a"),
      "SRR13884247" = c("1p-" = "1a", "2p+" = "2b", "6p+" = "6c", "10q+" = "10b", "11q-" = "11b", "13q-" = "13b", "17q+" = "17d", "20p-" = "20a"),
      "SRR13884248" = c("2p+" = "2b", "6p+" = "6a"),
      "SRR13884249" = c("1q+" = "1b", "2p+" = "2a", "13cnloh" = "13b", "16q-" = "16b"),
      "SRR14800534" = c("1q+" = "1b", "16q-" = "16c"), # 1q only GT 3; maybe reconsider
      "SRR14800535" = c("1q+" = "1b", "13cnloh" = "13a", "16q-" = "16b"), # 1q only GT 3; maybe reconsider
      "SRR14800536" = c("1q+" = "1b", "13-" = "13b", "15q+" = "15b", "16q-" = "16b", "19q-" = "19c"), # 16q
      "SRR14800540" = c("1q+" = "1d", "3q-" = "3b", "6q+" = "6b", "7q+" = "7b", "11q-" = "11b", "12p-" = "12f", "13-" = "13b", "16q-" = "16e", "20q+" = "20b"),
      "SRR14800541" = c("1q+" = "1f", "2p+" = "2a", "5p+" = "5b", "6p+" = "6a", "6q-" = "6f", "8p-" = "8a", "10p+" = "10a", "11q+" = "11c", "12p-" = "12b", "13cnloh" = "13e", "15cnloh" = "15e", "16q-" = "16d", "19p+" = "19a", "21+" = "21a", "22+" = "22b"),
      "SRR14800543" = c("1q+" = "1c", "12p+" = "12a", "13cnloh" = "13a", "16q-" = "16c", "17p-" = "17a", "18+" = "18b"), # 1q/16q only in GT 3
      "SRR17960481" = c("1q+" = "1b", "2p+" = "2a", "6p+" = "6a", "9q+" = "9b", "13cnloh" = "13b"), # 1q/2p/6p
      "SRR17960484" = c("1q+" = "1c", "2p+" = "2a", "5qcnloh" = "5e", "6p+" = "6a", "9qcnloh" = "9b", "10q+" = "10b", "11pcnloh" = "11a", "15+" = "15b", "16q-" = "16c"), # 16q only in GT 4; 2,3,4 have 1q; interesting sample
      "SRR27187899" = c("6p+" = "6b", "16q-" = "16e"), # 1q/2p/6p
      # "SRR27187900" = c("1q+" = "1d", "6p+" = "6c", "16q-" = "16b"), # 1q/2p/6p
      # "SRR27187901" = c("1q+" = "1b", "2p+" = "2a", "6p+" = "6a", "9q+" = "9b", "13cnloh" = "13b"), # 1q/2p/6p
      "SRR27187902" = c("1q+" = "1k", "2p+" = "2b", "16q-" = "16a") # 1q/2p/6p
    )
  ),
  tar_target(
    large_filter_expressions,
    list(
      "SRR13884242" = c(
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1b"',
        'clone_opt %in% c(2,3,4) & p_cnv <= 0.5 & seg == "1b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "16b"',
        'clone_opt %in% c(2,3,4) & p_cnv <= 0.5 & seg == "16b"'
      ),
      "SRR13884243" = c(
        'clone_opt %in% c(2,3) & p_cnv <= 0.5 & seg == "16b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "16b"',
        'clone_opt %in% c(2,3) & p_cnv <= 0.5 & seg == "1b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1b"'
      ),
      "SRR13884244" = c(
        'clone_opt %in% c(2,3,4,5) & p_cnv <= 0.5 & seg == "1b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1b"',
        'clone_opt %in% c(2,3,4,5) & p_cnv <= 0.5 & seg == "16a"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "16a"'
      ),
      "SRR13884245" = c(
        'clone_opt %in% c(2,3) & p_cnv <= 0.5 & seg == "1b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1b"',
        'clone_opt %in% c(2,3) & p_cnv <= 0.5 & seg == "16b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "16b"'
      ),
      "SRR13884246" = c(
        'clone_opt %in% c(2,3,4,5,6) & p_cnv <= 0.5 & seg == "1b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1b"',
        'clone_opt %in% c(5,6) & p_cnv <= 0.5 & seg == "2b"',
        'clone_opt %in% c(1,2,3,4) & p_cnv > 0.5 & seg == "2b"'
      ),
      "SRR13884247" = c(
        'clone_opt %in% c(2,3,4,5,6) & p_cnv <= 0.5 & seg == "6c"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "6c"',
        'clone_opt %in% c(1,2,3,4,5) & p_cnv > 0.5 & seg == "2b"',
        'clone_opt %in% c(6) & p_cnv <= 0.5 & seg == "2b"'
      ),
      "SRR13884248" = c(
        'clone_opt %in% c(3) & p_cnv <= 0.5 & seg == "2b"',
        'clone_opt %in% c(1,2) & p_cnv > 0.5 & seg == "2b"',
        'clone_opt %in% c(2) & p_cnv <= 0.5 & seg == "6a"'
      ),
      "SRR13884249" = c(
        'clone_opt %in% c(2,3) & p_cnv <= 0.5 & seg == "1b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1b"',
        'clone_opt %in% c(1,2) & p_cnv > 0.5 & seg == "2a"',
        'clone_opt %in% c(3) & p_cnv <= 0.5 & seg == "2a"'
      ),
      "SRR14800534" = c(
        'clone_opt %in% c(3) & p_cnv <= 0.5 & seg == "1b"',
        'clone_opt %in% c(1,2) & p_cnv > 0.5 & seg == "1b"',
        'clone_opt %in% c(2,3) & p_cnv <= 0.5 & seg == "16c"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "16c"'
      ),
      "SRR14800535" = c(
        'clone_opt %in% c(2,3) & p_cnv <= 0.5 & seg == "16b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "16b"',
        'clone_opt %in% c(3) & p_cnv <= 0.5 & seg == "1b"',
        'clone_opt %in% c(1,2) & p_cnv > 0.5 & seg == "1b"'
      ),
      "SRR14800536" = c(
        'clone_opt %in% c(3,4) & p_cnv <= 0.5 & seg == "1b"',
        'clone_opt %in% c(1,2) & p_cnv > 0.5 & seg == "1b"',
        'clone_opt %in% c(2,3,4) & p_cnv <= 0.5 & seg == "16b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "16b"'
      ),
      "SRR14800540" = c(
        'clone_opt %in% c(3) & p_cnv <= 0.5 & seg == "1d"',
        'clone_opt %in% c(1,2) & p_cnv > 0.5 & seg == "1d"',
        'clone_opt %in% c(3) & p_cnv <= 0.5 & seg == "1e"',
        'clone_opt %in% c(1,2) & p_cnv > 0.5 & seg == "1e"',
        'clone_opt %in% c(1,3) & p_cnv > 0.5 & seg == "6b"',
        'clone_opt %in% c(2) & p_cnv <= 0.5 & seg == "6b"',
        'clone_opt %in% c(1,3) & p_cnv > 0.5 & seg == "16e"',
        'clone_opt %in% c(2) & p_cnv <= 0.5 & seg == "16e"'
      ),
      "SRR14800541" = c(
        'clone_opt %in% c(2,3,4,5,6,7) & p_cnv <= 0.5 & seg == "1f"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1f"',
        'clone_opt %in% c(2,3,4,5,6,7) & p_cnv <= 0.5 & seg == "2a"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "2a"',
        'clone_opt %in% c(2,3,4,5,6,7) & p_cnv <= 0.5 & seg == "6a"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "6a"',
        'clone_opt %in% c(3,4) & p_cnv <= 0.5 & seg == "16d"',
        'clone_opt %in% c(1,2,5,6,7) & p_cnv > 0.5 & seg == "16d"'
      ),
      "SRR14800543" = c(
        'clone_opt %in% c(3) & p_cnv <= 0.5 & seg == "1c"',
        'clone_opt %in% c(1,2) & p_cnv > 0.5 & seg == "1c"',
        'clone_opt %in% c(1,2,4) & p_cnv > 0.5 & seg == "16c"',
        'clone_opt %in% c(3) & p_cnv <= 0.5 & seg == "16c"'
      ),
      "SRR17960481" = c(
        'clone_opt %in% c(2,3) & p_cnv <= 0.5 & seg == "1b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1b"'
      ),
      "SRR17960484" = c(
        'clone_opt %in% c(2,3,4,5) & p_cnv <= 0.5 & seg == "1c"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1c"',
        'clone_opt %in% c(4,5) & p_cnv <= 0.5 & seg == "2a"',
        'clone_opt %in% c(1,2,3) & p_cnv > 0.5 & seg == "2a"',
        'clone_opt %in% c(1,2) & p_cnv > 0.5 & seg == "6a"',
        'clone_opt %in% c(3,4,5) & p_cnv <= 0.5 & seg == "6a"'
      ),
      "SRR27187899" = c(
        'clone_opt %in% c(2,3,4) & p_cnv <= 0.5 & seg == "6b"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "6b"',
        'clone_opt %in% c(1,2,3) & p_cnv > 0.5 & seg == "16e"',
        'clone_opt %in% c(4) & p_cnv <= 0.5 & seg == "16e"'
      ),
      # "SRR27187900" = c(
      # 	'clone_opt %in% c(2,3) & p_cnv <= 0.5 & seg == "1d"',
      # 	'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1d"',
      # 	'clone_opt %in% c(2,3) & p_cnv <= 0.5 & seg == "6c"',
      # 	'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "6c"'
      # ),
      # "SRR27187901" = c(
      # 	'clone_opt %in% c(2,3,4,5) & p_cnv <= 0.5 & seg == "1c"',
      # 	'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1c"',
      # 	'clone_opt %in% c(4,5) & p_cnv <= 0.5 & seg == "2a"',
      # 	'clone_opt %in% c(1,2,3) & p_cnv > 0.5 & seg == "2a"',
      # 	'clone_opt %in% c(1,2) & p_cnv > 0.5 & seg == "6a"',
      # 	'clone_opt %in% c(3,4,5) & p_cnv <= 0.5 & seg == "6a"'
      # ),
      "SRR27187902" = c(
        'clone_opt %in% c(2,3,4) & p_cnv <= 0.5 & seg == "1k"',
        'clone_opt %in% c(1) & p_cnv > 0.5 & seg == "1k"',
        'clone_opt %in% c(1,2,3) & p_cnv > 0.5 & seg == "16a"',
        'clone_opt %in% c(4) & p_cnv <= 0.5 & seg == "16a"'
      )
    )
  ),
  tar_target(
    clone_comparison_table,
    tabulate_clone_comparisons(large_clone_comparisons)
  ),
  tar_target(whole_pseudobulks,
    score_whole_pseudobulks(numbat_rds_files, subtype_markers),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),
  tar_target(
    unfiltered_derived_pseudobulk_subtype_scores,
    derive_pseudobulk_subtype_scores(unfiltered_seus),
  ),
  tar_target(
    filtered_derived_pseudobulk_subtype_scores,
    derive_pseudobulk_subtype_scores(final_seus),
  ),

  # tar_target(
  #   current_params,
  #   retrieve_snakemake_params(numbat_rds_files),
  #   pattern = map(numbat_rds_files),
  #   iteration = "list"
  # ),
  #
  # tar_target(
  #   current_init_k,
  #   retrieve_current_param(current_params, "init_k")
  # ),

  # total metadata ------------------------------
  tar_target(study_cell_stats, collect_study_metadata()),
  tar_target(fig_s15, plot_study_metadata(study_cell_stats)),
  tar_target(table_s01, make_table_s01(study_cell_stats)), # umi, genes detected, mito %
  tar_target(table_s02, make_table_s02()), # detailed sample metadata
  tar_target(table_s03, make_table_s03(cluster_dictionary)), # tally clusters removed (filtered out) with marker genes
  tar_target(table_s04_fig_s16, make_fig_s16_table_s04()), # RB SCNA frequency in TCGA by cancer type (Taylor et al. 2018).
  tar_target(table_s07, make_table_s07()), # For each 1q+ sample, percentage of clone in cluster.
  tar_target(table_s09, make_table_s09(seu_path = "output/seurat/integrated_16q/integrated_seu_16q_complete.rds", table_path = "results/table_s09.csv")), # For each 16q- sample, percentage of clone in cluster.
  tar_target(table_s10, make_table_s10(seu_path = "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", table_path = "results/table_s10.csv")), # For each 16q- sample, percentage of clone in cluster.
  tar_target(total_metadata, read_tsv("data/metadata.tsv")),
  tar_target(
    excluded_samples,
    c(
      "SRR14800541" # no diploid ancestral clone; only diploid are B2M and APOE
    )
  ),
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
  tar_target(subtype_violins,
    score_and_vlnplot_seu(final_seus, numbat_rds_files, large_clone_simplifications, subtype_markers),
    pattern = map(final_seus, numbat_rds_files, large_clone_simplifications),
    iteration = "list"
  ),

  tar_target(
  	subtype_violin_files,
  	qpdf::pdf_combine(unlist(subtype_violins), "results/scna_vlns.pdf")
  ),
  #
  # tar_target(merged_subtype_violins,
  #            score_and_vlnplot_seu("output/seurat/SRR14800534_SRR14800535_SRR14800536_seu.rds", subtype_markers, leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv", y_lim = 1.8, step = 0.2),
  # ),

  tar_target(
    patchwork_phase_plot,
    plot_phase_distribution_of_all_samples_by_scna(final_seus)
  ),

  # tar_target(
  #   mp_violins,
  #   score_and_vlnplot_seu(final_seus, mps[["Cancer"]], leiden_cluster_file = "results/adata_filtered_metadata_0.2.csv"),
  #   pattern = map(final_seus),
  #   iteration = "list"
  # ),

  tar_target(mp_heatmaps,
    score_and_heatmap_seu(final_seus, mps[["Cancer"]], leiden_cluster_file = "results/adata_filtered_metadata_0.2.csv"),
    pattern = map(final_seus),
    iteration = "list"
  ),

  # tar_target(
  #   stemness_violins,
  #   score_and_vlnplot_seu(final_seus, stemness_markers[["smith"]]),
  #   pattern = map(final_seus),
  #   iteration = "list"
  # ),

  tar_target(
    marker_gene_vlnplots_by_clone,
    plot_putative_marker_across_samples(unlist(interesting_genes), scna_seus, plot_type = VlnPlot, group_by = "scna")
  ),
  
  tar_target(
  	cluster_vlns_all,
    plot_putative_marker_across_samples(unlist(interesting_genes), scna_seus, plot_type = VlnPlot, group_by = "clusters")
  ),
  
  tar_target(
  	cluster_vlns_1q,
  	plot_putative_marker_across_samples(interesting_genes[["1q"]], debranched_seus_1q, plot_type = VlnPlot, group_by = "clusters", extension = "1q")
  ),
  
  tar_target(
  	cluster_vlns_2p,
  	plot_putative_marker_across_samples(interesting_genes[["2p"]], debranched_seus_2p, plot_type = VlnPlot, group_by = "clusters", extension = "2p")
  ),
  
  tar_target(
  	cluster_vlns_6p,
  	plot_putative_marker_across_samples(interesting_genes[["6p"]], debranched_seus_6p, plot_type = VlnPlot, group_by = "clusters", extension = "6p")
  ),
  
  tar_target(
  	cluster_vlns_16q,
  	plot_putative_marker_across_samples(interesting_genes[["16q"]], debranched_seus_16q, plot_type = VlnPlot, group_by = "clusters", extension = "16q")
  ),
  
  tar_target(
  	cluster_vlns_s2,
  	plot_putative_marker_across_samples(interesting_genes[["S2"]], scna_seus, plot_type = VlnPlot, group_by = "clusters", extension = "s2")
  ),
  
  tar_target(
  	cluster_vlns_s1,
  	plot_putative_marker_across_samples(interesting_genes[["S1"]], scna_seus, plot_type = VlnPlot, group_by = "clusters", extension = "s1")
  ),
  
  tar_target(factor_heatmaps_1q,
  					 plot_seu_gene_heatmap(debranched_seus_1q, large_clone_comparisons, scna_of_interest = "1q", w = 8, h = 12),
  					 pattern = map(debranched_seus_1q),
  					 iteration = "list"
  ),
  
  tar_target(factor_heatmaps_1q_combined, qpdf::pdf_combine(factor_heatmaps_1q, "results/factor_heatmaps_1q.pdf")),
  
  tar_target(factor_heatmaps_2p,
  					 plot_seu_gene_heatmap(debranched_seus_2p, large_clone_comparisons, scna_of_interest = "2p", w = 8, h = 12),
  					 pattern = map(debranched_seus_2p),
  					 iteration = "list"
  ),
  
  tar_target(factor_heatmaps_2p_combined, qpdf::pdf_combine(factor_heatmaps_2p, "results/factor_heatmaps_2p.pdf")),
  
  tar_target(factor_heatmaps_6p,
  					 plot_seu_gene_heatmap(debranched_seus_6p, large_clone_comparisons, scna_of_interest = "6p", w = 8, h = 12),
  					 pattern = map(debranched_seus_6p),
  					 iteration = "list"
  ),
  
  tar_target(factor_heatmaps_6p_combined, qpdf::pdf_combine(factor_heatmaps_6p, "results/factor_heatmaps_6p.pdf")),
  
  tar_target(factor_heatmaps_16q,
  					 plot_seu_gene_heatmap(debranched_seus_16q, large_clone_comparisons, scna_of_interest = "16q", w = 8, h = 12),
  					 pattern = map(debranched_seus_16q),
  					 iteration = "list"
  ),
  
  tar_target(factor_heatmaps_16q_combined, qpdf::pdf_combine(factor_heatmaps_16q, "results/factor_heatmaps_16q.pdf")),
  
  #factorized 
  tar_target(all_factor_heatmaps, qpdf::pdf_combine(list(factor_heatmaps_1q_combined, factor_heatmaps_2p_combined, factor_heatmaps_6p_combined, factor_heatmaps_16q_combined), "results/fig_s16.pdf")),
  
  tar_target(
    marker_gene_featureplots_by_cluster,
    plot_putative_marker_across_samples(unlist(interesting_genes), scna_seus, plot_type = FeaturePlot, group_by = "clusters")
  ),
  
  tar_target(
    cluster_markers,
    collect_all_markers(numbat_rds_files, "results/clusters.xlsx")
  ),
  tar_target(
    large_clusters,
    collect_clusters_from_seus(final_seus)
  ),
  tar_target(large_numbat_pdfs,
    convert_numbat_pngs(numbat_rds_files),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),

  # tar_target(large_plot_files,
  #   make_numbat_plot_files(numbat_rds_files, cluster_dictionary, clone_simplifications = large_clone_simplifications),
  #   pattern = map(numbat_rds_files),
  #   iteration = "list"
  # ),

  tarchetypes::tar_file(
    branch_dictionary_file,
    "data/branch_dictionary.csv"
  ),
  tar_target(
    branch_dictionary,
    pull_branches(branch_dictionary_file)
  ),
  # tarchetypes::tar_file_fast(
  # 	cluster_file,
  # 	"data/raw_cluster_ids.csv"
  # ),
  tar_target(cluster_file,
  					 "data/scna_cluster_order.csv", 
  					 format = "file_fast"),
  
  tar_target(cluster_orders,
    pull_cluster_orders(cluster_file)
  ),
  
  tar_target(overall_seus,
  					 merge_orders(debranched_seus,
  					 						 final_seus)
  ),
  
  tarchetypes::tar_file(
  	cluster_comparisons_file,
  	"data/cluster_comparisons_by_phase_for_disctinct_clones.csv"
  ),
  
  tar_target(cluster_comparisons,
  					 pull_cluster_comparisons(cluster_comparisons_file)),
  
  tar_target(cluster_comparisons_by_phase_for_disctinct_clones,
  					 make_cluster_comparisons_by_phase_for_disctinct_clones(cluster_comparisons, overall_seus),
  					 pattern = map(cluster_comparisons),
  					 iteration = "list"
  ),
  
  tar_target(filtered_large_plot_files,
    make_numbat_plot_files(numbat_rds_files, final_seus, cluster_dictionary, large_filter_expressions, large_clone_simplifications, extension = "_filtered"),
    pattern = map(numbat_rds_files, final_seus, cluster_dictionary, large_filter_expressions, large_clone_simplifications),
    iteration = "list"
  ),
  
  tar_target(unfiltered_seus,
    prep_unfiltered_seu(numbat_rds_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions, extension = "_unfiltered"),
    pattern = map(numbat_rds_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions),
    iteration = "list"
  ),
  
  tar_target(original_seus,
  					 c("output/seurat/SRR13884242_unfiltered_seu.rds", 
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
  
  
  tar_target(filtered_seus,
    filter_cluster_save_seu(numbat_rds_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions, cells_to_remove, extension = "_filtered", leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv"),
    pattern = map(numbat_rds_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions),
    iteration = "list"
  ),

  tar_target(
    final_seus,
    set_final_seus(interesting_samples)
  ),
  
  # tar_target(cc_plots_wo_arms,
  # 					 filter_cluster_save_seu(numbat_rds_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions, cells_to_remove, extension = "_filtered", leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv"),
  # 					 pattern = map(numbat_rds_files, cluster_dictionary, large_clone_simplifications, large_filter_expressions),
  # 					 iteration = "list"
  # ),
  
  tar_target(
  	cc_plots_wo_arms,
  	plot_phase_wo_arm(final_seus)
  ),
  
  
  
  # tar_target(debranched_seus,
  #   debranch_seus(final_seus, branch_dictionary)
  # ),
  
  tar_file(debranched_seu_files,
  				 c(
  				 	"SRR13884242" = "output/seurat/SRR13884242_filtered_seu.rds", 
  				 	"SRR13884243" = "output/seurat/SRR13884243_filtered_seu.rds", 
  				 	"SRR13884246_branch_5" = "output/seurat/SRR13884246_branch_5_filtered_seu.rds", 
  				 	"SRR13884246_branch_6" = "output/seurat/SRR13884246_branch_6_filtered_seu.rds", 
  				 	"SRR13884247_branch_6" = "output/seurat/SRR13884247_branch_6_filtered_seu.rds", 
  				 	"SRR13884247_branch_4" = "output/seurat/SRR13884247_branch_4_filtered_seu.rds", 
  				 	"SRR13884247_branch_5" = "output/seurat/SRR13884247_branch_5_filtered_seu.rds", 
  				 	"SRR13884248" = "output/seurat/SRR13884248_filtered_seu.rds", 
  				 	"SRR13884249" = "output/seurat/SRR13884249_filtered_seu.rds", 
  				 	"SRR14800534" = "output/seurat/SRR14800534_filtered_seu.rds", 
  				 	"SRR14800535" = "output/seurat/SRR14800535_filtered_seu.rds", 
  				 	"SRR14800536" = "output/seurat/SRR14800536_filtered_seu.rds", 
  				 	"SRR14800540_branch_2" = "output/seurat/SRR14800540_branch_2_filtered_seu.rds", 
  				 	"SRR14800540_branch_3" = "output/seurat/SRR14800540_branch_3_filtered_seu.rds", 
  				 	"SRR14800541_branch_4" = "output/seurat/SRR14800541_branch_4_filtered_seu.rds", 
  				 	"SRR14800541_branch_7" = "output/seurat/SRR14800541_branch_7_filtered_seu.rds", 
  				 	"SRR14800543_branch_3" = "output/seurat/SRR14800543_branch_3_filtered_seu.rds", 
  				 	"SRR14800543_branch_4" = "output/seurat/SRR14800543_branch_4_filtered_seu.rds", 
  				 	"SRR17960481" = "output/seurat/SRR17960481_filtered_seu.rds", 
  				 	"SRR17960484" = "output/seurat/SRR17960484_filtered_seu.rds", 
  				 	"SRR27187899" = "output/seurat/SRR27187899_filtered_seu.rds", 
  				 	"SRR27187902_branch_3" = "output/seurat/SRR27187902_branch_3_filtered_seu.rds", 
  				 	"SRR27187902_branch_4" = "output/seurat/SRR27187902_branch_4_filtered_seu.rds"
  				 )
  				 ),
  
  tar_target(debranched_seus,
  					 set_names(debranched_seu_files, str_remove(fs::path_file(debranched_seu_files), "_filtered_seu.rds"))
  					 ),
  
  
  tar_target(
    filtered_recurrent_genes,
    pull_common_markers(final_seus, mps[["Cancer"]])
  ),
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
  tar_target(
    filtered_recurrent_heatmap_file,
    qpdf::pdf_combine(filtered_recurrent_heatmaps, "results/filtered_recurrent_heatmaps.pdf")
  ),
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
  tar_target(
    regressed_recurrent_genes,
    pull_common_markers(regressed_seus, mps[["Cancer"]])
  ),
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
  
  tar_target(regressed_recurrent_heatmap_file,
    qpdf::pdf_combine(regressed_recurrent_heatmaps, "results/regressed_recurrent_heatmaps.pdf")
  ),

  # tar_target(filtered_scanpys,
  #            convert_seu_to_scanpy(final_seus),
  #            pattern = map(final_seus),
  #            iteration = "list"
  # ),

  # tar_target(regressed_scanpys,
  #            convert_seu_to_scanpy(regressed_seus),
  #            pattern = map(regressed_seus),
  #            iteration = "list"
  # ),

  tar_target(regressed_seus,
    regress_filtered_seu(final_seus),
    pattern = map(final_seus),
    iteration = "list"
  ),

  # tar_target(clustify_plots,
  #            plot_celltype_predictions(regressed_seus, interesting_samples, plae_ref, group.by = "SCT_snn_res.0.6"),
  #            pattern = map(regressed_seus, interesting_samples),
  #            iteration = "list"
  # ),

  tar_target(effect_of_filtering,
    plot_effect_of_filtering(unfiltered_seus, final_seus),
    pattern = map(unfiltered_seus, final_seus),
    iteration = "list"
  ),

  tar_target(regression_effect_plots,
    plot_effect_of_regression(final_seus, regressed_seus, w = 18, h = 12),
    pattern = map(final_seus, regressed_seus),
    iteration = "list"
  ),
  
  # sdfg
  tar_target(fig_s02_04, qpdf::pdf_combine(unlist(regression_effect_plots), "results/fig_s02_04.pdf")),
  
  tar_target(regression_ora_plots,
    ora_effect_of_regression(final_seus, regressed_seus),
    pattern = map(final_seus, regressed_seus),
    iteration = "list"
  ),
  tar_target(cin_score_plots,
    score_chrom_instability(unfiltered_seus),
    pattern = map(unfiltered_seus),
    iteration = "list"
  ),

  # tar_target(large_heatmaps,
  #   make_numbat_heatmaps_old(numbat_rds_files, large_filter_expressions, cluster_dictionary, p_min = 0.5, line_width = 0.1, extension = "_filtered"),
  #   pattern = map(numbat_rds_files, large_filter_expressions, cluster_dictionary),
  #   iteration = "list"
  # ),

  # unfiltered_numbat_heatmaps
  tar_target(fig_s03a_plots,
    make_numbat_heatmaps(original_seus, numbat_rds_files, p_min = 0.9, line_width = 0.1, extension = "_unfiltered"),
    pattern = map(original_seus),
    iteration = "list"
  ),
  
  # rb scna ideograms
  # tar_target(fig_s03c,
  # 					 plot_fig_s03c()),
  
  tar_target(fig_s03a, qpdf::pdf_combine(unlist(fig_s03a_plots), "results/unfiltered_heatmaps.pdf")),
  
  tar_target(
    debranched_samples,
    str_remove(fs::path_file(debranched_seus), "_filtered_seu.rds")
  ),
  
  tar_target(clustrees,
    make_clustrees_for_sample(debranched_seus, mylabel = debranched_samples, assay = "SCT"),
    pattern = map(debranched_seus, debranched_samples),
    iteration = "list"
  ),
  
  # tar_target(clustrees,
  # 					 make_clustrees_for_sample("output/seurat/integrated_1q/integrated_seu_1q_complete.rds", mylabel = "sdfg", assay = "SCT")
  # ),
  
  tarchetypes::tar_files(clustree_compilation, qpdf::pdf_combine(unlist(clustrees), "results/clustrees.pdf"), format = "file"),
  
  # filtered_numbat_heatmaps
  tar_target(fig_s13,
    make_numbat_heatmaps(final_seus, numbat_rds_files, p_min = 0.5, line_width = 0.1, extension = "_filtered"),
    pattern = map(final_seus, numbat_rds_files),
    iteration = "list",
    cue = tar_cue("always")
  ),
  
  tar_target(filtered_numbat_heatmaps_file, qpdf::pdf_combine(map_chr(fig_s13, 1), "results/filtered_heatmaps.pdf")),
  
  tar_target(filtered_large_scna_prob_file, qpdf::pdf_combine(map_chr(fig_s13, 2), "results/filtered_scna_probabilities.pdf")),
  
  tar_target(
    large_montage_pdfs,
    make_pdf_montages(filtered_large_plot_files, large_heatmaps)
  ),
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
  tar_target(all_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(debranched_seus, numbat_rds_files, large_clone_comparisons, cluster_orders, location = "all"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),

  # tar_target(cis_louvain_comparisons,
  #            find_diffex_bw_clones_for_each_louvain(final_seus, numbat_rds_files, large_clone_comparisons, location = "in_segment"),
  #            pattern = map(final_seus, numbat_rds_files, large_clone_comparisons),
  #            iteration = "list"
  # ),
  #
  # tar_target(trans_louvain_comparisons,
  #            find_diffex_bw_clones_for_each_louvain(final_seus, numbat_rds_files, large_clone_comparisons, location = "out_of_segment"),
  #            pattern = map(final_seus, numbat_rds_files, large_clone_comparisons),
  #            iteration = "list"
  # ),
  #
  # tar_target(all_louvain_comparisons,
  #            find_diffex_bw_clones_for_each_louvain(final_seus, numbat_rds_files, large_clone_comparisons, location = "all"),
  #            pattern = map(final_seus, numbat_rds_files, large_clone_comparisons),
  #            iteration = "list"
  # ),

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
  tar_target(
    volcano_all,
    make_volcano_diffex_clones(
      all_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_all.pdf",
      all_diffex_clones,
      "results/diffex_bw_clones_all.pdf"
    ),
  ),
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
    find_diffex_clones(debranched_seus, numbat_rds_files, large_clone_comparisons, location = "in_segment"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  
  tar_target(cis_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(debranched_seus, numbat_rds_files, large_clone_comparisons, cluster_orders, location = "in_segment"),
    pattern = map(debranched_seus),
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
    volcano_cis,
    make_volcano_diffex_clones(
      cis_diffex_clones_for_each_cluster,
      "results/diffex_bw_clones_per_cluster_large_in_segment.pdf",
      cis_diffex_clones,
      "results/diffex_bw_clones_large_in_segment.pdf"
    ),
  ),
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
  tar_target(trans_diffex_clones_for_each_cluster,
    find_diffex_bw_clones_for_each_cluster(debranched_seus, numbat_rds_files, large_clone_comparisons, cluster_orders, location = "out_of_segment"),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
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
  #   gse_plot_from_cluster_diffex(cis_diffex_clones_for_each_cluster),
  #   pattern = map(cis_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),

  # tar_target(
  #   large_out_of_segment_cluster_gse_plots,
  #   gse_plot_from_cluster_diffex(trans_diffex_clones_for_each_cluster),
  #   pattern = map(trans_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),

  # tar_target(
  #   large_in_segment_cluster_gse_plots,
  #   gse_plot_from_clone_diffex(cis_diffex_clones_for_each_cluster),
  #   pattern = map(cis_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),
  #

  # tar_target(
  #   large_out_of_segment_cluster_gse_plots,
  #   gse_plot_from_clone_diffex(trans_diffex_clones_for_each_cluster),
  #   pattern = map(trans_diffex_clones_for_each_cluster),
  #   iteration = "list"
  # ),

  # tar_target(large_diffex_bw_clusters_for_each_clone,
  #   find_diffex_bw_clusters_for_each_clone(final_seus, cluster_dictionary),
  #   pattern = map(final_seus, cluster_dictionary),
  #   iteration = "list"
  # ),

  tar_target(
    large_numbat_expression,
    retrieve_numbat_plot_type(large_numbat_pdfs, "exp_roll_clust.pdf")
  ),
  tar_target(
    large_numbat_expression_file,
    qpdf::pdf_combine(large_numbat_expression, "results/numbat_expression.pdf")
  ),
  tar_target(
    clone_trees,
    retrieve_numbat_plot_type(filtered_large_plot_files, "tree_filtered.pdf")
  ),
  
  tar_target(paired_clone_distribution_plots,
    calculate_clone_distribution(scna_seus, cluster_orders, pairwise = TRUE),
    pattern = map(scna_seus, cluster_orders),
    iteration = "list"
  ),
  
  tar_target(clone_pearls_1q,
  					 plot_clone_pearls(debranched_seus_1q, var_y = "phase_level"),
  					 pattern = map(debranched_seus_1q),
  					 iteration = "list"
  ),
  
  tar_target(clone_pearls_2p,
  					 plot_clone_pearls(debranched_seus_2p, var_y = "clusters"),
  					 pattern = map(debranched_seus_2p),
  					 iteration = "list"
  ),
  
  tar_target(cluster_terms,
  					 enrich_by_cluster(scna_seus),
  					 pattern = map(scna_seus),
  					 iteration = "list"
  ),
  
  # enrich_by_cluster
  
  tar_target(clone_pearls_6p,
  					 plot_clone_pearls(debranched_seus_6p, var_y = "clusters"),
  					 pattern = map(debranched_seus_6p),
  					 iteration = "list"
  ),
  
  tar_target(clone_pearls_16q,
  					 plot_clone_pearls(debranched_seus_16q, var_y = "phase_level"),
  					 pattern = map(debranched_seus_16q),
  					 iteration = "list"
  ),
  
  tar_target(unpaired_clone_distribution_plots,
    calculate_clone_distribution(final_seus, cluster_orders, pairwise = FALSE),
    pattern = map(final_seus, cluster_orders),
    iteration = "list"
  ),
  
  tar_target(
    cluster_scorecard,
    score_clusters_up_down(paired_clone_distribution_plots, interesting_samples)
  ),
  
  tar_target(clusters_and_markers,
    plot_seu_clusters_and_markers(debranched_seus, cluster_orders),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  
  tar_target(debranched_clone_tree_files,
    save_clone_tree_from_path(debranched_seus, numbat_rds_files, large_clone_simplifications, label = "_debranched_clone_tree", legend = FALSE, horizontal = FALSE),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  
  tar_target(debranched_clone_trees,
  					 plot_clone_tree_from_path(debranched_seus, numbat_rds_files, large_clone_simplifications, label = "_debranched_clone_tree", legend = FALSE, horizontal = FALSE),
  					 pattern = map(debranched_seus),
  					 iteration = "list"
  ),
  
  tar_target(
  	debranched_clone_trees_file,
  	qpdf::pdf_combine(debranched_clone_tree_files, "results/debranched_clone_trees.pdf")
  ),
  
  tar_target(final_clone_tree_files,
    save_clone_tree_from_path(final_seus, numbat_rds_files, large_clone_simplifications, label = "_clone_tree", legend = FALSE, horizontal = FALSE),
    pattern = map(final_seus),
    iteration = "list"
  ),
  
  tar_target(
  	final_clone_trees_file,
  	qpdf::pdf_combine(final_clone_trees_files, "results/final_clone_trees.pdf")
  ),
  
  tar_target(heatmap_collages,
    plot_seu_marker_heatmap_all_resolutions(debranched_seus, numbat_rds_files, large_clone_simplifications),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  
  tar_target(heatmap_collages_6p,
  					 plot_seu_marker_heatmap_all_resolutions(debranched_seus_6p, numbat_rds_files, large_clone_simplifications),
  					 pattern = map(debranched_seus_6p),
  					 iteration = "list"
  ),
  
  # # check effect of regression
  # tar_target(heatmap_collages,
  # 					 plot_seu_marker_heatmap_all_resolutions(regressed_seus, numbat_rds_files, large_clone_simplifications),
  # 					 pattern = map(regressed_seus),
  # 					 iteration = "list"
  # ),
  
  tar_target(annotated_heatmap_collages,
    plot_seu_marker_heatmap(debranched_seus, cluster_orders, numbat_rds_files, large_clone_simplifications),
    pattern = map(debranched_seus),
    iteration = "list"
  ),
  
  tar_target(silhouette_plots,
  					 calc_silhouette(debranched_seus),
  					 pattern = map(debranched_seus),
  					 iteration = "list"
  					 ),
  
  tar_target(debranched_seus_1q,
  	c(
  		"SRR13884246" = "output/seurat/SRR13884246_branch_5_filtered_seu_1q.rds", 
  		"SRR13884249" = "output/seurat/SRR13884249_filtered_seu_1q.rds", 
  		"SRR14800534" = "output/seurat/SRR14800534_filtered_seu_1q.rds", 
  		"SRR14800535" = "output/seurat/SRR14800535_filtered_seu_1q.rds", 
  		"SRR14800536" = "output/seurat/SRR14800536_filtered_seu_1q.rds"
  		)
  	),
  
  tar_target(debranched_seus_2p,
  		c(
  			# "SRR13884246" = "output/seurat/SRR13884246_branch_5_filtered_seu_2p.rds", 
  			"SRR13884247" = "output/seurat/SRR13884247_branch_6_filtered_seu.rds", 
  			"SRR13884248" = "output/seurat/SRR13884248_filtered_seu_2p.rds", 
  			"SRR13884249" = "output/seurat/SRR13884249_filtered_seu_2p.rds", 
  			"SRR17960481" = "output/seurat/SRR17960481_filtered_seu.rds", 
  			"SRR17960484" = "output/seurat/SRR17960484_filtered_seu_2p.rds" # too many changes at once to conclude anything
  		)
  		),
  
  tar_target(debranched_seus_6p,  
  		c(
  			"SRR13884247" = "output/seurat/SRR13884247_filtered_seu.rds",  
  			"SRR13884248" = "output/seurat/SRR13884248_filtered_seu_6p.rds", 
  			"SRR17960484" = "output/seurat/SRR17960484_filtered_seu_6p.rds",
  			"SRR27187899"= "output/seurat/SRR27187899_filtered_seu.rds"
  		)
  		),
  
  tar_target(debranched_seus_16q,  
  		c(
  			"SRR14800534" = "output/seurat/SRR14800534_filtered_seu_16q.rds", 
  			"SRR14800535" = "output/seurat/SRR14800535_filtered_seu_16q.rds", 
  			"SRR14800536" = "output/seurat/SRR14800536_filtered_seu_16q.rds"
  		)
  ),
  
  tar_target(scna_seus,  
  					 unlist(list(
  					 	"1q" = debranched_seus_1q, 
  					 	"2p" = debranched_seus_2p, 
  					 	"6p" = debranched_seus_6p,
  					 	"16q" = debranched_seus_16q
  					 ))
  ),
  
  
# medium resolution ------------------------------
  
  tar_target(collages_1q,
  					 plot_seu_marker_heatmap_by_scna(unlist(debranched_seus_1q), cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "1q"),
  					 pattern = map(debranched_seus_1q),
  					 iteration = "list"
  ),

  tar_target(collages_2p,
  					 plot_seu_marker_heatmap_by_scna(unlist(debranched_seus_2p), cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p"),
  					 pattern = map(debranched_seus_2p),
  					 iteration = "list"
  ),

	tar_target(fig_02,
						 plot_fig_02(cluster_orders),
						 ),
	
	tar_target(fig_s07,
						 plot_fig_s04_06(integrated_seus_1q, cluster_orders, "results/fig_s07.pdf")
	),

tar_target(fig_s08,
					 plot_fig_s04_06(integrated_seus_16q, cluster_orders, "results/fig_s08.pdf")
),

tar_target(fig_s10,
					 plot_fig_s07_08(integrated_seus_2p, cluster_orders, "results/fig_s10.pdf")
),

tar_target(fig_s12,
					 plot_fig_s07_08(integrated_seus_6p, cluster_orders, "results/fig_s12.pdf")
),

tar_target(fig_s0x,
					 plot_fig_s0x()
),

tar_target(fig_03,
					 plot_fig_03(cluster_orders),
),

tar_target(fig_04,
					 plot_fig_04_05(c("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", "output/seurat/SRR13884247_branch_6_filtered_seu.rds"), corresponding_clusters_enrichments[[6]], plot_path = "results/fig_04.pdf", integrated_seu_paths = integrated_seus_2p, cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "2p", width = 14, height = 10)
),

tar_target(fig_05,
					 plot_fig_04_05("output/seurat/integrated_6p/integrated_seu_6p_duo.rds", corresponding_clusters_enrichments[[7]], integrated_seu_paths = integrated_seus_6p, plot_path = "results/fig_05.pdf", cluster_order = cluster_orders, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "6p", width = 18, height = 10)
),

tar_target(fig_s06,
					 plot_fig_s06()), # Fig. S06: Differential expression comparisons between integrated 1q clusters of interest…

tar_target(integrated_seus_1q,
					 c(
					 	"SRR13884249_filtered_seu.rds" = "output/seurat/integrated_1q/SRR13884249_integrated_1q_filtered_seu.rds", 
					 	"SRR14800534_filtered_seu.rds" = "output/seurat/integrated_1q/SRR14800534_integrated_1q_filtered_seu.rds", 
					 	"SRR14800535_filtered_seu.rds" = "output/seurat/integrated_1q/SRR14800535_integrated_1q_filtered_seu.rds", 
					 	"SRR14800536_filtered_seu.rds" = "output/seurat/integrated_1q/SRR14800536_integrated_1q_filtered_seu.rds"
					 )
),

tar_target(integrated_seus_16q,
					 c(
					 	"SRR14800534_filtered_seu.rds" = "output/seurat/integrated_16q/SRR14800534_integrated_16q_filtered_seu.rds",
					 	"SRR14800535_filtered_seu.rds" = "output/seurat/integrated_16q/SRR14800535_integrated_16q_filtered_seu.rds",
					 	"SRR14800536_filtered_seu.rds" = "output/seurat/integrated_16q/SRR14800536_integrated_16q_filtered_seu.rds"
					 )
),

tar_target(integrated_seus_2p,
					 c(
					 	"SRR13884248_integrated_2p_filtered_seu.rds" = "output/seurat/integrated_2p/SRR13884248_integrated_2p_filtered_seu.rds", 
					 	"SRR17960484_integrated_2p_filtered_seu.rds" = "output/seurat/integrated_2p/SRR17960484_integrated_2p_filtered_seu.rds"
					 )
),

tar_target(integrated_seus_6p,
					 c(
					 	"SRR13884248_integrated_6p_filtered_seu.rds" = "output/seurat/integrated_6p/SRR13884248_integrated_6p_filtered_seu.rds", 
					 	"SRR17960484_integrated_6p_filtered_seu.rds" = "output/seurat/integrated_6p/SRR17960484_integrated_6p_filtered_seu.rds"
					 )
),

tar_target(fig_07a_input,
					 find_diffex_bw_clones_for_each_cluster(integrated_seus_1q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "1q+"),
					 pattern = map(integrated_seus_1q),
					 iteration = "list"
					 ),

tar_target(fig_07a,
					 plot_fig_07_08(fig_07a_input, plot_path = "results/fig_07a.pdf", plot_title = "fig_07a: 1q+ cluster diffex after integration", height = 5, width = 4),
					 iteration = "list"
),

tar_target(fig_07b_input,
					 find_diffex_bw_clones_for_each_cluster(debranched_seus_1q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "1q+"),
					 pattern = map(debranched_seus_1q),
					 iteration = "list"
),

tar_target(fig_07b,
					 plot_fig_07_08(fig_07b_input, plot_path = "results/fig_07b.pdf", plot_title = "fig_07b: 1q+ cluster diffex without integration", height = 5, width = 4),
					 iteration = "list"
),

tar_target(fig_07c,
					 plot_fig_07c() # Fig. S06: Differential expression comparisons between integrated 1q clusters of interest…
),

tar_target(fig_08a_input,
					 find_diffex_bw_clones_for_each_cluster(integrated_seus_16q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "16q-"),
					 pattern = map(integrated_seus_16q),
					 iteration = "list"
),

tar_target(fig_08a,
					 plot_fig_07_08(fig_08a_input, plot_title = "fig_08a: 16q- cluster diffex after integration", plot_path = "results/fig_08a.pdf", height = 5, width = 6),
					 iteration = "list"
),

tar_target(fig_08b_input,
					 find_diffex_bw_clones_for_each_cluster(debranched_seus_16q, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "all", scna_of_interest = "16q-"),
					 pattern = map(debranched_seus_16q),
					 iteration = "list"
),

tar_target(fig_08b,
					 plot_fig_07_08(fig_08b_input, plot_title = "fig_08b: 16q- cluster diffex without integration", plot_path = "results/fig_08b.pdf", p_adj_threshold = 0.05, height = 5, width = 6),
					 iteration = "list"
),


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

tar_target(fig_09,
					 plot_fig_09_10(corresponding_seus_2p, corresponding_seus, corresponding_clusters_diffex, corresponding_clusters_enrichments, recurrence_threshold = 3, plot_path = "results/fig_09.pdf", widths = rep(4, 3), heights = rep(8,3), common_seus = c("SRR13884248_filtered_seu_2p.rds", "SRR17960484_filtered_seu_2p.rds"))
),

tar_target(
	corresponding_seus_6p,
	c(
		"SRR13884248_filtered_seu_6p.rds" = "output/seurat/SRR13884248_filtered_seu_6p.rds", 
		"SRR17960484_filtered_seu_6p.rds" = "output/seurat/SRR17960484_filtered_seu_6p.rds"
	)
),

tar_target(fig_10,
					 plot_fig_09_10(corresponding_seus_6p, corresponding_seus, corresponding_clusters_diffex, corresponding_clusters_enrichments, recurrence_threshold = 2, plot_path = "results/fig_10.pdf", widths = rep(4,3), heights = c(12, 4, 12), common_seus = c("SRR13884248_filtered_seu_6p.rds", "SRR17960484_filtered_seu_6p.rds"))
),

  tar_target(collages_6p,
  					 plot_seu_marker_heatmap_by_scna(unlist(debranched_seus_6p), cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "6p"),
  					 pattern = map(debranched_seus_6p),
  					 iteration = "list"
  ),

  tar_target(collages_16q,
  					 plot_seu_marker_heatmap_by_scna(unlist(debranched_seus_16q), cluster_orders, numbat_rds_files, large_clone_simplifications, rb_scna_samples = rb_scna_samples, large_clone_comparisons = large_clone_comparisons, scna_of_interest = "16q"),
  					 pattern = map(debranched_seus_16q),
  					 iteration = "list"
  ),

  tar_target(collected_scna_collages,
  					 qpdf::pdf_combine(
  					 	unlist(list(collages_1q, collages_2p, collages_6p, collages_16q)), 
  					 	"results/fig_s17.pdf")
  					 ),

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

# 2p ------------------------------

tar_target(states_dictionary_2p,
					 make_corresponding_states_dictionary(tibble::tribble(
					 	~file_name, ~w_scna, ~wo_scna, ~scna_of_interest,
					 	"seurat_2p_integrated_duo.rds",     "g1_1-g1_4",      "g1_9", "2p",
					 	"seurat_2p_integrated_duo.rds",     "g1_4",      "g1_9", "2p",
					 	"seurat_2p_integrated_duo.rds",     "g1_1",      "g1_9", "2p",
					 	"seurat_2p_integrated_duo.rds",     "g1_1",      "g1_4", "2p"
					 ))
),

tar_target(corresponding_clusters_diffex_2p,
					 find_diffex_clusters_between_corresponding_states("output/seurat/integrated_2p/seurat_2p_integrated_duo.rds", states_dictionary_2p, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all")
),

tar_target(corresponding_clusters_volcanos_2p,
					 plot_corresponding_clusters_diffex_volcanos(corresponding_clusters_diffex_2p, "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds")
),

tar_target(corresponding_clusters_enrichments_2p,
					 plot_corresponding_enrichment(corresponding_clusters_diffex_2p, "output/seurat/integrated_2p/seurat_2p_integrated_duo.rds"),
),

# 6p ------------------------------

tar_target(states_dictionary_6p,
					 make_corresponding_states_dictionary(tibble::tribble(
					 	~file_name, ~w_scna, ~wo_scna, ~scna_of_interest,
					 	"SRR13884247_filtered_seu.rds",   "g1_0-g1_1",    "g1_4-g1_5", "6p",
					 	"SRR17960484_filtered_seu_6p.rds",     "g1_0",      "g1_1", "6p"
					 ))
),

tar_target(corresponding_state_6p_seus,
					 list(
					 	"output/seurat/SRR13884247_filtered_seu.rds",
					 	"output/seurat/SRR17960484_filtered_seu_6p.rds"
					 )),

tar_target(corresponding_clusters_diffex_6p,
					 find_diffex_clusters_between_corresponding_states(unlist(corresponding_state_6p_seus), states_dictionary_6p, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all"),
					 pattern = map(corresponding_state_6p_seus, states_dictionary_6p),
					 iteration = "list"
),

tar_target(corresponding_clusters_volcanos_6p,
					 plot_corresponding_clusters_diffex_volcanos(corresponding_clusters_diffex_6p, unlist(corresponding_state_6p_seus)),
					 pattern = map(corresponding_clusters_diffex_6p, corresponding_state_6p_seus),
					 iteration = "list"
),

tar_target(corresponding_clusters_enrichments_6p,
					 plot_corresponding_enrichment(corresponding_clusters_diffex_6p, unlist(corresponding_state_6p_seus)),
					 pattern = map(corresponding_clusters_diffex_6p, corresponding_state_6p_seus),
					 iteration = "list"
),

# rest ------------------------------

	tar_target(corresponding_clusters_diffex,
						 find_diffex_clusters_between_corresponding_states(unlist(corresponding_seus), corresponding_states_dictionary, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all"),
						 pattern = map(corresponding_seus, corresponding_states_dictionary),
						 iteration = "list"
	),

tar_target(corresponding_clusters_volcanos,
					 plot_corresponding_clusters_diffex_volcanos(corresponding_clusters_diffex, unlist(corresponding_seus)),
					 pattern = map(corresponding_clusters_diffex, corresponding_seus),
					 iteration = "list"
),

tar_target(corresponding_clusters_heatmaps,
					 plot_corresponding_clusters_diffex_heatmaps(corresponding_clusters_diffex, unlist(corresponding_seus), corresponding_states_dictionary, large_clone_comparisons, numbat_rds_files = numbat_rds_files, location = "all"),
					 pattern = map(corresponding_clusters_diffex, corresponding_seus, corresponding_states_dictionary),
					 iteration = "list"
),

tar_target(corresponding_clusters_enrichments,
					 plot_corresponding_enrichment(corresponding_clusters_diffex, unlist(corresponding_seus)),
					 pattern = map(corresponding_clusters_diffex, corresponding_seus),
					 iteration = "list"
),

tar_target(clone_cc_plots_by_scna_1q,
					 clone_cc_plots_by_scna(debranched_seus_1q, scna_of_interest = "1q", large_clone_comparisons = large_clone_comparisons),
),

tar_target(clone_cc_plots_by_scna_16q,
					 clone_cc_plots_by_scna(debranched_seus_16q, scna_of_interest = "16q", large_clone_comparisons = large_clone_comparisons),
),

  
  tar_target(diffex_2p_g1, # 2p+ diffex; diffex 2p+
  					 find_diffex_clones_in_phase(corresponding_seus_2p, phase = "g1", scna_of_interest = "2p", numbat_rds_files, large_clone_comparisons, location = "all"),
  					 pattern = map(corresponding_seus_2p),
  					 iteration = "list"
  ),

tar_target(enrichment_2p_g1,
					 compare_enrichment(diffex_2p_g1)
),
  
  tar_target(diffex_6p_g1,
  					 find_diffex_clones_in_phase(debranched_seus_6p, phase = "g1", scna_of_interest = "6p", numbat_rds_files, large_clone_comparisons, location = "all"),
  					 pattern = map(debranched_seus_6p),
  					 iteration = "list"
  ),

tar_target(enrichment_6p_g1,
					 compare_enrichment(diffex_6p_g1)
),
  
  tar_target(integrated_seus,
  					 list(
  					 	"diploid_v_16q" = "output/seurat/SRR1_diploid_v_16q_filtered_seu.rds",
  					 	"16q_v_16q-1q" = "output/seurat/SRR1_16q_v_16q-1q_filtered_seu.rds",
  					 "diploid_v_16q-1q" = "output/seurat/SRR1_diploid_v_16q-1q_filtered_seu.rds"
  					 )
  					 ),
  
  tar_target(annotated_integrated_heatmap_collages,
  					 plot_seu_marker_heatmap_integrated(integrated_seus),
  					 pattern = map(integrated_seus),
  					 iteration = "list"
  ),

  tar_target(
    collage_compilation,
    qpdf::pdf_combine(
    	annotated_heatmap_collages, 
    	"results/heatmap_collages.pdf")
  ),
  
  tar_target(
  	collage_compilation_all_resolutions,
  	qpdf::pdf_combine(heatmap_collages, "results/heatmap_collages_all_resolutions.pdf")
  ),
  
  tarchetypes::tar_file(divergent_cluster_file, "data/clustree_divergent_clusters.csv"),
  tar_target(
    clustree_tables,
    pull_clustree_tables(clustrees, divergent_cluster_file)
  ),
  tar_target(clustree_diffexes,
    find_all_diffex_from_clustree(clustree_tables, debranched_seus, clone_comparisons = large_clone_comparisons),
    pattern = map(clustree_tables),
  ),
  tar_target(
    clustree_cis_changes,
    find_candidate_cis_in_clustree_diffexes(clustree_diffexes)
    # pattern = map(clustree_diffexes),
  ),
  tar_target(
    clustree_trans_changes,
    find_candidate_trans_in_clustree_diffexes(clustree_diffexes)
    # pattern = map(clustree_diffexes),
  ),
  tar_target(
    clustree_all_changes,
    find_candidate_all_in_clustree_diffexes(clustree_diffexes)
    # pattern = map(clustree_diffexes),
  ),

  # tar_target(clustree_diffex_test,
  # 					 find_all_diffex_from_clustree(clustree_tables[[1]], final_seus, clone_comparisons = large_clone_comparisons),
  # ),

  # tar_target(collage_compilation,
  #          qpdf::pdf_combine(heatmap_collages, "results/heatmap_collages.pdf")),

  # tar_target(large_clone_dist, retrieve_numbat_plot_type(large_plot_files, "exp_roll_clust.pdf")),

  # tar_target(large_numbat_expression_pdf, qpdf::pdf_combine(large_numbat_expression, "results/numbat_sridhar_large/large_numbat_expression.pdf")),
  #
  tar_target(large_numbat_heatmap_files, retrieve_numbat_plot_type(large_numbat_pdfs, "panel_2.pdf")),

  tar_target(large_numbat_heatmap_pdf, qpdf::pdf_combine(large_numbat_heatmap_files, "results/large_numbat_heatmaps.pdf")),
  #
  # tar_target(large_numbat_bulk_clones_final, retrieve_numbat_plot_type(large_numbat_pdfs, "bulk_clones_final.pdf")),
  #
  # tar_target(large_numbat_bulk_clones_final_pdf, qpdf::pdf_combine(large_numbat_bulk_clones_final, "results/numbat_sridhar_large/large_bulk_clones_final.pdf")),
  #
  # # large_numbat_sample_pdfs ------------------------------
  # tar_target(large_numbat_sample_pdfs,
  #   reroute_done_to_results_pdf(numbat_rds_files, "_large"),
  #   pattern = map(numbat_rds_files),
  #   iteration = "list",
  #   format = "file"
  # ),

  tar_target(large_numbat_sample_pdfs,
    reroute_done_to_results_pdf(numbat_rds_files, "_large"),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),
  tar_target(
    filtering_timelines,
    qpdf::pdf_combine(dir_ls("results", regexp = ".*SRR[0-9]*_filtering_timeline.pdf"), "results/filtering_timelines.pdf")
  ),

  # tar_target(
  #   oncoprint_input_by_scna_for_each_cluster_unfiltered,
  #   make_oncoprint_diffex_unfiltered(large_filter_expressions, cluster_dictionary, interesting_samples, cis_diffex_clones_for_each_cluster, trans_diffex_clones_for_each_cluster, large_clone_comparisons, n_slice = 20)
  # ),

  tar_target(
    num_diffex_clone_scna_tally,
    tally_num_diffex(oncoprint_input_by_scna_unfiltered)
  ),

  # tar_target(
  #   num_diffex_cluster_scna_tally,
  #   tally_num_diffex(oncoprint_input_by_scna_for_each_cluster_unfiltered)
  # ),
  
  tar_target(
  	rb_scna_samples,
  	list(
  	"1q" = c("SRR13884246", "SRR13884249", "SRR14800534", "SRR14800535", "SRR14800536"),
  	"2p"= c("SRR13884246", "SRR13884247", "SRR13884248", "SRR13884249", "SRR17960481", "SRR17960484"),
  	"6p"= c("SRR13884247", "SRR13884248", "SRR17960484"),
  	"16q" = c("SRR14800534", "SRR14800535", "SRR14800536")
  )),
  
  tar_target(
  	unfiltered_oncoprint_input_by_scna,
  	make_oncoprint_diffex(large_filter_expressions, cluster_dictionary, debranched_ids, cis_diffex_clones, trans_diffex_clones, all_diffex_clones, large_clone_comparisons, rb_scna_samples, n_slice = 20)
  ),

  tar_target(
  	oncoprint_input_by_scna,
    filter_oncoprint_diffex(unfiltered_oncoprint_input_by_scna, oncoprint_settings)
  ),
  
  tar_file(oncoprint_settings_file, "data/oncoprint_settings.tsv"),
  
  tar_file(oncoprint_settings_per_cluster_file, "data/oncoprint_settings_by_cluster.tsv"),
  
  tar_target(oncoprint_settings, read_tsv(oncoprint_settings_file)),
  
  tar_target(oncoprint_settings_per_cluster, read_tsv(oncoprint_settings_per_cluster_file)),

  tar_target(
    unfiltered_oncoprint_input_by_scna_for_each_cluster,
    make_oncoprint_diffex(large_filter_expressions, cluster_dictionary, debranched_ids, cis_diffex_clones_for_each_cluster, trans_diffex_clones_for_each_cluster, all_diffex_clones_for_each_cluster, large_clone_comparisons, rb_scna_samples, by_cluster = TRUE, n_slice = 20)
  ),
  
  tar_target(
  	oncoprint_input_by_scna_for_each_cluster,
  	filter_oncoprint_diffex(unfiltered_oncoprint_input_by_scna_for_each_cluster, oncoprint_settings_by_cluster)
  ),
  
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
  tar_target(rod_rich_samples,
    score_samples_for_rod_enrichment(numbat_rds_files),
    pattern = map(numbat_rds_files),
    iteration = "list"
  ),

tar_target(table_06,
					 rod_rich_samples),


  tar_target(celltype_rich_samples,
    score_samples_for_celltype_enrichment(unfiltered_seus, final_seus, celltype_markers),
    pattern = map(unfiltered_seus, final_seus),
    iteration = "list"
  ),
  tar_target(
    oncoprint_enrich_clones_plots_gobp,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_gobp,
      cis_plot_file = "results/enrichment/cis_enrichment_plots_by_clone_gobp.pdf",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_clone_gobp.pdf",
      cis_table_file = "results/enrichment/cis_enrichment_tables_by_clone_gobp.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_clone_gobp.xlsx"
    )
  ),
  tar_target(
    oncoprint_enrich_clones_plots_hallmark,
    compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_hallmark,
      cis_plot_file = "results/enrichment/cis_enrichment_plots_by_clone_hallmark.pdf",
      trans_plot_file = "results/enrichment/trans_enrichment_plots_by_clone_hallmark.pdf",
      cis_table_file = "results/enrichment/cis_enrichment_tables_by_clone_hallmark.xlsx",
      trans_table_file = "results/enrichment/trans_enrichment_tables_by_clone_hallmark.xlsx"
    )
  ),
  tar_target(
    subtype_markers,
    pull_subtype_genes(supp_excel = "data/Liu_Radvanyi_2022_supp_data/41467_2021_25792_MOESM6_ESM.xlsx")
  ),
  tar_target(
    liu_lu_supp_data,
    read_liu_lu_supp_tables()
  ),
  tar_target(
    mps,
    read_mps("/dataVolume/storage/Homo_sapiens/3ca/ITH_hallmarks/MPs_distribution/MP_list.RDS")
  ),
  tar_target(
    stemness_markers,
    pull_stem_cell_markers()
  ),

  # clusters enrichment ------------------------------
  # # enrichment by cluster msigdb gobp
  # tar_target(
  #   oncoprint_enrich_clusters_gobp,
  #   enrich_oncoprints_clusters(large_filter_expressions,
  #     cluster_dictionary,
  #     debranched_ids,
  #     cis_diffex_clones_for_each_cluster,
  #     trans_diffex_clones_for_each_cluster,
  #     all_diffex_clones_for_each_cluster,
  #     large_clone_comparisons,
  #     gene_set = "gobp"
  #   )
  # ),

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

  # # enrichment plots and tables by cluster gobp
  # tar_target(
  #   oncoprint_enrich_clusters_plots_gobp,
  #   compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clusters_gobp,
  #     cis_plot_file = "results/enrichment/cis_enrichment_plots_by_cluster_gobp",
  #     trans_plot_file = "results/enrichment/trans_enrichment_plots_by_cluster_gobp",
  #     cis_table_file = "results/enrichment/cis_enrichment_tables_by_cluster_gobp.xlsx",
  #     trans_table_file = "results/enrichment/trans_enrichment_tables_by_cluster_gobp.xlsx",
  #     by_cluster = TRUE
  #   )
  # ),

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
  
  tar_target(oncoprint_input_by_region,
    inspect_oncoprints(oncoprint_input_by_scna$cis, oncoprint_input_by_scna$trans, oncoprint_input_by_scna$all)
  ),

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

	tar_target(chosen_resolution_seus,
						 assign_designated_phase_clusters(scna_seus, cluster_orders, resolution_dictionary),
						 pattern = map(scna_seus),
						 iteration = "list"
	),
  
  tar_target(oncoprint_plots,
    make_oncoprint_plots(oncoprint_input_by_scna, debranched_clone_trees, oncoprint_settings, label = "_by_clone")
  ),
  
  tar_target(
  	oncoprint_plots_by_cluster,
  	make_oncoprint_plots(oncoprint_input_by_scna_for_each_cluster$cis, oncoprint_input_by_scna_for_each_cluster$trans, oncoprint_input_by_scna_for_each_cluster$all, debranched_clone_trees, oncoprint_settings_by_cluster, label = "_by_cluster", p_val_threshold = 1)
  ),
  
  tar_target(stachelek_score_plots,
    score_stachelek(final_seus, oncoprint_input_by_scna),
    pattern = map(final_seus),
    iteration = "list"
  ),

  #
  # tar_target(
  #   large_cluster_markers,
  #   collect_all_markers(numbat_rds_files, "results/clusters.xlsx")
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
