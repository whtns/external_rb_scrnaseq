# File-tracking targets, config loading, sample definitions, and reference data.
# Defines: pipeline_targets_inputs (spliced into tar_plan in _targets.R)

pipeline_targets_inputs <- c(

  # --- numbat / seurat file tracking ---

  tarchetypes::tar_files(numbat_rds_all, retrieve_numbat_rds_files("output/numbat_sridhar/"), format = "file"),
  tarchetypes::tar_files(seus_all, retrieve_seus("output/seurat/"), format = "file"),
  tarchetypes::tar_files(numbat_rds_files, retrieve_numbat_rds_files("output/numbat_sridhar/", interesting_samples), format = "file"),
  tarchetypes::tar_files(numbat_rds_filtered_files, retrieve_numbat_rds_files("output/numbat_sridhar_filtered/", interesting_samples), format = "file"),
  # Subset of Seurat files for interesting_samples; entries are named by SRR ID.
  tarchetypes::tar_files(seus_interesting, retrieve_seus("output/seurat/", interesting_samples), format = "file"),

  list(

    # --- sample definitions ---
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
        "SRR27187902",
        # "SRR14800538",  # no numbat output; needs full pipeline re-run
        "SRR14800539",
        "SRR17960482",
        "SRR17960483",
        "SRR17960480",
        # "SRR14800542",  # numbat stopped: no CNV after entropy filter; needs re-run with higher max_entropy
        # "SRR13633760",  # numbat stopped: no CNV after entropy filter; needs re-run with higher max_entropy
        "SRR13633762",
        "SRR13884240",
        "SRR13884241"
      )
    ),

    # tar_target(interesting_samples,
    #   c(
    #     "SRR13884242", "SRR13884243", "SRR13884244", "SRR13884245",
    #     "SRR13884246", "SRR13884247", "SRR13884248", "SRR13884249",
    #     "SRR14800534", "SRR14800535", "SRR14800536",
    #     "SRR14800540", "SRR14800541", "SRR14800543",
    #     "SRR17960481", "SRR17960484",
    #     "SRR27187899", "SRR27187900", "SRR27187901", "SRR27187902"
    #   )
    # ),

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
    tar_target(excluded_samples,
      c("SRR14800541")  # no diploid ancestral clone; only diploid are B2M and APOE
    ),

    # IDs for per-branch debranched Seurat objects.
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

    # Which samples carry each SCNA.
    tar_target(rb_scna_samples,
      list(
        "1q"  = c("SRR13884246", "SRR13884249", "SRR14800534", "SRR14800535", "SRR14800536"),
        "2p"  = c("SRR13884246", "SRR13884247", "SRR13884248", "SRR13884249", "SRR17960481", "SRR17960484"),
        "6p"  = c("SRR13884247", "SRR13884248", "SRR17960484"),
        "16q" = c("SRR14800534", "SRR14800535", "SRR14800536")
      )
    ),

    # Durable, LLM-readable mapping of branch IDs to SRR IDs.
    tarchetypes::tar_file(branch_registry_file, "config/branch_registry.csv"),
    tar_target(branch_registry,
      readr::read_csv(branch_registry_file, show_col_types = FALSE)
    ),
    tar_target(branch_registry_validated,
      {
        missing_srr <- setdiff(interesting_samples, unique(branch_registry$srr_id))
        stopifnot(
          !anyNA(branch_registry$branch_id),
          !anyNA(branch_registry$srr_id),
          !anyNA(branch_registry$seu_path),
          !anyDuplicated(branch_registry$branch_id),
          !anyDuplicated(branch_registry$seu_path),
          all(startsWith(branch_registry$branch_id, branch_registry$srr_id)),
          length(missing_srr) == 0
        )
        branch_registry
      }
    ),
    tar_target(branch_registry_manifest_csv,
      {
        out <- "results/branch_registry_manifest.csv"
        readr::write_csv(branch_registry_validated, out)
        out
      },
      format = "file"
    ),
    tar_target(branch_registry_manifest_jsonl,
      {
        out <- "results/branch_registry_manifest.jsonl"
        lines <- apply(branch_registry_validated, 1, function(row) {
          jsonlite::toJSON(as.list(row), auto_unbox = TRUE)
        })
        readr::write_lines(lines, out)
        out
      },
      format = "file"
    ),

    # --- clone comparison configs ---

    tar_target(large_clone_comparisons,
      yaml::read_yaml(here::here("config/large_clone_comparisons.yaml"))
    ),
    tar_target(large_clone_simplifications,
      yaml::read_yaml(here::here("config/large_clone_simplifications.yaml"))
    ),
    tar_target(large_filter_expressions,
      yaml::read_yaml(here::here("config/large_filter_expression.yaml"))
    ),
    tar_target(clone_comparison_table,
      tabulate_clone_comparisons(large_clone_comparisons)
    ),

    # --- cluster annotation configs ---

    tar_target(cluster_dictionary, read_cluster_dictionary("data/cluster_dictionary.tsv")),

    tarchetypes::tar_file(branch_dictionary_file, "data/branch_dictionary.csv"),
    tar_target(branch_dictionary, pull_branches(branch_dictionary_file)),

    tar_target(cluster_file, "data/scna_cluster_order.csv", format = "file"),
    tar_target(cluster_orders, pull_cluster_orders(cluster_file)),
    tar_target(cluster_orders_sqlite,
      encode_cluster_order_to_hash_table(cluster_orders),
      cue = tar_cue("always")
    ),

    tarchetypes::tar_file(cluster_comparisons_file, "data/cluster_comparisons_by_phase_for_disctinct_clones.csv"),
    tar_target(cluster_comparisons, pull_cluster_comparisons(cluster_comparisons_file)),

    tarchetypes::tar_file(divergent_cluster_file, "data/clustree_divergent_clusters.csv"),

    # --- oncoprint configs ---

    tar_file(oncoprint_settings_file, "data/oncoprint_settings.tsv"),
    tar_file(oncoprint_settings_per_cluster_file, "data/oncoprint_settings_by_cluster.tsv"),
    tar_target(oncoprint_settings, read_tsv(oncoprint_settings_file)),
    tar_target(oncoprint_settings_per_cluster, read_tsv(oncoprint_settings_per_cluster_file)),

    # --- QC and metadata ---

    tar_target(cells_to_remove, read_cells_to_remove("data/cells_to_remove_final.xlsx")),
    tar_target(total_metadata, read_tsv("data/metadata.tsv")),
    tar_target(hallmark_gene_sets, msigdbr(species = "Homo sapiens", category = "H")),

    # --- gene lists ---

    tar_target(celltype_markers,
      list(
        "cones"                = c("PDE6H", "ARR3", "MYL4", "RCVRN", "GNAT2"),
        "rods"                 = c("GNGT1", "RHO", "CNGA1", "SAG", "PDE6G"),
        "Müller glia"          = c("TF", "NPVF", "RLBP1", "RRH", "RGR"),
        "retinal astrocytes"   = c("TRH", "MLC1", "NDP", "PAX5", "PAX2"),
        "microglia"            = c("RGS1", "C1QB", "C1QA", "C1QC", "CCL3"),
        "bipolar cells"        = c("CABP5", "TMEM215", "GABRA1", "CA10", "PCP2"),
        "retinal ganglion cells" = c("NEFL", "NEFM", "PVALB", "SNCG", "STMN2"),
        "amacrine cells"       = c("C1QL2", "GAD2", "CARTPT", "GAD1", "SLC32A1"),
        "RPE"                  = c("TTR", "RPE65", "RGR", "RLBP1", "BEST1")
      )
    ),

    tar_target(interesting_genes,
      list(
        "1q"         = c("ANP32E", "ASH1L", "ASPM", "CENPF", "CKS1B", "CNIH4", "CRABP2", "ENAH", "ENO1", "IPO9", "KIF14", "NEK2", "NENF", "NUF2", "PSMD4", "RXRG", "UBE2T"),
        "2p"         = c("CEP68", "MEIS1", "MDH1", "OST4", "PDIA6", "POLE4", "RAB1A", "SNRPG", "SOX11"),
        "6p"         = c("ARID1B", "CASP8AP2", "CLIC1", "CUTA", "DDX39B", "DEK", "DST", "GLO1", "HDAC2", "HDDC2", "HMGA1", "LIN28B", "MCM3", "MNDA", "PRDM1", "SOX4"),
        "16q"        = c("CHD9", "CNGB1", "CYBA", "MT1E", "MT1X", "MT2A"),
        "S1"         = c("ARR3", "CRX", "PDC"),
        "S2"         = c("EBF3", "DCX", "ROBO1", "SOX11", "GAP43", "PCDHB10", "STMN2", "NEFM", "POU4F2", "TFF1", "CD24"),
        "interesting" = c("SETD6", "RBL2", "DDX19A", "EZH2", "RB1", "MDM4", "ESRRG"),
        "cell_cycle"  = c("CCNA2", "CCNB1")
      )
    ),

    tar_target(subtype_markers,
      pull_subtype_genes(supp_excel = "data/Liu_Radvanyi_2022_supp_data/41467_2021_25792_MOESM6_ESM.xlsx")
    ),
    tar_target(liu_lu_supp_data, read_liu_lu_supp_tables()),
    tar_target(mps, read_mps("/dataVolume/storage/Homo_sapiens/3ca/ITH_hallmarks/MPs_distribution/MP_list.RDS")),
    tar_target(stemness_markers, pull_stem_cell_markers())

  )
)
