# Smoke test for the lowest-resolution-first hypoxia split.
#
#   SRX10264518 — flags at low res AND carries residual hypoxia that only
#                 separates at higher resolution, so it exercises the
#                 confirmatory drop.
#   SRX10264520 — control.
#
# Verifies: round=1 sweep rows over the full grid + round=2 confirmatory rows,
# r_flag is the lowest flagging resolution, the three stage PDFs are emitted, and
# the persisted low object still carries SCT_snn_res.0.6 (the downstream contract).
options(warn = 1)  # print warnings as they happen; Rscript otherwise defers them to
                   # the end and collapses them into "There were N warnings", which
                   # hid the cause of the skipped collages across three runs
suppressMessages(library(tidyverse))  # the plot fns use unqualified tidyr/dplyr verbs
library(seuratTools)   # subset_seu_by_expression() calls find_all_markers()
devtools::load_all("/project2/cobrinik_1090/rpkgs/numbat_helpers")

samples <- c("SRX10264518", "SRX22868102", "SRX10264520")

for (s in samples) {
  in_path <- glue::glue("output/seurat/{s}_seu_hypoxia.rds")
  message("\n===== ", s, " =====")

  paths <- split_hypoxia_by_clusters(
    in_path,
    resolutions       = seq(0.2, 1.2, by = 0.2),
    recluster_step    = 0.4,
    confirmatory_drop = TRUE,
    write_diagnostics = TRUE,
    min_hyp_markers   = 2,
    max_dom_frac      = 0.70,
    split_assay = "gene", low_assay = "gene", high_assay = "SCT"
  )
  print(paths)

  # --- split log: the sweep + the confirmatory pass ---------------------------
  log_path <- glue::glue("results/hypoxia_cluster_split/{s}_hypoxia_split_log.csv")
  if (file.exists(log_path)) {
    lg <- readr::read_csv(log_path, show_col_types = FALSE)
    message("log rows by round x resolution:")
    print(dplyr::summarise(dplyr::group_by(lg, round, resolution),
                           n_clusters = dplyr::n(),
                           n_flagged  = sum(flagged),
                           .groups = "drop"))
    flagged <- dplyr::filter(lg, flagged)
    if (nrow(flagged) > 0) {
      message("r_flag (lowest flagging resolution, round 1): ",
              min(flagged$resolution[flagged$round == 1]))
      message("confirmatory (round 2) flagged clusters: ",
              sum(flagged$round == 2))
    } else {
      message("nothing flagged -> r_flag = NA, all cells low")
    }
  } else {
    message("NO SPLIT LOG at ", log_path)
  }

  # --- the three diagnostic stage PDFs ---------------------------------------
  pdfs <- fs::dir_ls("results", glob = glue::glue("*{s}_hypoxia_*heatmap_phase_scatter_patchwork.pdf"))
  message("stage PDFs (", length(pdfs), "):")
  print(as.character(pdfs))

  # --- downstream contract: low object must still carry SCT_snn_res.0.6 -------
  low <- paths[["low"]]
  if (!is.na(low) && file.exists(low)) {
    lseu <- readRDS(low)
    message("low object: ", ncol(lseu), " cells; SCT_snn_res.0.6 present: ",
            "SCT_snn_res.0.6" %in% colnames(lseu@meta.data))
    message("res columns: ",
            paste(grep("_res\\.", colnames(lseu@meta.data), value = TRUE),
                  collapse = ", "))
    rm(lseu); gc()
  }
}

message("\nsmoke test done")
