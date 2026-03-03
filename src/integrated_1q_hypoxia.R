#!/usr/bin/env Rscript

source("packages.R")
lapply(list.files("./R", full.names = TRUE), source)

integration_by_scna_clones <- function(seu_paths, scna_of_interest = "1q", filter_expression = , clone_comparisons){ {
	seus <- seu_paths |> 
		set_names()  |> 
		set_names(str_extract, "SRR[0-9]*")
	
	seus <- seus |> 
		map(readRDS)
	
	subset_seu_by_clones <- function(seu, sample_id, scna_of_interest = "1q", clone_comparisons){
		clone_comparisons <- names(clone_comparisons[[sample_id]])
		retained_clones <- clone_comparisons %>%
			str_extract("[0-9]_v_[0-9]") %>%
			str_split("_v_", simplify = TRUE)

		mode(retained_clones) <- "integer"

		retained_clones <- retained_clones[which.min(rowSums(retained_clones)),]

		seu <- seu[, seu$clone_opt %in% retained_clones]

		return(seu)
	}

	seus <- imap(seus,
		subset_seu_by_clones, scna_of_interest = scna_of_interest, clone_comparisons = clone_comparisons)
	
	integrated_seu <- seuratTools::integration_workflow(seus)

	seu_path <- tempfile(pattern = paste0("integrated_", scna_of_interest, "_"), tmpdir = "output/seurat/", fileext = "_filtered_seu.rds")

	saveRDS(integrated_seu, seu_path)

	integrated_seu <- ScaleData(integrated_seu)
	integrated_seu <- seurat_reduce_dimensions(integrated_seu)
	integrated_seu <- seurat_cluster(integrated_seu, resolution = seq(0.2, 2.0, by = 0.2), seurat_assay = "integrated")
	saveRDS(integrated_seu, seu_path)
	return(seu_path)
}
