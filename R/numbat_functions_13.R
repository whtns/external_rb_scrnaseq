# Numbat Functions (13)

#' Perform make table s02 operation
#'
#' @param meta_path File path
#' @return Function result
#' @export
# Performance optimizations applied:
# - long_pipe_chain: Consider breaking into intermediate variables for readability and debugging

make_table_s02 <- function(meta_path = "data/metadata.tsv") {
	
	project_ids <- c(
		"liu" = "PRJNA1051579",
		"collin" = "PRJNA699542", 
		"wu" = "PRJNA707413", 
		"yang" = "PRJNA737188", 
		"field" = "PRJNA804875"
	)
	
	retained_columns <- c("Run",  
												"sample_id", "Age", "Assay Type", "AvgSpotLen", "Bases", 
												"BIOMATERIAL_PROVIDER", "BioSampleModel", "Bytes", "Center Name", 
												"DATASTORE region", "Instrument", "Isolate", "Library Name", 
												"Organism", "Replicate", "Sample Name", "sex", "SRA Study", "Tissue", 
												"Treatment", "age_mo", "Developmental_stage", "GEO_Accession (exp)", 
												"source_name", "tissue", "familial", "iirc_tumor_class", "Experiment", "LibrarySelection", "LibrarySource", 
												"LibraryLayout", "Platform", "BioProject", "BioSample", "num_clone")
	
	field_bulk_samples <- glue("SRR17960{seq(497,512,1)}")
	collin_bad_samples <- c("SRR13633761", "SRR13633762")
	
	collated_metadata <-
		meta_path |> 
		read_tsv() |> 
		dplyr::filter(BioProject %in% project_ids) |> 
		dplyr::filter(is.na(`Assay Type`) | `Assay Type` == "RNA-Seq") |>
		dplyr::filter(is.na(LibrarySource) | LibrarySource != "Genomic") |>
		dplyr::filter(is.na(Tissue) | Tissue != "Low-passaged primary retinoblastoma cell line") |> 
		dplyr::filter(!Run %in% field_bulk_samples) |> 
		dplyr::filter(!Run %in% collin_bad_samples) |> 
		dplyr::left_join(tibble::enframe(project_ids, "author", "BioProject"), by = "BioProject") |> 
		dplyr::select(author, any_of(retained_columns)) |> 
		dplyr::arrange(author, Run) |> 
		identity()
	
	table_path <- "results/table_s02.csv"
	write_csv(collated_metadata, table_path)
	
	return(table_path)
}

#' Perform clustering analysis
#'
#' @param seu_path File path
#' @param cluster_dictionary Cluster information
#' @return Data frame
#' @export
check_cluster_marker_gene <- function(seu_path, cluster_dictionary){
	
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	seu <- readRDS(seu_path)
	
	df0 <- cluster_dictionary[[tumor_id]] |> 
		dplyr::filter(remove == "1")
	
	if(nrow(df0) > 0){
		cluster_sizes <-
			seu@meta.data |>
			dplyr::group_by(`gene_snn_res.0.2`) |> 
			dplyr::summarize(n_cells = dplyr::n()) |> 
			dplyr::mutate(Cluster = as.numeric(gene_snn_res.0.2))
		
		df1 <- seu@misc$markers$gene_snn_res.0.2$presto |> 
			dplyr::mutate(Cluster = as.numeric(Cluster)) |> 
			dplyr::inner_join(df0, by = c("Cluster" = "gene_snn_res.0.2")) |> 
			dplyr::group_by(Cluster) |> 
			dplyr::slice_head(n = 5) |> 
			dplyr::inner_join(cluster_sizes, by = c("Cluster"))
		
	}
	
	return(df1)
	
}

#' Perform make table s03 operation
#'
#' @param cluster_dictionary Cluster information
#' @param table_path File path
#' @return Function result
#' @export
make_table_s03 <- function(cluster_dictionary, table_path = "results/table_s03.csv"){
	myseunames <-c("SRR13884242", "SRR13884243", "SRR13884246", "SRR13884247",
		"SRR13884248", "SRR13884249", "SRR14800534", "SRR14800535", "SRR14800536",
		"SRR14800540", "SRR14800541", "SRR14800543", "SRR17960481", "SRR17960484",
		"SRR27187899", "SRR27187902")
	
	seu_paths <- fs::path("output/seurat/", glue("{myseunames}_seu.rds")) |> 
		set_names(myseunames)
	
	possibly_check_cluster_marker_gene <- possibly(check_cluster_marker_gene)
	
	df2 <- map(seu_paths, possibly_check_cluster_marker_gene, cluster_dictionary) |> 
		purrr::compact() |> 
		dplyr::bind_rows() |> 
		dplyr::group_by(sample_id, abbreviation, gene_snn_res.0.2, n_cells) |> 
		summarise(marker_genes = paste(Gene.Name, collapse = ", ")) |> 
		identity()
		
		test0 <- 
			map_int(seu_paths, ~ncol(readRDS(.x))) |>
			tibble::enframe("sample_id", "total_cells") |> 
			dplyr::full_join(df2, by = "sample_id") |> 
			dplyr::group_by(sample_id) |> 
			dplyr::mutate(percent_filtered = sum(n_cells)/total_cells) |> 
			dplyr::mutate(percent_filtered = replace_na(percent_filtered, 0)) |>
			identity()
		
		# test0 |> 
		# 	dplyr::ungroup() |> 
		# 	dplyr::distinct(sample_id, .keep_all = TRUE) |> 
		# 	dplyr::summarise(mean(percent_filtered)) |>
		# 	identity()
	
	write_csv(test0, table_path)
	
	return(table_path)
	
}

#' Perform retrieve genes in cis operation
#'
#' @param nb_path File path
#' @param tumor_id Character string (default: "asdf")
#' @return Function result
#' @export
retrieve_genes_in_cis <- function(nb_path, tumor_id = "asdf"){
	
	mynb <- readRDS(nb_path)
	
	segments <- 
		mynb$segs_consensus |> 
		dplyr::filter(!is.na(sample)) |> 
		dplyr::select(seqnames = CHROM, start = seg_start, end = seg_end, seg, cnv_state) |> 
		plyranges::as_granges() %>%
		identity()
	
	genes_in_cis <- 
		annotables::grch38 %>%
		dplyr::distinct(ensgene, .keep_all = TRUE) %>%
		dplyr::mutate(seqnames = chr) %>%
		dplyr::filter(!is.na(start), !is.na(end)) %>%
		plyranges::as_granges() %>%
		plyranges::join_overlap_intersect(segments) %>%
		as_tibble() |> 
		dplyr::mutate(sample_id = tumor_id) |> 
		dplyr::arrange(seg, start, end) |> 
		dplyr::select(sample_id, seg, everything()) |> 
		dplyr::filter(symbol != "", !is.na(entrez)) |> 
		dplyr::distinct(symbol, .keep_all = TRUE)
	
	return(genes_in_cis)
	
}

#' Perform make table s08 operation
#'
#' @param table_path File path
#' @return Function result
#' @export
make_table_s08 <- function(table_path = "results/table_s04.csv"){
	numbat_names <- c("SRR13884242", "SRR13884243", "SRR13884246", "SRR13884247", 
		"SRR13884248", "SRR13884249", "SRR14800534", "SRR14800535", "SRR14800536", 
		"SRR14800540", "SRR14800541", "SRR14800543", "SRR17960481", "SRR17960484", 
		"SRR27187899", "SRR27187902")
	
	nb_paths <- fs::path("output/numbat_sridhar/", glue("{numbat_names}_numbat.rds")) |> 
		set_names(numbat_names)
	
	test0 <- imap(nb_paths, retrieve_genes_in_cis) |> 
		dplyr::bind_rows() |> 
		write_csv(table_path)
	
	return(table_path)
}

