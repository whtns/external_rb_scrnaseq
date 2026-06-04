copy_paper_data <- function(clustree_compilation, collage_compilation, table_all_diffex_clones, oncoprint_enrich_clones_plots_gobp, oncoprint_enrich_clusters_plots_gobp, oncoprint_plots){
	fs::dir_create("results/.figures_and_tables")
	
	# figures ------------------------------
	
	my_figures <- list(
		# fig 01
		fs::file_copy("results/unfiltered_study_stats.pdf", "results/.figures_and_tables/fig_02.pdf", overwrite = TRUE),
		
		fs::file_copy("results/bad_qc_study_stats.pdf", "results/.figures_and_tables/fig_02b.pdf", overwrite = TRUE),
		fs::file_copy("results/bad_scna_study_stats.pdf", "results/.figures_and_tables/fig_02c.pdf", overwrite = TRUE),
		
		# fig 03
		fs::file_copy("results/subtype_scores_by_clone.pdf", "results/.figures_and_tables/fig_03.pdf", overwrite = TRUE),
		
		# fig 04
		fs::file_copy("results/SRR14800534_SRR14800535_SRR14800536_filtered_heatmap.pdf", "results/.figures_and_tables/fig_04.pdf", overwrite = TRUE),
		
		# fig 05
		qpdf::pdf_combine(c("results/diffex_oncoprints_in_segment.pdf", "results/diffex_oncoprints_out_of_segment.pdf"), output = "results/.figures_and_tables/fig_05.pdf"),
		
		# fig s07
		fs::file_copy(clustree_compilation, "results/.figures_and_tables/fig_s07.pdf", overwrite = TRUE),
		
		# fig s09
		fs::file_copy(collage_compilation, "results/.figures_and_tables/fig_s09.pdf", overwrite = TRUE)
	)
	
	# tables ------------------------------
	my_tables <- list(
		# table S01
		fs::file_copy("results/study_cell_stats.csv", "results/.figures_and_tables/table_s01.csv", overwrite = TRUE),
		
		# table S02
		fs::file_copy(table_all_diffex_clones$total[[2]], "results/.figures_and_tables/table_s02a.xlsx", overwrite = TRUE),
		fs::file_copy(table_all_diffex_clones$cluster[[2]], "results/.figures_and_tables/table_s02b.xlsx", overwrite = TRUE),
		
		# table S03
		fs::file_copy(oncoprint_enrich_clones_plots_gobp[[4]], "results/.figures_and_tables/table_s03a.xlsx", overwrite = TRUE),
		fs::file_copy(oncoprint_enrich_clusters_plots_gobp[[4]], "results/.figures_and_tables/table_s03b.xlsx", overwrite = TRUE),
		
		# table S08
		fs::file_copy(oncoprint_plots$cis$table, "results/.figures_and_tables/table_s08a.xlsx", overwrite = TRUE),
		fs::file_copy(oncoprint_plots$trans$table, "results/.figures_and_tables/table_s08b.xlsx", overwrite = TRUE),
		fs::file_copy(oncoprint_plots$all$table, "results/.figures_and_tables/table_s08c.xlsx", overwrite = TRUE)
	) %>% 
		purrr::set_names(fs::path_file(fs::path_ext_remove(.)))
	
	return(list("figures" = my_figures, "tables" = my_tables))
	
	# fig fe
	heatmaps <- dir_ls("output/numbat/", glob = "*phylo_heatmap.png", recurse = TRUE) %>%
		purrr::set_names(path_file(path_dir(.)))
	dir_create("results/.figures_and_tables/fig_fe/")
	imap(heatmaps, ~file_copy(.x, glue("results/.figures_and_tables/fig_fe/{.y}.png"), overwrite = TRUE))
	
	# fig S01
	exp_plots <- dir_ls("output/numbat/", glob = "*exp_roll_clust.png", recurse = TRUE) %>%
		purrr::set_names(path_file(path_dir(.)))
	dir_create("results/.figures_and_tables/fig_s01")
	imap(exp_plots, ~file_copy(.x, glue("results/.figures_and_tables/fig_s01/{.y}.png"), overwrite = TRUE))
	
	
}