# 2p ------------------------------

tar_load("enrichment_2p_g1")

names(enrichment_2p_g1) <- names(debranched_seus_2p)

test0 <-
	enrichment_2p_g1[c(1:2, 6)] |>
	# g1_2p_effects |> 
	purrr::list_flatten()

merged_test0 <-
	test0 |> 
	clusterProfiler::merge_result()

merged_test0@compareClusterResult <-
	merged_test0@compareClusterResult |> 
	dplyr::group_by(Description) |> 
	dplyr::mutate(n_recur = dplyr::n()) |> 
	dplyr::mutate(match_sign = all_same_sign(NES)) |> 
	dplyr::filter(match_sign) |>
	dplyr::filter(n_recur > 1) |>
	identity()

plot_enrichment(merged_test0, p_val_cutoff = 1, result_slot = "compareClusterResult") + 
	scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "_", " "), width = 10)) +
	labs(title = "expression of 2p+ enriched G1 states")
scna_of_interest = "2p"
phase = "g1"
plot_path <- ggsave(glue("results/{scna_of_interest}_{phase}_alternative_states_enrichment.pdf"), w = 10, h = 10)

browseURL(plot_path)

# 6p ------------------------------

tar_load(c("debranched_seus_6p", "enrichment_6p_g1"))

names(enrichment_6p_g1) <- names(debranched_seus_6p)

test0 <-
	enrichment_6p_g1[c(1,3)] |>
	purrr::list_flatten() |> 
	clusterProfiler::merge_result()

test0@compareClusterResult <-
	test0@compareClusterResult |> 
	dplyr::group_by(Description) |> 
	dplyr::mutate(n_recur = dplyr::n()) |> 
	dplyr::mutate(match_sign = all_same_sign(NES)) |> 
	dplyr::filter(match_sign) |>
	dplyr::filter(n_recur > 1) |>
	identity()

plot_enrichment(test0, p_val_cutoff = 1, result_slot = "compareClusterResult") + 
	scale_x_discrete(labels = function(x) str_wrap(str_replace_all(x, "_", " "), width = 10)) +
	labs(title = "expression of 6p+ enriched G1 states")
scna_of_interest = "6p"
phase = "g1"
plot_path <- ggsave(glue("results/{scna_of_interest}_{phase}_alternative_states_enrichment.pdf"), w = 10, h = 14)

browseURL(plot_path)
