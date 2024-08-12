library(targets)
source("packages.R")
source("functions.R")


tar_load(c("debranched_seus_2p", "debranched_seus_6p"))

myseu_paths <- 
	c(debranched_seus_2p, debranched_seus_6p) |> 
	sort() |> 
	identity()

myseu_paths <- set_names(myseu_paths, fs::path_file(myseu_paths))

cluster_orders <- pull_cluster_orders2("data/scna_cluster_order2.csv")

# seu <- readRDS(myseu_paths[10])
# 
# # debug(assign_phase_clusters)
# 
# assign_phase_clusters(seu, "SRR27187899_filtered_seu.rds", cluster_orders, resolution = "0", group.by = "SCT_snn_res.0.6")

seus <- 
	myseu_paths |> 
	map(readRDS) |> 
	imap(assign_phase_clusters, cluster_orders, resolution = "0", group.by = "SCT_snn_res.0.6") |>
	identity()


map2(seus, myseu_paths, saveRDS)


# debug(make_clone_distribution_figure)

# debug(plot_distribution_of_clones_across_clusters)

map(debranched_seus_2p, make_clone_distribution_figure, cluster_order = cluster_orders, scna_of_interest = "2p", width =14, height = 10) |> 
	map(browseURL)

map(debranched_seus_6p, make_clone_distribution_figure, cluster_order = cluster_orders, scna_of_interest = "6p", width =14, height = 10) |> 
	map(browseURL)

