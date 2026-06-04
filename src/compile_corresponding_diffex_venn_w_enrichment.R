library(targets)
library(ggVennDiagram)
library(ggplot2)
library(tidyverse)
source("packages.R")
source("functions.R")



tar_load(c("corresponding_seus", "corresponding_clusters_diffex"))

corresponding_clusters_diffex

set.seed(20231214)

diffex_genes_per_sample <- function(sample_diffex_list){
	# browser()
	diffex_out <- map(sample_diffex_list, "all") |> 
		map(dplyr::filter, p_val_adj < 0.05) |> 
		map(~dplyr::filter(.x, abs(avg_log2FC) > 0.5)) |> 
		dplyr::bind_rows(.id = "cluster_comparison") |> 
		dplyr::group_by(symbol) |> 
		dplyr::filter(all_same_sign(avg_log2FC)) |> 
		dplyr::distinct(symbol, .keep_all = TRUE) |> 
		identity()
}

diffex_by_sample <- map(corresponding_clusters_diffex, diffex_genes_per_sample) |> 
	set_names(fs::path_file(unlist(corresponding_seus)))

possible_plot_venn_w_genes <- possibly(plot_venn_w_genes)

# 2p g1 ------------------------------

g1_2p <- plot_venn_w_genes(diffex_by_sample[1:4], "2p g1", "g1_")

# 2p s ------------------------------

s_2p <- plot_venn_w_genes(diffex_by_sample[1:4], "2p s", "s_")

# # 2p all ------------------------------
# 
# all_2p <- plot_venn_w_genes(diffex_by_sample[1:4], "2p all")

# 2p g1 ------------------------------

g1_6p <- plot_venn_w_genes(diffex_by_sample[5:6], "6p g1", "g1_")

# 2p s ------------------------------

s_6p <- possible_plot_venn_w_genes(diffex_by_sample[5:6], "6p s", "s_")

# # 2p all ------------------------------
# 
# all_6p <- plot_venn_w_genes(diffex_by_sample[5:6], "6p all")

results <- 
list(
	"g1_2p"  = g1_2p, 
	"s_2p"  = s_2p, 
	# "all_2p"  = all_2p, 
	"g1_6p"  = g1_6p
	# "s_6p"  = s_6p
	# "all_6p" = all_6p
	) 

pdf("results/venn_diagrams_2p_6p.pdf")
print(map(results, "plot"))
dev.off()

browseURL("results/venn_diagrams_2p_6p.pdf")

tables <- 
	results |> 
	purrr::compact() |> 
	map("genes") |>
	# map(tail, n = 1) |>
	# map(~dplyr::select(.x, all_of(c("name", "item")))) |>
	# map(tidyr::unnest, item) |>
	# dplyr::bind_rows(.id = "comparison") |> 
	# dplyr::rename(symbol = item) |>
	# dplyr::group_by(symbol) |> 
	# tidyr::separate(symbol, c("symbol", "sign"), "_") |> 
	identity()

tables_2p_g1 <-
	tables$g1_2p |> 
		dplyr::group_by(id) |>
		dplyr::mutate(num_comparisons = length(str_split_1(id, "\\/"))) |> 
		dplyr::filter(num_comparisons > 2) |> 
		dplyr::select(all_of(c("name", "item"))) |>
		tidyr::unnest(item) |>
		# dplyr::bind_rows(.id = "comparison") |> 
		dplyr::rename(symbol = item) |>
		dplyr::group_by(symbol) |>
		tidyr::separate(symbol, c("symbol", "sign"), "_") |>
	dplyr::arrange(sign, symbol) |>
	identity()

down_csv(tables_2p_g1)

tables_6p_g1 <-
	tables$g1_6p |> 
	dplyr::group_by(id) |>
	dplyr::mutate(num_comparisons = length(str_split_1(id, "\\/"))) |> 
	dplyr::filter(num_comparisons > 1) |> 
	dplyr::select(all_of(c("name", "item"))) |>
	tidyr::unnest(item) |>
	# dplyr::bind_rows(.id = "comparison") |> 
	dplyr::rename(symbol = item) |>
	dplyr::group_by(symbol) |>
	tidyr::separate(symbol, c("symbol", "sign"), "_") |>
	dplyr::arrange(sign, symbol) |>
	identity()

down_csv(tables_6p_g1)

# ora_plot_and_tables <- function(mycomparison, tables) {
# 	# browser() 
#   compiled_enrichment <- 
#   	tables |> 
#   	dplyr::filter(comparison == mycomparison) |> 
#   	dplyr::left_join(annotables::grch38, by = "symbol") |> 
#   	dplyr::distinct(entrez) |>
#   	dplyr::pull(entrez) |> 
#   	enrichGO(keyType = "ENTREZID", 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
#   
#   compiled_enrichment |> 
#   	# dropGO(level = 1) |> 
#   	clusterProfiler::simplify() |> 
#   	clusterProfiler::dotplot() |>
#   	identity()
# }
# 
# comparisons <- unique(all_tables$comparison)
# 
# names(comparisons) <- comparisons
# 
# possible_ora_plot_and_tables <- possibly(ora_plot_and_tables)
# 
# # test5 <- possible_ora_plot_and_tables("g1_2p", all_tables)
# 
# ora_plots_and_tables <- map(c(g1_6p = "g1_6p", g1_2p = "g1_2p"), possible_ora_plot_and_tables, all_tables)
# 
# pdf("results/venn_enrichment_2p_6p.pdf")
# ora_plots_and_tables |> 
# 	imap(~{.x + labs(title = .y)})
# dev.off()
# 
# browseURL("results/venn_enrichment_2p_6p.pdf")

save_plot_set <- function(path, plot_list){
	pdf(path, h = 4)
	print(plot_list)
	dev.off()
	return(path)
}

# heatmaps ------------------------------

mygenes <- split(tables_2p_g1, tables_2p_g1$sign) |> 
	map(pull, symbol)

seus_2p <- corresponding_seus[1:4] |> map(readRDS) |> 
	map(annotate_cluster_scna_percentage)

names(seus_2p) <- fs::path_file(unlist(corresponding_seus[1:4]))

pdf("results/g1_2p_venn_feature_heatmaps.pdf")
for(seu_name in names(seus_2p)){
	for(direction in names(mygenes)){
		
		scna_genes <- select_genes_to_plot(seus_2p[[seu_name]], mygenes, direction = direction, in_scna = TRUE)
		
		scna_heatmap <- seu_gene_heatmap(seus_2p[[seu_name]], features = scna_genes, group.by = c("clusters", "scna", "nCount_gene"), col_arrangement = c("clusters", "scna"), column_split = "clusters") + 
			labs(title = seu_name, subtitle = direction, cluster_rows = TRUE)
		print(scna_heatmap)
		
		non_scna_genes <- select_genes_to_plot(seus_2p[[seu_name]], mygenes, direction = direction, in_scna = FALSE)
		
		non_scna_heatmap <- seu_gene_heatmap(seus_2p[[seu_name]], features = non_scna_genes, group.by = c("clusters", "scna", "nCount_gene"), col_arrangement = c("clusters", "scna"), column_split = "clusters") + 
			labs(title = seu_name, subtitle = direction, cluster_rows = TRUE)
		print(scna_heatmap)
		
		
	}
}
dev.off()

browseURL("results/g1_2p_venn_feature_heatmaps.pdf")


# violins ------------------------------

plot_venn_violins <- function(seu, mygenes, sample_id = "asdf", direction = "-1") {
	# browser()
	
		gene_data <- FetchData(seu, mygenes) |> 
			tibble::rownames_to_column("cell")
		
		# seu <- annotate_cluster_scna_percentage(seu)
		
		seu_meta <- seu@meta.data |> 
			tibble::rownames_to_column("cell") |> 
			dplyr::left_join(gene_data, by = "cell") |> 
			tibble::column_to_rownames("cell") |> 
			identity()
		
		mygenes <- str_subset(mygenes, "^MT", negate = TRUE)
  
		myvlns <- map(mygenes, ~{
			ggpubr::ggviolin(seu_meta, x = "clusters", y = .x, fill = "percent_scna", add = "boxplot", add.params = list(fill = "white")) +
				gradient_fill(c("blue", "white", "red")) +
											 	labs(title=.x) + 
				theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
											 })
		
		mypatch <- wrap_plots(myvlns) + 
			plot_annotation(title = sample_id, subtitle = direction) +
			plot_layout(guides = "collect", ncol = 5) & theme(legend.position = 'bottom')
  
  return(mypatch)
}

possible_plot_venn_violins <- possibly(plot_venn_violins)


pdf("results/g1_2p_venn_feature_violin.pdf", h =8, w = 12)


for(seu_name in names(seus_2p)){
	for(direction in names(mygenes)){
		
		plotted_genes <- select_genes_to_plot(seus_2p[[seu_name]], mygenes, direction = direction, max_genes = 100, in_scna = TRUE)
		
		myvln <- plot_venn_violins(seus_2p[[seu_name]], plotted_genes, seu_name, direction) + 
			# labs(title = seu_name, subtitle = direction) + 
			NULL 
		
		print(myvln)
	}
}

dev.off()

browseURL("results/g1_2p_venn_feature_violin.pdf")


# 6p ------------------------------

# heatmaps ------------------------------

mygenes <- split(tables_6p_g1, tables_6p_g1$sign) |> 
	map(pull, symbol)

seus_6p <- corresponding_seus[2:4] |> map(readRDS) |> 
	map(annotate_cluster_scna_percentage)

names(seus_6p) <- fs::path_file(unlist(corresponding_seus[2:4]))

pdf("results/g1_6p_venn_feature_heatmaps.pdf")
for(seu_name in names(seus_6p)){
	for(direction in names(mygenes)){
		
		scna_genes <- select_genes_to_plot(seus_6p[[seu_name]], mygenes, direction = direction, in_scna = TRUE, only_scna = TRUE, scna_of_interest = "06p")
		
		scna_heatmap <- seu_gene_heatmap(seus_6p[[seu_name]], features = scna_genes, group.by = c("clusters", "scna", "nCount_gene"), col_arrangement = c("clusters", "scna"), column_split = "clusters") + 
			labs(title = seu_name, subtitle = direction, cluster_rows = TRUE)
		print(scna_heatmap)
		
		non_scna_genes <- select_genes_to_plot(seus_6p[[seu_name]], mygenes, direction = direction, in_scna = FALSE, only_scna = TRUE, scna_of_interest = "06p")
		
		non_scna_heatmap <- seu_gene_heatmap(seus_6p[[seu_name]], features = non_scna_genes, group.by = c("clusters", "scna", "nCount_gene"), col_arrangement = c("clusters", "scna"), column_split = "clusters") + 
			labs(title = seu_name, subtitle = direction, cluster_rows = TRUE)
		print(non_scna_heatmap)
		
		
	}
}
dev.off()

browseURL("results/g1_6p_venn_feature_heatmaps.pdf")


# violins ------------------------------

plot_venn_violins <- function(seu, mygenes, sample_id = "asdf", direction = "-1") {
	# browser()
	
	gene_data <- FetchData(seu, mygenes) |> 
		tibble::rownames_to_column("cell")
	
	# seu <- annotate_cluster_scna_percentage(seu)
	
	seu_meta <- seu@meta.data |> 
		tibble::rownames_to_column("cell") |> 
		dplyr::left_join(gene_data, by = "cell") |> 
		tibble::column_to_rownames("cell") |> 
		identity()
	
	mygenes <- str_subset(mygenes, "^MT", negate = TRUE)
	
	myvlns <- map(mygenes, ~{
		ggpubr::ggviolin(seu_meta, x = "clusters", y = .x, fill = "percent_scna", add = "boxplot", add.params = list(fill = "white")) +
			gradient_fill(c("blue", "white", "red")) +
			labs(title=.x) + 
			theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
	})
	
	mypatch <- wrap_plots(myvlns) + 
		plot_annotation(title = sample_id, subtitle = direction) +
		plot_layout(guides = "collect", ncol = 5) & theme(legend.position = 'bottom')
	
	return(mypatch)
}

possible_plot_venn_violins <- possibly(plot_venn_violins)


pdf("results/g1_2p_venn_feature_violin.pdf", h =8, w = 12)

for(seu_name in names(seus_2p)){
	for(direction in names(mygenes)){
		
		plotted_genes <- select_genes_to_plot(seus_2p[[seu_name]], mygenes, direction = direction, max_genes = 100, in_scna = TRUE)
		
		myvln <- plot_venn_violins(seus_2p[[seu_name]], plotted_genes, seu_name, direction) + 
			# labs(title = seu_name, subtitle = direction) + 
			NULL 
		
		print(myvln)
	}
}

dev.off()

browseURL("results/g1_2p_venn_feature_violin.pdf")




debug(seu_gene_heatmap)




seu$scna_status

seu <- ScaleData(seu)


my_heatmaps <- seus_2p |> 
	set_names(fs::path_file(unlist(corresponding_seus[1:4]))) |> 
	imap(~{
	ggplotify::as.ggplot(seu_gene_heatmap(.x, features = mygenes[[2]], group.by = c("clusters", "scna"), col_arrangement = c("clusters", "scna"))) + labs(title = .y) 
}
)

# section ------------------------------

plot_venn_violin_multi_seu <- function(seu, gene) {
	
	myvlns <- split(mygenes, mygenes$sign) |> 
		map(pull, symbol) |> 
		map(~VlnPlot(seu, features = .x, group.by = "clusters", combine = FALSE, split.by = "scna_status"))
	map_depth(2, ~{.x + facet_wrap(~sample_id)})
	# imap(~save_plot_set(glue("results/{mycomparison}_{.y}.pdf"), .x))
	
	return(myvlns)
}