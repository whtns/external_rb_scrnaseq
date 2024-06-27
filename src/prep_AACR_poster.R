source("packages.R")
source("functions.R")
library(targets)


# cell stats plot ------------------------------

tar_load("study_cell_stats")
# 
# debug(plot_study_metadata)
debug(plot_study_cell_stats)

test0 <- plot_study_metadata(study_cell_stats, plot_height = 6, width = 12)




#  cell type filtering plot ------------------------------

tar_load("unfiltered_seus")

# 4249

# unfiltered 

unfiltered_seu <- 
	unfiltered_seus[[6]] |> 
	readRDS() |>
	identity()

unfiltered_seu_wo_rods <- unfiltered_seu[,!unfiltered_seu$type == "Rods"]

DimPlot(unfiltered_seu, group.by = "gene_snn_res.0.2")

marker_plot <-
	poster_plot_markers(unfiltered_seu, metavar = "type", num_markers = 6) +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank()
	) + 
	# coord_flip() + 
	NULL

mypalette <- scales::hue_pal()(3) |> 
	set_names(unique(unfiltered_seu$type))

unfiltered_plot <-
	DimPlot(unfiltered_seu, group.by = "type") +
	labs(title = "unfiltered") + 
	scale_color_discrete(labels = label_wrap_gen(10)) +
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		legend.position = "top"
	) +
	scale_color_manual(values=mypalette)

filtered_plot <-
	DimPlot(unfiltered_seu_wo_rods, group.by = "type") +
	labs(title = "filtered") + 
	scale_color_discrete(labels = label_wrap_gen(10)) + 
	theme(
		axis.title.x = element_blank(),
		axis.title.y = element_blank(),
		legend.position = "top"
	) + 
	scale_color_manual(values=mypalette)

design <- 
	"AACC
  AACC
  BBCC
  BBCC"

design <- 
	"AABB
  AABB
  CCCC
  CCCC"

mypatch <- wrap_plots(unfiltered_plot, filtered_plot, marker_plot) + 
	plot_layout(design = design)


patch_path <- ggsave("results/cell_type_filter_poster.pdf", mypatch, width = 7, height = 10)

browseURL(patch_path)


# trinity plots (umap and cc) ------------------------------


tar_load("filtered_seus")

filtered_seu <- 
	filtered_seus[[7]] |> 
	readRDS() |>
	identity()

DimPlot(filtered_seu, group.by = "clusters") + 
	DimPlot(filtered_seu, group.by = "Phase") + 
	DimPlot(filtered_seu, group.by = "scna")


patch_path <- ggsave("results/trinity_umap_plots.pdf", width = 12, height = 3)

browseURL(patch_path)


plot_ks_cell_cycle <- function(seu, color_by = "scna", whether_facet = TRUE) {
	
	# browser()
	
	labels <- data.frame(clusters=unique(seu[[]][["clusters"]]), label =unique(seu[[]][["clusters"]])) %>%
		# dplyr::rename({{group.by}} := cluster) %>%
		identity()
	
  cc_data <- FetchData(seu, unique(c("clusters", "G2M.Score", "S.Score", "Phase", "scna", color_by)))
  
  centroid_data <-
  	cc_data %>%
  	dplyr::group_by(clusters) %>%
  	dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
  	dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
  	dplyr::mutate(centroid = "centroids") %>%
  	identity()
  
  centroid_plot <-
  	cc_data %>%
  	ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[[color_by]])) +
  	geom_point(size = 0.1) +
  	theme_light() +
  	theme(
  		strip.background = element_blank(),
  		strip.text.x = element_blank()
  	) +
  	geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
  	NULL
  
  
  cell_cycle_plot <-
  	cc_data %>%
  	ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[[color_by]])) +
  	geom_point(size = 0.1) +
  	geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
  	theme_light() +
  	geom_label(data = labels,
  						 aes(label = label),
  						 # x = Inf,
  						 # y = -Inf,
  						 x = max(cc_data$S.Score)+0.05,
  						 y = max(cc_data$G2M.Score)-0.1,
  						 hjust=1,
  						 vjust=1,
  						 inherit.aes = FALSE) +
  	theme(
  		strip.background = element_blank(),
  		strip.text.x = element_blank(),
  		axis.title.y = element_blank(),
  		plot.title = element_text(size = 16)
  	) +
  	guides(
  				 fill = "none",
  				 color = guide_legend(override.aes = list(size = 4))) +
  	NULL
  
  facet_cell_cycle_plot <- 
  	cell_cycle_plot + 
  	facet_wrap(~.data[["clusters"]], ncol = 2)
  
  if(whether_facet){
  	return(facet_cell_cycle_plot)
  } else {
  	return(cell_cycle_plot)	
  }
}

## faceted 
plot_ks_cell_cycle(filtered_seu, color_by = "clusters") + 
	guides(color = "none") + 
	labs(title = "Cluster") +
	plot_ks_cell_cycle(filtered_seu, color_by = "Phase") + 
	theme(
		axis.title.y = element_blank()
	) +
	labs(title = "Phase") +
	plot_ks_cell_cycle(filtered_seu, color_by = "scna") +
	theme(
		axis.title.y = element_blank()
	) +
	labs(title = "Clone")


patch_path <- ggsave("results/trinity_cc_plots.pdf", width = 12, height = 7)

browseURL(patch_path)

## not faceted
plot_ks_cell_cycle(filtered_seu, color_by = "clusters", whether_facet = FALSE) + 
	guides(color = "none") + 
	labs(title = "Cluster") +
	plot_ks_cell_cycle(filtered_seu, color_by = "Phase", whether_facet = FALSE) + 
	theme(
		axis.title.y = element_blank()
	) +
	labs(title = "Phase") +
	plot_ks_cell_cycle(filtered_seu, color_by = "scna", whether_facet = FALSE) +
	theme(
		axis.title.y = element_blank()
	) +
	labs(title = "Clone")


patch_path <- ggsave("results/trinity_cc_all_plots.pdf", width = 12, height = 3)

browseURL(patch_path)


# clustrees ------------------------------

tar_load(c("clustrees", "debranched_samples"))

names(clustrees) <- debranched_samples

# clone pearls ------------------------------
tar_load("clone_pearls")

clone_pearls |> unlist() |> qpdf::pdf_combine("results/clone_pearls.pdf") |> browseURL()

tar_load("debranched_clone_pearls")

debranched_clone_pearls |> unlist() |> qpdf::pdf_combine("results/debranched_clone_pearls.pdf") |> browseURL()

# trans enrichment ------------------------------


# compile_cis_trans_enrichment_recurrence ------------------------------

tar_load("oncoprint_enrich_clones_hallmark")

# debug(compile_cis_trans_enrichment_recurrence)
debug(plot_enrichment_recurrence)

test1 <- 
	compile_cis_trans_enrichment_recurrence(oncoprint_enrich_clones_hallmark,
																					cis_plot_file = "results/enrichment/cis_enrichment_plots_by_clone_hallmark.pdf",
																					trans_plot_file = "results/enrichment/trans_enrichment_plots_by_clone_hallmark.pdf",
																					cis_table_file = "results/enrichment/cis_enrichment_tables_by_clone_hallmark.xlsx",
																					trans_table_file = "results/enrichment/trans_enrichment_tables_by_clone_hallmark.xlsx"
	)

