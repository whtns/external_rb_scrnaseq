library(Seurat)
library(tidyverse)

source("packages.R")
source("functions.R")

seu <- readRDS("output/seurat/SRR14800534_filtered_seu.rds")

# seu <- seu[,!seu$clusters == "other_10"]

seu$scna <- factor(seu$scna, levels = c("", "16q-", "16q- 1q+"))
seu$scna_status <- factor(ifelse(str_detect(seu$scna, "1q"), "w/ 1q+", "w/o 1q+"), levels = c("w/o 1q+", "w/ 1q+"))

seu$clusters <-
	factor(seu$clusters)

labels <- data.frame(clusters = unique(seu[[]][["clusters"]]), label = unique(seu[[]][["clusters"]])) %>%
	# dplyr::rename({{group.by}} := cluster) %>%
	identity()

cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))

centroid_data <-
	cc_data %>%
	dplyr::group_by(clusters) %>%
	dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
	dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
	dplyr::mutate(centroid = "centroids") %>%
	identity()

dot_colors <- c("clusters", "Phase", "scna")

umap_plots <- list()
centroid_plots <- list()
ccplots <- list()
for(dot_color in dot_colors){
	ccplots[[dot_color]] <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[[dot_color]])) +
		geom_point(size = 0.1) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
		facet_wrap(~ .data[["clusters"]], nrow = 2) +
		theme_light() +
		geom_label(
			data = labels,
			aes(label = label),
			# x = Inf,
			# y = -Inf,
			x = max(cc_data$S.Score) + 0.05,
			y = max(cc_data$G2M.Score) - 0.1,
			hjust = 1,
			vjust = 1,
			inherit.aes = FALSE
		) +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank()
		) +
		# guides(color = "none") +
		NULL
	
	centroid_plots[[dot_color]] <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[[dot_color]])) +
		geom_point(size = 0.1) +
		theme_light() +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank()
		) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23, colour = "black") +
		guides(fill = "none", color = "none") +
		NULL
	
	umap_plots[[dot_color]] <- 
		DimPlot(seu, group.by = {{dot_color}})
	
}

clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
	seu,
	seu_name = glue("asdf"), var_x = "scna", var_y = "clusters", signif = FALSE, plot_type = "clone"
)
ggsave("results/fig_01f.pdf", w = 3.5, h = 3)
browseURL("results/fig_01f.pdf")

pdf("results/fig_01d.pdf", h = 3, w = 4)
print(centroid_plots)
dev.off()
browseURL("results/fig_01d.pdf")

pdf("results/fig_01e.pdf", h = 4, w = 12)
print(ccplots)
dev.off()

qpdf::pdf_combine(
	list(
		"results/fig_01d.pdf",
		"results/fig_01e.pdf",
		"results/fig_01f.pdf"
	), 
	"results/fig_01.pdf"
)

browseURL("results/fig_01.pdf")

browseURL("results/fig_01.pdf")
