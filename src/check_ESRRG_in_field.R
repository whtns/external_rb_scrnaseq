seu <- readRDS("output/seurat/integrated_1q/integrated_seu_1q_complete.rds")


library(seuratTools)
library(fs)
library(tidyverse)
library(glue)

VlnPlot(seu, features="ESRRG", group.by = "scna")


field_seus <-
	dir_ls("output/seurat/", regexp = ".*SRR179604[0-9]*_seu.rds") |> 
	set_names(str_extract, "SRR[0-9]*") |> 
	map(readRDS) |>
	identity()

filtered_field_seus <-
	dir_ls("output/seurat/", regexp = ".*SRR179604[0-9]*_filtered_seu.rds") |> 
	set_names(str_extract, "SRR[0-9]*") |> 
	map(readRDS) |>
	identity()

plot_path = "results/ESSRG_correlation_w_1q.pdf"
pdf(plot_path)

map(field_seus, FeaturePlot, features = "ESRRG") |> 
	imap(~{.x + labs(subtitle = .y)})

map(filtered_field_seus, FeaturePlot, features = "ESRRG", split.by = "scna") |> 
	imap(~{.x + labs(subtitle = glue("{.y} filtered"))})

map(filtered_field_seus, plot_violin, features = "ESRRG", plot_var = "scna") |> 
	imap(~{.x + labs(subtitle = glue("{.y} filtered"))})

dev.off()
browseURL(plot_path)


# my_comparisons <- list( c("0.5", "1"), c("1", "2"), c("0.5", "2") )
# pbox + stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
# 	stat_compare_means(label.y = 50)     