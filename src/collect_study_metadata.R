#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(seuratTools)
library(glue)
library(ggridges)
library(patchwork)

retrieve_cell_stats <- function(seu_path){

  seu <- readRDS(seu_path)

  stats = seu@meta.data[c("nCount_gene", "nFeature_gene", "percent.mt", "gene_snn_res.0.2")] %>%
    tibble::rownames_to_column("cell")

  return(stats)
}

seus <-
  dir_ls("output/seurat/", regexp = ".*[0-9]+_seu.rds") %>%
  set_names(str_extract(., "SRR[0-9]*"))


# collin ------------------------------

collin_cell_stats <- seus[c("SRX10031191", "SRX10031192", "SRX10031193", "SRX10031194")] %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# field ------------------------------

field_cell_stats <- seus[c("SRX14116948", "SRX14116947", "SRX14116946", "SRX14116945", "SRX14116944")] %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# wu ------------------------------

wu_cell_stats <- seus[c("SRX10264517", "SRX10264518", "SRX10264519", "SRX10264520", "SRX10264521", "SRX10264522", "SRX10264523", "SRX10264524", "SRX10264525", "SRX10264526")] %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# yang ------------------------------

yang_cell_stats <- seus[c("SRX11133594", "SRX11133593", "SRX11133592", "SRX11133591", "SRX11133590", "SRX11133589", "SRX11133588", "SRX11133587", "SRX11133586", "SRX11133585")] %>%
  map_dfr(retrieve_cell_stats, .id = "sample_id")

# combined ------------------------------

study_cell_stats <- dplyr::bind_rows(list("collin" = collin_cell_stats, "field" = field_cell_stats, "wu" = wu_cell_stats, "yang" = yang_cell_stats), .id = "study")

write_csv(study_cell_stats, "results/study_cell_stats.")

plot_study_cell_stats <- function(study_cell_stats, cell_stats_plot_file) {
  umis_per_cell <- ggplot(study_cell_stats, aes(x = nCount_gene, y = sample_id, fill = study, group = sample_id)) +
    geom_density_ridges() +
    scale_x_log10() +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values = scales::hue_pal()(4)) +
    theme(
      axis.title.y=element_blank(),  #remove y axis labels,
      # axis.text.y=element_blank(),  #remove y axis labels
      # axis.ticks.y=element_blank()  #remove y axis ticks
    ) +
    labs(title = "UMIs per cell")

  genes_per_cell <- ggplot(study_cell_stats, aes(x = nFeature_gene, y = sample_id, fill = study, group = sample_id)) +
    geom_density_ridges() +
    scale_x_log10() +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values = scales::hue_pal()(4)) +
    theme(
      axis.title.y=element_blank(),  #remove y axis labels,
      axis.text.y=element_blank(),  #remove y axis labels
      axis.ticks.y=element_blank()  #remove y axis ticks
    ) +
    labs(title = "genes per cell")

  percent_mito_per_cell <- ggplot(study_cell_stats, aes(x = `percent.mt`, y = sample_id, fill = study, group = sample_id)) +
    geom_density_ridges() +
    scale_y_discrete(limits=rev) +
    scale_fill_manual(values = scales::hue_pal()(4)) +
    theme(
      axis.title.y=element_blank(),  #remove y axis labels,
      axis.text.y=element_blank(),  #remove y axis labels
      axis.ticks.y=element_blank()  #remove y axis ticks
    ) +
    labs(title = "percent mito per cell") +
    xlim(0, 25)

  study_stats_patchwork <- umis_per_cell + genes_per_cell + percent_mito_per_cell + plot_layout(guides = 'collect')

  ggsave(cell_stats_plot_file, study_stats_patchwork)
}

plot_study_cell_stats(study_cell_stats, "results/study_stats.pdf")

dropped_samples <- c("SRX10031191", "SRX10031192", "SRX10031193", "SRX10031194", "SRX14116945", "SRX11133590")

filtered_study_cell_stats <-
  study_cell_stats %>%
  dplyr::filter(!sample_id %in% dropped_samples)

plot_study_cell_stats(filtered_study_cell_stats, "results/filtered_study_stats.pdf")


# plot cell numbers ------------------------------

all_study_metadata <- read_csv("~/single_cell_projects/resources/external_rb_scrnaseq_proj/results/metadata_by_study.csv")

cells_per_sample <-
  all_study_metadata %>%
  group_by(study) %>%
  dplyr::count(sample_id) %>%
  identity()

write_csv(cells_per_sample, "results/cells_per_sample.csv")




