library(ggpubr)
library(ggplot2)
library(tidyverse)

meta <- read_csv("results/tmp.csv") %>%
  dplyr::rename(cell = `...1`) %>%
  dplyr::filter(clone_opt %in% c("1", "2", "3")) %>%
  dplyr::group_by(leiden, scna)  %>% 
  dplyr::mutate(leiden = factor(leiden))  %>% 
  identity()


# boxplots -------------------------------

ggplot(meta, aes(x = reorder(scna, subtype2, FUN = mean, na.rm = TRUE), y = subtype2, fill = scna)) +
  geom_violin() +
  # geom_jitter() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "subtype2 score by scna") +
  NULL

ggsave("results/subtype2_scores_anndata_aggregate.pdf", height = 8)
browseURL("results/subtype2_scores_anndata_aggregate.pdf")


ggplot(meta, aes(x = reorder(scna, subtype1, FUN = mean, na.rm = TRUE), y = subtype1, fill = scna)) +
  geom_violin() +
  # geom_jitter() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "subtype1 score by scna") +
  NULL

ggsave("results/subtype1_scores_anndata_aggregate.pdf", height = 8)
browseURL("results/subtype1_scores_anndata_aggregate.pdf")

# boxplots facet scna -------------------------------

ggplot(meta, aes(x = reorder(scna, subtype2, FUN = mean, na.rm = TRUE), y = subtype2, fill = scna)) +
  geom_violin() +
  # geom_jitter() + 
  facet_wrap(~sample_id, ncol = 1, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "subtype2 score by scna") +
  NULL

ggsave("results/subtype2_scores_anndata_facet.pdf", height = 36)
browseURL("results/subtype2_scores_anndata_facet.pdf")


ggplot(meta, aes(x = reorder(scna, subtype1, FUN = mean, na.rm = TRUE), y = subtype1, fill = scna)) +
  geom_violin() +
  # geom_jitter() + 
  facet_wrap(~sample_id, ncol = 1, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "subtype1 score by scna") +
  NULL

ggsave("results/subtype1_scores_anndata_facet.pdf", height = 36)
browseURL("results/subtype1_scores_anndata_facet.pdf")


# boxplots facet scna -------------------------------

ggplot(meta, aes(x = reorder(leiden, subtype2, FUN = mean, na.rm = TRUE), y = subtype2, fill = leiden)) +
  geom_violin() +
  # geom_jitter() + 
  facet_wrap(~sample_id, ncol = 1, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "subtype2 score by scna") +
  NULL

ggsave("results/subtype2_scores_anndata_facet.pdf", height = 36)
browseURL("results/subtype2_scores_anndata_facet.pdf")


ggplot(meta, aes(x = reorder(leiden, subtype1, FUN = mean, na.rm = TRUE), y = subtype1, fill = leiden)) +
  geom_violin() +
  # geom_jitter() + 
  facet_wrap(~sample_id, ncol = 1, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(title = "subtype1 score by scna") +
  NULL

ggsave("results/subtype1_scores_anndata_facet.pdf", height = 36)
browseURL("results/subtype1_scores_anndata_facet.pdf")

