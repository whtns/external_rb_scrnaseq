#!/usr/bin/env Rscript
# Boxplot of per-cluster mean hypoxia score for each sample at each resolution.
# Each boxplot summarises the cluster-level mean_score distribution.
#
# Points are coloured by the THREE states of the current two-gate rule, because
# being a score outlier is necessary but no longer sufficient to be filtered
# (see doc/hypoxia_splitting_strategy.md):
#   flagged  - is_outlier AND >=2 hypoxia markers AND phase-mixed (dom_frac<0.70)
#              -> these are the clusters actually moved to the high subset.
#   spared   - is_outlier but blocked by the marker and/or phase gate -> stays low.
#              Mostly phase-restricted clusters we deliberately keep.
#   cluster  - not a score outlier.
# Outlier points (flagged or spared) are labelled with cluster id and cell count.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

in_csv  <- "results/hypoxia_cluster_split/hypoxia_split_log_all.csv"
out_pdf <- "results/hypoxia_cluster_split/hypoxia_mean_score_boxplots.pdf"

d <- read_csv(in_csv, show_col_types = FALSE) %>%
  mutate(
    resolution = factor(resolution, levels = sort(unique(resolution))),
    res_lab    = paste0("res ", as.character(resolution)),
    is_outlier = as.logical(is_outlier),
    flagged    = as.logical(flagged),
    status = factor(
      dplyr::case_when(
        flagged             ~ "flagged",
        is_outlier          ~ "spared",
        TRUE                ~ "cluster"
      ),
      levels = c("cluster", "spared", "flagged")
    )
  )
d$res_lab <- factor(d$res_lab, levels = paste0("res ", levels(d$resolution)))

# Label every score outlier (flagged or spared) with cluster id and cell count.
d <- d %>%
  mutate(
    lab = ifelse(is_outlier,
                 paste0("c", cluster, " (n=", n_cells, ")"),
                 NA_character_)
  )

samples   <- sort(unique(d$sample_id))
n_out_tot <- sum(d$is_outlier, na.rm = TRUE)
n_flag_tot <- sum(d$flagged, na.rm = TRUE)
message("plotting ", length(samples), " samples on one page; ",
        n_out_tot, " score outliers, ", n_flag_tot, " flagged")

# Per-sample flagged / outlier counts in the facet strip label.
n_by <- d %>%
  group_by(sample_id) %>%
  summarise(n_out = sum(is_outlier, na.rm = TRUE),
            n_flag = sum(flagged, na.rm = TRUE), .groups = "drop")
d <- d %>%
  left_join(n_by, by = "sample_id") %>%
  mutate(strip = paste0(sample_id, " (", n_flag, "/", n_out, " flagged)"))

# All samples side-by-side, one page, shared (fixed) y axis across facets.
p <- ggplot(d, aes(x = res_lab, y = mean_score)) +
  geom_boxplot(outlier.shape = NA, width = 0.6,
               fill = "grey92", colour = "grey45") +
  geom_jitter(aes(colour = status),
              width = 0.12, height = 0, size = 1.1, alpha = 0.85) +
  geom_text_repel(aes(label = lab, colour = status), size = 2.0,
                  min.segment.length = 0, box.padding = 0.25,
                  max.overlaps = Inf,
                  na.rm = TRUE, seed = 1, show.legend = FALSE) +
  scale_colour_manual(
    values = c(cluster = "grey35", spared = "#e08214", flagged = "firebrick"),
    labels = c(cluster = "not an outlier",
               spared  = "outlier, spared by gates (stays low)",
               flagged = "flagged → high subset"),
    name = NULL, drop = FALSE) +
  facet_wrap(~ strip, ncol = 6, scales = "fixed") +
  labs(
    title    = "Cluster mean hypoxia score by resolution — all samples (shared y axis)",
    subtitle = paste0(
      n_flag_tot, " of ", n_out_tot,
      " score-outlier clusters flagged across ", length(samples),
      " samples. Outlier = median + 3*MAD fence; flagged also requires ",
      ">=2 hypoxia markers and phase-mixed (dom_frac < 0.70)."),
    x = "clustering resolution", y = "cluster mean_score"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold"),
        strip.text = element_text(size = 8),
        axis.text.x = element_text(size = 8))

pdf(out_pdf, width = 20, height = 24)
print(p)
invisible(dev.off())

message("wrote ", out_pdf)
