#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(targets)

source("packages.R")
source("functions.R")

tar_load(c("debranched_seus_1q", "debranched_seus_2p", "debranched_seus_6p", "debranched_seus_16q"))

plotted_seus <- c(debranched_seus_1q, debranched_seus_2p, debranched_seus_6p, debranched_seus_16q)

plotted_seus <- unique(plotted_seus)

names(plotted_seus) <- str_remove(fs::path_file(plotted_seus), "_filtered_seu.rds")

# undebug(plot_cc_space_plot)

test0 <- map(plotted_seus, plot_cc_space_plot, color_by = "scna")

map(test0, browseURL)
