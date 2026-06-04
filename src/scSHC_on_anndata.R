#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')
library(scSHC)

mymat2 <- data.table::fread("results/scSHC/X.csv", header = FALSE) %>%
  as.matrix() %>%
  t()

vars <- read.csv("results/scSHC/var.csv")
obs <- read.csv("results/scSHC/obs.csv")

colnames(mymat2) <- obs[["X"]]
rownames(mymat2) <- vars[["X"]]

# scSHC(mymat2)

scSHC::testClusters(mymat2, obs$leiden)
