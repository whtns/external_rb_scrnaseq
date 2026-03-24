
# Project-local R library — takes priority over site-library.
# Install modified packages here with: devtools::install_local("path/to/pkg", lib = ".R_libs/")
local({
  proj_lib <- "/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/.R_libs"
  dir.create(proj_lib, showWarnings = FALSE, recursive = TRUE)
  .libPaths(c(proj_lib, .libPaths()))
})

source("~/.Rprofile")
source("/home/skevin/single_cell_projects/resources/external_rb_scrnaseq_proj/R/targets_utils.R")

# Sys.setenv(RETICULATE_PYTHON = "/opt/miniconda3/envs/scvi-env/bin/python")
Sys.setenv(RETICULATE_PYTHON = "/dataVolume/miniconda3/envs/scvelo/bin/python")
# Sys.setenv(RETICULATE_PYTHON = "/opt/miniconda3/envs/scvelo/bin/python")
# Sys.setenv(RETICULATE_PYTHON = "/opt/miniconda3/envs/velocyto/bin/python")
