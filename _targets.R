
# All R scripts in ./R/ are sourced below, including packages.R and functions.R if present.
## Load your packages, e.g. library(targets).
suppressPackageStartupMessages(source("./packages.R"))
## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

tar_option_set(
  memory = "transient",
  garbage_collection = TRUE,
  error = "continue",
  workspace_on_error = TRUE
  # controller = crew_controller_local(workers = 4)  # temporarily disabled to debug
)

## _targets.R must return a list of tar_target objects.
## Each pipeline_targets_* variable is a list of targets; c() flattens them.
c(
  pipeline_targets_inputs,
  pipeline_targets_seurat,
  pipeline_targets_diffex,
  pipeline_targets_integration,
  pipeline_targets_figures
)

# =============================================================================
# Sample-level analysis notes (2023-03-07)
# =============================================================================
#
# INCLUDED SAMPLES (interesting_samples):
#   wu  SRR13884242: likely 16q GT; f31 sample; cluster 1 (prolif.) enriched GT 2;
#                    can attribute proliferation to 16q- acquisition in GT 2
#   wu  SRR13884243: likely 16q GT; biological replicate of SRR13884242;
#                    cluster 4 analogue of SRR13884242 c3
#   wu  SRR13884247: likely 16q GT; no GSEA diffex output; c5 marker DLK1 in
#                    miRNA cluster with MEG3
#   wu  SRR13884249: possible 1q/2p/16q GT; GT 1 (only 16q-) does not contribute
#                    to proliferating c2/c4 (TOP2A high); split by phase
#   yang SRR14800534: GT1 lacks 16q-; does not contribute to proliferating c2/c4;
#                     GT 2 lacks 1q—can't identify tx distinction b/w GT2/3
#   yang SRR14800535: 16q GT; GT 1 decreased in c1 (high TOP2A)
#   yang SRR14800536: possible 16q GT; GT1 decreased in c2/c3 (high G2M markers)
#   yang SRR14800540: clear 16q GT; three GTs with SCNAs; each progressively more
#                     proliferative; GT5/c4 markers: C1QA/B, CD74, HLA genes
#   yang SRR14800541: clear 1q/6p/16q GTs; GT1 not proliferating (no c2 G2M contribution)
#   yang SRR14800543: possible minor 16q/1q GT; dual/individual 1q+16q contribution
#                     with 13q CNLOH; MYCN marks clusters
#   field SRR17960481: 6p interesting; c2/c4 notable; c2 has mito genes; likely
#                      clonal 6p with PRs and stressed cells in GT 1
#   field SRR17960484: c1 enriched for wt GT 1; also has high Xist expression
#
# MAYBE (not in interesting_samples):
#   SRR13884240: possible 2p GT
#   SRR13884241: possible 1q GT
#   SRR13884244: possible 1q GT
#   SRR13884245: possible 1q GT
#   SRR13884246: possible 16q GT
#   SRR17960480: possible minor 16q- GT
#   SRR14800539: possible 16q GT
#
# EXCLUDED:
#   SRR14800537: possible 16q GT — excluded
#   SRR17960482: too complicated
#   SRR13884248: clear 6p (missing possible 2p in expression)
#   SRR17960483: cluster 6 maybe interesting
