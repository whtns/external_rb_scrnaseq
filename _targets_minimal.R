## Load your packages, e.g. library(targets).
source("./packages.R")
source("./functions.R")
# plan(callr)
# debug(make_numbat_heatmaps)
# debug(make_numbat_plot_files)
# debug(diffex_cells)
# debug(filter_numbat_cells)
# debug(diffex_by_cluster)
# debug(enrich_diffex_by_cluster)
# debug(enrich_by_cluster)
# debug(find_diffex_bw_clusters_for_each_clone)

tar_option_set(memory = "transient", garbage_collection = TRUE)

output_plot_extensions <- c(
  "dimplot.pdf",
  "merged_marker.pdf",
  "combined_marker.pdf",
  "phylo_probability.pdf",
  "clone_distribution.pdf"
)

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

## tar_plan supports drake-style targets and also tar_target()
tar_plan(
  tar_target(
    analyzed_samples, c(
      "SRR13884242",
      "SRR13884243",
      "SRR13884247",
      "SRR13884249",
      "SRR14800534",
      "SRR14800535",
      "SRR14800536",
      "SRR14800540",
      "SRR14800541",
      "SRR14800543",
      "SRR17960481",
      "SRR17960484"
    )
  ),

  # end of plan------------------------------
)

# notes 2023-03-07 ------------------------------

# yes ------------------------------
# wu SRR13884242 (RB2 rep 1): likely 16q GT; f31 sample; cluster 1 (prolif.) has higher proportion of GT 2; can attribute this proliferation to acquisition of 16q- in GT 2
# wu SRR13884243 (RB2 rep 2): likely 16q GT; biological replication of SRR13884242; cluster 4 is an analogue of SRR13884242 c3 # nolint

# wu SRR13884247: likely 16q GT; no GSEA diffex output; cluster 5 is interesting; c5 marker DLK1 in amicroRNA cluster with MEG3
# wu SRR13884249: possible 1q/2p/16q GT; evidence that GT 1 (only 16q-) does not contribute to proliferating clusters c2 and c4 (TOP2A high); split by phase
# yang SRR14800534: GT1 lacks 16q-; does not contribute to proliferating clusters c2 and c4; GT 2 lacks 1q--can't identify tx distinction b/w gt2/3
# yang SRR14800535: 16q GT; GT 1 decreased in c1; c1 high TOP2A
# yang SRR14800536: possible 16q GT; GT1 decreased in c2/3 high G2M markers
# yang SRR14800540: clear 16q GT; very complicated tumor; three GTs with SCNAs; each is proliferating to a greater degree; confused about GT5 with high c4; markers are C1QA/B, CD74, HLA genes
# yang SRR14800541: clear 1q 6p and 16q GTs; GT1 not proliferating--no contribution to c2 with G2M markers
# yang SRR14800543: possible minor 16q/1q GT; can identify dual or individual contribution of 1q and 16q with 13q CNLOH; strange MYCN marker of clusters
# field SRR17960481: 6p interesting; cluster 2 and 4; cluster 2 has mito genes; not promising; likely clonal 6p with PRs and stressed cells composing GT 1
# field SRR17960484: cluster 1 enriched for wt GT 1; also has high Xist expression

# maybe ------------------------------
# SRR13884240: possible 2p GT
# SRR13884241: possible 1q GT
# SRR13884244: possible 1q GT
# SRR13884245: possible 1q GT
# SRR13884246: possible 16q GT
# SRR17960480: possible minor 16q- GT
# SRR14800539: possible 16q GT

# no ------------------------------
# SRR14800537: possible 16q GT
# SRR17960482: too complicated
# SRR13884248: clear 6p (missing possible 2p in expression)
# SRR17960483: cluster 6 maybe interesting
