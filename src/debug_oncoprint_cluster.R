source("packages.R")
source("functions.R")
library(targets)
tar_load("oncoprint_input_by_scna")

tar_load(c("large_filter_expressions", "cluster_dictionary", "interesting_samples", "large_in_segment_diffex_clones_for_each_cluster", "large_out_of_segment_diffex_clones_for_each_cluster", "large_clone_comparisons")
)

debug(enrich_oncoprints_clusters)

large_in_segment_diffex_clones_for_each_cluster <- map(large_in_segment_diffex_clones_for_each_cluster, read_csv)

large_out_of_segment_diffex_clones_for_each_cluster <- map(large_out_of_segment_diffex_clones_for_each_cluster, read_csv)

oncoprint_enrich_plots <- enrich_oncoprints_clusters(large_filter_expressions,
                                                     cluster_dictionary,
                                                     interesting_samples,
                                                     large_in_segment_diffex_clones_for_each_cluster,
                                                     large_out_of_segment_diffex_clones_for_each_cluster,
                                                     large_clone_comparisons)
