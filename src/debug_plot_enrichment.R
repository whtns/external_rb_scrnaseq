source("packages.R")
source("functions.R")
library(targets)
tar_load(c("oncoprint_enrich_clusters"))

debug(plot_cis_trans_enrichment_recurrence)

plot_cis_trans_enrichment_recurrence(oncoprint_enrich_clusters)

tar_load(c("oncoprint_enrich_clones"))

plot_cis_trans_enrichment_recurrence(oncoprint_enrich_clones)

