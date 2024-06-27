source("packages.R")
source("functions.R")
library(targets)
tar_load("oncoprint_input_by_scna")

tar_load(c("large_filter_expressions", "cluster_dictionary", "interesting_samples", "large_in_segment_diffex_clones", "large_out_of_segment_diffex_clones", "large_clone_comparisons")
)

all_comps_cis <- oncoprint_input_by_scna$cis

all_comps_trans <- oncoprint_input_by_scna$trans

# debug(enrich_oncoprints)
# debug(enrichment_analysis)

oncoprint_enrich_plots <- enrich_oncoprints(large_filter_expressions, cluster_dictionary, interesting_samples, large_in_segment_diffex_clones, large_out_of_segment_diffex_clones, large_clone_comparisons)


# plotting ------------------------------


cis_plots <-
  oncoprint_enrich_plots$cis %>%
  # map(~map(.x, list_flatten)) %>%
  purrr::list_flatten() %>%
  purrr::list_flatten() %>%
  imap(~{.x + labs(title = .y)}) %>%
  identity()

trans_plots <-
  oncoprint_enrich_plots$trans %>%
  # map(~map(.x, list_flatten)) %>%
  purrr::list_flatten() %>%
  purrr::list_flatten() %>%
  imap(~{.x + labs(title = .y)}) %>%
  identity()


pdf("~/tmp/cis.pdf")
cis_plots
dev.off()

browseURL("~/tmp/cis.pdf")

pdf("~/tmp/trans.pdf")
trans_plots
dev.off()

browseURL("~/tmp/trans.pdf")
