zinovyev_genes <- read_zinovyev_genes() %>%
  tibble::enframe("zphase", "symbol") %>%
  tidyr::unnest(symbol) %>%
  identity()

summary_of_cluster_diffex <-
  table_large_cis_diffex_clones[["cluster"]][[1]] %>%
  myreadxl() %>%
  dplyr::bind_rows(.id = "sample_id") %>%
  dplyr::filter(chr %in% c(1, 2, 6, 13, 16)) %>%
  dplyr::left_join(zinovyev_genes, by = "symbol") %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(recurrence = dplyr::n()) %>%
  dplyr::arrange(desc(recurrence), symbol) %>%
  identity()

summary_of_total_diffex <-
  table_large_cis_diffex_clones[["total"]][[1]] %>%
  myreadxl() %>%
  dplyr::bind_rows(.id = "sample_id") %>%
  dplyr::filter(chr %in% c(1, 2, 6, 13, 16)) %>%
  dplyr::left_join(zinovyev_genes, by = "symbol") %>%
  dplyr::group_by(symbol) %>%
  dplyr::mutate(recurrence = dplyr::n()) %>%
  dplyr::arrange(desc(recurrence), symbol) %>%
  identity()


