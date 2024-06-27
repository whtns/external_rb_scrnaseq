#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

plot_cis_trans_enrichment_recurrence <- function(enrichment_tables){

  plot_enrichment_recurrence <- function(enrichment_tables, num_recur = 2, mytitle = "", n_slice = 10){
    # browser()

    df <-
      enrichment_tables %>%
      purrr::list_flatten() %>%
      purrr::discard(is.na) %>%
      map(~{.x@result}) %>%
      dplyr::bind_rows(.id = "comparison") %>%
      identity()


    test0 <-
      df %>%
      # dplyr::arrange(symbol, sample_id) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::group_by(Description) %>%
      dplyr::mutate(neg_log10_p_val_adj = -log(p.adjust, base = 10)) %>%
      dplyr::filter(p.adjust < 0.05) %>%
      dplyr::filter(n_distinct(comparison) >= num_recur) %>%
      identity()

    test1 <-
      test0 %>%
      dplyr::summarize(mean_NES = mean(NES)) %>%
      dplyr::slice_max(abs(mean_NES), n = n_slice) %>%
      dplyr::inner_join(test0, by = "Description") %>%
      dplyr::mutate(comparison = str_replace(comparison, "_", "\n"))

    enrich_plot <-
      test1 %>%
      ggplot(aes(x = comparison, y = Description, size = neg_log10_p_val_adj, color = NES)) +
      scale_color_gradient2() +
      geom_point() +
      labs(title = mytitle, size  = "-log10 p.adj")

    return(enrich_plot)
  }

  # cis ------------------------------
  recurrences <- c(3,1,1,1)

  titles = c(
    "1q+ cis enrichment",
    "2p+ cis enrichment",
    "6p+ cis enrichment",
    "16q- cis enrichment"
             )

  cis_enrich_plots <- pmap(list(enrichment_tables$cis, recurrences, titles), plot_enrichment_recurrence)

  # trans ------------------------------

  recurrences <- c(2,2,2,3)

  titles = c(
    "1q+ trans enrichment",
    "2p+ trans enrichment",
    "6p+ trans enrichment",
    "16q- trans enrichment"
  )

  trans_enrich_plots <- pmap(list(enrichment_tables$trans, recurrences, titles), plot_enrichment_recurrence)

}
