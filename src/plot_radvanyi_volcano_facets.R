
genes_diff_expressed <-
  "data/Liu_Radvanyi_2022_supp_data/41467_2021_25792_MOESM6_ESM.xlsx" %>%
  readxl::read_excel(sheet = 1, skip = 2) %>%
  janitor::clean_names() %>%
  dplyr::rename(symbol = gene) %>%
  # dplyr::filter(gene_cluster %in% c("1.2", "2")) %>%
  dplyr::mutate(gene_cluster = paste0("subtype", str_remove(gene_cluster, "\\..*"))) %>%
  dplyr::left_join(annotables::grch38, by  = "symbol") %>%
  dplyr::arrange(chr) %>%
  dplyr::distinct(symbol, .keep_all = TRUE) %>%
  tibble::column_to_rownames("symbol") %>%
  dplyr::mutate(region = case_when(
    (chr == "1" & start >= "123400000") ~ "1q+",
    (chr == "2" & start < "93900000") ~ "2p+",
    (chr == "6" & start < "59800000") ~ "6p+",
    chr == "16" ~ "16q-",
    .default = "other")) %>%
  # dplyr::filter(chr %in% c("1", "6", "16")) %>%
  identity()

regions <- factor(genes_diff_expressed$region)
mycols <- scales::hue_pal()(nlevels(mypal))
mypal <- regions
levels(mypal) <- mycols
names(mypal) <- regions


labels <- data.frame(region = unique(genes_diff_expressed$region),
                     label = unique(genes_diff_expressed$region)
                     )

EnhancedVolcano::EnhancedVolcano(genes_diff_expressed,
                                 x = "log_fc_subtype_2_vs_subtype_1",
                                 y = "adjusted_p_value",
                                 lab = NA,
                                 pointSize = 0.5,
                                 # lab = rownames(genes_diff_expressed),
                                 colCustom = mypal) +
  facet_wrap(~region, ncol = 2) +
  # facet_wrap(~chr) +
  theme_light() +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_text(
      size = 20)
  ) +
  guides(color = "none") +
  xlim(-6,6) +
  ylim(0,14) +
  geom_text(data = labels, aes(label = label),
             x = Inf, y = -Inf, hjust=1, vjust=-1,
             inherit.aes = FALSE, size = 10) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  NULL


  ggsave("results/radvanyi_chr_distribution.pdf", height = 6, width = 6)

browseURL("results/radvanyi_chr_distribution.pdf")
