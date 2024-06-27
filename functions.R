## functions

#' Plot numbat
#'
#' @param nb
#' @param myseu
#' @param myannot
#' @param mytitle
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_numbat <- function(nb, myseu, myannot, mytitle, sort_by = "scna", ...) {
  # celltypes <-
  #   myseu@meta.data["type"] %>%
  #   tibble::rownames_to_column("cell") %>%
  #   dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
  #   identity()
  #
  # myannot <- dplyr::left_join(myannot, celltypes, by = "cell")

  # mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

  num_cols <-  length(unique(nb$clone_post$clone_opt))
  mypal <- scales::hue_pal()(num_cols) %>%
    set_names(seq(num_cols))

  myheatmap <- nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = FALSE,
    sort_by = sort_by,
    annot_bar_width = 1,
    raster = FALSE,
    ...
  ) +
    labs(title = mytitle)
  
  return(myheatmap)
}

#' Title
#'
#' @param nb
#' @param myseu
#' @param myannot
#' @param mytitle
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_numbat_w_phylo <- function(nb, myseu, myannot, mytitle, ...) {
  # browser()
  celltypes <-
    myseu@meta.data["type"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    identity()

  myannot <- dplyr::left_join(myannot, celltypes, by = "cell")

  mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

  nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = TRUE,
    annot_bar_width = 1,
    ...
  ) +
    labs(title = mytitle)
}

safe_plot_numbat <- safely(plot_numbat, otherwise = NA_real_)

safe_plot_numbat_w_phylo <- safely(plot_numbat_w_phylo, otherwise = NA_real_)

# plot variability at SCNA
#' Title
#'
#' @param phylo_plot_output
#' @param chrom
#' @param p_min
#'
#' @return
#' @export
#'
#' @examples
plot_variability_at_SCNA <- function(phylo_plot_output, chrom = "1", p_min = 0.9){
  # browser()
  test0 <-
    phylo_plot_output %>%
    dplyr::mutate(seg = factor(seg, levels = str_sort(unique(seg), numeric = TRUE))) %>%
    # dplyr::filter(CHROM == chrom) %>%
    identity()

  p_cnv_plot <- ggplot(test0, aes(x = cell_index, y = p_cnv, color = cnv_state)) +
    geom_point(size = 0.1, alpha = 0.1) +
    # geom_hline(aes(yintercept = p_min)) +
    scale_x_reverse() +
    scale_color_manual(values = c("amp" = "#7f180f",
                                  "bamp" = "pink",
                                  "del" = "#010185",
                                  "loh" = "#387229")) +
    facet_wrap(~seg)

  return(p_cnv_plot)
}

#' Title
#'
#' @param myseus
#'
#' @return
#' @export
#'
#' @examples
read_regress_save <- function(myseus){
  seu <- readRDS(myseus[[sample_id]])

  seu <- ScaleData(seu, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seu))

  # cell cycle effects strongly mitigated in PCA
  seu <- RunPCA(seu, features = VariableFeatures(seu), nfeatures.print = 10)

  seu <- FindNeighbors(seu, dims = 1:10)
  seu <- FindClusters(seu, resolution = 0.15)

  seu <- RunUMAP(seu, dims = 1:10, min.dist = 0.01)

  saveRDS(seu, myseus)
}

#' Title
#'
#' @param myseus
#'
#' @return
#' @export
#'
#' @examples
read_unregress_cc_save <- function(myseus){
  seu <- readRDS(myseus[[sample_id]])

  seu <- seuratTools::clustering_workflow(seu, resolution = c(0.2, 0.4))

  # cell cycle effects strongly mitigated in PCA
  # seu <- seuratTools::seurat_reduce_dimensions(seu)

  # seu <- RunPCA(seu, features = VariableFeatures(seu), nfeatures.print = 10)
  #
  # seu <- RunUMAP(seu, dims = 1:30)

  # saveRDS(seu, myseus)
}

annotate_seu_with_rb_subtype_gene_expression <- function(seu){

  hallmark_gene_sets = msigdbr(species = "Homo sapiens", category = "H")

  subtype_genes <- list(
    gp1 = c("EGF", "TPBG", "GUCA1C", "GUCA1B", "GUCA1A", "GNAT2", "GNGT2", "ARR3", "PDE6C", "PDE6H", "OPN1SW"),
    gp2 = c("TFF1", "CD24", "EBF3", "GAP43", "STMN2", "POU4F2", "SOX11", "EBF1", "DCX", "ROBO1", "PCDHB10")
  )

  seu <- AddModuleScore(
    object = seu,
    features = subtype_genes,
    ctrl = 5,
    name = 'exprs_gp'
  )

  subtype_hallmarks <-
    list(gp1 = c(
      "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_INFLAMMATORY_RESPONSE",
      "HALLMARK_COMPLEMENT", "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
      "HALLMARK_FATTY_ACID_METABOLISM", "HALLMARK_PEROXISOME", "HALLMARK_BILE_ACID_METABOLISM",
      "HALLMARK_PROTEIN_SECRETION"
    ), gp2 = c(
      "HALLMARK_G2M_CHECKPOINT",
      "HALLMARK_E2F_TARGETS", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MITOTIC_SPINDLE",
      "HALLMARK_MYC_TARGETS_V2"
    ))

  pull_hallmark_genes <- function(hallmark_gene_set){
    hallmark_gene_sets %>%
      dplyr::filter(gs_name == hallmark_gene_set) %>%
      dplyr::pull(gene_symbol)
  }

  subtype_hallmarks <- unlist(subtype_hallmarks) %>%
    set_names(.)

  hallmark_genes <- map(subtype_hallmarks, pull_hallmark_genes)

  for (i in seq_along(hallmark_genes)){
    hallmark_genes[[i]] <- hallmark_genes[[i]][hallmark_genes[[i]] %in% rownames(seu)]
  }

  seu <- AddModuleScore(
    object = seu,
    features = hallmark_genes,
    ctrl = 5,
    name = 'hallmark'
  )

  names(seu@meta.data)[which(names(seu@meta.data) %in% paste0("hallmark", seq(1, length(hallmark_genes))))] <- names(hallmark_genes)

  return(seu)
}

#' plot distribution of clones across clusters
#'
#' @param seu
#' @param seu_name
#' @param clusters
#'
#' @return
#' @export
#'
#' @examples
plot_distribution_of_clones_across_clusters <- function(seu, seu_name, var_x = "scna", var_y = "SCT_snn_res.0.6", plot_type = c("both", "clone", "cluster"), avg_line = NULL, signif = FALSE){
	
	plot_type = match.arg(plot_type)
	
  seu_meta <- seu@meta.data %>%
    identity()

  summarized_clones <-
    seu_meta %>%
    dplyr::select(.data[[var_x]], .data[[var_y]]) %>%
    dplyr::mutate(scna ="all")

  cluster_plot <- ggplot(seu_meta) +
    geom_bar(position = "fill", aes(x = .data[[var_x]], fill = .data[[var_y]])) +
    geom_bar(data = summarized_clones, position = "fill", aes(x = .data[[var_x]], fill = .data[[var_y]])) +
    # scale_x_discrete(limits = rev) +
  	scale_x_discrete(labels = function(x) str_wrap(x, width = 20), limits = rev) + 
    coord_flip()

  summarized_clusters <-
    seu_meta %>%
    dplyr::select(.data[[var_x]], .data[[var_y]]) %>%
    dplyr::mutate(clusters ="all")
  
  # browser()

  cluster_levels <- c("all", levels(seu_meta$clusters))
  
  test0 <- dplyr::bind_rows(seu_meta, summarized_clones) %>% 
  	dplyr::mutate(clusters = factor(clusters, levels = cluster_levels))
  
  if(signif){
  	
  	cluster_levels <- c("all", levels(seu_meta$clusters))
  	
  	summarized_clusters_tabyl <-
  		summarized_clusters |> 
  		dplyr::mutate(scna = as.character(scna)) |> 
  		janitor::tabyl(clusters, scna) |>
  		identity()
  	
  	fisher_results <-
  		seu_meta |> 
  		dplyr::mutate(scna = as.character(scna)) |> 
  		janitor::tabyl(clusters, scna) %>%
  		split(.$clusters) |> 
  		map(bind_rows, summarized_clusters_tabyl) |> 
  		map(tibble::column_to_rownames, "clusters") |> 
  		# map(fisher.test, simulate.p.value=TRUE, B=1e5) %>%
  		map(fisher.test) %>%
  		map(broom::tidy) %>%
  		dplyr::bind_rows(.id = "clusters") %>%
  		dplyr::mutate(p.adjust = p.adjust(p.value)) |> 
  		dplyr::mutate(signif = 
  										symnum(p.adjust, corr = FALSE,
  													 cutpoints = c(0,  .001,.01,.05, .1, 1),
  													 symbols = c("***","**","*","."," "))
  		) |> 
  		dplyr::select(clusters, signif) |> 
  		dplyr::mutate(clusters = factor(clusters, levels = cluster_levels)) |> 
  		identity()
  	
  	cluster_levels = tidyr::unite(fisher_results, "clusters", clusters, signif, sep = " ") |> 
  		dplyr::pull(clusters)
  	
  	
  	clone_input <- 
  		seu_meta |> 
  		dplyr::bind_rows(summarized_clusters) |> 
  		dplyr::left_join(fisher_results, by = "clusters") |>
  		tidyr::unite("clusters", clusters, signif, sep = " ") |> 
  		dplyr::mutate(clusters = str_replace(clusters, "all NA", "all")) |> 
  		# mutate(clusters = factor(paste(clusters, signif, sep="_"))) |>
  		dplyr::mutate(clusters = factor(clusters, levels = c("all", cluster_levels))) |>
  		identity()
  	
  	clone_plot <- 
  		ggplot(clone_input) +
  		geom_bar(position = "fill", aes(x = .data[[var_y]], fill = .data[[var_x]])) +
  		# geom_bar(data = summarized_clusters, position = "fill", aes(x = .data[[var_y]], fill = .data[[var_x]])) +
  		scale_x_discrete(limits = rev) +
  		coord_flip() +
  		NULL
  	
  } else {
  	clone_plot <- ggplot(test0) +
  		geom_bar(position = "fill", aes(x = .data[[var_y]], fill = .data[[var_x]])) +
  		geom_bar(data = summarized_clusters, position = "fill", aes(x = .data[[var_y]], fill = .data[[var_x]])) +
  		scale_x_discrete(limits = rev) +
  		coord_flip() +
  		NULL
  }

  if(!is.null(avg_line)){
  	
    clone_plot <-
      clone_plot +
      geom_hline(yintercept = avg_line)
  }
  
  plot_return <- switch(plot_type,
  											clone = clone_plot,
  											cluster = cluster_plot,
  											both = (clone_plot / cluster_plot) +
  												plot_layout(ncol = 1) +
  												plot_annotation(title = seu_name))

  return(plot_return)
}


#' plot distribution of clones across clusters
#'
#' @param seu
#' @param seu_name
#' @param clusters
#'
#' @return
#' @export
#'
#' @examples
plot_distribution_of_clones_pearls <- function(seu, seu_name, var_x = "scna", var_y = "SCT_snn_res.0.6", avg_line = NULL){
	seu_meta <- seu@meta.data %>%
		identity()
	
	summarized_clones <-
		seu_meta %>%
		dplyr::select(.data[[var_x]], .data[[var_y]]) %>%
		dplyr::mutate(scna ="all")
	
	phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star")
	
	df <-
		dplyr::bind_rows(seu_meta, summarized_clones) %>% 
		# dplyr::mutate(clusters = fct_rev(clusters)) %>%
		group_by(.data[[var_x]], .data[[var_y]]) %>%
		dplyr::count() %>% 
		group_by(.data[[var_x]]) %>% 
		mutate(per = prop.table(n) * 100) %>%
		mutate(label = paste0(.data[[var_y]], " Fraction: \n", 
													round(per, 2), "%")) |> 
		mutate(label_lag = dplyr::lag(per, default = 0)) %>%
		mutate(label_position = cumsum(per)-0.5*per) %>%
		mutate(phase = str_remove(clusters, "_[0-9]*$")) |> 
		mutate(phase = factor(phase, levels = phase_levels)) |>
		mutate(scna = factor(scna, levels = c("all", levels(seu_meta$scna))))
		
	make_pearls_plot <- function(df, y_setting = 4){
		browser()
		
		shifter <- function(x, n = 1) {
			if (n == 0) x else c(tail(x, -n), head(x, n))
		}
		
		only_phases_df <- 
			df %>% 
			mutate(cluster_number = str_extract(clusters, "[0-9]*$")) |> 
			dplyr::filter(!phase %in% c("hsp", "hypoxia", "other", "s_star")) |> 
			group_by(scna, phase) |> 
			mutate(cluster_sequence = dplyr::row_number()) |> 
			mutate(cluster_sequence = scale(cluster_sequence, scale = FALSE) + y_setting) |> 
			dplyr::ungroup() |>
			tidyr::nest(data = -phase) |>
			dplyr::mutate(to_x = dplyr::lead(phase)) |>
			dplyr::mutate(to_x = dplyr::coalesce(to_x, phase)) |>
			# dplyr::mutate(to_x = shifter(phase)) |>
			tidyr::unnest(cols = c(data)) |>
			arrange(scna, clusters) |>
			identity()
		
		non_phases_df <- 
			df |> 
			dplyr::filter(phase %in% c("hsp", "hypoxia", "other", "s_star")) |> 
			mutate(cluster_number = str_extract(clusters, "[0-9]*$")) |>
			group_by(scna, phase) |> 
			# mutate(cluster_sequence = runif(dplyr::n())) |>
			mutate(cluster_sequence = dplyr::row_number()) |> 
			mutate(cluster_sequence = scale(cluster_sequence, scale = FALSE) + 1) |> 
			group_by(scna) |> 
			mutate(phase = sample(unique(only_phases_df$phase), dplyr::n(), replace = TRUE)) |>  
			mutate(cluster_number = clusters) |> 
			identity()
		
		n_groups <- n_distinct(only_phases_df$phase)
		
		ggplot(only_phases_df, aes(x = phase, y=cluster_sequence, size = per)) + 
			geom_rect(aes(xmin = 0, xmax = n_groups+0.5, ymin = 0, ymax = 1.9), fill = "white") +
			geom_segment(data=only_phases_df,aes(x=phase,y=y_setting, xend=to_x, yend=y_setting), size = 1, arrow=arrow(angle=15,ends='last',length=unit(0.03,'npc'), type='closed'), alpha = 0.5) +
			geom_point(aes(color = cluster_number)) +
			geom_text_repel(aes(x = phase, y=cluster_sequence, label = cluster_number), size = 4, color = "black") +
			geom_point(data = non_phases_df, aes(x = phase, y = cluster_sequence, color = cluster_number)) +
			geom_text_repel(data = non_phases_df, aes(x = phase, y=cluster_sequence, label = cluster_number), size = 4, color = "black") +
			scale_size(range = c(0, 10)) +
			coord_polar() + 
			facet_wrap(~scna, nrow = 1) + 
			ylim(0, y_setting+2) + # adjust as you like
			theme_minimal() +
			guides(
				color = "none",
				size=guide_legend(override.aes=list(alpha=0.2), theme = theme(legend.key = element_rect(colour = NA, fill = NA) ))) +
			theme(
				axis.title.x = element_blank(),
				axis.title.y = element_blank(),
				axis.text.y = element_blank(),
				text = element_text(size = 16),
				axis.text = element_text(size = 20),
				axis.title = element_text(size = 20),
				plot.title = element_text(size = 20),
				legend.text = element_text(size = 12),
				legend.title = element_text(size = 20),
				strip.text = element_text(size = 20),
				legend.spacing.y = unit(0.03, 'npc'),
				legend.position = "top"
			) +
			labs(size = "%") + 
			NULL
		
		# ggsave("~/tmp/test.pdf", height = 8, width = 8) |> 
		# 	browseURL()
	}
	
	y_setting <- 
		df |> 
		group_by(scna, phase) |> 
		summarize(n = n()) |> 
		dplyr::ungroup()
	
	pearls_plots <- 
		make_pearls_plot(df, y_setting = max(y_setting$n)+2)
	
	return(pearls_plots)
}

#' Title
#'
#' @param seu
#' @param seu_name
#' @param clusters
#'
#' @return
#' @export
#'
#' @examples
table_distribution_of_clones_across_clusters <- function(seu, seu_name, clusters = "SCT_snn_res.0.4"){
  meta <- seu@meta.data %>%
    dplyr::mutate(cluster = factor(.data[[clusters]])) %>%
    identity()

  cluster_per_scna <-
    janitor::tabyl(meta, cluster, scna) %>%
    janitor::adorn_percentages("col") %>%
    # tidyr::pivot_longer(-cluster, names_to = "scna", values_to = "percent_scna") %>%
    dplyr::mutate(sample_id = seu_name) %>%
    identity()

  scna_per_cluster <-
    janitor::tabyl(meta, scna, cluster) %>%
    janitor::adorn_percentages("col") %>%
    # tidyr::pivot_longer(-scna, names_to = "cluster", values_to = "percent_cluster") %>%
    dplyr::mutate(sample_id = seu_name) %>%
    identity()

  return(list("cluster_per_scna" = cluster_per_scna, "scna_per_cluster" = scna_per_cluster))
}

table_cluster_markers  <- function(seu, assay = "SCT"){
  cluster_names <- glue("{assay}_snn_res.{seq(0.2, 2.0, by = 0.2)}") %>%
    set_names(.)

  map(cluster_names, ~(seu@misc$markers[[.x]]$presto))

}

#' Title
#'
#' @param numbat_rds_file
#' @param cluster_dictionary
#' @param filter_expressions
#' @param clone_simplifications
#' @param extension
#'
#' @return
#' @export
#'
#' @examples
make_numbat_plot_files <- function(numbat_rds_file, seu_path, cluster_dictionary, filter_expressions = NULL, clone_simplifications = NULL, extension = ""){
  # browser()

  output_plots <- list()

  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(seu_path)

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[,!is.na(seu$clone_opt)]

  phylo_heatmap_data <- mynb$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::left_join(mynb$joint_post, by = "cell")

  plot_markers(seu, metavar = "abbreviation", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", seurat_assay = "SCT") +
    ggplot2::scale_y_discrete(position = "left") +
    labs(title = sample_id)

  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_sample_marker{extension}.pdf"))

  # output_plots[["merged_marker"]] + output_plots[["sample_marker"]]
  # ggsave(glue("results/{numbat_dir}/{sample_id}_combined_marker{extension}.pdf"))

  # seu <- annotate_seu_with_rb_subtype_gene_expression(seu)

  DimPlot(seu, group.by = c("abbreviation", "clone_opt", "Phase")) +
    plot_annotation(title = sample_id)
  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_dimplot{extension}.pdf"), width = 10)


  ## clone distribution ------------------------------
  plot_distribution_of_clones_across_clusters(seu, sample_id)
  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_clone_distribution{extension}.pdf"), width = 8, height = 4)

  ## clone tree ------------------------------

  nclones <- length(unique(mynb$clone_post$clone_opt))

  mypal <- scales::hue_pal()(nclones) %>%
    set_names(1:nclones)

  rb_scnas = clone_simplifications[[sample_id]]

  mynb <- simplify_gt(mynb, rb_scnas)

  mynb$plot_mut_history(pal = mypal, legend = FALSE, horizontal = FALSE) +
    labs(title = sample_id)

  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_tree{extension}.pdf"), width = 2, height = 5)

  # plot types ------------------------------
  plot_types <- c("dimplot", "sample_marker", "clone_distribution", "tree")

  plot_files <- glue("results/{numbat_dir}/{sample_id}/{sample_id}_{plot_types}{extension}.pdf") %>%
    set_names(plot_types)

  return(plot_files)
}

#' filter_cluster_save_seu
#'
#' @param numbat_rds_file
#' @param cluster_dictionary
#' @param filter_expressions
#' @param extension
#'
#' @return
#' @export
#'
#' @examples
filter_cluster_save_seu <- function(numbat_rds_file, cluster_dictionary, large_clone_simplifications, filter_expressions = NULL, cells_to_remove, extension = "", leiden_cluster_file = "results/adata_filtered_metadata_0.25.csv"){

  output_plots <- list()

  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  # clone post
  nb_clone_post <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_clone_post)

  seu <- seu[,!is.na(seu$clone_opt)]

  all_cells_meta <- seu@meta.data

  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")

  seu <- AddMetaData(seu, test0)

  phylo_heatmap_data <- mynb$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::left_join(mynb$joint_post, by = "cell")

  # filter out cells
  low_numbat_prob_cells <- map(filter_expressions[[sample_id]], pull_cells_matching_expression, phylo_heatmap_data) %>%
    unlist()

  # seu <- seu[,!colnames(seu) %in% low_numbat_prob_cells]

  # label SCNAs
  large_clone_simplifications <-
    tibble::enframe(large_clone_simplifications[[sample_id]], "scna", "seg")

  simplify_gt_col <- function(gt_val, scna_key){
    # browser()
    if(gt_val != ""){
      gt_vals =
        gt_val %>%
        str_split(pattern = ",") %>%
        tibble::enframe("name", "seg") %>%
        tidyr::unnest(seg) %>%
        dplyr::left_join(scna_key, by = "seg") %>%
        dplyr::filter(!is.na(scna)) %>%
        dplyr::mutate(scna = paste(scna, collapse = ",")) %>%
        dplyr::pull(scna) %>%
        unique()
    } else {
      gt_vals = ""
    }
    return(gt_vals)
  }

  scna_metadata <-
    nb_clone_post %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(scna = simplify_gt_col(GT_opt, large_clone_simplifications)) %>%
    dplyr::mutate(scna = wrap_scna_labels(scna)) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, scna_metadata)

  scna_meta <- seu@meta.data

  seu <- filter_sample_qc(seu,
                          mito_threshold = 10, nCount_threshold = 1000, nFeature_threshold = 1000)

  qc_meta <- seu@meta.data

  # remove clusters of non tumor diploid cells
  clusters_to_remove <-
    cluster_dictionary[[sample_id]] %>%
    dplyr::filter(remove == "1") %>%
    dplyr::pull(`gene_snn_res.0.2`)

  # seu <- seu[,!((seu$gene_snn_res.0.2 %in% clusters_to_remove) & seu$scna == "")]
  seu <- seu[,!(seu$gene_snn_res.0.2 %in% clusters_to_remove)]

  mysample = sample_id

  # leiden_clusters <-
  #   leiden_cluster_file %>%
  #   read_csv() %>%
  #   dplyr::filter(sample_id == mysample) %>%
  #   select(cell = `...1`, leiden) %>%
  #   dplyr::filter(cell %in% colnames(seu)) %>%
  #   tibble::column_to_rownames("cell")
  # 
  # if(nrow(leiden_clusters)>0){
  #   seu <- Seurat::AddMetaData(seu, leiden_clusters)
  # }

  # seu <- seu[,!is.na(seu$leiden)]

  # exclude low read count cells (MALAT1) ------------------------------
  if("MALAT1" %in% unique(seu$abbreviation)){
    seu <- seu[,seu$abbreviation != "MALAT1"]
  }

  seu <- seu[,!colnames(seu) %in% cells_to_remove[[sample_id]][["cell"]]]

  # run sctransform
  seu <- SCTransform(seu, assay = "gene", verbose = FALSE)

  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)

  seu <- seurat_cluster(seu = seu, resolution = seq(0.2, 2.0, by = 0.2),
                        reduction = "pca")

  seu <- find_all_markers(seu, seurat_assay = "SCT")

  filtered_seu_path <- glue("output/seurat/{sample_id}_filtered_seu.rds")
  saveRDS(seu, filtered_seu_path)

  cell_type_meta <- seu@meta.data

  plot_filtering_timeline(all_cells_meta, scna_meta, qc_meta, cell_type_meta, sample_id)

  ggsave(glue("results/{sample_id}_filtering_timeline_new.pdf"), width = 8, height = 4)

  return(filtered_seu_path)

}


prep_unfiltered_seu <- function(numbat_rds_file, cluster_dictionary, large_clone_simplifications, filter_expressions = NULL, extension = ""){
  # browser()

  output_plots <- list()

  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  # clone post
  nb_clone_post <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_clone_post)
  seu <- seu[,!is.na(seu$clone_opt)]

  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")

  seu <- AddMetaData(seu, test0)

  phylo_heatmap_data <- mynb$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::left_join(mynb$joint_post, by = "cell")

  large_clone_simplifications <-
    tibble::enframe(large_clone_simplifications[[sample_id]], "scna", "seg")

  simplify_gt_col <- function(gt_val, scna_key){
    # browser()
    if(gt_val != ""){
      gt_vals =
        gt_val %>%
        str_split(pattern = ",") %>%
        tibble::enframe("name", "seg") %>%
        tidyr::unnest(seg) %>%
        dplyr::left_join(scna_key, by = "seg") %>%
        dplyr::filter(!is.na(scna)) %>%
        dplyr::mutate(scna = paste(scna, collapse = ",")) %>%
        dplyr::pull(scna) %>%
        unique()
    } else {
      gt_vals = ""
    }
    return(gt_vals)
  }

  scna_metadata <-
    nb_clone_post %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(scna = simplify_gt_col(GT_opt, large_clone_simplifications)) %>%
    dplyr::mutate(scna = wrap_scna_labels(scna)) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, scna_metadata)

  mysample = sample_id

  seu <- SCTransform(seu, assay = "gene", verbose = FALSE)

  seu <- RunPCA(seu, verbose = FALSE)
  seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)

  seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)

  seu <- seurat_cluster(seu = seu, resolution = c(0.2, 0.4, 0.6),
                                  reduction = "pca")

  seu <- find_all_markers(seu, seurat_assay = "SCT")

  unfiltered_seu_path <- glue("output/seurat/{sample_id}_unfiltered_seu.rds")

  saveRDS(seu, unfiltered_seu_path)

  return(unfiltered_seu_path)

}

regress_filtered_seu <- function(filtered_seu_path){
  regressed_seu_path <- str_replace(filtered_seu_path, "_filtered", "_regressed")

  filtered_seu <- readRDS(filtered_seu_path)

  regressed_seu <- filtered_seu

  # store mitochondrial percentage in object meta data
  regressed_seu <- PercentageFeatureSet(regressed_seu, pattern = "^MT-", col.name = "percent.mt")

  # run sctransform
  regressed_seu <- SCTransform(regressed_seu, assay = "gene", vars.to.regress = c("percent.mt", "S.Score", "G2M.Score"), verbose = FALSE)

  regressed_seu <- RunPCA(regressed_seu, verbose = FALSE)
  regressed_seu <- RunUMAP(regressed_seu, dims = 1:30, verbose = FALSE)

  regressed_seu <- FindNeighbors(regressed_seu, dims = 1:30, verbose = FALSE)

  regressed_seu <- seurat_cluster(seu = regressed_seu, resolution = seq(0.2, 2.0, by = 0.2),
                                  reduction = "pca")

  regressed_seu <- find_all_markers(regressed_seu, seurat_assay = "SCT")

  saveRDS(regressed_seu, regressed_seu_path)

  return(regressed_seu_path)
}

#' Title
#'
#' @param numbat_rds_file
#' @param p_min
#' @param line_width
#' @param filter_expressions
#' @param cluster_dictionary
#' @param extension
#'
#' @return
#' @export
#'
#' @examples
make_numbat_heatmaps_old <- function(numbat_rds_file, filter_expressions = NULL, cluster_dictionary, p_min = 0.9, line_width = 0.1, extension = ""){
  # browser()

  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}/"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot <- mynb$clone_post[,c("cell")]

  if(!is.null(filter_expressions)){
    # filter cells
    phylo_heatmap_data <- mynb$clone_post %>%
      dplyr::select(cell, clone_opt) %>%
      dplyr::left_join(mynb$joint_post, by = "cell") %>%
      dplyr::left_join(myannot, by = "cell")

    excluded_cells <- map(filter_expressions[[sample_id]], pull_cells_matching_expression, phylo_heatmap_data) %>%
      unlist()

    seu <- seu[,!colnames(seu) %in% excluded_cells]
  }

  seu <- seu[,colnames(seu) %in% myannot$cell]

  clusters_to_remove <-
    cluster_dictionary[[sample_id]] %>%
    dplyr::filter(remove == "1") %>%
    dplyr::pull(`gene_snn_res.0.2`)

  seu <- seu[,!seu$gene_snn_res.0.2 %in% clusters_to_remove]

  myannot <-
    seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, clone_opt, nCount_gene) %>%
    identity()

  ## numbat ------------------------------
  numbat_heatmap <- safe_plot_numbat(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = p_min, line_width = line_width)[["result"]]
  # numbat_heatmap <- plot_numbat(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = p_min, line_width = line_width)[["result"]]

  scna_variability_plot <- plot_variability_at_SCNA(numbat_heatmap[[3]][["data"]])
  # patchwork::wrap_plots(numbat_heatmap, scna_variability_plot, ncol = 1)
  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_numbat_heatmap{extension}.pdf"), numbat_heatmap)

  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_numbat_probability{extension}.pdf"), scna_variability_plot)

  # ## numbat phylo ------------------------------
  # numbat_heatmap_w_phylo <- safe_plot_numbat_w_phylo(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9)[["result"]]
  #
  # scna_variability_plot_w_phylo <- plot_variability_at_SCNA(numbat_heatmap_w_phylo[[3]][["data"]])
  #
  # patchwork::wrap_plots(numbat_heatmap_w_phylo, scna_variability_plot_w_phylo, ncol = 1)
  # ggsave(glue("results/{numbat_dir}/{sample_id}_numbat_phylo_probability{extension}.pdf"))

  plot_types <- c("numbat_heatmap", "numbat_probability")

  plot_files <- glue("results/{numbat_dir}/{sample_id}/{sample_id}_{plot_types}{extension}.pdf") %>%
    set_names(plot_types)

  return(plot_files)
  # return(numbat_heatmap)
}

#' Title
#'
#' @param numbat_rds_file
#' @param p_min
#' @param line_width
#' @param filter_expressions
#' @param cluster_dictionary
#' @param extension
#'
#' @return
#' @export
#'
#' @examples
make_numbat_heatmaps <- function(seu_path, numbat_rds_file, p_min = 0.9, line_width = 0.1, extension = ""){
  # browser()

  sample_id <- str_extract(seu_path, "SRR[0-9]*")

  numbat_dir = "numbat_sridhar"

  dir_create(glue("results/{numbat_dir}/"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(seu_path)

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  myannot <- mynb$clone_post[,c("cell")]

  seu <- seu[,colnames(seu) %in% myannot$cell]

  myannot <-
    seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, scna) %>%
    identity()

  ## numbat ------------------------------
  numbat_heatmap <- safe_plot_numbat(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = p_min, line_width = line_width)[["result"]]
  # numbat_heatmap <- plot_numbat(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = p_min, line_width = line_width)

  scna_variability_plot <- plot_variability_at_SCNA(numbat_heatmap[[3]][["data"]])
  # patchwork::wrap_plots(numbat_heatmap, scna_variability_plot, ncol = 1)
  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_numbat_heatmap{extension}.pdf"), numbat_heatmap)

  ggsave(glue("results/{numbat_dir}/{sample_id}/{sample_id}_numbat_probability{extension}.pdf"), scna_variability_plot)

  # ## numbat phylo ------------------------------
  # numbat_heatmap_w_phylo <- safe_plot_numbat_w_phylo(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9)[["result"]]
  #
  # scna_variability_plot_w_phylo <- plot_variability_at_SCNA(numbat_heatmap_w_phylo[[3]][["data"]])
  #
  # patchwork::wrap_plots(numbat_heatmap_w_phylo, scna_variability_plot_w_phylo, ncol = 1)
  # ggsave(glue("results/{numbat_dir}/{sample_id}_numbat_phylo_probability{extension}.pdf"))

  plot_types <- c("numbat_heatmap", "numbat_probability")

  plot_files <- glue("results/{numbat_dir}/{sample_id}/{sample_id}_{plot_types}{extension}.pdf") %>%
    set_names(plot_types)

  return(plot_files)
  # return(numbat_heatmap)
}

#' Title
#'
#' @param output_file
#' @param sample_id
#' @param myseus
#' @param mynbs
#' @param merged_metadata
#' @param myexpressions
#'
#' @return
#' @export
#'
#' @examples
make_filtered_numbat_plots <- function(output_file, sample_id, myseus, mynbs, merged_metadata, myexpressions){

  output_plots <- list()
  seu <- readRDS(myseus[[sample_id]])

  merged_metadata_transfer <-
    merged_metadata %>%
    dplyr::filter(sample_id == {{sample_id}}) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

  mynb <- readRDS(mynbs[[sample_id]])

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot = mynb$clone_post[,c("cell", "GT_opt", "clone_opt")]

  if(!'' %in% unlist(myexpressions)){
    final_phylo_heatmap <- filter_phylo_plot(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9, expressions = myexpressions)

    test0 <- final_phylo_heatmap[[3]][["data"]] %>%
      plot_variability_at_SCNA()

    patchwork::wrap_plots(final_phylo_heatmap / test0)

    ggsave(output_file)

  }

  return(output_file)
}

#' Title
#'
#' @param sample_id
#' @param myseus
#' @param mynbs
#' @param merged_metadata
#' @param myexpressions
#'
#' @return
#' @export
#'
#' @examples
filter_numbat_cells <- function(sample_id, myseus, mynbs, merged_metadata, myexpressions){

  seu <- readRDS(myseus[[sample_id]])

  seu_meta <- seu@meta.data %>%
    tibble::rownames_to_column("cell")

  merged_metadata_transfer <-
    merged_metadata %>%
    dplyr::filter(sample_id == {{sample_id}}) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)

  mynb <- readRDS(mynbs[[sample_id]])

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "-", "\\.")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  myannot = seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    select(cell, GT_opt, clone_opt, nCount_gene, nFeature_gene) %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-"))

  if(!'' %in% unlist(myexpressions)){
    filtered_nb <- filter_phylo_plot(mynb, seu, myannot, sample_id, clone_bar = FALSE, p_min = 0.9, expressions = myexpressions)
  }

  returned_meta <- dplyr::left_join(
    filtered_nb[[3]][["data"]],
    myannot,
    by = "cell") %>%
    dplyr::mutate(cell = str_replace(cell, "-", ".")) %>%
    dplyr::filter(!is.na(cell)) %>%
    identity()

  return(returned_meta)
}

#' Title
#'
#' @param seu
#' @param seu_name
#' @param subtype_hallmarks
#'
#' @return
#' @export
#'
#' @examples
plot_rb_subtype_expression <- function(seu, seu_name, subtype_hallmarks){

  myfeatures <- c("exprs_gp1", "exprs_gp2", c(subtype_hallmarks))
  featureplots <- FeaturePlot(seu, myfeatures, combine = FALSE)

  featureplots <- map(featureplots, ~(.x + labs(subtitle = seu_name))) %>%
    set_names(myfeatures)

}

#' Title
#'
#' @param seu
#' @param checked_cluster_markers
#'
#' @return
#' @export
#'
#' @examples
plot_markers_by_cell_cycle <- plot_cluster_markers_by_cell_type <- function(seu, checked_cluster_markers){
  cluster_plots <- map(checked_cluster_markers, ~VlnPlot(seu, features = .x, group.by = "Phase"))

  cluster_plots <- map2(cluster_plots, names(checked_cluster_markers), ~(.x + labs(subtitle = .y)))

  return(cluster_plots)

}

#' Title
#'
#' @param seu
#' @param seu_title
#' @param myvar
#'
#' @return
#' @export
#'
#' @examples
plot_feature_across_seus <- function(seu, seu_title, myvar){

  varplot <- seu %>%
    FeaturePlot(features = myvar)

  dimplot <- DimPlot(seu, group.by = "gene_snn_res.0.2")

  phase_plot <- DimPlot(seu, group.by = "Phase")

  (varplot | (phase_plot / dimplot)) +
    plot_annotation(title = seu_title)

  # mypatch = wrap_plots(varplot, dimplot, phase_plot, nrow = 1) +
  # 	plot_annotation(title = seu_title)
}

compplot_feature_and_clusters <- function(seu, feature){
  fp <- FeaturePlot(seu, feature)

  cp <- DimPlot(seu, group.by = "Phase")

  dp1 <- DimPlot(seu, group.by = c("gene_snn_res.0.15", "merged_leiden"))

  mypatch <- wrap_plots(fp, cp, nrow = 1) / dp1

  return(mypatch)
}

#' Title
#'
#' @param meta_path
#'
#' @return
#' @export
#'
#' @examples
get_merged_metadata <- function(meta_path){
  metadata = read_csv(meta_path) %>%
    set_names(c("cell", "sample_id", "merged_leiden")) %>%
    dplyr::mutate(sample_id = str_remove(sample_id, ".h5ad")) %>%
    dplyr::mutate(cell = str_replace(cell, "-", ".")) %>%
    identity()

  return(metadata)
}

#' Title
#'
#' @param nb
#' @param myseu
#' @param myannot
#' @param mytitle
#' @param expressions
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
filter_phylo_plot <- function(nb, myseu, myannot, mytitle, expressions, ...) {
  # browser()
  celltypes <-
    myseu@meta.data["type"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    identity()

  myannot <- dplyr::left_join(myannot, celltypes, by = "cell")

  mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")

  initial_phylo_heatmap <- nb$plot_phylo_heatmap(
    pal_clone = mypal,
    annot = myannot,
    show_phylo = FALSE,
    sort_by = "GT_opt",
    ...
  ) +
    labs(title = mytitle)

  phylo_heatmap_data <- initial_phylo_heatmap$data %>%
    dplyr::left_join(myannot, by = "cell")

}

diffex_groups_old <- function(sample_id, myseus, celldf, ...){

  seu <- readRDS(myseus[[sample_id]])

  celldf <-
    celldf %>%
    dplyr::distinct(cell, .keep_all = TRUE) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, celldf)

  Seurat::FindMarkers(seu, ...)

}


#' Title
#'
#' @param numbat_rds_file
#' @param cluster_dictionary
#'
#' @return
#' @export
#'
#' @examples
retrieve_numbat_seurat <- function(numbat_rds_file, cluster_dictionary){
  # browser()
  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)


  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")

  seu <- AddMetaData(seu, test0)

  return(seu)

}


#' Title
#'
#' @param numbat_rds_file
#' @param filter_expressions
#' @param idents
#'
#' @return
#' @export
#'
#' @examples
diffex_groups <- function(numbat_rds_file, filter_expressions, idents = NULL){

  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  ident.1 = idents[[sample_id]]

  filter_expressions <- filter_expressions[[sample_id]]

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  joined_nb_meta <-
    mynb$clone_post %>%
    dplyr::select(-p_cnv) %>%
    dplyr::left_join(mynb$joint_post, by = "cell")

  excluded_cells <- map(filter_expressions, pull_cells_matching_expression, joined_nb_meta) %>%
    unlist()

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    dplyr::filter(!cell %in% excluded_cells) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[,!is.na(seu$clone_opt)]

  diffex <- Seurat::FindMarkers(seu, ident.1)

  gse_plot_path <- glue("results/{numbat_dir}/{sample_id}_gsea.pdf")
  gse_plot <- enrichment_analysis(diffex) +
    labs(title = glue("{sample_id}"))
  ggsave(gse_plot_path)

  diffex_path <- glue("results/{numbat_dir}/{sample_id}_diffex.csv")
  diffex <-
    diffex %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::distinct(ensgene, .keep_all = TRUE)

  write_csv(diffex, diffex_path)

  return(list(sample_id, diffex_path, gse_plot_path))

}

#' Title
#'
#' @param myseus
#' @param cells.1
#' @param cells.2
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
diffex_cells <- function(myseus, cells.1, cells.2, ...){

  seu <- readRDS(myseus[[sample_id]])

  Seurat::FindMarkers(seu$gene, cells.1 = cells.1, cells.2 = cells.2)

}

#' Title
#'
#' @param returned_meta
#' @param myseg
#'
#' @return
#' @export
#'
#' @examples
split_by_cnv <- function(returned_meta, myseg = "2a"){
  test0 <-
    returned_meta %>%
    dplyr::filter(!is.na(cell)) %>%
    dplyr::filter(seg == myseg) %>%
    dplyr::mutate(cnv_status = ifelse(p_cnv > 0.5, "present", "absent")) %>%
    dplyr::group_by(cnv_status) %>%
    dplyr::group_split() %>%
    map(pull, cell) %>%
    identity()
}

#' Title
#'
#' @param tbl
#' @param sample_id
#'
#' @return
#' @export
#'
#' @examples
plot_pcnv_by_nsnp <- function(tbl, sample_id){
  ggplot(tbl, aes(x = p_cnv, y = n_snp)) +
    geom_point(size = 0.1, alpha = 0.1) +
    facet_wrap(~seg, ncol = 1)

  ggsave(glue("results/{sample_id}_pcnv_by_nsnp.pdf"))

  return(glue("results/{sample_id}_pcnv_by_nsnp.pdf"))

}

ora_analysis <- function(seu, clusters  = "gene_snn_res.0.2", only_unique_terms = FALSE){
  # browser()

  df_list <- seu@misc$markers[[clusters]][["presto"]] %>%
    dplyr::rename(symbol = `Gene.Name`) %>%
    dplyr::left_join(annotables::grch38, by = "symbol", relationship = "many-to-many") %>%
    dplyr::distinct(Cluster, entrez, .keep_all = TRUE) %>%
    split(.[["Cluster"]])

  run_enrich_go <- function(df, fold_change_col = "Average.Log.Fold.Change"){
    # we want the log2 fold change
    original_gene_list <- df[[fold_change_col]]

    # name the vector
    names(original_gene_list) <- df$entrez

    # omit any NA values
    gene_list<-na.omit(original_gene_list)

    # sort the list in decreasing order (required for clusterProfiler)
    gene_list = sort(gene_list, decreasing = TRUE)

    clusterProfiler::enrichGO(gene = names(gene_list),
                              OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                              ont = "BP",
                              readable = TRUE)

  }

  safe_enrich_go <- purrr::safely(run_enrich_go, otherwise = NA_real_)

  ora_output <- map(df_list, safe_enrich_go) %>%
    map("result") %>%
    identity()

  ora_tables <-
    ora_output %>%
    map(tibble::as_tibble) %>%
    dplyr::bind_rows(.id = "cluster") %>%
    dplyr::group_by(Description) %>%
    # dplyr::filter(dplyr::n() == 1) %>%
    split(.$cluster) %>%
    identity()

  unique_terms <-
    ora_tables %>%
    map(dplyr::pull, ID) %>%
    identity()

  filter_ora_result_by_terms <- function(cluster_name, ora_results, terms_list){
    ora_result <- ora_results[[cluster_name]]

    myterms <- terms_list[[cluster_name]]

    ora_result@result <- ora_result@result[rownames(ora_result@result) %in% myterms,]

    return(ora_result)
  }

  if(only_unique_terms){
    ora_output <-
      names(unique_terms) %>%
      set_names(.) %>%
      map(filter_ora_result_by_terms, ora_output, unique_terms)
  }

  ora_plots <-
    ora_output %>%
    imap(~clusterProfiler::dotplot(.x, title = .y, showCategory = 20))

  return(list("result" = ora_output, "tables" = ora_tables, "plots" = ora_plots))

}

enrichment_analysis <- function(df, fold_change_col = "avg_log2FC", analysis_method = "gsea", gene_set = "hallmark", pvalueCutoff = 0.5){
  # browser()

  df <-
    df %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::distinct(entrez, .keep_all = TRUE)

  # we want the log2 fold change
  original_gene_list <- df[[fold_change_col]]

  # name the vector
  names(original_gene_list) <- df$entrez

  # omit any NA values
  gene_list<-na.omit(original_gene_list)

  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)

  if(analysis_method == "gsea"){
    if(gene_set == "hallmark"){
      gse <- clusterProfiler::GSEA(geneList=gene_list,
                                   minGSSize = 3,
                                   maxGSSize = 800,
                                   pvalueCutoff = pvalueCutoff,
                                   verbose = TRUE,
                                   TERM2GENE = msig_h,
                                   pAdjustMethod = "BH")

    } else if(gene_set == "gobp"){
      gse <- clusterProfiler::gseGO(geneList=gene_list,
                                    ont = "BP",
                                    OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                    keyType = "ENTREZID",
                                   minGSSize = 3,
                                   maxGSSize = 800,
                                   pvalueCutoff = pvalueCutoff,
                                   verbose = TRUE,
                                   pAdjustMethod = "BH") %>%
        clusterProfiler::simplify() # for GSEGO

    }

    return(gse)
  }

}

plot_enrichment <- function(gse, p_val_cutoff = 0.1, signed = TRUE, result_slot = "result"){
	
	if(signed){
		make_signed_dotplot <- function(gse, showCategory=10){
			clusterProfiler::dotplot(gse, showCategory=showCategory, split=".sign", font.size = 12) + facet_grid(.sign~., scales = "free_y")
		}
		
		possible_dotplot <- purrr::possibly(make_signed_dotplot, otherwise = ggplot())
		
	} else {
		make_dotplot <- function(gse, showCategory=10){
			clusterProfiler::dotplot(gse, showCategory=showCategory, font.size = 12)
		}
		
		possible_dotplot <- purrr::possibly(make_dotplot, otherwise = ggplot())
	}
  showCategories <- slot(gse, result_slot) |> 
  	dplyr::filter(p.adjust <= p_val_cutoff) |> 
  	pull(ID)
  
  # mydotplot <- possible_dotplot(gse, showCategory = showCategories)
  mydotplot <- possible_dotplot(gse)

  return(mydotplot)
}

make_cell_cycle_plot <- function(sample_id, myseus){

  output_plots <- list()
  seu <- readRDS(myseus[[sample_id]])

  DimPlot(seu, group.by = c("Phase")) +
    plot_annotation(title = sample_id)
  ggsave(glue("results/{sample_id}_cell_cycle.pdf"), width = 7, height = 7)

  return(glue("results/{sample_id}_cell_cycle.pdf"))
}

diffex_by_cluster <- function(sample_id, myseus, celldf, ...){
  seu <- readRDS(myseus[[sample_id]])

  celldf <-
    celldf %>%
    dplyr::distinct(cell, .keep_all = TRUE) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, celldf)

  clusters <- unique(seu$gene_snn_res.0.2) %>%
    set_names(.)

  filter_seu_to_cluster <- function(cluster, seu){
    seu[,(seu@meta.data$gene_snn_res.0.2 == cluster)]
  }

  split_seu <- map(clusters, filter_seu_to_cluster, seu)

  safe_FindMarkers <- purrr::safely(FindMarkers, otherwise = NA_real_)

  cluster_diffex <- map(split_seu, safe_FindMarkers) %>%
    map("result") %>%
    identity()

  cluster_diffex <- cluster_diffex[!is.na(cluster_diffex)] %>%
    map(tibble::rownames_to_column, "symbol")

  write_xlsx(cluster_diffex, glue("results/{sample_id}_cluster_diffex.xlsx"))

  return(glue("results/{sample_id}_cluster_diffex.xlsx"))

}

enrich_diffex_by_cluster <- function(sample_id, myseus, celldf, ...){
  seu <- readRDS(myseus[[sample_id]])

  celldf <-
    celldf %>%
    dplyr::distinct(cell, .keep_all = TRUE) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- Seurat::AddMetaData(seu, celldf)

  clusters <- unique(seu$gene_snn_res.0.2) %>%
    set_names(.)

  filter_seu_to_cluster <- function(cluster, seu){
    seu[,(seu@meta.data$gene_snn_res.0.2 == cluster)]
  }

  clusters <- janitor::tabyl(seu@meta.data, gene_snn_res.0.2, clone_opt) %>%
    rowwise() %>%
    mutate(empty_clone = min(c_across(!any_of("gene_snn_res.0.2")))) %>%
    dplyr::filter(!empty_clone < 3) %>%
    dplyr::pull(`gene_snn_res.0.2`) %>%
    set_names(.) %>%
    identity()

  split_seu <- map(clusters, filter_seu_to_cluster, seu)

  # cluster_diffex <- map(split_seu, Seurat::FindMarkers, ...)

  possible_FindMarkers <- purrr::possibly(FindMarkers, otherwise = NA_real_)

  cluster_diffex <- map(split_seu, possible_FindMarkers, ...) %>%
    # map("result") %>%
    identity()

  cluster_diffex <- cluster_diffex[!is.na(cluster_diffex)]

  safe_enrichment_analysis <- purrr::safely(enrichment_analysis, otherwise = NA_real_)

  enrich_plots <- map(cluster_diffex, safe_enrichment_analysis) %>%
    map("result") %>%
    identity()

  enrich_plots <- enrich_plots[!is.na(enrich_plots)] %>%
    compact()

  enrich_plots <- compact(enrich_plots) %>%
    imap(~(.x + labs(title = .y))) %>%
    identity()

  pdf(glue("results/{sample_id}_diffex_cluster_enrichment.pdf"), width = 10)
  print(enrich_plots)
  dev.off()

  return(glue("results/{sample_id}_diffex_cluster_enrichment.pdf"))


}

enrich_by_cluster  <- function(seu_path){
  # browser()

	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
  
  seu <- readRDS(seu_path)

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  cluster_diffex <- seu@misc$markers$clusters$presto %>%
    split(.$Cluster) %>%
    map(tibble::column_to_rownames, "Gene.Name")

  drop_cc_genes <- function(df, cc.genes = Seurat::cc.genes){
    cc_genes <- unlist(cc.genes)

    dplyr::filter(df, !rownames(df) %in% cc_genes)
  }

  cluster_diffex <-
    cluster_diffex %>%
    # map(drop_cc_genes) |> 
  	identity()

  safe_enrichment_analysis <- purrr::safely(enrichment_analysis, otherwise = NA_real_)

  enrich_results <- map(cluster_diffex, safe_enrichment_analysis, fold_change_col = "Average.Log.Fold.Change") %>%
    map("result") %>%
    identity()
  
  enrich_plots <- enrich_results |> 
  	map(plot_enrichment, p_val_cutoff = 1, signed = FALSE) |> 
  	compact() %>%
    imap(~(.x + labs(title = sample_id, subtitle = .y))) %>%
    identity()

  gse_plot_path <- glue("results/{sample_id}_cluster_gsea.pdf")
  
  pdf(gse_plot_path, width = 8, height = 10)
  print(enrich_plots)
  dev.off()

  return(gse_plot_path)

}

plot_pcnv_by_reads <- function(tbl, sample_id){
  ggplot(tbl, aes(x = p_cnv, y = nCount_gene)) +
    geom_point(size = 0.1, alpha = 0.1) +
    facet_wrap(~seg, ncol = 1)

  ggsave(glue("results/{sample_id}_pcnv_by_reads.pdf"))

  return(glue("results/{sample_id}_pcnv_by_reads.pdf"))

}

convert_numbat_pngs <- function(numbat_rds_file){

  numbat_output_dir <- str_remove(numbat_rds_file, "_numbat.*")

  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_pngs <- dir_ls(numbat_output_dir, glob = "*.png") %>%
    set_names(path_file(.))

  numbat_pdfs <- stringr::str_replace(path_file(numbat_pngs), ".png", ".pdf")

  numbat_pdfs <- glue("results/{sample_id}/{numbat_pdfs}")
  dir_create(glue("results/{sample_id}"))

  numbat_images <- purrr::map(numbat_pngs,image_read) %>%
    imap(~image_annotate(.x, sample_id, size = 50))

  map2(numbat_images, numbat_pdfs, ~image_write(.x, format = "pdf", .y))

  return(numbat_pdfs)
}

compare_infercnv <- function(myseus, sample_id){
  seu <- readRDS(myseus[[sample_id]])
}

make_all_numbat_plots <- function(numbat_dir, num_iter = 2, min_LLR = 2, genome = "hg38", init_k = 3, gtf = gtf_hg38, overwrite = FALSE){

  sample_id = path_file(numbat_dir)

  print(numbat_dir)

  for (i in seq(num_iter)){

    bulk_clone_path = glue("{numbat_dir}/bulk_clones_{i}.png")
    if(!file.exists(bulk_clone_path) | (overwrite = TRUE)){
      # Plot bulk clones
      bulk_clones <- read_tsv(glue("{numbat_dir}/bulk_clones_{i}.tsv.gz"), col_types = cols())
      p = plot_bulks(bulk_clones, min_LLR = min_LLR, use_pos = TRUE,
                     genome = genome) +
        labs(title = sample_id)

      ggsave(bulk_clone_path, p,
             width = 13, height = 2 * length(unique(bulk_clones$sample)),
             dpi = 250)
      print(glue("plotted {numbat_dir}/bulk_clones_{i}.png"))
    }


    bulk_subtrees_path = glue("{numbat_dir}/bulk_subtrees_{i}.png")
    if(!file.exists(bulk_subtrees_path) | (overwrite = TRUE)){

      # Plot bulk subtrees
      bulk_subtrees <- read_tsv(glue("{numbat_dir}/bulk_subtrees_{i}.tsv.gz"), col_types = cols())
      p = plot_bulks(bulk_subtrees, min_LLR = min_LLR,
                     use_pos = TRUE, genome = genome) +
        labs(title = sample_id)

      ggsave(glue("{numbat_dir}/bulk_subtrees_{i}.png"),
             p, width = 13, height = 2 * length(unique(bulk_subtrees$sample)),
             dpi = 250)
    }

  }

  final_bulk_clones_path = glue("{numbat_dir}/bulk_clones_final.png")
  if(!file.exists(final_bulk_clones_path) | (overwrite = TRUE)){

    bulk_clones <- read_tsv(glue("{numbat_dir}/bulk_clones_final.tsv.gz"), col_types = cols())
    p = plot_bulks(bulk_clones, min_LLR = min_LLR, use_pos = TRUE,
                   genome = genome) +
      labs(title = sample_id)
    ggsave(final_bulk_clones_path, p, width = 13,
           height = 2 * length(unique(bulk_clones$sample)),
           dpi = 250)

  }

  # exp_clust_path = glue("{numbat_dir}/exp_roll_clust.png")
  # if(!file.exists(exp_clust_path)){
  #
  #   # # Plot single-cell smoothed expression magnitude heatmap
  #   gexp_roll_wide <- read_tsv(glue("{numbat_dir}/gexp_roll_wide.tsv.gz"), col_types = cols()) %>%
  #     column_to_rownames("cell")
  #   hc <- readRDS(glue("{numbat_dir}/hc.rds"))
  #   p = plot_exp_roll(gexp_roll_wide = gexp_roll_wide,
  #                     hc = hc, k = init_k, gtf = gtf, n_sample = 10000)
  #   labs(title = sample_id)
  #   ggsave(exp_clust_path, p,
  #          width = 8, height = 4, dpi = 200)
  #
  # }

  # phylo_heatmap_path = glue("{numbat_dir}/phylo_heatmap.png")
#
#   if(!file.exists(phylo_heatmap_path)){
#
#     # # Plot single-cell CNV calls along with the clonal phylogeny
#     nb <- readRDS(glue("{numbat_dir}_numbat.rds"))
#     mypal = c('1' = 'gray', '2' = "#377EB8", '3' = "#4DAF4A", '4' = "#984EA3")
#
#     phylo_heatmap <- nb$plot_phylo_heatmap(
#       pal_clone = mypal,
#       show_phylo = TRUE
#     ) +
#       labs(title = sample_id)
#     ggsave(phylo_heatmap_path, width = 13,
#            height = 10,
#            dpi = 250)
#
#   }


    return("success!")

}

retrieve_numbat_rds_files <- function(numbat_dir, kept_samples = NULL){

  numbat_rds_files <- fs::dir_ls(numbat_dir, regexp = ".*SRR[0-9]*_numbat.rds", recurse = TRUE) %>%
    set_names(str_extract(., "SRR[0-9]*"))

  if(!is.null(kept_samples)){
    numbat_rds_files <- numbat_rds_files[names(numbat_rds_files) %in% kept_samples]
  }

  return(numbat_rds_files)
}

retrieve_seus <- function(seu_dir, kept_samples = NULL){

  seus <- fs::dir_ls(seu_dir, regexp = "./SRR[0-9]*_seu.rds", recurse = TRUE) %>%
    set_names(str_extract(., "SRR[0-9]*"))

  if(!is.null(kept_samples)){
    seus <- seus[names(seus) %in% kept_samples]
  }

  return(seus)
}


retrieve_numbat_plot_type <- function(numbat_plots, plot_type = "exp_roll_clust.png"){
  browser()
  retrieved_plot_types <- map(numbat_plots, ~set_names(.x, fs::path_file(.x))) %>%
    map(str_detect, plot_type)

  retrieved_plots <- map2(numbat_plots, retrieved_plot_types, ~{.x[.y]}) %>%
    unlist()


}

reroute_done_to_results_pdf <- function(numbat_rds_file, label = ""){

  sample_id = str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  numbat_pdfs <- dir_ls(glue("results/{numbat_dir}/{sample_id}/"), glob = "*.pdf")

  results_file <- glue("results/{numbat_dir}/{sample_id}{label}.pdf")

  qpdf::pdf_combine(numbat_pdfs, results_file)

  return(results_file)

}

retrieve_snakemake_params <- function(numbat_rds_file){

  str_extract(numbat_rds_file, "SRR[0-9]*")

  log_file <- fs::path(path_dir(numbat_rds_file), "log.txt")

  log <- read_lines(log_file)[3:26] %>%
    str_split(" = ") %>%
    transpose() %>%
    identity()

  params <- log[[2]] %>%
    set_names(log[[1]])

  return(list(sample_id, params))

}

retrieve_current_param <- function(current_params, myparam){

  sample_ids <- map(current_params, 1)

  param_values <- map(current_params, 2) %>%
    map(myparam) %>%
    set_names(sample_ids)
}

plot_putative_marker_across_samples <- function(mymarkers, seu_paths, plot_type = FeaturePlot, group_by = "gene_snn_res.0.2", cluster_dictionary, extension = "_filtered"){
  print(mymarkers)
  plot_markers_in_sample <- function(seu_path, mymarkers, plot_type = FeaturePlot, group_by = group_by, cluster_dictionary){
    # browser()
    sample_id <- str_extract(seu_path, "SRR[0-9]*")

    numbat_dir = "numbat_sridhar"

    dir_create(glue("results/{numbat_dir}"))
    dir_create(glue("results/{numbat_dir}/{sample_id}"))

    seu <- readRDS(seu_path)

    mymarkers <- mymarkers[mymarkers %in% rownames(seu)]

    if(identical(plot_type, VlnPlot)){
      feature_plots_first <- plot_type(seu, features = mymarkers, group.by = group_by, combine = FALSE, pt.size  = 0) %>%
        set_names(mymarkers)

      max_ys = map(feature_plots_first, ~layer_scales(.x)$y$get_limits()) %>%
        map(2) %>%
        identity()

      feature_plots <- map2(feature_plots_first, max_ys, ~{
        .x +
          # expand_limits(y = c(0, .y*2.5)) +
          stat_compare_means(comparisons = list(c(1, 2)), method = "t.test", label.y = .y*0.9) +
          # stat_compare_means() +
          geom_boxplot(width = 0.2) +
          theme(legend.position="none",
                axis.title.x=element_blank(),
                axis.title.y=element_blank()) +
        # stat_compare_means(method = "anova", label.y= 0.4) +
        NULL
      })

    } else if(identical(plot_type, FeaturePlot)){
      feature_plots <- plot_type(seu, features = mymarkers, combine = FALSE) %>%
        set_names(mymarkers)
    }


    feature_plots <- map(feature_plots, ~(.x + labs(title = sample_id)))

    return(feature_plots)

  }

  sample_ids <- str_extract(seu_paths, "SRR[0-9]*")

  myplots <- map(seu_paths, plot_markers_in_sample, mymarkers = mymarkers, plot_type = plot_type, group_by = group_by, cluster_dictionary) %>%
    set_names(sample_ids)


  myplots0 <-
    myplots %>%
    transpose() %>%
    imap(~{
      patchwork::wrap_plots(.x) +
        plot_annotation(
          title = .y
        )

    })

  if(identical(plot_type, VlnPlot)){
    plot_type_label = "VlnPlot"
  } else if(identical(plot_type, FeaturePlot)){
    plot_type_label = "FeaturePlot"
  }


  plot_paths <- glue("results/numbat_sridhar/gene_plots/{names(myplots0)}_{plot_type_label}_{group_by}_{extension}.pdf")

  map2(plot_paths, myplots0, ~ggsave(.x, .y, height = 8, width = 14))

  return(plot_paths)

}

find_diffex_bw_clusters_for_each_clone <- function(numbat_rds_file, cluster_dictionary, ident.1 = "G2M", ident.2 = "cone"){
  # browser()
  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")


  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[,!is.na(seu$clone_opt)]

  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")

  seu <- AddMetaData(seu, test0)

  cluster_diff_per_clone <- function(clone_for_diffex, seu){
    # browser()
    seu0 <- seu[,seu$clone_opt == clone_for_diffex]

    Idents(seu0) <- seu0$abbreviation

    diffex <- FindAllMarkers(seu0, group.by = "abbreviation")

  }

  myclones <- sort(unique(seu$clone_opt)) %>%
    set_names(.)

  possible_cluster_diff_per_clone <- possibly(cluster_diff_per_clone)

  diffex <- map(myclones, possible_cluster_diff_per_clone, seu)

  diffex0 <- map(diffex, compact) %>%
    map(tibble::rownames_to_column, "symbol") %>%
    dplyr::bind_rows(.id = "clone") %>%
    dplyr::select(-symbol) %>%
    dplyr::rename(symbol = gene) %>%
    dplyr::mutate(sample_id = sample_id) %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::select(symbol, description, everything()) %>%
    dplyr::distinct(symbol, clone, .keep_all = TRUE) %>%
    dplyr::arrange(cluster, clone, p_val_adj) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    identity()

  diffex_path <- glue("results/{numbat_dir}/{sample_id}_clone_diffex.csv")
  write_csv(diffex0, diffex_path)


  compare_diffex_cluster_by_clone <- function(diffex_bw_clusters_for_each_clone){
    test0 <-
      diffex_bw_clusters_for_each_clone %>%
      dplyr::group_by(clone, cluster) %>%
      dplyr::slice_head() %>%
      dplyr::arrange(cluster, clone)

  }

  diffex1 <- compare_diffex_cluster_by_clone(diffex0)

  cluster_clone <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, clone_opt, abbreviation) %>%
    tidyr::unite(cluster_clone, abbreviation, clone_opt) %>%
    dplyr::select(cell, cluster_clone) %>%
    tibble::column_to_rownames("cell") %>%
    identity()

  seu <- AddMetaData(seu, cluster_clone)

  diffex_cluster_by_clone_dotplot <-
    DotPlot(seu, features = unique(diffex1$symbol), group.by = "cluster_clone") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    coord_flip() +
    labs(title = sample_id)

  dotplot_path <- glue("results/{numbat_dir}/{sample_id}_cluster_clone_diffex_dotplot.pdf")
  ggsave(dotplot_path, diffex_cluster_by_clone_dotplot, width = 10, height = 12)

  trend_genes <- diffex1 %>%
    slice_head(n = 1)

  diffex_cluster_by_clone_trendplot <-
    plot_gene_clone_trend(seu, trend_genes$symbol) +
    labs(title = sample_id)

  trendplot_path <- glue("results/{numbat_dir}/{sample_id}_gene_trend.pdf")
  ggsave(trendplot_path, diffex_cluster_by_clone_trendplot, width = 10, height = 30)

  return(list("diffex" = diffex_path, "plot" = dotplot_path))


}

find_diffex_bw_divergent_clusters <- function(sample_id, tumor_id, seu, mynb, to_SCT_snn_res. = "SCT_snn_res.1", to_clust = c("1", "10"), clone_comparisons){
	
	possible_make_cluster_comparison <- possibly(make_cluster_comparison)
	
	clone_comparisons <- clone_comparisons[[tumor_id]] %>% 
		tibble::enframe("comparison", "segment") %>% 
		dplyr::mutate(clones = str_extract(comparison, "[0-9]_v_[0-9]")) %>%
		dplyr::mutate(clones = str_split(clones, "_v_")) %>%
		dplyr::rowwise() %>% 
		dplyr::filter(all(clones %in% seu$clone_opt)) %>%
		dplyr::select(comparison, segment) %>% 
		tibble::deframe() %>% 
		identity()
		
	message(clone_comparisons)
	
	diffex <- imap(clone_comparisons, make_cluster_comparison, seu, mynb, to_SCT_snn_res., to_clust) %>% 
		identity()
	
	return(diffex)
}

#' find_diffex_bw_clones_for_each_cluster
#'
#' @param numbat_rds_file
#' @param clone_comparisons
#' @param cluster_dictionary annotated louvain clusters
#' @param location in_segment or out_of_segment
#'
#' @return
#' @export
#'
#' @examples
find_diffex_bw_clones_for_each_cluster <- function(seu_path, numbat_rds_files, large_clone_comparisons, cluster_dictionary, location = "in_segment"){
	# browser()
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	numbat_rds_files <- numbat_rds_files %>% 
		set_names(str_extract(., "SRR[0-9]*"))
	
	mynb <- readRDS(numbat_rds_files[[tumor_id]])
	
	seu <- readRDS(seu_path)
	
	seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))
	
	clone_diff_per_cluster <- function(cluster_for_diffex, seu, group.by = "SCT_snn_res.0.6"){
		# browser()
		seu <- seu[,seu[[]][["clusters"]] == cluster_for_diffex]
		print_table_tally(seu$clone_opt)
		
		Idents(seu) <- seu$clone_opt
		
		# diffex <- FindAllMarkers(seu) %>%
		#   dplyr::rename(clone_opt = cluster)
		
		diffex <- imap(large_clone_comparisons[[sample_id]], make_clone_comparison, seu, mynb, location = location)
		
		return(diffex)
		
	}
	
	myclusters <- sort(unique(seu@meta.data[["clusters"]])) %>%
		set_names(.)
	
	possible_clone_diff_per_cluster <- possibly(clone_diff_per_cluster)
	
	diffex <- map(myclusters, possible_clone_diff_per_cluster, seu) %>%
		compact() %>%
		map(bind_rows, .id = "clone_comparison") %>%
		bind_rows(.id = "cluster") %>%
		# dplyr::arrange(cluster, p_val_adj) %>%
		identity()
	
	diffex_path <- glue("results/{sample_id}_cluster_clone_comparison_diffex_{location}.csv")
	write_csv(diffex, diffex_path)
	
	return(diffex_path)
	
}

find_diffex_bw_clones_for_each_cluster_integrated <- function(seu_path, kept_samples = c('SRR14800534', 'SRR14800535', 'SRR14800536'), clone_comparisons = list("2_v_1_16q-" = c("16c", "16b"), "3_v_2_1q+" = c("1b")), location = "in_segment"){
  # browser()

  seu <- readRDS(seu_path)

  numbat_dir = "numbat_sridhar"

  mynbs <- glue("output/{numbat_dir}/{kept_samples}_numbat.rds") %>%
    map(readRDS)

  location  = "in_segment"

  seu <- seu[,!is.na(seu$clone_opt)]

  clone_diff_per_cluster <- function(cluster_for_diffex, seu, group.by){
    # browser()
    seu0 <- seu[,seu[["clusters"]] == cluster_for_diffex]

    Idents(seu0) <- seu0$clone_opt

    # diffex <- FindAllMarkers(seu0) %>%
    #   dplyr::rename(clone_opt = cluster)

    diffex <- imap(clone_comparisons, make_clone_comparison_integrated, seu0, mynbs, location = location)

    return(diffex)

  }

  myclusters <- sort(unique(seu$gene_snn_res.0.2)) %>%
    set_names(.)

  possible_clone_diff_per_cluster <- possibly(clone_diff_per_cluster)

  diffex <- map(myclusters, possible_clone_diff_per_cluster, seu) %>%
    compact() %>%
    map(bind_rows, .id = "clone_comparison") %>%
    bind_rows(.id = "cluster") %>%
    # dplyr::arrange(cluster, p_val_adj) %>%
    identity()

  kept_samples_slug = paste(kept_samples, collapse = "_")

  diffex_path <- glue("{numbat_dir}/{kept_samples_slug}_cluster_clone_comparison_diffex_{location}.csv")
  write_csv(diffex, diffex_path)

  return(diffex_path)


}

gse_plot_from_cluster_diffex <- function(diffex_path){
  browser()

  sample_id <- str_extract(diffex_path, "SRR[0-9]*")

  numbat_dir <- path_split(diffex_path)[[1]][[2]]

  location <- str_extract(diffex_path, "(?<=diffex_).*_segment")

  diffex <-
    diffex_path %>%
    read_csv() %>%
    split(.$clone_comparison) %>%
    map(~split(.x, .x[["cluster"]])) %>%
    identity()

  add_cluster_label <- function(plot_list, mylabel){
    map(plot_list, ~{.x + labs(subtitle = mylabel)})
  }

  annotable_cols <- colnames(annotables::grch38)
  annotable_cols <- annotable_cols[!annotable_cols == "symbol"]

  make_gse_plot <- function(diffex_list, sample_id){
    # browser()
    diffex_list <-
      diffex_list %>%
      map(dplyr::select, -annotable_cols) %>%
      map(dplyr::distinct, symbol, .keep_all = TRUE) %>%
      map(tibble::column_to_rownames, "symbol") %>%
      # map(dplyr::filter, p_val_adj < 0.05) %>%
      identity()

    diffex_list <- diffex_list[map_int(diffex_list, ~dim(.x)[1]) > 0]

    gse_plots <-
      diffex_list %>%
      map(enrichment_analysis) %>%
      imap(~{.x + labs(title = sample_id, subtitle = .y)})

    return(gse_plots)
  }

  gse_plots <-
    map(diffex, make_gse_plot, sample_id) %>%
    purrr::compact()

  drop_empty_plots <- function(plot_list){
    # browser()
    plot_content <-
      plot_list %>%
      map(~{dim(.x[["data"]])}) %>%
      map_lgl(is.null) %>%
      identity()

    plot_list <- plot_list[!plot_content]
  }

  gse_plots <-
    gse_plots %>%
    map(drop_empty_plots) %>%
    purrr::compact()

  gse_plot_path <- glue("results/{numbat_dir}/{sample_id}_cluster_clone_comparison_diffex_{location}.pdf")

  if(length(gse_plots) > 0){
    pdf(gse_plot_path)
    gse_plots
    dev.off()
  }

  return(gse_plots)

}

gse_plot_from_clone_diffex <- function(diffex_path){
  browser()

  sample_id <- str_extract(diffex_path, "SRR[0-9]*")

  numbat_dir <- path_split(diffex_path)[[1]][[2]]

  location <- str_extract(diffex_path, "(?<=diffex_).*_segment")

  diffex <-
    diffex_path %>%
    read_csv() %>%
    split(.$clone_comparison) %>%
    identity()

  annotable_cols <- colnames(annotables::grch38)
  annotable_cols <- annotable_cols[!annotable_cols == "symbol"]

  make_gse_plot <- function(diffex_list, sample_id){
    browser()
    diffex_list <-
      diffex_list %>%
      map(dplyr::select, -any_of(annotable_cols)) %>%
      map(dplyr::distinct, symbol, .keep_all = TRUE) %>%
      map(tibble::column_to_rownames, "symbol") %>%
      map(dplyr::filter, p_val_adj < 0.05)

    diffex_list <- diffex_list[map_int(diffex_list, ~dim(.x)[1]) > 0]

    gse_plots <-
      diffex_list %>%
      map(enrichment_analysis) %>%
      imap(~{.x + labs(title = sample_id, subtitle = .y)})

    return(gse_plots)
  }

  gse_plots <-
    make_gse_plot(diffex, sample_id) %>%
    purrr::compact()

  drop_empty_plots <- function(plot_list){
    # browser()
    plot_content <-
      plot_list %>%
      map(~{dim(.x[["data"]])}) %>%
      map_lgl(is.null) %>%
      identity()

    plot_list <- plot_list[!plot_content]
  }

  gse_plots <-
    gse_plots %>%
    map(drop_empty_plots) %>%
    purrr::compact()

  gse_plot_path <- glue("results/{numbat_dir}/{sample_id}_clone_comparison_diffex_{location}.pdf")

  pdf(gse_plot_path)
  gse_plots
  dev.off()

  return(gse_plots)

}


plot_feature_in_seu <- function(numbat_rds_file, ...){
  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[,!is.na(seu$clone_opt)]

  Seurat::FeaturePlot(seu, ...)
}


collect_clusters_from_seus <- function(filtered_seus){
  # browser()

  resolutions = glue("SCT_snn_res.{seq(0.2, 2.0, by = 0.2)}")

  gather_clusters <- function(filtered_seu, resolutions){
    # browser()
    # sample_id <- str_extract(filtered_seu, "SRR[0-9]*")

    # numbat_dir = fs::path_split(filtered_seu)[[1]][[2]]

    seu <- readRDS(filtered_seu)

    clusters <- seu@meta.data[resolutions] %>%
      tidyr::pivot_longer(everything(), names_to = "resolution", values_to = "cluster") %>%
      dplyr::group_by(resolution, cluster) %>%
      dplyr::summarize(n_cells = dplyr::n()) %>%
      split(.$resolution) %>%
      map(ungroup) %>%
      imap(~dplyr::select(.x, {{.y}} := cluster, n_cells)) %>%
      identity()

    return(clusters)

  }

  sample_ids = str_extract(filtered_seus, "SRR[0-9]*")

  names(filtered_seus) <- sample_ids

  myclusters <- map(filtered_seus, gather_clusters, resolutions) %>%
    purrr::transpose() %>%
    purrr::map(dplyr::bind_rows, .id = "sample_id") %>%
    identity()

  return(myclusters)

}

read_cluster_dictionary <- function(cluster_dictionary_path = "data/cluster_dictionary.csv"){
  cluster_dictionary <- read_tsv(cluster_dictionary_path) %>%
    split(.$sample_id)

  return(cluster_dictionary)
}

read_cells_to_remove <- function(cell_remove_file = "data/cells_to_remove_final.xlsx"){

  mysheets <- excel_sheets(cell_remove_file) %>%
    set_names(.)

  cells_to_remove <-
    mysheets %>%
    map(~readxl::read_xlsx(cell_remove_file, .x))

  return(cells_to_remove)
}

make_pdf_montages <- function(plot_files, heatmaps, tile = '6'){

  numbat_dir <-
    path_split(plot_files[[1]][[1]])[[1]][[2]] %>%
    identity()

  plot_files <- map2(plot_files, heatmaps, c)

  # plot_files <- map2(plot_files, expression_files, c)

  sample_ids <-
    plot_files %>%
    map(1) %>%
    map(fs::path_split) %>%
    map(c(1,3))

  names(plot_files) <- sample_ids

  montage_paths <- imap(plot_files, montage_images, tile = tile)

  return(montage_paths)

}

make_expression_heatmap_comparison <- function(large_numbat_pdfs, heatmaps){
  browser()

  numbat_dir <-
    path_split(large_numbat_pdfs[[1]][[1]])[[1]][[2]] %>%
    identity()

  expression_files <- map(large_numbat_pdfs, 6)

  heatmaps <- map(heatmaps, 1)

  plot_files <- map2(expression_files, heatmaps, c)

  sample_ids <-
    plot_files %>%
    map(1) %>%
    map(fs::path_split) %>%
    map(c(1,3))

  names(plot_files) <- sample_ids

  montage_images <- function(plot_files, sample_id){
    # browser()
    plot_images <- magick::image_read(plot_files, density = 600)

    my_montage <- magick::image_montage(plot_images, tile = '2', geometry='800x', shadow = FALSE)

    # my_montage <- image_montage(plot_images, geometry = c('x200+10+10', 'x800+10+10', 'x100+10+10'), tile = '3x', shadow = FALSE)

    montage_path <- glue("results/{numbat_dir}/{sample_id}_exp_v_heatmap.pdf")

    image_write(my_montage, format = "pdf", montage_path)

    return(montage_path)
  }

  montage_paths <- imap(plot_files, montage_images)

  return(montage_paths)

}

browse_celltype_expression <- function(sridhar_seu, symbol){

  pdf(glue("~/tmp/{symbol}.pdf"))
  FeaturePlot(sridhar_seu, features = c(glue("{symbol}")), split.by = "CellType_predict", combine = FALSE)
  dev.off()
  browseURL(glue("~/tmp/{symbol}.pdf"))

  return(glue("~/tmp/{symbol}.pdf"))

}

collect_markers <- function(numbat_rds_file, metavar = "gene_snn_res.0.2", num_markers = 5, selected_values = NULL, return_plotly = TRUE, marker_method = "presto", seurat_assay = "gene", hide_technical = "all", unique_markers = FALSE, p_val_cutoff = 1, ...) {

  # # by default only resolution markers are calculated in pre-processing
  # seu <- find_all_markers(numbat_rds_file, metavar, seurat_assay = seurat_assay, p_val_cutoff = p_val_cutoff)

  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")
  
  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]
  
  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()
  
  marker_table <- seu@misc$markers[[metavar]][[marker_method]]
  
  markers <-
    marker_table %>%
    seuratTools:::enframe_markers() %>%
    dplyr::mutate(dplyr::across(.fns = as.character))
  
  if (!is.null(hide_technical)) {
    markers <- purrr::map(markers, c)

    if (hide_technical == "pseudo") {
      markers <- purrr::map(markers, ~ .x[!.x %in% pseudogenes[[seurat_assay]]])
    } else if (hide_technical == "mito_ribo") {
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^MT-")])
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^RPS")])
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^RPL")])
    } else if (hide_technical == "all") {
      markers <- purrr::map(markers, ~ .x[!.x %in% pseudogenes[[seurat_assay]]])
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^MT-")])
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^RPS")])
      markers <- purrr::map(markers, ~ .x[!stringr::str_detect(.x, "^RPL")])
    }

    min_length <- min(purrr::map_int(markers, length))

    markers <- purrr::map(markers, head, min_length) %>%
      dplyr::bind_cols()
  }

  colnames(markers) <- glue("{colnames(markers)}_{sample_id}")

  return(markers)
}

collect_all_markers <- function(numbat_rds_files, excel_output = "results/markers.xlsx"){



  names(numbat_rds_files) <- str_extract(numbat_rds_files, "SRR[0-9]*")

  marker_tables <- purrr::map(numbat_rds_files, collect_markers)

  writexl::write_xlsx(marker_tables, excel_output)

}

pull_cells_matching_expression <- function(myexpression, joint_post){
  # browser()
  excluded_cells <-
    joint_post %>%
    dplyr::filter(!!parse_expr(myexpression)) %>%
    dplyr::pull(cell) %>%
    identity()

  return(excluded_cells)

}


# merged_metadata_transfer <-
#   merged_metadata %>%
#   dplyr::filter(sample_id == {{sample_id}}) %>%
#   dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
#   # dplyr::filter(cell %in% colnames(seu)) %>%
#   tibble::column_to_rownames("cell") %>%
#   identity()
#
# seu <- Seurat::AddMetaData(seu, merged_metadata_transfer)
#
# plot_markers(seu, metavar = "merged_leiden", marker_method = "presto", return_plotly = FALSE) +
#   labs(title = sample_id)
# ggsave(glue("results/{numbat_dir}/{sample_id}_merged_marker.png"))


score_filtration <- function(numbat_rds_file, cluster_dictionary, filter_expressions = NULL){
  # browser()

  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  dir_create(glue("results/{numbat_dir}"))
  dir_create(glue("results/{numbat_dir}/{sample_id}"))

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  mynb <- readRDS(numbat_rds_file)

  nb_meta <- mynb[["clone_post"]][,c("cell", "clone_opt", "GT_opt")] %>%
    dplyr::mutate(cell = str_replace(cell, "\\.", "-")) %>%
    tibble::column_to_rownames("cell")

  seu <- Seurat::AddMetaData(seu, nb_meta)

  seu <- seu[,!is.na(seu$clone_opt)]

  test0 <- seu@meta.data["gene_snn_res.0.2"] %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(gene_snn_res.0.2 = as.numeric(gene_snn_res.0.2)) %>%
    dplyr::left_join(cluster_dictionary[[sample_id]], by = "gene_snn_res.0.2") %>%
    dplyr::select("cell", "abbreviation") %>%
    tibble::column_to_rownames("cell")

  seu <- AddMetaData(seu, test0)

  phylo_heatmap_data <- mynb$clone_post %>%
    dplyr::select(cell, clone_opt) %>%
    dplyr::left_join(mynb$joint_post, by = "cell")

  # filter out cells
  excluded_cells <- map(filter_expressions[[sample_id]], pull_cells_matching_expression, phylo_heatmap_data) %>%
    unlist()

  seu_filtered <- seu[,!colnames(seu) %in% excluded_cells]

  dist_list <- list(
    "unfiltered" = score_pca(seu),
    "filtered" = score_pca(seu_filtered)
  )

  return(dist_list)
}

score_pca <- function(seu){

  mygroups <- seu$clone_opt %>%
    tibble::enframe("cell", "group")

  mypca <- seu@reductions$pca@cell.embeddings %>%
    as.data.frame() %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(-c("cell"), names_to = "PC", values_to = "value") %>%
    dplyr::left_join(mygroups, by = "cell")

  test0 <-
    mypca %>%
    dplyr::group_by(PC, group) %>%
    dplyr::summarize(value = mean(value)) %>%
    tidyr::pivot_wider(names_from = "PC", values_from = "value") %>%
    tibble::column_to_rownames("group") %>%
    as.matrix() %>%
    t() %>%
    cor()  %>%
    # dist() %>%
    identity()

  test1 <- dist(1-test0)

  return(test1)

  # # aggregate values within categories using 'mean'
  # mean_df = rep_df.groupby(level=0).mean()
  #
  # import scipy.cluster.hierarchy as sch
  # from scipy.spatial import distance
  #
  # corr_matrix = mean_df.T.corr(method=cor_method)
  # corr_condensed = distance.squareform(1 - corr_matrix)
  # z_var = sch.linkage(
  #   corr_condensed, method=linkage_method, optimal_ordering=optimal_ordering
  # )
  # dendro_info = sch.dendrogram(z_var, labels=list(categories), no_plot=True)
}



filter_sample_qc <- function(seu, mito_threshold = 5, nCount_threshold = 500, nFeature_threshold = 500){
  # browser()
  seu <-
    seu %>%
    subset(subset = percent.mt < mito_threshold) %>%
    subset(subset = nFeature_gene > nFeature_threshold) %>%
    subset(subset = nCount_gene > nCount_threshold) %>%
    identity()

}


plot_gene_clone_trend <- function(seu, mygenes = c('CRABP2', 'MEG3')){
  gene_df <-
    FetchData(seu, vars = c(mygenes, "cluster_clone")) %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(-c("cell", "cluster_clone"), names_to = "gene", values_to = "counts")

  ggplot(gene_df, aes(cluster_clone, counts)) +
    geom_jitter(width = 0.1) +
    facet_wrap(~gene, ncol = 1) +
    NULL

}

#' Title
#'
#' @param numbat_rds_file
#' @param clone_comparisons
#' @param cluster_dictionary annotated louvain clusters
#' @param location in_segment or out_of_segment
#'
#' @return
#' @export
#'
#' @examples
find_diffex_clones <- function(seu_path, numbat_rds_files, large_clone_comparisons, location = "in_segment"){
  # browser()
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	numbat_rds_files <- numbat_rds_files %>% 
		set_names(str_extract(., "SRR[0-9]*"))
	
	mynb <- readRDS(numbat_rds_files[[tumor_id]])

  seu <- readRDS(glue("output/seurat/{sample_id}_filtered_seu.rds"))

  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  possible_clone_comparison <- possibly(make_clone_comparison)

  diffex <- imap(large_clone_comparisons[[sample_id]], possible_clone_comparison, seu, mynb, location = location) %>%
    purrr::compact()

  return(diffex)

}

#' Title
#'
#' @param numbat_rds_file
#' @param clone_comparisons
#' @param cluster_dictionary annotated louvain clusters
#' @param location in_segment or out_of_segment
#'
#' @return
#' @export
#'
#' @examples
find_diffex_clones_in_phase <- function(seu_path, phase = "g1", scna_of_interest = "2p", numbat_rds_files, large_clone_comparisons, location = "in_segment"){
	# browser()
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	numbat_rds_files <- numbat_rds_files %>% 
		set_names(str_extract(., "SRR[0-9]*"))
	
	mynb <- readRDS(numbat_rds_files[[tumor_id]])
	
	seu <- readRDS(glue("output/seurat/{sample_id}_filtered_seu.rds"))
	
	seu <- seu[,seu$phase_level %in% phase]
	
	seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))
	
	possible_clone_comparison <- possibly(make_clone_comparison)
	
	large_clone_comparisons <- large_clone_comparisons[[sample_id]]
	
	large_clone_comparisons <- large_clone_comparisons[str_detect(names(large_clone_comparisons), scna_of_interest)]
	
	diffex <- imap(large_clone_comparisons, possible_clone_comparison, seu, mynb, location = location) %>%
		purrr::compact()
	
	# browser()
	
	diffex <- diffex[str_detect(names(diffex), scna_of_interest)]
	
	enrich_diffex <- function(df){
		df |>
			tibble::column_to_rownames("symbol") |>
			dplyr::select(-any_of(colnames(annotables::grch38))) |>
			enrichment_analysis(gene_set = "hallmark")
	}
	
	enrichments <- diffex |> 
		purrr::list_flatten() |> 
		map(enrich_diffex)

	return(diffex)
	
}

#' Title
#'
#' @param numbat_rds_file
#' @param clone_comparisons
#' @param cluster_dictionary annotated louvain clusters
#' @param location in_segment or out_of_segment
#'
#' @return
#' @export
#'
#' @examples
compare_enrichment <- function(diffex, ...
){
	# browser()
	enrich_diffex <- function(df){
		df |>
			tibble::column_to_rownames("symbol") |>
			dplyr::select(-any_of(colnames(annotables::grch38))) |>
			enrichment_analysis(...)
	}
	
	enrichments <- diffex |> 
		purrr::list_flatten() |> 
		map(enrich_diffex)
	
	return(enrichments)
	
}

#' Title
#'
#' @param numbat_rds_file
#' @param clone_comparisons
#' @param cluster_dictionary annotated louvain clusters
#' @param location in_segment or out_of_segment
#'
#' @return
#' @export
#'
#' @examples
plot_several_diffex_clones_in_phase <- function(enrichments, scna_of_interest = "2p", phase = "g1"){
	
	enrichment_plot <-
		enrichments |>
		purrr::list_flatten() |> 
		clusterProfiler::merge_result() |> 
		plot_enrichment(p_val_cutoff = 1, result_slot = "compareClusterResult")

	pdf_path = ggsave(glue("results/{scna_of_interest}_{phase}_enrichment.pdf"), h = 10, w = 8)
	
}

find_diffex_clones_integrated <- function(seu_path, kept_samples, clone_comparisons, location = "in_segment"){
  # browser()

  seu <- readRDS(seu_path)

  seu <- seu[,!is.na(seu$clone_opt)]

  mynbs <- glue("output/numbat_sridhar/{kept_samples}_numbat.rds") %>%
    map(readRDS)


  diffex <- imap(clone_comparisons, make_clone_comparison_integrated, seu, mynbs, location = location)


  return(diffex)

}


tabulate_diffex_clones <- function(cluster_diffex_clones,
                                   cluster_xlsx = "results/diffex_bw_clones_per_cluster_large.xlsx",
                                   cluster_by_chr_xlsx = "results/diffex_bw_clones_per_cluster_large_by_chr.xlsx",
                                   total_diffex_clones,
                                   total_xlsx = "results/straight_diffex_bw_clones_large.xlsx",
                                   total_by_chr_xlsx = "results/straight_diffex_bw_clones_large_by_chr.xlsx"){

  kooi_candidates <- read_csv("data/kooi_candidates.csv")

  cc_genes <- Seurat::cc.genes %>%
    tibble::enframe("phase_of_gene", "symbol") %>%
    tidyr::unnest(symbol)

  sample_ids = str_extract(cluster_diffex_clones, "SRR[0-9]*")

  annotate_percent_segment_diffex <- function(diffex){

    if("genes_in_segment" %in% colnames(diffex)){
      diffex <-
        diffex %>%
        dplyr::group_by(seg, genes_in_segment) %>%
        dplyr::select(symbol, genes_in_segment) %>%
        dplyr::summarize(diffex_genes = list(symbol)) %>%
        dplyr::mutate(genes_in_segment = str_split(genes_in_segment, ", ")) %>%
        mutate(diffex_genes_in_segment=map2(diffex_genes, genes_in_segment, safely(base::intersect))) %>%
        mutate(diffex_genes_in_segment=map(diffex_genes_in_segment, "result")) %>%
        dplyr::mutate(ratio_genes_diffex_in_segment = map2_dbl(diffex_genes_in_segment, genes_in_segment, ~(length(.x)/length(.y)))) %>%
        dplyr::select(seg, ratio_genes_diffex_in_segment) %>%
        dplyr::left_join(diffex, by = "seg") %>%
        dplyr::select(-c("genes_in_segment")) %>%
        identity()
    }

    return(diffex)
  }

  total_diffex_clones <-
    total_diffex_clones %>%
    set_names(sample_ids) %>%
    map(bind_rows, .id = "clone_comparison") %>%
    map(dplyr::left_join, cc_genes, by = "symbol") %>%
    map(group_by, clone_comparison) %>%
    map(dplyr::filter, p_val_adj < 1) %>%
    purrr::keep(~nrow(.x) > 0) %>%
    map(dplyr::arrange, clone_comparison, p_val_adj) %>%
    map(dplyr::select, clone_comparison, chr, symbol, description, everything()) %>%
    map(dplyr::left_join, kooi_candidates, by = "symbol") %>%
    map(annotate_percent_segment_diffex) %>%
    identity()

  # cluster ------------------------------
  cluster_diffex_clones <-
    cluster_diffex_clones %>%
    set_names(str_extract(., "SRR[0-9]*")) %>%
    map(read_csv) %>%
    purrr::keep(~nrow(.x) > 0) %>%
    map(dplyr::filter, p_val_adj < 1) %>%
    purrr::keep(~nrow(.x) > 0) %>%
    map(dplyr::left_join, cc_genes, by = "symbol") %>%
    map(dplyr::arrange, clone_comparison, cluster, p_val_adj) %>%
    map(dplyr::select, clone_comparison, cluster, chr, symbol, description, everything()) %>%
    map(dplyr::left_join, kooi_candidates, by = "symbol") %>%
    map(annotate_percent_segment_diffex) %>%
    identity()

  cluster_diffex_clones <-
  cluster_diffex_clones %>%
    # map2(total_diffex_clones, annotate_cluster_membership, "is_total_clone_diffex") %>%
    map(dplyr::distinct, clone_comparison, cluster, chr, symbol, .keep_all = TRUE)

  write_xlsx(cluster_diffex_clones, cluster_xlsx)

  cluster_diffex_clones_by_chr <-
    cluster_diffex_clones %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::distinct(sample_id, clone_comparison, cluster, chr, symbol, .keep_all = TRUE) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(num_samples = length(unique(sample_id))) %>%
    dplyr::arrange(desc(num_samples), symbol) %>%
    dplyr::filter(num_samples > 1) %>%
    dplyr::filter(!str_detect(chr, "CHR_")) %>%
    split(.$chr)

  write_xlsx(cluster_diffex_clones_by_chr, cluster_by_chr_xlsx)

  # total ------------------------------

  total_diffex_clones <-
    total_diffex_clones %>%
    # map2(total_diffex_clones, cluster_diffex_clones, annotate_cluster_membership, "is_cluster_clone_diffex") %>%
    map(dplyr::distinct, clone_comparison, chr, symbol, .keep_all = TRUE)

  write_xlsx(total_diffex_clones, total_xlsx)

  total_diffex_clones_by_chr <-
    total_diffex_clones %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::distinct(sample_id, clone_comparison, chr, symbol, .keep_all = TRUE) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(num_samples = length(unique(sample_id))) %>%
    dplyr::arrange(desc(num_samples), symbol) %>%
    dplyr::filter(num_samples > 1) %>%
    dplyr::filter(!str_detect(chr, "CHR_")) %>%
    split(.$chr)

  write_xlsx(total_diffex_clones_by_chr, total_by_chr_xlsx)

  return(list("total" = c(total_xlsx, total_by_chr_xlsx), "cluster" = c(cluster_xlsx, cluster_by_chr_xlsx)))

}

make_volcano_diffex_clones <- function(cluster_diffex_clones,
                                   cluster_pdf = "results/diffex_bw_clones_per_cluster_large.pdf",
                                   total_diffex_clones,
                                   total_pdf = "results/straight_diffex_bw_clones_large.pdf"
                                   ){
  # browser()

  kooi_candidates <- read_csv("data/kooi_candidates.csv")

  cc_genes <- Seurat::cc.genes %>%
    tibble::enframe("phase_of_gene", "symbol") %>%
    tidyr::unnest(symbol)

  sample_ids = str_extract(cluster_diffex_clones, "SRR[0-9]*")

  # total ------------------------------
  total_diffex_clones <-
    total_diffex_clones %>%
    set_names(sample_ids) %>%
    map(bind_rows, .id = "clone_comparison") %>%
    purrr::discard(~nrow(.x) < 1) %>%
    map(dplyr::left_join, cc_genes, by = "symbol") %>%
    map(group_by, clone_comparison) %>%
    # map(dplyr::filter, p_val_adj < 0.05) %>%
    map(dplyr::arrange, clone_comparison, p_val_adj) %>%
    map(dplyr::select, clone_comparison, chr, symbol, description, everything()) %>%
    map(dplyr::left_join, kooi_candidates, by = "symbol") %>%
    identity()

  # cluster------------------------------
  cluster_diffex_clones <-
    cluster_diffex_clones %>%
    set_names(str_extract(., "SRR[0-9]*")) %>%
    map(read_csv) %>%
    purrr::keep(~nrow(.x) > 0) %>%
    map(dplyr::left_join, cc_genes, by = "symbol") %>%
    # map(dplyr::filter, p_val_adj < 0.05) %>%
    map(dplyr::arrange, clone_comparison, cluster, p_val_adj) %>%
    map(dplyr::select, clone_comparison, cluster, chr, symbol, description, everything()) %>%
    map(dplyr::left_join, kooi_candidates, by = "symbol") %>%
    identity()

  cluster_diffex_clones <-
  cluster_diffex_clones %>%
    # map2(total_diffex_clones, annotate_cluster_membership, "is_total_clone_diffex") %>%
    map(dplyr::distinct, clone_comparison, cluster, chr, symbol, .keep_all = TRUE)

  total_diffex_clones <-
  total_diffex_clones %>%
    # map2(cluster_diffex_clones, annotate_cluster_membership, "is_cluster_clone_diffex") %>%
    map(dplyr::distinct, clone_comparison, chr, symbol, .keep_all = TRUE)

  # cluster ------------------------------
  volcano_plot_clone_clusters <- function(clone_diffex, sample_id){
  	# browser()
    myres <- clone_diffex %>%
      group_by(clone_comparison, cluster)

    myres <-
      myres %>%
      group_split() %>%
      map(tibble::column_to_rownames, "symbol") %>%
      identity()

    names(myres) <- myres %>%
      map_chr(~glue("{unique(.x$clone_comparison)} {unique(.x$cluster)}"))

    volcano_plots <- myres %>% 
    	map(dplyr::mutate, diffex_comparison = str_replace(str_extract(clone_comparison, "[0-9]_v_[0-9]"), "_v_", "_")) %>% 
    	imap(make_volcano_plots, sample_id)

    return(volcano_plots)
  }

  clone_cluster_comparison_volcanos <- imap(cluster_diffex_clones, volcano_plot_clone_clusters)

  pdf(cluster_pdf)
  print(clone_cluster_comparison_volcanos)
  dev.off()

  # total ------------------------------
  volcano_plot_clones <- function(clone_diffex, sample_id){
    myres <- clone_diffex %>%
      group_by(clone_comparison)

    myres <-
      myres %>%
      group_split() %>%
      map(tibble::column_to_rownames, "symbol") %>%
      identity()

    names(myres) <-
      myres %>%
      map_chr(~glue("{unique(.x$clone_comparison)}"))

    volcano_plots <- 
    	myres %>%
    	map(dplyr::mutate, diffex_comparison = str_replace(str_extract(clone_comparison, "[0-9]_v_[0-9]"), "_v_", "_")) %>% 
    	imap(make_volcano_plots, sample_id)

    return(volcano_plots)
  }

  clone_comparison_volcanos <- imap(total_diffex_clones, volcano_plot_clones)

  pdf(total_pdf)
  print(clone_comparison_volcanos)
  dev.off()

  return(list("cluster" = cluster_pdf, "total" = total_pdf))

}

compare_per_cluster_and_total_clone_diffex <- function(per_cluster_clone_diffex, total_clone_diffex){
  per_cluster_clone_diffex <-
    per_cluster_clone_diffex %>%
    readxl()

  total_clone_diffex <-
    total_clone_diffex %>%
    readxl()
}


#' Make a differential expression comparison between numbat clones either cis, trans, or all
#'
#' @param mysegs
#' @param comparison
#' @param seu
#' @param mynb
#' @param location
#'
#' @return
#' @export
#'
#' @examples
make_clone_comparison <- function(mysegs, comparison, seu, mynb, location = "in_segment"){
  # browser()

  idents <-
    comparison %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split(pattern = "_v_") %>%
    unlist()

  segments <- mynb$clone_post %>%
    dplyr::left_join(mynb$joint_post, by = "cell") %>%
    dplyr::filter(clone_opt %in% idents) %>%
    dplyr::filter(seg %in% mysegs) %>%
    dplyr::distinct(CHROM, seg, seg_start, seg_end, cnv_state_map) %>%
    dplyr::mutate(seqnames = CHROM, start = seg_start, end = seg_end) %>%
    dplyr::filter(!cnv_state_map == "neu") %>%
    plyranges::as_granges() %>%
    identity()

  if(location=="in_segment"){
    diffex <- FindMarkers(seu, ident.1 = idents[[1]], ident.2 = idents[[2]], group.by = "clone_opt",  logfc.threshold = 0.1, test.use = "MAST", assay = "gene") %>%
      tibble::rownames_to_column("symbol") %>%
      dplyr::left_join(annotables::grch38, by = "symbol") %>%
      dplyr::distinct(ensgene, .keep_all = TRUE) %>%
      dplyr::mutate(seqnames = chr) %>%
      dplyr::filter(!is.na(start), !is.na(end)) %>%
      plyranges::as_granges() %>%
      plyranges::join_overlap_intersect(segments) %>%
      as_tibble() %>%
      dplyr::mutate(log2_sign = dplyr::case_when(cnv_state_map == "amp" ~ -1,
                                                 cnv_state_map == "del" ~ 1)) %>%
      # dplyr::filter(sign(log2_sign) == sign(avg_log2FC)) %>%
      dplyr::select(-c("CHROM", "seg_start",
                       "seg_end", "cnv_state_map", "log2_sign")) %>%
      dplyr::filter(!str_detect(chr, "CHR_")) %>%
      dplyr::distinct(symbol, .keep_all = TRUE)

    genes_in_segments <-
      annotables::grch38 %>%
      dplyr::rename(seqnames = chr) %>%
      plyranges::as_granges() %>%
      plyranges::join_overlap_intersect(segments) %>%
      as_tibble() %>%
      dplyr::select(seg, symbol) %>%
      dplyr::group_by(seg) %>%
      dplyr::summarize(genes_in_segment = paste(symbol, collapse = ", ")) %>%
      identity()

    diffex <-
      diffex %>%
      dplyr::left_join(genes_in_segments, by = "seg")


  } else if(location=="out_of_segment"){
    diffex <- FindMarkers(seu, ident.1 = idents[[1]], ident.2 = idents[[2]], group.by = "clone_opt", logfc.threshold = 0.1, test.use = "MAST", assay = "gene") %>%
      tibble::rownames_to_column("symbol") %>%
      dplyr::left_join(annotables::grch38, by = "symbol") %>%
      dplyr::distinct(ensgene, .keep_all = TRUE) %>%
      dplyr::mutate(seqnames = chr) %>%
      dplyr::filter(!is.na(start), !is.na(end)) %>%
      plyranges::as_granges()

    out_of_segment_ranges <-
      diffex %>%
      plyranges::setdiff_ranges(segments)

    diffex <-
      diffex %>%
      plyranges::join_overlap_intersect(out_of_segment_ranges) %>%
      as_tibble() %>%
      dplyr::filter(!str_detect(chr, "CHR_")) %>%
      dplyr::distinct(symbol, .keep_all = TRUE)

  } else if(location=="all"){
    diffex <- FindMarkers(seu, ident.1 = idents[[1]], ident.2 = idents[[2]], group.by = "clone_opt", logfc.threshold = 0.1, test.use = "MAST", assay = "gene") %>%
      tibble::rownames_to_column("symbol") %>%
      dplyr::left_join(annotables::grch38, by = "symbol") %>%
      dplyr::distinct(ensgene, .keep_all = TRUE) %>%
      dplyr::mutate(seqnames = chr) %>%
      dplyr::filter(!is.na(start), !is.na(end)) %>%
      dplyr::filter(!str_detect(chr, "CHR_")) %>%
      dplyr::distinct(symbol, .keep_all = TRUE)
  }

  diffex <-
    diffex %>%
    append_clone_nums(comparison, seu)

  return(diffex)

}

#' Make a differential expression clone_comparison between numbat clones either cis, trans, or all
#'
#' @param mysegs
#' @param clone_comparison
#' @param seu
#' @param mynb
#' @param location
#'
#' @return
#' @export
#'
#' @examples
make_cluster_comparison <- function(mysegs, clone_comparison, seu, mynb, to_SCT_snn_res. = "SCT_snn_res.1", to_clust = c("1", "10")){
	# browser()
	message(glue("{to_SCT_snn_res.} {paste(to_clust, collapse = '_')} {clone_comparison}"))
	
	idents <-
		clone_comparison %>%
		str_extract("[0-9]_v_[0-9]") %>%
		str_split(pattern = "_v_") %>%
		unlist()
	
	segments <- mynb$clone_post %>%
		dplyr::left_join(mynb$joint_post, by = "cell") %>%
		dplyr::filter(clone_opt %in% idents) %>%
		dplyr::filter(seg %in% mysegs) %>%
		dplyr::distinct(CHROM, seg, seg_start, seg_end, cnv_state_map) %>%
		dplyr::mutate(seqnames = CHROM, start = seg_start, end = seg_end) %>%
		dplyr::filter(!cnv_state_map == "neu") %>%
		plyranges::as_granges() %>%
		identity()
	
	diffex_all <- FindMarkers(seu, ident.1 = to_clust[[1]], ident.2 = to_clust[[2]], group.by = to_SCT_snn_res., logfc.threshold = 0.1, test.use = "MAST", assay = "gene") %>%
		tibble::rownames_to_column("symbol") %>%
		dplyr::left_join(annotables::grch38, by = "symbol") %>%
		dplyr::distinct(ensgene, .keep_all = TRUE) %>%
		dplyr::mutate(seqnames = chr) %>%
		dplyr::filter(!is.na(start), !is.na(end)) %>%
		dplyr::filter(!str_detect(chr, "CHR_")) %>%
		dplyr::distinct(symbol, .keep_all = TRUE) %>% 
		append_clone_nums(clone_comparison, seu) %>% 
		dplyr::mutate(to_SCT_snn_res. := {{to_SCT_snn_res.}},
					  to_clust := paste({{to_clust}}, collapse = "_")) %>% 
		dplyr::select(-strand)
	
	diffex_cis <- 
		diffex_all %>% 
		plyranges::as_granges() %>%
		plyranges::join_overlap_intersect(segments) %>%
		as_tibble() %>%
		dplyr::mutate(log2_sign = dplyr::case_when(cnv_state_map == "amp" ~ -1,
												   cnv_state_map == "del" ~ 1)) %>%
		# dplyr::filter(sign(log2_sign) == sign(avg_log2FC)) %>%
		dplyr::select(-c("CHROM", "seg_start",
						 "seg_end", "cnv_state_map", "log2_sign")) %>%
		dplyr::filter(!str_detect(chr, "CHR_")) %>%
		dplyr::distinct(symbol, .keep_all = TRUE) %>% 
		dplyr::select(-strand)
		
	trans_ranges <-
		diffex_all %>%
		plyranges::as_granges() %>%
		plyranges::setdiff_ranges(segments)
	
	diffex_trans <-
		diffex_all %>%
		plyranges::as_granges() %>%
		plyranges::join_overlap_intersect(trans_ranges) %>%
		as_tibble() %>%
		dplyr::filter(!str_detect(chr, "CHR_")) %>%
		dplyr::distinct(symbol, .keep_all = TRUE) %>% 
		dplyr::select(-strand)
	
	return(list("all" = diffex_all, "cis" = diffex_cis, "trans" = diffex_trans))
	
}

make_clone_comparison_integrated <- function(mysegs, comparison, seu, mynbs, location = "in_segment"){
  # browser()

  idents <-
    comparison %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split(pattern = "_v_") %>%
    unlist()

  pull_segment_from_nb <- function(mynb, idents){
    segments <- mynb$clone_post %>%
      dplyr::left_join(mynb$joint_post, by = "cell") %>%
      dplyr::filter(clone_opt %in% idents) %>%
      dplyr::filter(seg %in% mysegs) %>%
      dplyr::distinct(CHROM, seg, seg_start, seg_end, cnv_state_map) %>%
      dplyr::mutate(seqnames = CHROM, start = seg_start, end = seg_end) %>%
      dplyr::filter(!cnv_state_map == "neu") %>%
      plyranges::as_granges() %>%
      identity()

  }

  segments <- map(mynbs, pull_segment_from_nb, idents)

  segments <- unlist(as(segments, "GRangesList"))

  if(location=="in_segment"){
    diffex <- FindMarkers(seu, ident.1 = idents[[1]], ident.2 = idents[[2]], group.by = "clone_opt",  logfc.threshold = 0.1) %>%
      tibble::rownames_to_column("symbol") %>%
      dplyr::left_join(annotables::grch38, by = "symbol") %>%
      dplyr::distinct(ensgene, .keep_all = TRUE) %>%
      dplyr::mutate(seqnames = chr) %>%
      dplyr::filter(!is.na(start), !is.na(end)) %>%
      plyranges::as_granges() %>%
      plyranges::join_overlap_intersect(segments) %>%
      as_tibble() %>%
      dplyr::mutate(log2_sign = dplyr::case_when(cnv_state_map == "amp" ~ -1,
                                                 cnv_state_map == "del" ~ 1)) %>%
      # dplyr::filter(sign(log2_sign) == sign(avg_log2FC)) %>%
      dplyr::select(-c("CHROM", "seg_start",
                       "seg_end", "cnv_state_map", "log2_sign")) %>%
      dplyr::filter(!str_detect(chr, "CHR_")) %>%
      dplyr::distinct(symbol, .keep_all = TRUE)

  } else if(location=="out_of_segment"){
    diffex <- FindMarkers(seu, ident.1 = idents[[1]], ident.2 = idents[[2]], group.by = "clone_opt", logfc.threshold = 0.1) %>%
      tibble::rownames_to_column("symbol") %>%
      dplyr::left_join(annotables::grch38, by = "symbol") %>%
      dplyr::distinct(ensgene, .keep_all = TRUE) %>%
      dplyr::mutate(seqnames = chr) %>%
      dplyr::filter(!is.na(start), !is.na(end)) %>%
      plyranges::as_granges()

    out_of_segment_ranges <-
      diffex %>%
      plyranges::setdiff_ranges(segments)

    diffex <-
      diffex %>%
      plyranges::join_overlap_intersect(out_of_segment_ranges) %>%
      as_tibble() %>%
      dplyr::filter(!str_detect(chr, "CHR_")) %>%
      dplyr::distinct(symbol, .keep_all = TRUE)

  }

}


convert_volcano_to_plotly <- function(myplot, out_html){
  browser()

  myplotly <- ggplotly(myplot + aes(x= avg_log2FC, y= -log10(p_val_adj), symbol = symbol), tooltip = "symbol")

  saveWidget(myplotly, out_html, selfcontained = T, libdir = "lib")

}



tabulate_clone_comparisons <- function(large_clone_comparisons){

  clone_comparison_table <-
    large_clone_comparisons %>%
    map(~{.x %>%
        tibble::enframe("comparison", "segment") %>%
        tidyr::unnest(segment)
        }) %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    identity()

  return(clone_comparison_table)
}

remove_non_tumor_cells <- function(seu){
  # identify retinal cell type clusters
  # identify cells with no SCNAs

  # filter out cells that are
  # 1) not cones or RB cells
  # 2) have no SCNAs

}

seu_integrate_rbs <- function(numbat_dir = "output/numbat_sridhar", kept_samples = c('SRR14800534', 'SRR14800535', 'SRR14800536'), cluster_dictionary = cluster_dictionary){

  numbat_rds_files <- retrieve_numbat_rds_files(numbat_dir, kept_samples)

  seus <- map(numbat_rds_files, retrieve_numbat_seurat, cluster_dictionary)

  integrated_seu <- seuratTools::seurat_integration_pipeline(seus, resolution = c(0.2, 0.4))

  sample_slug = paste(kept_samples, collapse = "_")

  seu_path = glue("output/seurat/{sample_slug}_seu.rds")

  saveRDS(integrated_seu, seu_path)

  return(seu_path)

}

simplify_gt <- function(mynb, rb_scnas = c("1" = "1q", "2" = "2p", "6" = "6p", "8" = "8p", "11" = "11p", "13" = "13q", "16" = "16q")){

  mynb$mut_graph

  to_labels <- igraph::edge_attr(mynb$mut_graph, "to_label") %>%
    str_split(",") %>%
    map(~(paste(names(rb_scnas[rb_scnas %in% .x]), collapse = ","))) %>%
    map_chr(str_pad, side = "left", width = 2) %>%
    identity()

  from_labels <- igraph::edge_attr(mynb$mut_graph, "from_label") %>%
    str_split(",") %>%
    map(~(paste(names(rb_scnas[rb_scnas %in% .x]), collapse = ","))) %>%
    map_chr(str_pad, side = "left", width = 2) %>%
    identity()

  mynb$mut_graph <- mynb$mut_graph %>%
    igraph::set_edge_attr("to_label", value = to_labels) %>%
    igraph::set_edge_attr("from_label", value = from_labels) %>%
    identity()

  return(mynb)
}

filter_input_by_scna <- function(df, scna_of_interest, p_val_threshold, fc_threshold, n_recur){
	# browser()
	filtered_df <- 
		df %>% 
		dplyr::group_by(symbol) %>%
		dplyr::filter(p_val_adj <= p_val_threshold) %>%
		dplyr::filter(abs(avg_log2FC) >= fc_threshold) %>%
		dplyr::filter(n_distinct(sample_id) >= n_recur) %>%
		identity()
}

filter_input_by_region <- function(df, segment_region){
	# browser()
	
	region_settings <- 
		dplyr::filter(oncoprint_settings, region == segment_region)
	
	df <- df[[segment_region]]
	
	test0 <- pmap(list(df, region_settings$scna, region_settings$p_val, region_settings$fc, region_settings$recurrence), filter_input_by_scna)
	
}

filter_oncoprint_diffex <- function(unfiltered_oncoprint_input_by_scna, oncoprint_settings){
	
	oncoprint_input_by_scna <- map(c("cis" = "cis", "trans" = "trans", "all" = "all"), ~filter_input_by_region(unfiltered_oncoprint_input_by_scna, .x))
	
	return(oncoprint_input_by_scna)
	
}


#' Title
#'
#' @param large_filter_expressions
#' @param cluster_dictionary
#' @param interesting_samples
#' @param cis_diffex_clones
#' @param trans_diffex_clones
#' @param large_clone_comparisons
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
make_oncoprint_diffex <- function(large_filter_expressions, cluster_dictionary, interesting_samples, cis_diffex_clones, trans_diffex_clones, all_diffex_clones, large_clone_comparisons, rb_scna_samples, by_cluster = FALSE, ...) {

  filter_diffex <- function(df, n_slice = 10){
    # browser()

    test0 <-
      df %>%
      dplyr::arrange(p_val_adj) %>%
      dplyr::group_by(symbol) %>%
      dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
      dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
      identity()

    test1 <-
      test0 %>%
      dplyr::summarize(mean_FC = mean(abs(avg_log2FC))) %>%
      dplyr::inner_join(test0, by = "symbol")

    return(test1)
  }

  names(cis_diffex_clones) <- interesting_samples
  names(trans_diffex_clones) <- interesting_samples
  names(all_diffex_clones) <- interesting_samples
  
  select_scna_diffex <- function(diffex_clones, selected_scna = "1q+", rb_scna_samples, ...){
  	# browser()
  	
  	segment = str_extract(selected_scna, "[0-9]*[a-z]")
  	
  	sign = str_extract(selected_scna, "[+,-]")
  	
  	comparisons <- 
  		diffex_clones[str_extract(names(diffex_clones), "SRR[0-9]*") %in% rb_scna_samples[[segment]]] %>% 
  		map(~.x[str_detect(names(.x), glue("{segment}\\{sign}$"))]) %>%
  		compact() %>%
  		map(dplyr::bind_rows, .id = "clone_comparison") %>%
  		dplyr::bind_rows(.id = "sample_id") %>%
  		# dplyr::filter(seqnames == str_extract(segment, "[0-9]*")) %>%
  		filter_diffex(...) %>%
  		identity()
  	
  }
  
  if(by_cluster){
  	cis_diffex_clones <- map(cis_diffex_clones, read_csv) %>% 
  		purrr::compact() %>% 
  		map(~split(.x, .x$clone_comparison))
  	trans_diffex_clones <- map(trans_diffex_clones, read_csv) %>% 
  		purrr::compact() %>% 
  		map(~split(.x, .x$clone_comparison))
  	all_diffex_clones <- map(all_diffex_clones, read_csv) %>% 
  		purrr::compact() %>% 
  		map(~split(.x, .x$clone_comparison))
  }
  
  # cis ------------------------------
  
  comparisons_of_1q <- select_scna_diffex(cis_diffex_clones, "1q+", rb_scna_samples)
  
  comparisons_of_2p <- select_scna_diffex(cis_diffex_clones, "2p+", rb_scna_samples)
  
  comparisons_of_6p <- select_scna_diffex(cis_diffex_clones, "6p+", rb_scna_samples)
  
  comparisons_of_16q <- select_scna_diffex(cis_diffex_clones, "16q-", rb_scna_samples)

  cis_comps <- list(comparisons_of_1q,
                        comparisons_of_2p,
                        comparisons_of_6p,
                        comparisons_of_16q) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  inspect_cis_comps <-
    cis_comps %>%
    map(~dplyr::distinct(.x, symbol, description)) %>%
    dplyr::bind_rows(.id = "clone_comparison")

  # trans ------------------------------
  # browser()
  comparisons_of_1q <- select_scna_diffex(trans_diffex_clones, "1q+", rb_scna_samples)
  
  comparisons_of_2p <- select_scna_diffex(trans_diffex_clones, "2p+", rb_scna_samples)
  
  comparisons_of_6p <- select_scna_diffex(trans_diffex_clones, "6p+", rb_scna_samples)
  
  comparisons_of_16q <- select_scna_diffex(trans_diffex_clones, "16q-", rb_scna_samples)

  trans_comps <- list(comparisons_of_1q,
                          comparisons_of_2p,
                          comparisons_of_6p,
                          comparisons_of_16q) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # all ------------------------------

  comparisons_of_1q <- select_scna_diffex(all_diffex_clones, "1q+", rb_scna_samples)
  
  comparisons_of_2p <- select_scna_diffex(all_diffex_clones, "2p+", rb_scna_samples)
  
  comparisons_of_6p <- select_scna_diffex(all_diffex_clones, "6p+", rb_scna_samples)
  
  comparisons_of_16q <- select_scna_diffex(all_diffex_clones, "16q-", rb_scna_samples)

  all_comps <- list(comparisons_of_1q,
                      comparisons_of_2p,
                      comparisons_of_6p,
                      comparisons_of_16q) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()
  
  if(by_cluster){
  	cis_comps <- 
  		cis_comps %>% 
  		map(tidyr::unite, "sample_id", sample_id, cluster, sep = "-")
  	
  	trans_comps <- 
  		trans_comps %>% 
  		map(tidyr::unite, "sample_id", sample_id, cluster, sep = "-")
  	
  	all_comps <- 
  		all_comps %>% 
  		map(tidyr::unite, "sample_id", sample_id, cluster, sep = "-")
  }

  return(list("cis" = cis_comps, "trans" = trans_comps, "all" = all_comps))

}

#' Title
#'
#' @param large_filter_expressions
#' @param cluster_dictionary
#' @param interesting_samples
#' @param cis_diffex_clones
#' @param trans_diffex_clones
#' @param large_clone_comparisons
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
make_oncoprint_diffex_unfiltered <- function(large_filter_expressions, cluster_dictionary, interesting_samples, cis_diffex_clones, trans_diffex_clones, large_clone_comparisons, ...) {

  filter_diffex_for_recurrence <- function(df, num_recur = 2, n_slice = 10){
    # browser()

    test0 <-
      df %>%
      # dplyr::arrange(symbol, sample_id) %>%
      dplyr::arrange(p_val_adj) %>%
      dplyr::group_by(symbol) %>%
      dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
      dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      # dplyr::filter(n_distinct(sample_id) >= num_recur) %>%
      identity()

    test1 <-
      test0 %>%
      dplyr::summarize(mean_FC = mean(abs(avg_log2FC))) %>%
      # dplyr::slice_max(abs(mean_FC), n = n_slice) %>%
      dplyr::inner_join(test0, by = "symbol")

    return(test1)
  }

  # cis ------------------------------
  names(cis_diffex_clones) <- interesting_samples

  names(trans_diffex_clones) <- interesting_samples

  clone_comparisons <- map(cis_diffex_clones, names)


  comparisons_of_1q <- map(clone_comparisons, str_detect, "1q\\+$") %>%
    map2(cis_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(seqnames == "1") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_2p <- map(clone_comparisons, str_detect, "2p\\+$") %>%
    map2(cis_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(seqnames == "2") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_6p <- map(clone_comparisons, str_detect, "6p\\+$") %>%
    map2(cis_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(seqnames == "6") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_16q <- map(clone_comparisons, str_detect, "[0-9]_v_[0-9]_16q\\-$") %>%
    map2(cis_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(seqnames == "16") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  cis_comps <- list(comparisons_of_1q,
                        comparisons_of_2p,
                        comparisons_of_6p,
                        comparisons_of_16q) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  inspect_cis_comps <-
    cis_comps %>%
    map(~dplyr::distinct(.x, symbol, description)) %>%
    dplyr::bind_rows(.id = "clone_comparison")

  # trans ------------------------------
  names(trans_diffex_clones) <- interesting_samples

  names(trans_diffex_clones) <- interesting_samples

  clone_comparisons <- map(trans_diffex_clones, names)

  comparisons_of_1q <- map(clone_comparisons, str_detect, "1q\\+$") %>%
    map2(trans_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(abs(avg_log2FC) > 0.25) %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_2p <- map(clone_comparisons, str_detect, "2p\\+$") %>%
    map2(trans_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_6p <- map(clone_comparisons, str_detect, "6p\\+$") %>%
    map2(trans_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  comparisons_of_16q <- map(clone_comparisons, str_detect, "[0-9]_v_[0-9]_16q\\-$") %>%
    map2(trans_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(dplyr::bind_rows, .id = "clone_comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    filter_diffex_for_recurrence(num_recur = 0, ...) %>%
    identity()

  trans_comps <- list(comparisons_of_1q,
                          comparisons_of_2p,
                          comparisons_of_6p,
                          comparisons_of_16q) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  return(list("cis" = cis_comps, "trans" = trans_comps))

}

inspect_oncoprints <- function(cis_comps, trans_comps, all_comps){
  inspect_trans_comps <-
    trans_comps %>%
    map(~dplyr::distinct(.x, symbol, description, .keep_all = TRUE)) %>%
    dplyr::bind_rows(.id = "clone_comparison") %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(mean_FC = mean(avg_log2FC)) %>%
    dplyr::arrange(clone_comparison, desc(mean_FC)) %>%
    dplyr::select(symbol, description, everything())

  inspect_cis_comps <-
    cis_comps %>%
    map(~dplyr::distinct(.x, symbol, description, .keep_all = TRUE)) %>%
    dplyr::bind_rows(.id = "clone_comparison") %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(mean_FC = mean(avg_log2FC)) %>%
    dplyr::arrange(clone_comparison, desc(mean_FC)) %>%
    dplyr::select(symbol, description, everything())

  inspect_all_comps <-
    all_comps %>%
    map(~dplyr::distinct(.x, symbol, description, .keep_all = TRUE)) %>%
    dplyr::bind_rows(.id = "clone_comparison") %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(mean_FC = mean(avg_log2FC)) %>%
    dplyr::arrange(clone_comparison, desc(mean_FC)) %>%
    dplyr::select(symbol, description, everything())

  return(list("cis" = inspect_cis_comps, "trans" = inspect_trans_comps, "all" = inspect_all_comps))


}

plot_recurrence <- function(diffex_input, scna_of_interest, segment_region = "cis", oncoprint_settings, clone_trees, n_genes = 20){
  # browser()
	
	region_settings <- 
		dplyr::filter(oncoprint_settings, region == segment_region) |> 
		dplyr::filter(scna == scna_of_interest)
	
	clone_trees <- 
		clone_trees %>% 
		set_names(map(., c("labels", "title")))

  mytitle = scna_of_interest
  mysubtitle = glue("
  									recurrence: {region_settings$recurrence}
  									p_val <= {region_settings$p_val}
  									fold-change >= {region_settings$fc}
  									")

  n_recurrence = region_settings$recurrence

  diffex_table <-
    diffex_input %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(mean_FC = mean(avg_log2FC)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(mean_FC)) %>%
    dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) %>%
    dplyr::mutate(neg_log_p_val_adj = -log10(p_val_adj)) %>%
    dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
    dplyr::group_by(symbol) %>%
    dplyr::filter(n_distinct(sample_id) >= n_recurrence) %>%
    dplyr::ungroup() %>%
    # dplyr::slice_head(n= 10) %>%
    identity()

  plot_input <-
    diffex_table %>%
    dplyr::filter(p_val_adj <= region_settings$p_val) %>%
  	dplyr::filter(abs_log2FC >= region_settings$fc) %>%
  	# dplyr::filter(abs_log2FC >= 0.5) %>%
    dplyr::select(sample_id, symbol, description, neg_log_p_val_adj, abs_log2FC, avg_log2FC, mean_FC) %>%
    group_by(symbol) %>%
    dplyr::mutate(num_positive = sum(avg_log2FC > 0)) %>%
    dplyr::mutate(num_negative = sum(avg_log2FC < 0)) %>%
    dplyr::mutate(same_sign = all_same_sign(avg_log2FC)) %>%
    dplyr::mutate(major_sign = abs(num_positive - num_negative)) %>%
    dplyr::filter(same_sign > 0 | major_sign > 1) %>%
    dplyr::mutate(minor_sign = min(num_positive, num_negative)) %>%
    dplyr::filter(major_sign >= n_recurrence) %>%
    # dplyr::filter(minor_sign < 2) %>%
    identity()

  top_genes = unique(plot_input$symbol)[1:n_genes]
  top_genes = top_genes[!is.na(top_genes)]


  plot_input <-
    plot_input %>%
    # dplyr::filter(symbol %in% top_genes) %>%
    dplyr::filter(minor_sign < 2) %>%
    dplyr::mutate(comparison_sign = ifelse(num_positive > num_negative, 1, -1)) %>%
    dplyr::filter(sign(avg_log2FC) == sign(comparison_sign)) %>%
    dplyr::group_by(symbol) %>%
    dplyr::mutate(n_samples = n_distinct(sample_id)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(n_samples, desc(mean_FC)) %>%
    dplyr::mutate(symbol = factor(symbol, levels = unique(symbol))) %>%
    identity()
  
  # clone_trees

  diffex_plot <- ggplot(plot_input, aes(x=sample_id, y=symbol)) +
    geom_point(aes(color=neg_log_p_val_adj, size = abs_log2FC)) +
    # geom_tile(fill = NA, color = "black", linewidth = 0.5) +
    labs(title = mytitle, subtitle = mysubtitle, color = "-log10 \n p_adj") +
    theme_bw() +
    # scale_size_continuous(
    #   limits = c(0.1, 0.9)
    # ) +
    scale_color_continuous(
      limits = c(1, 150),
      oob = squish
    ) +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
      # legend.position="none"
    )
  
  sub_clone_trees <- 
  	clone_trees[names(clone_trees) %in% sort(unique(str_remove(plot_input$sample_id, "-.*")))] %>% 
  	map(~{.x + theme(plot.title = element_blank())}) %>% 
  	wrap_plots(nrow = 1)
  
  diffex_plot <- wrap_plots(diffex_plot, sub_clone_trees, ncol = 1)
  

  return(list("table" = diffex_table, "plot" = diffex_plot))


}


make_oncoprint_plots <- function(comps, clone_trees, oncoprint_settings, label = "_by_clone", ...) {
	# browser()

	comps_res <- list()
  # cis ------------------------------
	for(region in names(comps)){
		comps_res[[region]]  <- imap(comps[[region]], plot_recurrence, region, oncoprint_settings, clone_trees, ...)
	}
	

  cis_comps_plots <- map(comps_res[["cis"]], "plot")

  # wrap_plots(cis_comps_plots) +  plot_layout(guides = 'collect') +
  #   plot_annotation(title = "top 10 DE genes in segment")

  in_segment_plot_path = glue("results/diffex_oncoprints{label}_in_segment.pdf")
  
  pdf(in_segment_plot_path, height = 8, width = 6)
  cis_comps_plots %>% 
  	walk(print)
  dev.off()

  in_segment_table_path = glue("results/diffex_oncoprints{label}_in_segment.xlsx")

  cis_comps_tables <-
    comps_res[["cis"]] %>%
    map("table") %>%
    map(dplyr::select, -genes_in_segment) %>%
    writexl::write_xlsx(in_segment_table_path)

  # trans ------------------------------

  trans_comps_plots <- map(comps_res[["trans"]], "plot")

  # wrap_plots(trans_comps_plots) +  plot_layout(guides = 'collect') +
  #   plot_annotation(title = "top 10 DE genes out of segment")

  out_segment_plot_path = glue("results/diffex_oncoprints{label}_out_of_segment.pdf")
  
  pdf(out_segment_plot_path, height = 8, width = 6)
  trans_comps_plots %>% 
  	walk(print)
  dev.off()

  out_segment_table_path = glue("results/diffex_oncoprints{label}_out_of_segment.xlsx")

  comps_res[["trans"]] %>%
    map("table") %>%
    writexl::write_xlsx(out_segment_table_path)

  # all ------------------------------

  all_comps_plots <- map(comps_res[["all"]], "plot")

  all_plot_path = glue("results/diffex_oncoprints{label}_all.pdf")
  

  pdf(all_plot_path, height = 8, width = 6)
  all_comps_plots %>% 
  	walk(print)
  dev.off()

  all_table_path = glue("results/diffex_oncoprints{label}_all.xlsx")

  comps_res[["all"]] %>%
    map("table") %>%
    writexl::write_xlsx(all_table_path)

  return(
    list("cis" =
           list("plot" = in_segment_plot_path,
                "table" = in_segment_table_path),
         "trans" =
           list("plot" = out_segment_plot_path,
                "table" = out_segment_table_path),
         "all" =
           list("plot" = all_plot_path,
                "table" = all_table_path)
         ))

}

enrich_oncoprints <- function(large_filter_expressions, cluster_dictionary, interesting_samples, cis_diffex_clones, trans_diffex_clones, all_diffex_clones, large_clone_comparisons, ...) {

  filter_diffex_for_recurrence <- function(df, num_recur = 2, n_slice = 10){
    # browser()

    test0 <-
      df %>%
      # dplyr::arrange(symbol, sample_id) %>%
      dplyr::arrange(p_val_adj) %>%
      dplyr::group_by(symbol) %>%
      dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
      dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      dplyr::filter(n_distinct(sample_id) >= num_recur) %>%
      identity()

    test1 <-
      test0 %>%
      dplyr::summarize(mean_FC = mean(abs(avg_log2FC))) %>%
      # dplyr::slice_max(abs(mean_FC), n = n_slice) %>%
      dplyr::inner_join(test0, by = "symbol")

    return(test1)
  }

  # cis ------------------------------
  names(cis_diffex_clones) <- interesting_samples

  add_clone_comparison_column <- function(df, mycomparison){
    dplyr::mutate(df, clone_comparison = mycomparison)
  }

  cis_diffex_clones <- purrr::map(cis_diffex_clones, ~purrr::imap(.x, add_clone_comparison_column))

  clone_comparisons <- map(cis_diffex_clones, names)


  comparisons_of_1q <- map(clone_comparisons, str_detect, "1q\\+$") %>%
    map2(cis_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(~map(.x, dplyr::filter, seqnames == "1")) %>%
    identity()

  comparisons_of_2p <- map(clone_comparisons, str_detect, "2p\\+$") %>%
    map2(cis_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(~map(.x, dplyr::filter, seqnames == "2")) %>%
    identity()

  comparisons_of_6p <- map(clone_comparisons, str_detect, "6p\\+$") %>%
    map2(cis_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(~map(.x, dplyr::filter, seqnames == "6")) %>%
    identity()

  comparisons_of_16q <- map(clone_comparisons, str_detect, "[0-9]_v_[0-9]_16q\\-$") %>%
    map2(cis_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    map(~map(.x, dplyr::filter, seqnames == "16")) %>%
    identity()

  cis_comps <- list(comparisons_of_1q,
                        comparisons_of_2p,
                        comparisons_of_6p,
                        comparisons_of_16q) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

   # trans ------------------------------
  names(trans_diffex_clones) <- interesting_samples

  trans_diffex_clones <- purrr::map(trans_diffex_clones, ~purrr::imap(.x, add_clone_comparison_column))

  clone_comparisons <- map(trans_diffex_clones, names)

  comparisons_of_1q <- map(clone_comparisons, str_detect, "1q\\+$") %>%
    map2(trans_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    identity()

  comparisons_of_2p <- map(clone_comparisons, str_detect, "2p\\+$") %>%
    map2(trans_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    identity()

  comparisons_of_6p <- map(clone_comparisons, str_detect, "6p\\+$") %>%
    map2(trans_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    identity()

  comparisons_of_16q <- map(clone_comparisons, str_detect, "[0-9]_v_[0-9]_16q\\-$") %>%
    map2(trans_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    identity()

  trans_comps <- list(comparisons_of_1q,
                          comparisons_of_2p,
                          comparisons_of_6p,
                          comparisons_of_16q) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # all ------------------------------
  names(all_diffex_clones) <- interesting_samples

  all_diffex_clones <- purrr::map(all_diffex_clones, ~purrr::imap(.x, add_clone_comparison_column))

  clone_comparisons <- map(all_diffex_clones, names)

  comparisons_of_1q <- map(clone_comparisons, str_detect, "1q\\+$") %>%
    map2(all_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    identity()

  comparisons_of_2p <- map(clone_comparisons, str_detect, "2p\\+$") %>%
    map2(all_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    identity()

  comparisons_of_6p <- map(clone_comparisons, str_detect, "6p\\+$") %>%
    map2(all_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    identity()

  comparisons_of_16q <- map(clone_comparisons, str_detect, "[0-9]_v_[0-9]_16q\\-$") %>%
    map2(all_diffex_clones, ~{.y[.x]}) %>%
    compact() %>%
    identity()

  all_comps <- list(comparisons_of_1q,
                          comparisons_of_2p,
                          comparisons_of_6p,
                          comparisons_of_16q) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # proceed ------------------------------

  possible_prep_for_enrichment <- purrr::possibly(prep_for_enrichment, otherwise = NA_real_)

  cis_enrichment_tables <- modify_depth(cis_comps, 3, possible_prep_for_enrichment, .ragged = TRUE, ...)

  trans_enrichment_tables <- modify_depth(trans_comps, 3, possible_prep_for_enrichment, .ragged = TRUE, ...)

  all_enrichment_tables <- modify_depth(all_comps, 3, possible_prep_for_enrichment, .ragged = TRUE, ...)

  # cis_enrichment_plots <- modify_depth(cis_enrichment_tables, 3, plot_enrichment, .ragged = TRUE) %>%
  #   purrr::list_flatten() %>%
  #   purrr::list_flatten() %>%
  #   imap(~{.x + labs(title = .y)}) %>%
  #   identity()
  #
  # pdf("cis.pdf")
  # cis_enrichment_plots
  # dev.off()
  #
  # trans_enrichment_plots <- modify_depth(trans_enrichment_tables, 3, plot_enrichment, .ragged = TRUE) %>%
  #   purrr::list_flatten() %>%
  #   purrr::list_flatten() %>%
  #   imap(~{.x + labs(title = .y)}) %>%
  #   identity()
  #
  # pdf("trans.pdf")
  # trans_enrichment_plots
  # dev.off()

  return(list("cis" = cis_enrichment_tables, "trans" = trans_enrichment_tables, "all" = all_enrichment_tables))

}

enrich_oncoprints_clusters <- function(large_filter_expressions, cluster_dictionary, interesting_samples, cis_diffex_clones_for_each_cluster, trans_diffex_clones_for_each_cluster, large_clone_comparisons, ...) {

  cis_diffex_clones_for_each_cluster <- map(cis_diffex_clones_for_each_cluster, read_csv)

  trans_diffex_clones_for_each_cluster <- map(trans_diffex_clones_for_each_cluster, read_csv)

  filter_diffex_for_recurrence <- function(df, num_recur = 2, n_slice = 10){
    # browser()

    test0 <-
      df %>%
      # dplyr::arrange(symbol, sample_id) %>%
      dplyr::arrange(p_val_adj) %>%
      dplyr::group_by(symbol) %>%
      dplyr::mutate(neg_log_p_val_adj = -log(p_val_adj, base = 10)) %>%
      dplyr::mutate(abs_log2FC = abs(avg_log2FC)) %>%
      dplyr::filter(p_val_adj < 0.05) %>%
      dplyr::filter(n_distinct(sample_id) >= num_recur) %>%
      identity()

    test1 <-
      test0 %>%
      dplyr::summarize(mean_FC = mean(abs(avg_log2FC))) %>%
      # dplyr::slice_max(abs(mean_FC), n = n_slice) %>%
      dplyr::inner_join(test0, by = "symbol")

    return(test1)
  }

  # cis ------------------------------
  names(cis_diffex_clones_for_each_cluster) <- interesting_samples

  cis_diffex_clones_for_each_cluster <-
    cis_diffex_clones_for_each_cluster %>%
    keep(~nrow(.x) > 0)

  clone_comparisons <- map(cis_diffex_clones_for_each_cluster, ~unique(.x[["clone_comparison"]]))

  shape_cluster_comparison <- function(mylist){
    mylist %>%
      purrr::discard(~nrow(.x) == 0) %>%
      map(~split(.x, .x$cluster)) %>%
      identity()
  }

  comparisons_of_1q <- map(cis_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "1q\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_2p <- map(cis_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "2p\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_6p <- map(cis_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "6p\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_16q <- map(cis_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "[0-9]_v_[0-9]_16q\\-$")) %>%
    shape_cluster_comparison()

  cis_comps <- list(comparisons_of_1q,
                        comparisons_of_2p,
                        comparisons_of_6p,
                        comparisons_of_16q) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # trans ------------------------------

  names(trans_diffex_clones_for_each_cluster) <- interesting_samples

  trans_diffex_clones_for_each_cluster <-
    trans_diffex_clones_for_each_cluster %>%
    keep(~nrow(.x) > 0)

  clone_comparisons <- map(trans_diffex_clones_for_each_cluster, ~unique(.x[["clone_comparison"]]))

  comparisons_of_1q <- map(trans_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "1q\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_2p <- map(trans_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "2p\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_6p <- map(trans_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "6p\\+$")) %>%
    shape_cluster_comparison()

  comparisons_of_16q <- map(trans_diffex_clones_for_each_cluster, dplyr::filter, str_detect(clone_comparison, "[0-9]_v_[0-9]_16q\\-$")) %>%
    shape_cluster_comparison()

  trans_comps <- list(comparisons_of_1q,
                          comparisons_of_2p,
                          comparisons_of_6p,
                          comparisons_of_16q) %>%
    set_names(c("1q+", "2p+", "6p+", "16q-")) %>%
    identity()

  # proceed ------------------------------

  possible_prep_for_enrichment <- purrr::possibly(prep_for_enrichment, otherwise = NA_real_)

  cis_enrichment_tables <- modify_depth(cis_comps, 3, possible_prep_for_enrichment, .ragged = TRUE, ...)

  trans_enrichment_tables <- modify_depth(trans_comps, 3, possible_prep_for_enrichment, .ragged = TRUE, ...)

  # cis_enrichment_plots <- modify_depth(cis_enrichment_tables, 3, plot_enrichment, .ragged = TRUE) %>%
  #   purrr::list_flatten() %>%
  #   purrr::list_flatten() %>%
  #   imap(~{.x + labs(title = .y)}) %>%
  #   identity()
  #
  # pdf("cis_cluster.pdf")
  # cis_enrichment_plots
  # dev.off()
  #
  # trans_enrichment_plots <- modify_depth(trans_enrichment_tables, 3, plot_enrichment, .ragged = TRUE) %>%
  #   purrr::list_flatten() %>%
  #   purrr::list_flatten() %>%
  #   imap(~{.x + labs(title = .y)}) %>%
  #   identity()
  #
  # pdf("trans_cluster.pdf")
  # trans_enrichment_plots
  # dev.off()

  return(list("cis" = cis_enrichment_tables, "trans" = trans_enrichment_tables))

}

#' Make integrated numbat plots
#'
#' @param seu_path
#' @param extension
#'
#' @return
#' @export
#'
#' @examples
make_integrated_numbat_plots <- function(seu_path, extension = ""){
  # browser()

  sample_id = str_extract(seu_path, "SRR.*(?=_seu.rds)")

  dir_create(glue("results/numbat_sridhar/{sample_id}"))

  seu <- readRDS(seu_path) %>%
    filter_sample_qc()

  seu <- seu[,!seu$abbreviation %in% c("APOE", "MALAT1")]

  seu@meta.data$abbreviation <-
  seu@meta.data$abbreviation %>%
    str_replace_all("ARL1IP1", "ARL6IP1") %>%
    str_replace_all("PCLAF", "TFF1")


  seu <- Seurat::RenameCells(seu, new.names = str_replace(colnames(seu), "\\.", "-"))

  seu <- seu[,!is.na(seu$clone_opt)]

  # seu <- seurat_cluster(seu, resolution = 0.1)
  #
  # plot_markers(seu, metavar = "integrated_snn_res.0.1", marker_method = "presto", return_plotly = FALSE, hide_technical = "all") +
  #   ggplot2::scale_y_discrete(position = "left")

  markerplot <- plot_markers(seu, metavar = "abbreviation", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 3) +
    ggplot2::scale_y_discrete(position = "left")

  # ggsave(glue("results/numbat_sridhar/{sample_id}/{sample_id}_sample_marker{extension}.pdf"), height = 8, width = 6)

  dimplot <- DimPlot(seu, group.by = c("abbreviation", "Phase")) +
    plot_annotation(title = sample_id)
  # ggsave(glue("results/numbat_sridhar/{sample_id}/{sample_id}_dimplot{extension}.pdf"), width = 8, height = 4)


  ## clone distribution ------------------------------
  distplot <- plot_distribution_of_clones_across_clusters(seu, sample_id)

  # ggsave(glue("results/numbat_sridhar/{sample_id}/{sample_id}_clone_distribution{extension}.pdf"), width = 4, height = 4)

  seu_list <- list("markerplot" = markerplot, "dimplot" = dimplot, "distplot" = distplot)

  return(seu_list)
}

plot_scna_violins <- function(., readRDS, FeaturePlot, .x, .y, VlnPlot) {

  # plot drivers ------------------------------

  samples_1q <- c("SRR14800534", "SRR14800536", "SRR14800535", "SRR13884249",
                  "SRR17960484")

  samples_2p <- c("SRR13884249", "SRR17960484", "SRR17960481", "SRR13884247")

  samples_6p <- c("SRR17960484", "SRR13884247", "SRR17960481")

  samples_16q <- c("SRR14800540", "SRR13884242", "SRR14800535", "SRR14800534", "SRR14800543", "SRR14800536", "SRR13884243")

  # samples_16q <- c("SRR14800540", "SRR14800543", "SRR14800536")

  plot_markers_featureplot <- function(sample_ids, features){
    # browser()

    seus <- glue("output/seurat/{sample_ids}_filtered_seu.rds") %>%
      set_names(str_extract(., "SRR[0-9]*")) %>%
      map(readRDS) %>%
      identity()

    markerplots <- map(seus, FeaturePlot, features = features, combine = TRUE)

    markerplots0 <- imap(markerplots, ~{.x + labs(subtitle = .y)})

    vlnplots <- map(seus, VlnPlot, features = features, group.by = "GT_opt", combine = TRUE, pt.size = 0)

    vlnplots0 <- imap(vlnplots, ~{.x +
        labs(title = .y) +
        theme(legend.position="none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank())})

    # vln_wrapped = wrap_plots(vlnplots0, ncol = 2) + plot_annotation(title = features)

    return(vlnplots0)
  }

  vlns_1q <- plot_markers_featureplot(samples_1q, "CENPF")

  vlns_2p <- plot_markers_featureplot(samples_2p, "SOX11")

  vlns_6p <- plot_markers_featureplot(samples_6p, "DEK")

  vlns_16q <- plot_markers_featureplot(samples_16q, "CDT1")


  vlns_1q +
    wrap_plots(nrow = 1) + plot_annotation(title = "CENPF")
  ggsave("results/vln_plots_from_oncoprint_1q.pdf", width = 6, height = 8)
  vlns_2p +
    wrap_plots(nrow = 1) + plot_annotation(title = "SOX11")
  ggsave("results/vln_plots_from_oncoprint_2p.pdf", width = 6, height = 6)
  vlns_6p +
    wrap_plots(nrow = 1) + plot_annotation(title = "DEK")
  ggsave("results/vln_plots_from_oncoprint_6p.pdf", width = 6, height = 6)
  vlns_16q +
    wrap_plots(nrow = 1) + plot_annotation(title = "CDT1")
  ggsave("results/vln_plots_from_oncoprint_16q.pdf", width = 6, height = 6)

}


compile_cis_trans_enrichment_recurrence <- function(enrichment_tables,
																										cis_plot_file = "results/cis_enrichment_plots",
																										trans_plot_file = "results/trans_enrichment_plots",
																										cis_table_file = "results/cis_enrichment_tables.xlsx",
																										trans_table_file = "results/trans_enrichment_tables.xlsx",
																										...){
	
	# cis ------------------------------
	recurrences <- c(3,1,1,1)
	
	titles = c(
		"1q+ cis enrichment",
		"2p+ cis enrichment",
		"6p+ cis enrichment",
		"16q- cis enrichment"
	)
	
	cis_enrich_results <- pmap(list(enrichment_tables$cis, recurrences, titles), plot_enrichment_recurrence, ...)
	
	# trans ------------------------------
	
	recurrences <- c(2,2,2,3)
	
	titles = c(
		"1q+ trans enrichment",
		"2p+ trans enrichment",
		"6p+ trans enrichment",
		"16q- trans enrichment"
	)
	
	trans_enrich_results <- pmap(list(enrichment_tables$trans, recurrences, titles), plot_enrichment_recurrence, ...)
	
	# cis ------------------------------
	# enrich plot
	cis_enrich_plots <- map(cis_enrich_results, "enrich_plot") %>%
		purrr::discard(~nrow(.x$data) < 1)
	
	names(cis_enrich_plots) <- str_remove_all(names(cis_enrich_plots), "[+]|-")
	
	enrich_plot_widths <- map_int(cis_enrich_plots, ~n_distinct(.x$data$clone_comparison))*1+1
	
	fs::dir_create(fs::path_ext_remove(cis_plot_file))
	cis_pdfs <- map2(names(cis_enrich_plots), enrich_plot_widths, ~ggsave(glue("{fs::path_ext_remove(cis_plot_file)}/{.x}_enrich.pdf"), cis_enrich_plots[[.x]], width = .y, height = 7))
	
	qpdf::pdf_combine(cis_pdfs, cis_plot_file)
	
	per_scna_enrichment_cis <- imap(enrichment_tables$cis, plot_enrichment_per_scna, cis_plot_file)
	
	# trans ------------------------------
	# enrich plot
	trans_enrich_plots <- map(trans_enrich_results, "enrich_plot") %>%
		purrr::discard(~nrow(.x$data) < 1)
	
	names(trans_enrich_plots) <- str_remove_all(names(trans_enrich_plots), "[+]|-")
	
	enrich_plot_widths <- map_int(trans_enrich_plots, ~n_distinct(.x$data$clone_comparison))*1 + 1
	
	fs::dir_create(fs::path_ext_remove(trans_plot_file))
	trans_pdfs <- map2(names(trans_enrich_plots), enrich_plot_widths, ~ggsave(glue("{fs::path_ext_remove(trans_plot_file)}/{.x}_enrich.pdf"), trans_enrich_plots[[.x]], width = .y, height = 7, limitsize = FALSE))
	
	qpdf::pdf_combine(trans_pdfs, trans_plot_file)
	
	per_scna_enrichment_trans <- imap(enrichment_tables$trans, plot_enrichment_per_scna, trans_plot_file)
	
	# cis tables
	cis_tables <- map(cis_enrich_results, "table") %>%
		purrr::discard(~nrow(.x) < 1) %>%
		purrr::map(dplyr::mutate, comparison = purrr::pluck(str_split_1(clone_comparison, "_"), 2)) %>%
		purrr::map(dplyr::arrange, comparison, p.adjust)
	
	writexl::write_xlsx(cis_tables, cis_table_file)
	
	# trans tables
	trans_tables <- map(trans_enrich_results, "table") %>%
		purrr::discard(~nrow(.x) < 1) %>%
		purrr::map(dplyr::mutate, comparison = purrr::pluck(str_split_1(clone_comparison, "_"), 2)) %>%
		purrr::map(dplyr::arrange, comparison, p.adjust)
	
	writexl::write_xlsx(trans_tables, trans_table_file)
	
	return(
		list(
			"cis_enrich_compiled" = cis_plot_file,
			"trans_enrich_compiled" = trans_plot_file,
			"cis_table" = cis_table_file,
			"trans_table" = trans_table_file,
			"cis_per_sample" = per_scna_enrichment_cis,
			"trans_per_sample" = per_scna_enrichment_trans
		)
	)
	
}

compile_cis_trans_enrichment_recurrence_by_cluster <- function(enrichment_tables,
                                                    cis_plot_file = "results/cis_enrichment_plots",
                                                    trans_plot_file = "results/trans_enrichment_plots",
                                                    cis_table_file = "results/cis_enrichment_tables.xlsx",
                                                    trans_table_file = "results/trans_enrichment_tables.xlsx",
                                                    ...){

  # cis ------------------------------
  recurrences <- c(3,1,1,1)

  titles = c(
    "1q+ cis enrichment",
    "2p+ cis enrichment",
    "6p+ cis enrichment",
    "16q- cis enrichment"
  )

  cis_enrich_results <- pmap(list(enrichment_tables$cis, recurrences, titles), plot_enrichment_recurrence_by_cluster, ...)

  # trans ------------------------------

  recurrences <- c(2,2,2,3)

  titles = c(
    "1q+ trans enrichment",
    "2p+ trans enrichment",
    "6p+ trans enrichment",
    "16q- trans enrichment"
  )

  trans_enrich_results <- pmap(list(enrichment_tables$trans, recurrences, titles), plot_enrichment_recurrence_by_cluster, ...)

  # cis ------------------------------
  # enrich plot
  cis_enrich_plots <- map(cis_enrich_results, "enrich_plot_by_phase") %>%
    purrr::discard(~nrow(.x$data) < 1)

  names(cis_enrich_plots) <- str_remove_all(names(cis_enrich_plots), "[+]|-")

  enrich_plot_widths <- map_int(cis_enrich_plots, ~n_distinct(.x$data$clone_comparison))*0.3+1

  fs::dir_create(fs::path_ext_remove(cis_plot_file))
  cis_pdfs <- map2(names(cis_enrich_plots), enrich_plot_widths, ~ggsave(glue("{fs::path_ext_remove(cis_plot_file)}/{.x}_enrich.pdf"), cis_enrich_plots[[.x]], width = .y, height = 10))

  qpdf::pdf_combine(cis_pdfs, cis_plot_file)

  per_scna_enrichment_cis <- imap(enrichment_tables$cis, plot_enrichment_per_scna, cis_plot_file)

  # trans ------------------------------
  # enrich plot
  trans_enrich_plots <- map(trans_enrich_results, "enrich_plot_by_phase") %>%
    purrr::discard(~nrow(.x$data) < 1)

  names(trans_enrich_plots) <- str_remove_all(names(trans_enrich_plots), "[+]|-")

  enrich_plot_widths <- map_int(trans_enrich_plots, ~n_distinct(.x$data$clone_comparison))*0.3 + 1

  fs::dir_create(fs::path_ext_remove(trans_plot_file))
  trans_pdfs <- map2(names(trans_enrich_plots), enrich_plot_widths, ~ggsave(glue("{fs::path_ext_remove(trans_plot_file)}/{.x}_enrich.pdf"), trans_enrich_plots[[.x]], width = .y, height = 10, limitsize = FALSE))

  qpdf::pdf_combine(trans_pdfs, trans_plot_file)

  per_scna_enrichment_trans <- imap(enrichment_tables$trans, plot_enrichment_per_scna, trans_plot_file)

  # cis tables
  cis_tables <- map(cis_enrich_results, "table") %>%
    purrr::discard(~nrow(.x) < 1) %>%
    purrr::map(dplyr::mutate, comparison = purrr::pluck(str_split_1(clone_comparison, "_"), 2)) %>%
    purrr::map(dplyr::arrange, comparison, p.adjust)

  writexl::write_xlsx(cis_tables, cis_table_file)

  # trans tables
  trans_tables <- map(trans_enrich_results, "table") %>%
    purrr::discard(~nrow(.x) < 1) %>%
    purrr::map(dplyr::mutate, comparison = purrr::pluck(str_split_1(clone_comparison, "_"), 2)) %>%
    purrr::map(dplyr::arrange, comparison, p.adjust)

  writexl::write_xlsx(trans_tables, trans_table_file)

  return(
  	list(
  		"cis_enrich_compiled" = cis_plot_file,
  		"trans_enrich_compiled" = trans_plot_file,
  		"cis_table" = cis_table_file,
  		"trans_table" = trans_table_file,
  		"cis_per_sample" = per_scna_enrichment_cis,
  		"trans_per_sample" = per_scna_enrichment_trans
  		)
  	)

}

score_and_vlnplot_seu <- function(seu_path, nb_path, clone_simplifications, gene_lists, ...){
  # browser()

  mysample = str_extract(seu_path, "SRR[0-9]*")

  seu <- readRDS(seu_path)
  
  # gene_lists <- map(gene_lists, ~(.x[.x %in% rownames(seu$SCT)]))
  
  # DefaultAssay(seu) <- "gene"
  nbin <- floor(length(VariableFeatures(seu))/100)

  seu <- Seurat::AddModuleScore(seu, features = gene_lists, name = "subtype", nbin = nbin, ctrl = 100)

  module_names = paste0("subtype", seq(length(gene_lists)))

  seu@meta.data[names(gene_lists)] <- seu@meta.data[module_names]

  # cluster_seu <- seu[,!is.na(seu$leiden)]

  pub_violin <- function(score.by = "subtype2", group.by = "scna", seu, y_lim = 0.4, step = 0.1){

    mydf <- seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::mutate(scna = fct_reorder(scna, clone_opt)) %>%
      identity()

    if(is.factor(mydf[[group.by]])){
      mycomparisons = seq(1, length(levels(mydf[[group.by]])))
    } else{
      mycomparisons = rev(unique(mydf[[group.by]]))
    }

    mycomparisons = map2(head(mycomparisons, -1), mycomparisons[-1], c)

    ggpubr::ggviolin(mydf, x = group.by, y = score.by, fill = group.by, add = "boxplot", add.params = list(fill = "white")) +
    # stat_compare_means(label.y = y_lim) +
      stat_compare_means(comparisons = mycomparisons, label = "p.signif") +
      scale_y_continuous(limits = c(0, y_lim), breaks = seq(0, y_lim, step))
  }

  # cluster_vln = map(names(gene_lists), pub_violin, group.by = "abbreviation", seu = seu, ...)
  # 
  # cluster_vln <- map(cluster_vln, ~{.x + guides(fill="none") + scale_x_discrete(labels = function(x) wrap_scna_labels(x))})
  # 
  # names(cluster_vln) <- names(gene_lists)

  scna_vln = map(names(gene_lists), pub_violin, group.by = "clone_opt", seu = seu, ...)
  
  names(scna_vln) <- paste0(names(gene_lists), " genes")

  scna_vln <- imap(scna_vln, ~{.x + guides(fill="none") + scale_x_discrete(labels = function(x) wrap_scna_labels(x)) + labs(x = "clone", y = .y)})

  
  myclonetree <- plot_clone_tree(seu, mysample, nb_path, clone_simplifications)
  
  design <- 
 "AAAA
  AAAA
  BBCC
  BBCC
  BBCC"
  
  test0 <- list(myclonetree, scna_vln[[1]], scna_vln[[2]]) %>% 
  	wrap_plots() +
  	plot_layout(design = design)
  
  scna_vln_plot_path <- ggsave(glue("results/{mysample}_scna_vln.pdf"), height = 6, width = 6)

  return(scna_vln_plot_path)
}

score_and_heatmap_seu <- function(seu_path, gene_lists, group.by = "SCT_snn_res.0.4", leiden_cluster_file = NULL){
  # browser()

  mysample = str_extract(seu_path, "SRR[0-9]*")

  # for(geneset in names(gene_lists)){
  #   # seu <- Seurat::MetaFeature(seu, features = gene_lists[[geneset]], meta.name = geneset)
  # }

  seu <- readRDS(seu_path)

  seu <- Seurat::AddModuleScore(seu, features = gene_lists, name = "subtype")

  module_names = paste0("subtype", seq(length(gene_lists)))

  seu@meta.data[names(gene_lists)] <- seu@meta.data[module_names]

  cluster_mps <-
    seu@meta.data[,c(names(gene_lists), group.by)] %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(-c(group.by, "cell"), names_to = "mp", values_to = "score") %>%
    group_by(.data[[group.by]], mp) %>%
    dplyr::summarize(score = mean(score)) %>%
    dplyr::mutate({{group.by}} := as.factor(.data[[group.by]])) %>%
    dplyr::group_by(mp) %>%
    dplyr::filter(any(score > 0.1))

  cluster_heatmap <-
    cluster_mps %>%
    ggplot(aes(x = .data[[group.by]], y = mp, fill = score)) +
    geom_tile() +
    theme_minimal() +
    scale_fill_gradient(low = "red", high = "yellow", na.value = NA)

  clone_heatmap <- seu@meta.data[,c(names(gene_lists), "scna")] %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(-c("scna", "cell"), names_to = "mp", values_to = "score") %>%
    group_by(scna, mp) %>%
    dplyr::summarize(score = mean(score)) %>%
    dplyr::mutate(scna = as.factor(scna)) %>%
    dplyr::group_by(mp) %>%
    dplyr::filter(any(score > 0.1)) %>%
    ggplot(aes(x = scna, y = mp, fill = score)) +
    geom_tile() +
    theme_minimal()

  cluster_heatmap + clone_heatmap
  mp_score_heatmap_file = glue("results/{mysample}_mp_score_heatmaps.pdf")
  ggsave(mp_score_heatmap_file)


  cluster_mp_genes <- gene_lists[unique(cluster_mps$mp)] %>%
    tibble::enframe("mp", "symbol") %>%
    tidyr::unnest(symbol) %>%
    dplyr::filter(symbol %in% VariableFeatures(seu)) %>%
    dplyr::filter(symbol %in% rownames(seu$gene@scale.data)) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    dplyr::slice_sample(n = 50) %>%
    dplyr::arrange(mp) %>%
    # dplyr::bind_rows(.id = "mp") %>%
    # sample(100) %>%
    identity()

  row_ha = ComplexHeatmap::rowAnnotation(mp = rev(cluster_mp_genes$mp))

  ggplotify::as.ggplot(seu_complex_heatmap(seu, features = cluster_mp_genes$symbol, group.by = c("SCT_snn_res.0.4", "scna"), col_arrangement = c("SCT_snn_res.0.4", "scna"), cluster_rows = FALSE, right_annotation = row_ha)) +
    labs(title = mysample)

  cluster_mp_gene_heatmap_file = glue("results/{mysample}_mp_gene_heatmaps.pdf")
  ggsave(cluster_mp_gene_heatmap_file, height = 8, width = 10)

  return(list("mp_score_heatmap" = mp_score_heatmap_file, "mp_gene_heatmap" = cluster_mp_gene_heatmap_file))

}

pull_subtype_genes <- function(supp_excel = "data/Liu_Radvanyi_2022_supp_data/41467_2021_25792_MOESM6_ESM.xlsx", subtype_chr_preference_file = "results/chromosome_distribution_of_subtype_genes.xlsx"){
  genes_diff_expressed <-
    supp_excel %>%
    readxl::read_excel(sheet = 1, skip = 2) %>%
    janitor::clean_names() %>%
    dplyr::rename(symbol = gene) %>%
    dplyr::filter(gene_cluster %in% c("1.2", "2")) %>%
    dplyr::mutate(gene_cluster = paste0("subtype", str_remove(gene_cluster, "\\..*"))) %>%
    dplyr::left_join(annotables::grch38, by  = "symbol") %>%
    dplyr::arrange(chr) %>%
    dplyr::filter(abs(log_fc_subtype_2_vs_subtype_1) > 0.5) %>%
    dplyr::filter(adjusted_p_value < 0.05) %>%
    dplyr::filter(chr %in% c(1:22, "X")) %>%
    identity()

  list(
    "chr_distribution" = janitor::tabyl(genes_diff_expressed, chr, gene_cluster),
    "all_subtype_genes" = genes_diff_expressed) %>%
    writexl::write_xlsx(subtype_chr_preference_file)

  genes_diff_expressed <-
    genes_diff_expressed %>%
    dplyr::select(gene_cluster, symbol) %>%
    split(.$gene_cluster) %>%
    map(tibble::deframe)


}

compile_subtype_violins <- function(interesting_samples, subtype_violins, selected_samples = c("SRR13884242", "SRR13884243", "SRR13884248", "SRR13884249", "SRR14800534",
                                                                                               "SRR14800535", "SRR14800536", "SRR14800543", "SRR17960484")) {

  names(subtype_violins) <- interesting_samples

  if(!is.null(selected_samples)){
    subtype_violins <- subtype_violins[selected_samples]
  }

  all_subtype_scores <- subtype_violins %>%
    purrr::map(purrr::flatten) %>%
    purrr::imap(~{wrap_plots(.x) + plot_annotation(title = glue("{.y} subtype enrichment"))}) %>%
    identity()

  all_subtype_score_plot_file <- "results/all_subtype_scores.pdf"

  pdf(all_subtype_score_plot_file, height = 12, width = 20)
  print(all_subtype_scores)
  dev.off()

  subtype_scores_by_scna <-
    subtype_violins %>%
    purrr::map(2) %>%
    identity()

  rescale_and_clean_plots <- function(groblist){

    gglist <- map(groblist, ggplotify::as.ggplot)

    gglist <-
      gglist %>%
      map(~(.x + scale_fill_manual(values = scales::hue_pal()(7)))) %>%
      map(~(.x +
              theme(
        legend.position="none",
        axis.title.x = element_blank()) +
          NULL
          )) %>%
      identity()

    return(gglist)
  }

  test0 <-
    subtype_scores_by_scna %>%
    map(rescale_and_clean_plots) %>%
    purrr::imap(~{wrap_plots(.x) + plot_annotation(title = glue("{.y}"))}) %>%
    identity()

  plot_tags <- list(c(rbind(names(test0), rep("", length(names(test0))))))

  wrap_plots(test0, ncol = 2) +
    plot_annotation(tag_levels = plot_tags) &
    theme(
      plot.tag.position = c(1,0),
      plot.tag = element_text(size = 12, hjust = 0, vjust = 0)
      ) +
    NULL

  subtype_by_scna_plot_file = "results/subtype_scores_by_scna.pdf"

  ggsave(subtype_by_scna_plot_file, height = 25, width = 20, limitsize = FALSE)

  # pdf(subtype_by_scna_plot_file)
  # print(subtype_scores_by_scna)
  # dev.off()

  subtype_2_scores_by_scna <- subtype_violins %>%
    purrr::map(purrr::flatten) %>%
    purrr::map(4) %>%
    map(ggplotify::as.ggplot) %>%
    purrr::imap(~{.x + labs(title = glue("{.y} subtype 2 enrichment"))}) %>%
    identity()

  subtype_2_by_scna_plot_file = "results/subtype_2_scores_by_scna.pdf"


  wrap_plots(subtype_2_scores_by_scna, ncol = 3)

  ggsave(subtype_2_by_scna_plot_file, height = 12, width = 12, limitsize = FALSE)

  # pdf(subtype_2_by_scna_plot_file)
  # print(subtype_2_scores_by_scna)
  # dev.off()

  return(list("all" = all_subtype_score_plot_file, "patchwork" = subtype_by_scna_plot_file, "s2" = subtype_2_by_scna_plot_file))

}

#' Retrieve stem cell markers
#'
#' markers from
#'
#' @return
#' @export
#'
#' @examples
pull_stem_cell_markers <- function(){
  smith_markers <- list(
    adult = c(
      "VSNL1", "AKR1B1", "NOTCH4", "TMEM237", "SAMD5",
      "PKD2", "NAP1L1", "PTTG1", "CDK6", "CDCA7", "ACSL4", "HELLS",
      "IKBIP", "PLTP", "TMEM201", "CACHD1", "ILF3", "DNMT1", "USP31",
      "FAM216A", "SLC41A1", "PFKM", "KANK1", "SUPT16H", "ADCY3", "FGD1",
      "PTPN14", "C200rf27", "LGR6", "SLC16A7", "JAM3", "FBL", "NASP",
      "RANBP1", "PRNP", "DSE", "GPX7", "KDELC1", "FCHSD2", "SLCO3A1",
      "CONB11P1", "LOC284023", "NOL9", "NKRE", "NUP107", "RCC2", "ARHGAP25",
      "DDX46", "TCOF1", "GMPS"
    ),
    naive = c(
      "DNMT3L", "ALPPL2", "NLRP7",
      "SLC25A16", "DPPA5", "ATAD3B", "OLAH", "SAMHD1", "PYGB", "TFCP2L1",
      "CP1", "NEFH", "RAB15", "SUSD2", "VSIG10", "PRSS12", "PTPRU",
      "ASRGL1", "A4GALT", "DNAJC15", "CBFA2T2", "KHDC1L", "SLC8B1",
      "SLC35F6", "AARS2", "CDHR1", "SLC25A44", "SLC7A7", "REEP1", "PINK1",
      "GAS7", "KLHL18", "HYAL4", "DCAF4", "TGFBR3", "SLC23A2", "TUBB4A",
      "VAV2", "SLC16A10", "IL6R", "ARPC1B", "MYBL2", "TNS3", "CACNA2D2",
      "ITGAM", "GALNT6", "NDUFAB1", "KIE5", "UPP1", "DACT2"
    ),
    primed = c(
      "SALL2", "DUSP6", "LRRN1", "P3H2", "CDH2", "VRTN", "UCHL1", "STC1", "FGFBP3",
      "NELL2", "ANOS1", "FZD7", "THY1", "PHLDA1", "USP44", "NAP1L3",
      "SPRY4", "PTPRZ1", "EDNRB", "ADAMTS19", "FREM2", "PCYT1B", "TAGLN",
      "CAV1", "COLZA1", "CYP26A1", "MAP7", "PODXL", "GRPR", "NTS",
      "PLCH1", "COL18A1", "PCDH18", "CRABP1", "EPHA2", "VIM", "NECTIN3",
      "GI12", "FAM13A", "DPYSL2", "ATP8A1", "PMEL", "ZIC2", "GPC6",
      "FKBP10", "SEMA3A", "SALL1", "ROR1", "НЕРН", "МЕХЗ"
    )
  )

  return(list("smith" = smith_markers))
}

append_clone_nums <- function(diffex, clone_comparison, seu){
  # browser()
  idents <-
    clone_comparison %>%
    str_extract("[0-9]_v_[0-9]") %>%
    str_split(pattern = "_v_") %>%
    set_names(clone_comparison) %>%
    identity()

  clone_nums <- map_chr(idents, ~{paste(table(seu@meta.data[["clone_opt"]])[.x], collapse = "_v_")}) %>%
    identity()

  diffex0 <-
    diffex %>%
  dplyr::mutate(clone_nums = clone_nums)

  return(diffex0)
}

prep_for_enrichment <- function(df, ...){
  # browser()
  enrich_table <-
    df %>%
    dplyr::select(symbol, avg_log2FC, p_val, clone_nums, clone_comparison) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    tibble::column_to_rownames("symbol") %>%
    enrichment_analysis(...)

  return(enrich_table)

}

plot_enrichment_per_scna <- function(enrichment_table_list, scna, input_plot_file){
	
	scna_plot_file = str_replace(input_plot_file, ".pdf", glue("_{scna}.pdf"))
	
	pdf(scna_plot_file)
	enrichment_table_list %>% 
		purrr::list_flatten() %>% 
		purrr::compact() %>% 
		purrr::discard(~nrow(.x) == 0) %>% 
		map(clusterProfiler::dotplot) %>% 
		imap(~(.x + labs(title = .y))) %>% 
		map(print) %>% 
		identity()
	dev.off()
	
	return(scna_plot_file)
}

plot_enrichment_recurrence <- function(enrichment_tables, num_recur = 2, mytitle = "", n_slice = 10, by_cluster = FALSE, pvalueCutoff = 0.3){
	# browser()
	
	names(enrichment_tables) <- str_replace_all(names(enrichment_tables), "_", "-")
	
	for(sample in names(enrichment_tables)){
		names(enrichment_tables[[sample]]) <- str_replace_all(names(enrichment_tables[[sample]]), "_", "-")
	}
	
	df <-
		enrichment_tables %>%
		purrr::list_flatten() %>%
		purrr::discard(is.na) %>%
		map(~{.x@result}) %>% 
		dplyr::bind_rows(.id = "clone_comparison") %>%
		identity()
	
	test0 <-
		df %>%
		# dplyr::arrange(symbol, sample_id) %>%
		dplyr::arrange(p.adjust) %>%
		dplyr::filter(p.adjust <= pvalueCutoff) %>% 
		dplyr::group_by(Description) %>%
		dplyr::mutate(neg_log10_p_val_adj = -log(p.adjust, base = 10)) %>%
		dplyr::mutate(n_samples = n_distinct(clone_comparison)) %>%
		dplyr::arrange(desc(n_samples), Description)
	
	test0 <- 
		test0 %>% 
		dplyr::select(Description, core_enrichment) %>%
		dplyr::mutate(core_enrichment = str_split(core_enrichment, pattern = "/")) %>%
		tidyr::unnest(core_enrichment) %>%
		dplyr::mutate(core_enrichment = as.integer(core_enrichment)) %>%
		dplyr::left_join(annotables::grch38[c("symbol", "entrez")], by = c("core_enrichment" = "entrez")) %>%
		dplyr::select(Description, symbol) %>%
		dplyr::group_by(Description) %>%
		dplyr::summarize(genes = list(symbol), set_size = n_distinct(symbol)) %>%
		dplyr::rowwise() %>%
		dplyr::mutate(genes = paste(unique(genes), collapse = ",")) %>% 
		dplyr::left_join(test0, by = "Description") %>%
		dplyr::select(-c("set_size", "setSize")) %>%
		identity()
	
	
	test1 <-
		test0 %>%
		dplyr::filter(p.adjust < 0.5) %>%
		dplyr::group_by(Description) %>%
		# dplyr::filter(n_distinct(clone_comparison) >= num_recur) %>%
		dplyr::summarize(mean_NES = mean(NES)) %>%
		# dplyr::slice_max(abs(mean_NES), n = n_slice) %>%
		dplyr::inner_join(test0, by = c("Description")) %>%
		dplyr::mutate(comparison = str_remove(clone_comparison, "SRR[0-9]*_")) %>%
		dplyr::mutate(comparison = factor(str_remove(comparison, "_.*"))) %>%
		# dplyr::mutate(clone_comparison = str_replace_all(clone_comparison, "_", "\n")) %>%
		dplyr::mutate(clone_comparison = factor(clone_comparison)) %>%
		dplyr::mutate(clone_comparison = fct_reorder(clone_comparison, as.integer(comparison))) %>%
		dplyr::mutate(genes = str_split(genes, ",")) %>%
		# dplyr::mutate(phase = str_extract(clone_comparison, "(?<=_).*(?=-[0-9]*)")) %>% 
		# dplyr::mutate(sample_id = str_extract(clone_comparison, ".*(?=_)")) %>% 
		# dplyr::mutate(sample_id = factor(sample_id)) %>% 
		identity()
	
	color_scale_lim <- ceiling(max(abs(test1$NES)))
	
	color_scale_lim <- ifelse(is.infinite(color_scale_lim), 2, color_scale_lim)
	
	enrich_plot <-
		test1 %>% 
		ggplot(aes(x = clone_comparison, y = Description, size = neg_log10_p_val_adj, color = NES)) + # by_cluster
		scale_color_gradient2(breaks=seq((-color_scale_lim), color_scale_lim, 1), limits=c((-color_scale_lim-0.5),(color_scale_lim+0.5))) +
		# scale_color_gradient2() +
		geom_point() +
		labs(title = mytitle, size  = "-log10 p.adj") +
		theme(
			axis.title.x = element_blank(),
			axis.title.y = element_blank(),
			axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) 
		) +
		scale_y_discrete(labels = function(x) stringr::str_wrap(str_replace_all(x, "_", " "), width = 20)) +
		# scale_y_discrete(labels = function(x) str_wrap(x, width = 4, whitespace_only = FALSE)) +
		labs(color = "enrichment") +
		# facet_wrap(~comparison) +
		NULL
	
	# upset_plot <-
	#   test1 %>%
	#   ggplot(aes(x=genes, y = Description)) +
	#   geom_point() +
	#   scale_x_upset() +
	#   theme_minimal() +
	#   labs(title = glue("{mytitle} Gene set overlap"))
	
	return(list("table" = test0, "enrich_plot" = enrich_plot))
}

plot_enrichment_recurrence_by_cluster <- function(enrichment_tables, num_recur = 2, mytitle = "", n_slice = 10, by_cluster = TRUE, pvalueCutoff = 0.3){
  # browser()

	names(enrichment_tables) <- str_replace_all(names(enrichment_tables), "_", "-")
	
	for(sample in names(enrichment_tables)){
		names(enrichment_tables[[sample]]) <- str_replace_all(names(enrichment_tables[[sample]]), "_", "-")
	}
	
  df <-
    enrichment_tables %>%
    purrr::list_flatten() %>%
    purrr::discard(is.na) %>%
    map(~{.x@result}) %>% 
    dplyr::bind_rows(.id = "clone_comparison") %>%
    identity()

  test0 <-
    df %>%
    # dplyr::arrange(symbol, sample_id) %>%
    dplyr::arrange(p.adjust) %>%
  	dplyr::filter(p.adjust <= pvalueCutoff) %>% 
    dplyr::group_by(Description) %>%
    dplyr::mutate(neg_log10_p_val_adj = -log(p.adjust, base = 10)) %>%
    dplyr::mutate(n_samples = n_distinct(clone_comparison)) %>%
    dplyr::arrange(desc(n_samples), Description)
  
  test0 <- 
  	test0 %>% 
    dplyr::select(Description, core_enrichment) %>%
    dplyr::mutate(core_enrichment = str_split(core_enrichment, pattern = "/")) %>%
    tidyr::unnest(core_enrichment) %>%
    dplyr::mutate(core_enrichment = as.integer(core_enrichment)) %>%
    dplyr::left_join(annotables::grch38[c("symbol", "entrez")], by = c("core_enrichment" = "entrez")) %>%
    dplyr::select(Description, symbol) %>%
    dplyr::group_by(Description) %>%
    dplyr::summarize(genes = list(symbol), set_size = n_distinct(symbol)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(genes = paste(unique(genes), collapse = ",")) %>% 
    dplyr::left_join(test0, by = "Description") %>%
    dplyr::select(-c("set_size", "setSize")) %>%
    identity()

  
  test1 <-
    test0 %>%
    dplyr::filter(p.adjust < 0.5) %>%
    dplyr::group_by(Description) %>%
    # dplyr::filter(n_distinct(clone_comparison) >= num_recur) %>%
    dplyr::summarize(mean_NES = mean(NES)) %>%
    # dplyr::slice_max(abs(mean_NES), n = n_slice) %>%
    dplyr::inner_join(test0, by = c("Description")) %>%
    dplyr::mutate(comparison = str_remove(clone_comparison, "SRR[0-9]*_")) %>%
    dplyr::mutate(comparison = factor(str_remove(comparison, "_.*"))) %>%
    # dplyr::mutate(clone_comparison = str_replace_all(clone_comparison, "_", "\n")) %>%
    dplyr::mutate(clone_comparison = factor(clone_comparison)) %>%
    dplyr::mutate(clone_comparison = fct_reorder(clone_comparison, as.integer(comparison))) %>%
    dplyr::mutate(genes = str_split(genes, ",")) %>%
  	dplyr::mutate(phase = str_extract(clone_comparison, "(?<=_).*(?=-[0-9]*)")) %>% 
  	dplyr::mutate(sample_id = str_extract(clone_comparison, ".*(?=_)")) %>% 
  	dplyr::mutate(sample_id = factor(sample_id)) %>% 
    identity()
  
  phase_levels = c("g1", "g1-s", "s", "s-g2", "g2", "g2-m", "pm", "hsp", "hypoxia", "other", "s-2")
  
  phase_levels <- phase_levels[phase_levels %in% unique(test1$phase)]

  color_scale_lim <- ceiling(max(abs(test1$NES)))

  color_scale_lim <- ifelse(is.infinite(color_scale_lim), 2, color_scale_lim)

  enrich_plot_by_phase <-
  	test1 %>% 
  	dplyr::mutate(phase = factor(phase, levels = phase_levels)) %>% # by_cluster
  	dplyr::mutate(phase = as.numeric(phase)) %>% # by_cluster
    ggplot(aes(x = fct_reorder(clone_comparison, phase), y = Description, size = neg_log10_p_val_adj, color = NES)) + # by_cluster
    scale_color_gradient2(breaks=seq((-color_scale_lim), color_scale_lim, 1), limits=c((-color_scale_lim-0.5),(color_scale_lim+0.5))) +
    # scale_color_gradient2() +
    geom_point() +
    labs(title = mytitle, size  = "-log10 p.adj") +
    theme(
    	axis.title.x = element_blank(),
    	axis.title.y = element_blank(),
    	axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) 
    ) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
    labs(color = "enrichment") +
    # facet_wrap(~comparison) +
    NULL
  
  enrich_plot_by_sample <-
  	test1 %>% 
  	dplyr::mutate(sample_id = factor(sample_id)) %>% # by_cluster
  	dplyr::mutate(sample_id = as.numeric(sample_id)) %>% # by_cluster
  	ggplot(aes(x = fct_reorder(clone_comparison, sample_id), y = Description, size = neg_log10_p_val_adj, color = NES)) + # by_cluster
  	scale_color_gradient2(breaks=seq((-color_scale_lim), color_scale_lim, 1), limits=c((-color_scale_lim-0.5),(color_scale_lim+0.5))) +
  	# scale_color_gradient2() +
  	geom_point() +
  	labs(title = mytitle, size  = "-log10 p.adj") +
  	theme(
  		axis.title.x = element_blank(),
  		axis.title.y = element_blank(),
  		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1) 
  	) +
  	scale_y_discrete(labels = function(x) str_wrap(x, width = 40)) +
  	labs(color = "enrichment") +
  	# facet_wrap(~comparison) +
  	NULL

  # upset_plot <-
  #   test1 %>%
  #   ggplot(aes(x=genes, y = Description)) +
  #   geom_point() +
  #   scale_x_upset() +
  #   theme_minimal() +
  #   labs(title = glue("{mytitle} Gene set overlap"))

  return(list("table" = test0, "enrich_plot_by_phase" = enrich_plot_by_phase, "enrich_plot_by_sample" = enrich_plot_by_sample))
}

read_mps <- function(mp_file = "/dataVolume/storage/Homo_sapiens/3ca/ITH_hallmarks/MPs_distribution/MP_list.RDS"){
  mps <- readRDS(mp_file)

  mps <- purrr::map(mps, ~purrr::set_names(.x, janitor::make_clean_names(names(.x))))
}

tally_num_diffex <- function(oncoprint_input_by_scna_unfiltered){
  symbol_tally <- purrr::map_depth(oncoprint_input_by_scna_unfiltered, 2, ~{table(.x$symbol)})
  sample_tally <- purrr::map_depth(oncoprint_input_by_scna_unfiltered, 2, ~{table(.x$sample_id)})
  return(list("symbol" = symbol_tally, "sample" = sample_tally))
}

annotate_cluster_membership <- function(diffex_1, diffex_2, new_col_name){

  diffex_2 <-
    diffex_2 %>%
    dplyr::ungroup() %>%
    dplyr::select(clone_comparison, symbol) %>%
    dplyr::mutate({{new_col_name}} := TRUE) %>%
    identity()

  diffex_1 <-
    diffex_1 %>%
    dplyr::ungroup() %>%
    dplyr::left_join(diffex_2, by = c("clone_comparison", "symbol"))

  return(diffex_1)
}


tally_kooi_candidates <- function(cis_diffex_clones = "results/diffex_bw_clones_large_in_segment_by_chr.xlsx", trans_diffex_clones = "results/diffex_bw_clones_trans_by_chr.xlsx"){
  # browser()

  cis_diffex_clones <-
    cis_diffex_clones %>%
    excel_sheets() %>%
    set_names() %>%
    map(read_excel, path = cis_diffex_clones) %>%
    map(dplyr::filter, !is.na(kooi_region))

  cis_diffex_clones <-
  cis_diffex_clones[map_lgl(cis_diffex_clones, ~(nrow(.x) > 0))] %>%
    dplyr::bind_rows(.id = "chr")

  trans_diffex_clones <-
    trans_diffex_clones %>%
    excel_sheets() %>%
    set_names() %>%
    map(read_excel, path = trans_diffex_clones) %>%
    map(dplyr::filter, !is.na(kooi_region))

  trans_diffex_clones <-
    trans_diffex_clones[map_lgl(trans_diffex_clones, ~(nrow(.x) > 0))] %>%
    dplyr::bind_rows(.id = "chr")

  return(list("cis" = cis_diffex_clones, "trans" = trans_diffex_clones))

}


collect_study_metadata <- function() {
  # browser()
  retrieve_cell_stats <- function(seu_path){

    seu <- readRDS(seu_path)

    stats = seu@meta.data[c("nCount_gene", "nFeature_gene", "percent.mt")] %>%
      tibble::rownames_to_column("cell")

    return(stats)
  }

  seus <-
    dir_ls("output/seurat/", regexp = "\\/SRR[0-9]*_seu.rds") %>%
    set_names(str_extract(., "SRR[0-9]*"))

# collin ------------------------------

  collin_cell_stats <- seus[c("SRR13633759", "SRR13633760", "SRR13633761", "SRR13633762")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # field ------------------------------

  field_cell_stats <- seus[c("SRR17960480", "SRR17960481", "SRR17960482", "SRR17960483", "SRR17960484")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # wu ------------------------------

  wu_cell_stats <- seus[c("SRR13884240", "SRR13884241", "SRR13884242", "SRR13884243", "SRR13884244", "SRR13884245", "SRR13884246", "SRR13884247", "SRR13884248", "SRR13884249")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")

  # yang ------------------------------

  yang_cell_stats <- seus[c("SRR14800534", "SRR14800535", "SRR14800536", "SRR14800537", "SRR14800538", "SRR14800539", "SRR14800540", "SRR14800541", "SRR14800542", "SRR14800543")] %>%
    map_dfr(retrieve_cell_stats, .id = "sample_id")
  
  liu_cell_stats <- seus[c("SRR27187899", "SRR27187900", "SRR27187901", "SRR27187902")]  %>%
  	map_dfr(retrieve_cell_stats, .id = "sample_id")

  # combined ------------------------------

  study_cell_stats <- dplyr::bind_rows(list("collin" = collin_cell_stats, "field" = field_cell_stats, "wu" = wu_cell_stats, "yang" = yang_cell_stats, "liu" = liu_cell_stats), .id = "study")

  write_csv(study_cell_stats, "results/study_cell_stats.csv")

  return(study_cell_stats)

}

plot_study_cell_stats <- function(study_cell_stats, cell_stats_plot_file, umi_threshold = 1e3, genes_threshold = 1e3, mito_threshold = 5, mito_expansion = 0.8, bandwidth = 0.35, plot_height = NULL, ...) {
	# browser()
	
	study_cell_stats <- 
		study_cell_stats |> 
		dplyr::mutate(study = dplyr::case_when(
			study == "collin" ~ "Collin et al. 2021", 
			study == "wu" ~ "Wu et al. 2022", 
			study == "yang" ~ "Yang et al. 2021", 
			study == "field" ~ "Field et al. 2022", 
			study == "liu" ~ "Liu et al. 2024"
		)) |> 
		dplyr::mutate(study = factor(study, levels = c("Collin et al. 2021", "Wu et al. 2022", "Yang et al. 2021", "Field et al. 2022", "Liu et al. 2024")))
	
	mypal <- scales::hue_pal()(5) %>%
		set_names(levels(study_cell_stats$study))
	
	umis_per_cell <- ggplot(study_cell_stats, aes(x = nCount_gene, y = sample_id, fill = study, group = sample_id)) +
		geom_density_ridges(scale=1, panel_scaling = FALSE, bandwidth = bandwidth) +
		scale_x_log10() +
		scale_y_discrete(limits=rev) +
		scale_fill_manual(values = mypal) +
		theme(
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),  #remove y axis labels,
			# axis.text.y=element_blank(),  #remove y axis labels
			# axis.ticks.y=element_blank()  #remove y axis ticks
		) +
		geom_vline(xintercept = umi_threshold, linetype="dotted") +
		labs(title = "UMIs/cell")
	
	
	genes_per_cell <- ggplot(study_cell_stats, aes(x = nFeature_gene, y = sample_id, fill = study, group = sample_id)) +
		geom_density_ridges(scale=1, panel_scaling = FALSE, bandwidth = bandwidth) +
		scale_x_log10() +
		scale_y_discrete(limits=rev) +
		scale_fill_manual(values = mypal) +
		theme(
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),  #remove y axis labels,
			# axis.text.y=element_blank(),  #remove y axis labels
			# axis.ticks.y=element_blank()  #remove y axis ticks
		) +
		geom_vline(xintercept = genes_threshold, linetype="dotted") +
		labs(title = "genes/cell")
	
	percent_mito_per_cell <- ggplot(study_cell_stats, aes(x = `percent.mt`, y = sample_id, fill = study, group = sample_id)) +
		geom_density_ridges(scale=1, panel_scaling = FALSE, bandwidth = bandwidth) +
		scale_y_discrete(limits=rev, expand = expansion(add = c(0.55, mito_expansion))) +
		# scale_y_discrete(limits=rev) +
		scale_fill_manual(values = mypal) +
		theme(
			axis.title.x=element_blank(),
			axis.title.y=element_blank(),  #remove y axis labels,
			# axis.text.y=element_blank(),  #remove y axis labels
			# axis.ticks.y=element_blank()  #remove y axis ticks
		) +
		geom_vline(xintercept = mito_threshold, linetype="dotted") +
		labs(title = "% mito/cell") +
		xlim(0, 25) +
		NULL
	
	retention_labels <- study_cell_stats %>% 
		dplyr::select(study, sample_id, exclusion_criteria) %>% 
		dplyr::distinct() %>% 
		dplyr::mutate(exclusion_criteria = replace_na(exclusion_criteria, "retained")) %>% 
		dplyr::mutate(sample_id = factor(sample_id)) %>% 
		identity()
	
	retention_labels_plot <-
		retention_labels %>% 
		ggplot(aes(y = sample_id, x = 1, fill = exclusion_criteria)) + 
		geom_tile(color = "black") +
		scale_y_discrete(limits=rev) +
		theme_void() 
	
	list(umis_per_cell, genes_per_cell, percent_mito_per_cell, retention_labels_plot) %>%
		wrap_plots() +
		plot_layout(axes = "collect", guides = 'collect', widths = c(3,3,3,0.2))
	
	plot_height <- ifelse(is.null(plot_height), 0.33*n_distinct(study_cell_stats$sample_id), plot_height)
	
	ggsave(cell_stats_plot_file, height = plot_height, ...)
}


plot_study_metadata <- function(study_cell_stats, ...) {
	
	# study_cell_stats <- read_csv("results/study_cell_stats.csv")
	
	plot_mt_v_nUMI <- function(study_cell_stats, cell_stats_plot_file){
		# percent_mito_per_read <-
		study_cell_stats %>%
			dplyr::filter(nCount_gene > 1000) %>%
			# dplyr::mutate(mito_per_read = nCount_gene/`percent.mt`) %>%
			ggplot(aes(x = `percent.mt`, y = nCount_gene, group = sample_id)) +
			geom_hex(bins=70) +
			facet_wrap(~sample_id) +
			geom_vline(xintercept = 5, linetype="dotted") +
			geom_hline(yintercept = 1e3, linetype="dotted") +
			# scale_y_discrete(limits=rev, expand = expansion(add = c(0.55, mito_expansion))) +
			labs(title = "percent mito per cell") +
			scale_y_continuous(limits = c(0,1e5)) +
			xlim(0, 50) +
			ylim(0, 2e5) +
			annotate("rect", xmin = 0, xmax = 5, ymin = 1e3, ymax = 2e5,
							 alpha = .2, color = "yellow") +
			NULL
		
		ggsave(cell_stats_plot_file, height = 6, width = 8)
	}
	
	normal_ctrl_samples <- unlist(list(
		"collin" = c("SRR13633761", "SRR13633762")
	))
	
	bad_qc_sample_ids <- list(
		"collin" = c("SRR13633759", "SRR13633760"),
		"yang" = c("SRR14800538", "SRR14800539", "SRR14800540", "SRR14800541", "SRR14800542", "SRR14800543"),
		"field" = c("SRR17960480", "SRR17960482", "SRR17960484")
	) %>% 
		enframe("study", "sample_id") %>%
		unnest(sample_id) %>%
		identity()

	bad_scna_sample_ids <- list(
		"collin" = c(),
		"yang" = c("SRR14800537"),
		"field" = c("SRR17960483"),
		"wu" = c("SRR13884240", "SRR13884241", "SRR13884244", "SRR13884245")
	) %>% 
		enframe("study", "sample_id") %>% 
		unnest(sample_id) %>% 
		identity()
	
	excluded_sample_ids <- 
		list("bad_qc" = bad_qc_sample_ids, 
				 "bad_scna" = bad_scna_sample_ids) %>% 
		dplyr::bind_rows(.id = "exclusion_criteria")
		
	unfiltered_cell_stats_plot_file <- "results/unfiltered_study_stats_mt_v_nUMI.pdf"
	
	study_cell_stats <- 
		study_cell_stats %>% 
		dplyr::full_join(excluded_sample_ids, by = c("study", "sample_id")) %>% 
		identity()
	
	study_cell_stats %>%
		dplyr::filter(!sample_id %in% normal_ctrl_samples) %>%
		plot_mt_v_nUMI(unfiltered_cell_stats_plot_file)
	
	unfiltered_cell_stats_plot_file <- "results/unfiltered_study_stats.pdf"
	
	study_cell_stats %>%
		dplyr::filter(!sample_id %in% normal_ctrl_samples) %>%
		plot_study_cell_stats(unfiltered_cell_stats_plot_file, mito_expansion = 0.8, ...)
		
	
	# good_qc_study_cell_stats <-
	# 	study_cell_stats %>%
	# 	dplyr::filter(!sample_id %in% bad_qc_samples)
	# 
	# good_qc_cell_stats_plot_file <- "results/good_qc_study_stats.pdf"
	# plot_study_cell_stats(good_qc_study_cell_stats, good_qc_cell_stats_plot_file, mito_expansion = 1)
	# 
	# 
	# # retained ------------------------------
	# retained_study_cell_stats <-
	# 	study_cell_stats %>%
	# 	dplyr::filter(!sample_id %in% c(bad_scna_samples, bad_qc_samples))
	# 
	# retained_cell_stats_plot_file <- "results/retained_study_stats.pdf"
	# plot_study_cell_stats(retained_study_cell_stats, retained_cell_stats_plot_file, mito_expansion = 0.55)
	
	return(unfiltered_cell_stats_plot_file)
	
}

montage_images <- function(plot_files, sample_id, numbat_dir = "numbat_sridhar", tile = '6', label = "montage"){
  # browser()
  plot_images <- magick::image_read(plot_files, density = 600)

  my_montage <- magick::image_montage(plot_images, tile = tile, geometry='800x', shadow = FALSE)

  # my_montage <- image_montage(plot_images, geometry = c('x200+10+10', 'x800+10+10', 'x100+10+10'), tile = '3x', shadow = FALSE)

  montage_path <- glue("results/{numbat_dir}/{sample_id}_{label}.pdf")

  image_write(my_montage, format = "pdf", montage_path)

  return(montage_path)
}

make_volcano_plots <- function(myres, mysubtitle, sample_id){
	
	diffex_comparison = str_split(unique(myres[["diffex_comparison"]]), "_", simplify = TRUE)
	
	right_label = diffex_comparison[[1]] %||% "right"
	left_label = diffex_comparison[[2]] %||% "left"
	
  myres <-
    myres %>%
    dplyr::mutate(chr = case_when(chr == "X" ~ "23",
                                  chr == "Y" ~ "24",
                                  TRUE ~ as.character(chr))) %>%
    dplyr::mutate(chr = str_pad(chr, side = "left", pad = "0", width = 2)) %>%
    dplyr::mutate(clone_comparison = str_replace_all(clone_comparison, "_", " ")) %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::mutate(rownames = symbol) %>%
    tibble::column_to_rownames("rownames")

  mytitle = sample_id
  mysubtitle = mysubtitle

  ref_var <-
    myres$chr %>%
    set_names(.)

  chrs <- str_pad(as.character(1:24), side = "left", pad = "0", width = 2)

  mypal <- scales::hue_pal()(length(chrs))
  names(mypal) <- chrs

  custom_cols <- mypal[ref_var]

  FCcutoff = summary(abs(myres$avg_log2FC))[[5]]

  selected_genes <-
    myres %>%
    dplyr::filter(abs(avg_log2FC) > 0.05, p_val_adj < 0.1) %>%
    dplyr::pull(symbol)


  myplot <- EnhancedVolcano(myres,
                            lab = rownames(myres),
                            selectLab = selected_genes,
                            labSize = 4,
                            x = 'avg_log2FC',
                            y = 'p_val_adj',
                            FCcutoff = FCcutoff,
                            pCutoff = 5e-2,
                            colCustom = custom_cols,
                            max.overlaps = 25) +
    aes(color = chr) +
    # facet_wrap(~chr) +
    labs(title = mytitle, subtitle = mysubtitle)
  
  layer_scales(myplot)$y$range$range
  
  # plot_ymax = max(-log(myplot$data$p_val_adj, base = 10))
  plot_ymax = max(ggplot_build(myplot)$layout$panel_params[[1]]$y.range)
  
  myplot <- 
  myplot +
  	annotation_custom(
  		text_grob(
  			left_label,
  			size= 13,
  			color = "red",
  			face = "bold"),
  		xmin=-Inf,
  		xmax=-Inf,
  		ymin=plot_ymax+plot_ymax/10, 
  		ymax=plot_ymax+plot_ymax/10) +
  	annotation_custom(
  		text_grob(
  			right_label,
  			size= 13,
  			color = "red",
  			face = "bold"),
  		xmin=Inf,
  		xmax=Inf,
  		ymin=plot_ymax+plot_ymax/10, 
  		ymax=plot_ymax+plot_ymax/10) +
  	theme(plot.margin = unit(c(1,3,1,1), "lines")) +
  	coord_cartesian(clip = "off") + 
  	NULL
  	

  return(myplot)
}

score_samples_for_rod_enrichment <- function(numbat_rds_file){

  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  # count number of rods by RHO, ROM1, GNAT1, NR2E3?------------------------------

  seu <- AddModuleScore(seu, list("rod" = c("RHO", "ROM1", "GNAT1", "NR2E3")), name = "rod")

  rod_meta <-
    seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(cell, cluster = gene_snn_res.0.2, rod1) %>%
    dplyr::mutate(rod_identity = ifelse(rod1 > 1.5, "rod", "cell")) %>%
    identity()

  seu@meta.data["rod_identity"] <-
    rod_meta %>%
    dplyr::pull(rod_identity)

  # Seurat::FeaturePlot(seu, features = "rod1", split.by = "gene_snn_res.0.2")
  rod_score_plot <- Seurat::FeaturePlot(seu, features = "rod1")

  rod_id_plot <- Seurat::DimPlot(seu, group.by = "rod_identity")

  percent_rod_cells <-
    janitor::tabyl(rod_meta, rod_identity) %>%
    dplyr::filter(rod_identity == "rod") %>%
    dplyr::pull(percent) %>%
    round(2) %>%
    identity()

  rod_patch <- rod_score_plot + rod_id_plot +
    plot_annotation(title = sample_id, subtitle = scales::label_percent()(percent_rod_cells))

  dir_create("results/rod_plots")
  rod_plot_path <- glue("results/rod_plots/{sample_id}.pdf")

  ggsave(rod_plot_path)

  return(list("sample_id" = sample_id, "percent_rod_cells" = (scales::label_percent()(percent_rod_cells)), "rod_rich" = ifelse(percent_rod_cells > 0.05, 1, 0)))

}

score_whole_pseudobulks <- function(numbat_rds_file, subtype_markers){
  sample_id <- str_extract(numbat_rds_file, "SRR[0-9]*")

  numbat_dir = fs::path_split(numbat_rds_file)[[1]][[2]]

  seu <- readRDS(glue("output/seurat/{sample_id}_seu.rds")) %>%
    filter_sample_qc()

  bulk_expression <- GetAssayData(seu, slot = "data") %>%
    rowSums()

  subtype1_expression <- bulk_expression[names(bulk_expression) %in% subtype_markers$subtype1]

  subtype2_expression <- bulk_expression[names(bulk_expression) %in% subtype_markers$subtype2]

  return(list("sample_id" = sample_id, "s1" = mean(subtype1_expression), "s2" = mean(subtype2_expression)))


}

derive_pseudobulk_subtype_scores  <- function(seu_paths){
  # browser()
  pull_assay_data <- function(seu_path){

    seu <- readRDS(seu_path) %>%
      filter_sample_qc()

    bulk_data <- GetAssayData(seu, slot = "data") %>%
      rowSums()

    bulk_counts <- GetAssayData(seu, slot = "counts") %>%
      rowSums()

    return(list("counts" = bulk_counts, "data" = bulk_data))
  }

  bulk_assays <-
    seu_paths %>%
    set_names(str_extract(., "SRR[0-9]*")) %>%
    map(pull_assay_data)

  bulk_assay_counts <-
    bulk_assays %>%
    map("counts") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    tibble::column_to_rownames("sample_id") %>%
    as.matrix()

  bulk_assay_counts[is.na(bulk_assay_counts)] <- 0

  scaled_counts_datas <- prop.table(t(bulk_assay_counts), 2)

  mydend <- as.dendrogram(hclust(dist(t(scaled_counts_datas)), method = "ward.D2"))

  dend_plot <- ggplot(dendextend::as.ggdend(mydend))

  dend_groups <- dendextend:::cutree.dendrogram(mydend, 2) %>%
    tibble::enframe("sample_id", "group") %>%
    tibble::column_to_rownames("sample_id") %>%
    dplyr::mutate(group = factor(group))

  dds <- DESeqDataSetFromMatrix(t(bulk_assay_counts), dend_groups, ~group)

  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep,]

  dds <- DESeq(dds)
  res <- results(dds)
  res

  res_table <-
    res %>%
    as.data.frame() %>%
    tibble::rownames_to_column("symbol") %>%
    dplyr::inner_join(annotables::grch38, by = "symbol") %>%
    dplyr::arrange(log2FoldChange) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    dplyr::select(symbol, description, everything()) %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::distinct(entrez, .keep_all = TRUE)

  # we want the log2 fold change
  original_gene_list <- res_table[["log2FoldChange"]]

  # name the vector
  names(original_gene_list) <- res_table$entrez

  # omit any NA values
  gene_list<-na.omit(original_gene_list)

  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)

  gse <- clusterProfiler::gseGO(geneList=gene_list,
                                ont = "BP",
                                OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                keyType = "ENTREZID",
                                minGSSize = 3,
                                maxGSSize = 800,
                                pvalueCutoff = 0.05,
                                verbose = TRUE,
                                pAdjustMethod = "BH") %>%
    clusterProfiler::simplify() # for GSEGO

  gse_table <-
    dplyr::arrange(gse@result, NES)

  return(list("dend" = dend_plot, "diffex" = res_table, "enrichment" = gse_table))


}

plot_effect_of_filtering <- function(unfiltered_seu_path, filtered_seu_path, group.by = "SCT_snn_res.0.6"){
  # browser()

  plot_list = list()

  sample_id <- str_extract(unfiltered_seu_path, "SRR[0-9]*")

  unfiltered_seu <- readRDS(unfiltered_seu_path)

  filtered_seu <- readRDS(filtered_seu_path)

  fs::dir_create("results/effect_of_filtering")

  # distribution ------------------------------
  dir_create("results/effect_of_filtering/distribution")
  plot_distribution_of_clones_across_clusters(filtered_seu, sample_id, var_x = "scna", var_y = group.by)
  fs::dir_create(glue("results/effect_of_filtering/distribution/filtered/"))
  plot_path = glue("results/effect_of_filtering/distribution/filtered/{sample_id}_filtered_distribution.pdf")
  plot_list["filtered_distribution"] = plot_path
  ggsave(plot_path, height = 4, width = 8)

  filtered_dist_tables <- table_distribution_of_clones_across_clusters(filtered_seu, sample_id, clusters = group.by)

  table_path = glue("results/effect_of_filtering/distribution/filtered/{sample_id}_filtered_distribution.xlsx")
  plot_list["filtered_distribution_tables"] = table_path
  writexl::write_xlsx(filtered_dist_tables, table_path)

  plot_distribution_of_clones_across_clusters(unfiltered_seu, sample_id, var_y = group.by)
  fs::dir_create(glue("results/effect_of_filtering/distribution/filtered"))
  plot_path = glue("results/effect_of_filtering/distribution/filtered/{sample_id}_unfiltered_distribution.pdf")
  plot_list["unfiltered_distribution"] = plot_path
  ggsave(plot_path, height = 4, width = 8)

  unfiltered_dist_tables <- table_distribution_of_clones_across_clusters(unfiltered_seu, sample_id, clusters = group.by)

  table_path = glue("results/effect_of_filtering/distribution/filtered/{sample_id}_unfiltered_distribution.xlsx")
  plot_list["unfiltered_distribution_tables"] = table_path
  writexl::write_xlsx(unfiltered_dist_tables, table_path)

  # abbreviation markers ------------------------------
  dir_create("results/effect_of_filtering/abbreviation")
  (plot_markers(filtered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
     labs(title = "filtered")) +
    (plot_markers(unfiltered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
       labs(title = "unfiltered")) +
    plot_annotation(title = sample_id)

  plot_path = glue("results/effect_of_filtering/abbreviation/{sample_id}_abbreviation_markers.pdf")
  plot_list["abbreviation_markers"] = plot_path
  ggsave(plot_path, height = 12, width = 15)

  filtered_marker_tables <- table_cluster_markers(filtered_seu)

  table_path = glue("results/effect_of_filtering/abbreviation/{sample_id}_filtered_markers.xlsx")
  plot_list["filtered_marker_tables"] = table_path
  writexl::write_xlsx(filtered_marker_tables, table_path)

  unfiltered_marker_tables <- table_cluster_markers(unfiltered_seu) %>%
    purrr::compact()

  table_path = glue("results/effect_of_filtering/abbreviation/{sample_id}_unfiltered_markers.xlsx")
  plot_list["unfiltered_marker_tables"] = table_path
  writexl::write_xlsx(unfiltered_marker_tables, table_path)

  heatmap_features  <-
    table_cluster_markers(unfiltered_seu) %>%
    pluck(group.by) %>%
    group_by(Cluster) %>%
    slice_head(n=10) %>%
    dplyr::pull(Gene.Name) %>%
    identity()

  ggplotify::as.ggplot(
      seu_complex_heatmap(unfiltered_seu,
                          features = heatmap_features,
                          group.by = c(group.by, "Phase", "scna"),
                          col_arrangement = c(group.by, "Phase", "scna"),
                          cluster_rows = FALSE)) +
    labs(title = sample_id)

  plot_path = glue("results/effect_of_filtering/abbreviation/{sample_id}_abbreviation_heatmap.pdf")
  plot_list["abbreviation_heatmap"] = plot_path
  ggsave(plot_path, height = 8, width = 8)

  # nCount_gene umaps ------------------------------

  unfiltered_seu$log_nCount_gene <- log1p(unfiltered_seu$nCount_gene)
  filtered_seu$log_nCount_gene <- log1p(filtered_seu$nCount_gene)

  (FeaturePlot(unfiltered_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
     labs(title = "unfiltered")) +
  (FeaturePlot(filtered_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
     labs(title = "filtered")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_filtering/nCount_gene"))
  plot_path = glue("results/effect_of_filtering/nCount_gene/{sample_id}_nCount_gene_umaps.pdf")
  plot_list["nCount_gene_umaps"] = plot_path
  ggsave(plot_path, height = 8, width = 10)

  # abbreviation umaps ------------------------------
  mycols = scales::hue_pal()(length(unique(unfiltered_seu@meta.data[["abbreviation"]])))

  (DimPlot(unfiltered_seu, group.by = "abbreviation", cols = mycols) +
      labs(title = "unfiltered")) +
    (DimPlot(filtered_seu, group.by = "abbreviation", cols = mycols) +
       labs(title = "filtered")) +
    # (DimPlot(regressed_seu, group.by = "gene_snn_res.0.2") +
    #   labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_filtering/abbreviation"))
  plot_path = glue("results/effect_of_filtering/abbreviation/{sample_id}_abbreviation_umaps.pdf")
  plot_list["abbreviation_umaps"] = plot_path
  ggsave(plot_path, height = 8, width = 10)

  # scna umaps ------------------------------
  mycols = scales::hue_pal()(length(unique(unfiltered_seu@meta.data[["scna"]])))

  unfiltered_seu@meta.data$scna <- vec_split_label_line(unfiltered_seu@meta.data$scna, 3)
  filtered_seu@meta.data$scna <- vec_split_label_line(filtered_seu@meta.data$scna, 3)

  (DimPlot(unfiltered_seu, group.by = "scna", cols = mycols) +
      labs(title = "unfiltered")) +
  (DimPlot(filtered_seu, group.by = "scna", cols = mycols) +
      labs(title = "filtered")) +
  # (DimPlot(regressed_seu, group.by = "gene_snn_res.0.2") +
  #     labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_filtering/scna"))
  plot_path = glue("results/effect_of_filtering/scna/{sample_id}_scna_umaps.pdf")
  plot_list["scna_umaps"] = plot_path
  ggsave(plot_path, height = 8, width = 10)

  return(plot_list)

}

plot_effect_of_regression <- function(filtered_seu_path, regressed_seu_path, resolution = 0.4){
	# browser()
	
	sample_id <- str_extract(filtered_seu_path, "SRR[0-9]*")
	
	regressed_seu <- readRDS(regressed_seu_path)
	
	regressed_meta <- regressed_seu@meta.data %>% 
		tibble::rownames_to_column("cell") %>% 
		dplyr::select(cell, starts_with("SCT_snn_res")) %>% 
		tibble::column_to_rownames("cell") %>% 
		dplyr::rename_with(~ str_replace(.x, "SCT_snn_res.", "regressed."),
											 starts_with("SCT_")) %>% 
		identity()
	
	filtered_seu <- readRDS(filtered_seu_path)
	
	filtered_seu <- AddMetaData(filtered_seu, regressed_meta)
	
	filtered_seu <- find_all_markers(filtered_seu, colnames(regressed_meta))
	
	regressed_cluster = glue("regressed.{resolution}")
	
	filtered_cluster = glue("SCT_snn_res.{resolution}")
	
	regressed_features  <-
		filtered_seu@misc$markers[[regressed_cluster]]$presto %>%
		group_by(Cluster) %>%
		slice_head(n=10) %>%
		dplyr::select(regressed_cluster = Cluster, Gene.Name) %>%
		dplyr::ungroup() %>% 
		dplyr::distinct(Gene.Name, .keep_all = TRUE) %>% 
		dplyr::filter(Gene.Name %in% rownames(GetAssayData(filtered_seu, "SCT", "scale.data"))) %>% 
		identity()
	
	filtered_features  <-
		filtered_seu@misc$markers[[filtered_cluster]]$presto %>%
		group_by(Cluster) %>%
		slice_head(n=10) %>%
		dplyr::select(Cluster, Gene.Name) %>%
		dplyr::mutate(filtered_cluster = Cluster) %>% 
		dplyr::ungroup() %>% 
		dplyr::distinct(Gene.Name, .keep_all = TRUE) %>% 
		dplyr::filter(Gene.Name %in% rownames(GetAssayData(filtered_seu, "SCT", "scale.data"))) %>% 
		identity()

	# filtered ------------------------------
	
	heatmap_features <- filtered_features %>% 
		dplyr::left_join(regressed_features, by = "Gene.Name") %>% 
		dplyr::mutate(regressed_cluster = replace_na(regressed_cluster, "NA")) %>% 
		identity()
	
	filtered_heatmap <- ggplotify::as.ggplot(
		seu_complex_heatmap(filtered_seu,
												features = heatmap_features$Gene.Name,
												group.by = c("G2M.Score", "S.Score", "scna", regressed_cluster),
												col_arrangement = c(filtered_cluster, "scna"),
												cluster_rows = FALSE,
												column_split =  sort(filtered_seu@meta.data[[filtered_cluster]]),
												row_split = rev(heatmap_features$filtered_cluster),
												row_title_rot = 0
												# right_annotation = row_ha
		)) +
		labs(title = "filtered") +
		theme()
	
	cc_data <- FetchData(filtered_seu, c(filtered_cluster, regressed_cluster, "G2M.Score", "S.Score", "Phase", "scna"))
	
	centroid_data <-
		cc_data %>%
		dplyr::group_by(.data[[filtered_cluster]]) %>%
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
		# dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
		dplyr::mutate({{filtered_cluster}} := as.factor(.data[[filtered_cluster]])) %>% 
		dplyr::mutate(centroid = "centroids") %>%
		identity()
	
	filtered_facet_cell_cycle_plot <-
		cc_data %>% 
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[filtered_cluster]], color = .data[["scna"]])) +
		geom_point(size = 0.1) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[filtered_cluster]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
		facet_wrap(~.data[[filtered_cluster]], ncol = 2) +
		theme_light() +
		# geom_label(data = labels,
		# 					 aes(label = label),
		# 					 # x = Inf,
		# 					 # y = -Inf,
		# 					 x = max(cc_data$S.Score)+0.05,
		# 					 y = max(cc_data$G2M.Score)-0.1,
		# 					 hjust=1,
		# 					 vjust=1,
		# 					 inherit.aes = FALSE) +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank()
		) +
		labs(title = "filtered") + 
		# guides(color = "none") +
		NULL
	
	# regressed ------------------------------
	
	heatmap_features <- regressed_features %>% 
		dplyr::left_join(filtered_features, by = "Gene.Name") %>% 
		dplyr::mutate(filtered_cluster = replace_na(filtered_cluster, "NA")) %>% 
		identity()
	
	regressed_heatmap <- ggplotify::as.ggplot(
		seu_complex_heatmap(filtered_seu,
												features = heatmap_features$Gene.Name,
												group.by = c("G2M.Score", "S.Score", "scna", filtered_cluster),
												col_arrangement = c(regressed_cluster, "scna"),
												cluster_rows = FALSE,
												column_split =  sort(filtered_seu@meta.data[[regressed_cluster]]),
												row_split = rev(heatmap_features$regressed_cluster),
												row_title_rot = 0
												# right_annotation = row_ha
		)) +
		labs(title = "regressed") +
		theme()
	
	cc_data <- FetchData(filtered_seu, c(filtered_cluster, regressed_cluster, "G2M.Score", "S.Score", "Phase", "scna"))
	
	centroid_data <-
		cc_data %>%
		dplyr::group_by(.data[[regressed_cluster]]) %>%
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
		# dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
		dplyr::mutate({{regressed_cluster}} := as.factor(.data[[regressed_cluster]])) %>% 
		dplyr::mutate(centroid = "centroids") %>%
		identity()
	
	regressed_facet_cell_cycle_plot <-
		cc_data %>% 
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[regressed_cluster]], color = .data[["scna"]])) +
		geom_point(size = 0.1) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[regressed_cluster]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
		facet_wrap(~.data[[regressed_cluster]], ncol = 2) +
		theme_light() +
		# geom_label(data = labels,
		# 					 aes(label = label),
		# 					 # x = Inf,
		# 					 # y = -Inf,
		# 					 x = max(cc_data$S.Score)+0.05,
		# 					 y = max(cc_data$G2M.Score)-0.1,
		# 					 hjust=1,
		# 					 vjust=1,
		# 					 inherit.aes = FALSE) +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank()
		) +
		labs(title = "regressed") + 
		# guides(color = "none") +
		NULL

	
	# patchworks ------------------------------
	wrap_plots(filtered_heatmap, 
						 regressed_heatmap,
						 filtered_facet_cell_cycle_plot, 
						 regressed_facet_cell_cycle_plot) + 
		plot_annotation(title = sample_id)
	
	mypatch <- ggsave(glue("results/{sample_id}_regression_effects.pdf"), w = 24, h = 20)
	
	# markerplot <- plot_seu_marker_heatmap(filtered_seu_path, nb_path = numbat_rds_file, clone_simplifications = large_clone_simplifications, ...)
	return(mypatch)
	
}



plot_effect_of_regression_old <- function(filtered_seu_path, regressed_seu_path, resolution = "0.4", group.by = "SCT_snn_res.0.6"){
  # browser()

  plot_list = list()

  sample_id <- str_extract(filtered_seu_path, "SRR[0-9]*")

  filtered_seu <- readRDS(filtered_seu_path)

  regressed_seu <- readRDS(regressed_seu_path)

  fs::dir_create("results/effect_of_regression")

  # percent.mt umaps ------------------------------

  (FeaturePlot(filtered_seu, features = "percent.mt") +
     labs(title = "filtered")) +
    (FeaturePlot(regressed_seu, features = "percent.mt") +
       labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/percent_mt"))
  plot_path = glue("results/effect_of_regression/percent_mt/{sample_id}_percent_mt_umaps.pdf")
  plot_list["percent_mt_umaps"] = plot_path
  ggsave(plot_path, height = 8, width = 10)

  # distribution ------------------------------

  plot_distribution_of_clones_across_clusters(filtered_seu, sample_id, var_y = group.by)
  fs::dir_create(glue("results/effect_of_regression/distribution/filtered/"))
  plot_path = glue("results/effect_of_regression/distribution/filtered/{sample_id}_filtered_distribution.pdf")
  plot_list["filtered_distribution"] = plot_path
  ggsave(plot_path, height = 4, width = 8)

  filtered_dist_tables <- table_distribution_of_clones_across_clusters(filtered_seu, sample_id, clusters = group.by)

  table_path = glue("results/effect_of_regression/distribution/filtered/{sample_id}_filtered_distribution.xlsx")
  plot_list["filtered_distribution_tables"] = table_path
  writexl::write_xlsx(filtered_dist_tables, table_path)

  plot_distribution_of_clones_across_clusters(regressed_seu, sample_id, var_y = group.by)
  fs::dir_create(glue("results/effect_of_regression/distribution/regression"))
  plot_path = glue("results/effect_of_regression/distribution/regression/{sample_id}_regressed_distribution.pdf")
  plot_list["regressed_distribution"] = plot_path
  ggsave(plot_path, height = 4, width = 8)

  regressed_dist_tables <- table_distribution_of_clones_across_clusters(regressed_seu, sample_id, clusters = group.by)

  table_path = glue("results/effect_of_regression/distribution/regression/{sample_id}_regressed_distribution.xlsx")
  plot_list["regressed_distribution_tables"] = table_path
  writexl::write_xlsx(regressed_dist_tables, table_path)

  # abbreviation markers ------------------------------
  (plot_markers(filtered_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
     labs(title = "filtered")) +
    (plot_markers(regressed_seu, group.by, marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10) +
       labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  plot_path = glue("results/effect_of_regression/abbreviation/{sample_id}_abbreviation_markers.pdf")
  plot_list["abbreviation_markers"] = plot_path
  ggsave(plot_path, height = 12, width = 15)

  filtered_marker_tables <- table_cluster_markers(filtered_seu)

  table_path = glue("results/effect_of_regression/abbreviation/{sample_id}_filtered_markers.xlsx")
  plot_list["filtered_marker_tables"] = table_path
  writexl::write_xlsx(filtered_marker_tables, table_path)

  regressed_marker_tables <- table_cluster_markers(regressed_seu) %>%
    purrr::compact()

  table_path = glue("results/effect_of_regression/abbreviation/{sample_id}_regressed_markers.xlsx")
  plot_list["regressed_marker_tables"] = table_path
  writexl::write_xlsx(regressed_marker_tables, table_path)

  heatmap_features  <-
    table_cluster_markers(regressed_seu) %>%
    pluck(group.by) %>%
    group_by(Cluster) %>%
    slice_head(n=10) %>%
    dplyr::pull(Gene.Name) %>%
    identity()

  ggplotify::as.ggplot(
    seu_complex_heatmap(regressed_seu,
                        features = heatmap_features,
                        group.by = c(group.by, "Phase", "scna"),
                        col_arrangement = c(group.by, "Phase", "scna"),
                        cluster_rows = FALSE)) +
    labs(title = sample_id)

  plot_path = glue("results/effect_of_regression/abbreviation/{sample_id}_abbreviation_heatmap.pdf")
  plot_list["abbreviation_heatmap"] = plot_path
  ggsave(plot_path, height = 8, width = 8)

  # # nCount_gene umaps ------------------------------
  #
  # filtered_seu$log_nCount_gene <- log1p(filtered_seu$nCount_gene)
  # regressed_seu$log_nCount_gene <- log1p(regressed_seu$nCount_gene)
  #
  # (FeaturePlot(filtered_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
  #     labs(title = "filtered")) +
  #   (FeaturePlot(regressed_seu, features = "log_nCount_gene", cols = c("blue", "lightgrey")) +
  #      labs(title = "regressed")) +
  #   plot_annotation(title = sample_id)
  #
  # fs::dir_create(glue("results/effect_of_regression/nCount_gene"))
  # plot_path = glue("results/effect_of_regression/nCount_gene/{sample_id}_nCount_gene_umaps.pdf")
  # plot_list["nCount_gene_umaps"] = plot_path
  # ggsave(plot_path, height = 8, width = 10)

  # # cluster umaps ------------------------------
  # mycols = scales::hue_pal()(length(unique(filtered_seu@meta.data[[group.by]])))
  #
  # (DimPlot(filtered_seu, group.by = group.by, cols = mycols) +
  #   labs(title = "filtered")) +
  # (DimPlot(regressed_seu, group.by = group.by, cols = mycols) +
  #   labs(title = "regressed")) +
  #   plot_annotation(title = sample_id)
  #
  # fs::dir_create(glue("results/effect_of_regression/cluster"))
  # plot_path = glue("results/effect_of_regression/cluster/{sample_id}_cluster_umaps.pdf")
  # plot_list["cluster_umaps"] = plot_path
  # ggsave(plot_path, height = 8, width = 10)
  #
  # abbreviation umaps ------------------------------
  mycols = scales::hue_pal()(length(unique(filtered_seu@meta.data[[group.by]]))+3)

  (DimPlot(filtered_seu, group.by = group.by, cols = mycols) +
      labs(title = "filtered")) +
    (DimPlot(regressed_seu, group.by = group.by, cols = mycols) +
       labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/abbreviation"))
  plot_path = glue("results/effect_of_regression/abbreviation/{sample_id}_abbreviation_umaps.pdf")
  plot_list["abbreviation_umaps"] = plot_path
  ggsave(plot_path, height = 8, width = 10)

  # scna umaps ------------------------------
  mycols = scales::hue_pal()(length(unique(filtered_seu@meta.data[["scna"]]))+3)

  # filtered_seu@meta.data$scna <- vec_split_label_line(filtered_seu@meta.data$scna, 3)
  # regressed_seu@meta.data$scna <- vec_split_label_line(regressed_seu@meta.data$scna, 3)

  (DimPlot(filtered_seu, group.by = "scna", cols = mycols) +
      labs(title = "filtered")) +
    (DimPlot(regressed_seu, group.by = "scna", cols = mycols) +
       labs(title = "regressed")) +
    # (DimPlot(regressed_filtered_seu, group.by = group.by) +
    #     labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/scna"))
  plot_path = glue("results/effect_of_regression/scna/{sample_id}_scna_umaps.pdf")
  plot_list["scna_umaps"] = plot_path
  ggsave(plot_path, height = 8, width = 10)

  # # scna markers ------------------------------
  # regressed_seu$scna[regressed_seu$scna == ""] <- "none"
  #
  # plot_markers(regressed_seu, "scna", marker_method = "presto", return_plotly = FALSE, hide_technical = "all", num_markers = 10, unique_markers = TRUE) +
  #   labs(title = sample_id)
  #
  # plot_path = glue("results/effect_of_regression/{sample_id}_scna_markers.pdf")
  # plot_list["scna_markers"] = plot_path
  # ggsave(plot_path, height = 8, width = 6)

  # browseURL(glue("results/effect_of_regression/{sample_id}_markers_by_scna.pdf"))
  #
  # #   original_mt_plot <-
  # #     FeaturePlot(regressed_seu, features = "percent.mt") +
  # #     labs(title = "filtered")
  # #
  # #   regressed_mt_plot <- FeaturePlot(regressed_seu0, features = "percent.mt") +
  # #   labs(title = "regressed")
  # #
  # #   wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  # #     plot_annotation(title = sample_id)
  # #
  # # ggsave(glue("results/effect_of_regression/{sample_id}_percent_mt.pdf"), heigh = 6, width = 10)
  # # browseURL(glue("results/effect_of_regression/{sample_id}_percent_mt.pdf"))
  #
  # regressed_seu <- AddModuleScore(regressed_seu, subtype_markers)
  # regressed_seu0 <- AddModuleScore(regressed_seu0, subtype_markers)
  #
  #   original_mt_plot <-
  #     FeaturePlot(regressed_seu, features = "Cluster1") +
  #     labs(title = "filtered")
  #
  #   regressed_mt_plot <- FeaturePlot(regressed_seu0, features = "Cluster1") +
  #   labs(title = "regressed")
  #
  #   wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  #     plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_subtype1.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_subtype1.pdf"))
  #
  # original_mt_plot <-
  #   FeaturePlot(regressed_seu, features = "Cluster2") +
  #   labs(title = "filtered")
  #
  # regressed_mt_plot <- FeaturePlot(regressed_seu0, features = "Cluster2") +
  #   labs(title = "regressed")
  #
  # wrap_plots(original_mt_plot, regressed_mt_plot, nrow = 1) +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_subtype2.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_subtype2.pdf"))

  # DimPlot(regressed_seu, group.by = group.by) +
  #   labs(title = "filtered") +
  #   DimPlot(regressed_seu0, group.by = group.by) +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_louvain.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_louvain.pdf"))
  #
  # DimPlot(regressed_seu, group.by = "scna") +
  #   labs(title = "filtered") +
  #   DimPlot(regressed_seu0, group.by = "scna") +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_scna.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_scna.pdf"))
  #
  # DimPlot(regressed_seu, group.by = "Phase") +
  #   labs(title = "filtered") +
  #   DimPlot(regressed_seu0, group.by = "Phase") +
  #   labs(title = "regressed") +
  #   plot_annotation(title = sample_id)
  #
  # ggsave(glue("results/effect_of_regression/{sample_id}_phase.pdf"), heigh = 6, width = 10)
  # browseURL(glue("results/effect_of_regression/{sample_id}_phase.pdf"))

  # phase umaps ------------------------------

  (DimPlot(filtered_seu, group.by = "Phase") +
      labs(title = "filtered")) +
    (DimPlot(regressed_seu, group.by = "Phase") +
       labs(title = "regressed")) +
    plot_annotation(title = sample_id)

  fs::dir_create(glue("results/effect_of_regression/phase"))
  plot_path = glue("results/effect_of_regression/phase/{sample_id}_phase_umaps.pdf")
  plot_list["phase_umaps"] = plot_path
  ggsave(plot_path, height = 8, width = 10)


  return(plot_list)

}

ora_effect_of_regression <- function(filtered_seu_path, regressed_seu_path, resolution = "0.4"){
  # browser()

  plot_list = list()

  sample_id <- str_extract(filtered_seu_path, "SRR[0-9]*")

  filtered_seu <- readRDS(filtered_seu_path)

  regressed_seu <- readRDS(regressed_seu_path)

  fs::dir_create("results/effect_of_regression")

  # cluster ora_analysis ------------------------------
  fs::dir_create(glue("results/effect_of_regression/enrichment"))
  fs::dir_create(glue("results/effect_of_regression/enrichment/filtered"))

  filtered_ora_output <- ora_analysis(filtered_seu, "SCT_snn_res.0.4")

  filtered_ora_tables <- filtered_ora_output[["tables"]]

  table_path = glue("results/effect_of_regression/enrichment/filtered/{sample_id}_filtered_cluster_enrichment.xlsx")
  plot_list["filtered_enrichment_tables"] = table_path
  writexl::write_xlsx(filtered_ora_tables, table_path)

  filtered_ora_plots <- filtered_ora_output[["plots"]] %>%
    wrap_plots() +
    plot_annotation(title = sample_id)

  plot_path = glue("results/effect_of_regression/enrichment/filtered/{sample_id}_filtered_cluster_enrichment.pdf")
  plot_list["filtered_enrichment"] = plot_path
  ggsave(plot_path, height = 32, width = 40)

  fs::dir_create(glue("results/effect_of_regression/enrichment"))
  fs::dir_create(glue("results/effect_of_regression/enrichment/regressed"))

  regressed_ora_output <- ora_analysis(regressed_seu, "SCT_snn_res.0.4")
  regressed_ora_tables <- regressed_ora_output[["tables"]]

  table_path = glue("results/effect_of_regression/enrichment/regressed/{sample_id}_regressed_cluster_enrichment.xlsx")
  plot_list["regressed_enrichment_tables"] = table_path
  writexl::write_xlsx(regressed_ora_tables, table_path)

  regressed_ora_plots <- regressed_ora_output[["plots"]] %>%
    wrap_plots() +
    plot_annotation(title = sample_id)

  plot_path = glue("results/effect_of_regression/enrichment/regressed/{sample_id}_regressed_cluster_enrichment.pdf")
  plot_list["regressed_enrichment"] = plot_path
  ggsave(plot_path, height = 32, width = 40)

  # diffex ora_analysis ------------------------------
  fs::dir_create(glue("results/effect_of_regression/enrichment"))
  fs::dir_create(glue("results/effect_of_regression/enrichment/filtered"))

  filtered_ora_output <- ora_analysis(filtered_seu, "SCT_snn_res.0.4")

  filtered_ora_tables <- filtered_ora_output[["tables"]]

  table_path = glue("results/effect_of_regression/enrichment/filtered/{sample_id}_filtered_diffex_enrichment.xlsx")
  plot_list["filtered_enrichment_tables"] = table_path
  writexl::write_xlsx(filtered_ora_tables, table_path)

  filtered_ora_plots <- filtered_ora_output[["plots"]] %>%
    wrap_plots() +
    plot_annotation(title = sample_id)

  plot_path = glue("results/effect_of_regression/enrichment/filtered/{sample_id}_filtered_diffex_enrichment.pdf")
  plot_list["filtered_enrichment"] = plot_path
  ggsave(plot_path, height = 32, width = 40)

  fs::dir_create(glue("results/effect_of_regression/enrichment"))
  fs::dir_create(glue("results/effect_of_regression/enrichment/regressed"))

  regressed_ora_output <- ora_analysis(regressed_seu, "SCT_snn_res.0.4")
  regressed_ora_tables <- regressed_ora_output[["tables"]]

  table_path = glue("results/effect_of_regression/enrichment/regressed/{sample_id}_regressed_diffex_enrichment.xlsx")
  plot_list["regressed_enrichment_tables"] = table_path
  writexl::write_xlsx(regressed_ora_tables, table_path)

  regressed_ora_plots <- regressed_ora_output[["plots"]] %>%
    wrap_plots() +
    plot_annotation(title = sample_id)

  plot_path = glue("results/effect_of_regression/enrichment/regressed/{sample_id}_regressed_diffex_enrichment.pdf")
  plot_list["regressed_enrichment"] = plot_path
  ggsave(plot_path, height = 32, width = 40)


  return(plot_list)

}

split_label_line = function(label, n_comma_values = 2){
  label_vec <- label %>% stringr::str_split_1(",")

  label_groups <- ceiling(seq_along(label_vec)/n_comma_values)

  split_label <- split(label_vec, label_groups) %>%
    purrr::map(paste, collapse = ", ") %>%
    paste(collapse = "\n")
}

vec_split_label_line <- Vectorize(split_label_line)

score_samples_for_celltype_enrichment <- function(unfiltered_seu_path, filtered_seu_path, celltype_markers){
  # browser()
  sample_id <- str_extract(unfiltered_seu_path, "SRR[0-9]*")

  unfiltered_seu <- readRDS(unfiltered_seu_path) %>%
    Seurat::AddModuleScore(features = celltype_markers, name = "celltype")

  filtered_seu <- readRDS(filtered_seu_path) %>%
    Seurat::AddModuleScore(features = celltype_markers, name = "celltype")

  module_names = paste0("celltype", seq(length(celltype_markers)))

  unfiltered_seu@meta.data[names(celltype_markers)] <- unfiltered_seu@meta.data[module_names]
  unfiltered_seu@meta.data[module_names] <- NULL

  filtered_seu@meta.data[names(celltype_markers)] <- filtered_seu@meta.data[module_names]
  filtered_seu@meta.data[module_names] <- NULL

  # unfiltered featureplot ------------------------------
  (FeaturePlot(unfiltered_seu, names(celltype_markers)) +
    plot_annotation(subtitle = "unfiltered"))

  file_tag <- str_extract(unfiltered_seu_path, "SRR[0-9]*_[a-z]*")
  dir_create("results/celltype_plots")
  unfiltered_celltype_plot_path <- glue("results/celltype_plots/{file_tag}_featureplot.pdf")
  ggsave(unfiltered_celltype_plot_path, height = 8, width = 8)

  # filtered featureplot ------------------------------
  (FeaturePlot(filtered_seu, names(celltype_markers)) +
      plot_annotation(subtitle = "filtered"))

  file_tag <- str_extract(filtered_seu_path, "SRR[0-9]*_[a-z]*")
  dir_create("results/celltype_plots")
  filtered_celltype_plot_path <- glue("results/celltype_plots/{file_tag}_featureplot.pdf")
  ggsave(filtered_celltype_plot_path, height = 8, width = 8)

  unfiltered_seu <- score_binary_celltype_markers(unfiltered_seu, celltype_markers)
  filtered_seu <- score_binary_celltype_markers(filtered_seu, celltype_markers)

  # unfiltered dimplot ------------------------------
  (DimPlot(unfiltered_seu, group.by = paste0(names(celltype_markers), "_id")) +
     plot_annotation(subtitle = "unfiltered"))

  file_tag <- str_extract(unfiltered_seu_path, "SRR[0-9]*_[a-z]*")
  dir_create("results/celltype_plots")
  unfiltered_celltype_plot_path <- glue("results/celltype_plots/{file_tag}_dimplot.pdf")
  ggsave(unfiltered_celltype_plot_path, height = 8, width = 8)

  # filtered dimplot ------------------------------
  (DimPlot(filtered_seu, group.by = paste0(names(celltype_markers), "_id")) +
     plot_annotation(subtitle = "filtered"))

  file_tag <- str_extract(filtered_seu_path, "SRR[0-9]*_[a-z]*")
  dir_create("results/celltype_plots")
  filtered_celltype_plot_path <- glue("results/celltype_plots/{file_tag}_dimplot.pdf")
  ggsave(filtered_celltype_plot_path, height = 8, width = 8)

  return(list(filtered_celltype_plot_path, unfiltered_celltype_plot_path))

}

score_binary_celltype_markers <- function(seu, celltype_markers){
  # browser()
  upper_quartiles <-
    names(celltype_markers) %>%
    set_names(.) %>%
    map(~quantile(seu@meta.data[[.x]], 0.75))

  for(celltype in names(celltype_markers)){
    seu@meta.data[[paste0(celltype, "_id")]] <- seu@meta.data[[celltype]] > upper_quartiles[[celltype]]
  }

  return(seu)

}

score_binary_celltype_clusters <- function(seu, celltype_markers){
  browser()
  upper_quartiles <-
    names(celltype_markers) %>%
    set_names(.) %>%
    map(~quantile(seu@meta.data[[.x]], 0.95))

  seu <- seurat_cluster(seu, resolution = 2.0)

  pick_max_cluster <- function(upper_quartile, cluster_name, seu){


    test0 <-
      seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::group_by(`gene_snn_res.2`) %>%
      dplyr::filter(dplyr::n() < 30) %>%
      dplyr::summarize(score = mean(.data[[cluster_name]])) %>%
      # dplyr::filter(score > upper_quartile) %>%
      # dplyr::slice_max(score) %>%
      # dplyr::pull(`gene_snn_res.2`) %>%
      identity()

  }

  test0 <- imap(upper_quartiles, pick_max_cluster, seu)

  for(celltype in names(celltype_markers)){
    seu@meta.data[[paste0(celltype, "_id")]] <- seu@meta.data[["gene_snn_res.2"]] %in% test0[[celltype]]
  }

  return(seu)

}

convert_seu_to_scanpy <- function(seu_path){

  filtered_seu <- readRDS(seu_path)

  # filtered_seu <- DietSeurat(filtered_seu, misc = FALSE)

  filtered_seu = DietSeurat(
    filtered_seu,
    counts = TRUE, # so, raw counts save to adata.layers['counts']
    data = TRUE, # so, log1p counts save to adata.X when scale.data = False, else adata.layers['data']
    scale.data = FALSE, # if only scaled highly variable gene, the export to h5ad would fail. set to false
    features = rownames(filtered_seu), # export all genes, not just top highly variable genes
    assays = "gene",
    dimreducs = c("pca","umap"),
    graphs = c("gene_nn", "gene_snn"), # to RNA_nn -> distances, RNA_snn -> connectivities
    misc = FALSE
  )

  filtered_scanpy_path <- str_replace(seu_path, "seu.rds", "scanpy.h5ad") %>%
    str_replace("seurat", "scanpy")

  MuDataSeurat::WriteH5AD(filtered_seu, filtered_scanpy_path, assay="gene")

  return(filtered_scanpy_path)

}


score_chrom_instability <- function(seu_path) {

  sample_id <- str_extract(seu_path, "SRR[0-9]*")

  seu <- readRDS(seu_path)

  cin_scores <- list(
    "cin70" = c("TPX2", "PRC1", "FOXM1", "CDC2", "TGIF2", "MCM2", "H2AFZ", "TOP2A", "PCNA", "UBE2C", "MELK", "TRIP13", "CNAP1", "MCM7", "RNASEH2A", "RAD51AP1", "KIF20A", "CDC45L", "MAD2L1", "ESPL1", "CCNB2", "FEN1", "TTK", "CCT5", "RFC4", "ATAD2", "ch-TOG", "NUP205", "CDC20", "CKS2", "RRM2", "ELAVL1", "CCNB1", "RRM1", "AURKB", "MSH6", "EZH2", "CTPS", "DKC1", "OIP5", "CDCA8", "PTTG1", "CEP55", "H2AFX", "CMAS", "BRRN1", "MCM10", "LSM4", "MTB", "ASF1B", "ZWINT", "TOPK", "FLJ10036", "CDCA3", "ECT2", "CDC6", "UNG", "MTCH2", "RAD21", "ACTL6A", "GPIandMGC13096", "SFRS2", "HDGF", "NXT1", "NEK2", "DHCR7", "STK6", "NDUFAB1", "KIAA0286", "KIF4A"
    ),
    "pos_tri70" = c("ANXA7", "ATG7", "BDNF", "CDKN2B", "CHFR", "CNDP2", "ENO3", "F3", "GIPC2", "GJB5", "HSPB7", "IMPACT", "P2RY14", "PCDH7", "PKD1", "PLCG2", "SNCA", "SNCG", "TMEM140", "TMEM40"),
    "neg_tri70" = c("AURKA", "BCL11B", "BIRC5", "BLMH", "BRD8", "BUB1B", "CCNA2", "CDC5L", "CDK1", "CENPE", "CENPN", "CTH", "DLGAP5", "HMGB2", "IDH2", "ISOC1", "KIAA0101", "KIF22", "LIG1", "LSM2", "MCM2", "MCM5", "MCM7", "MYBL2", "NASP", "NCAPD2", "NCAPH", "NMI", "NUDT21", "PCNA", "PIGO", "PLK1", "PLK4", "POLD2", "RACGAP1", "RAD51", "RFC2", "RFC3", "RFC5", "RPA2", "SMAD4", "SMC4", "SSRP1", "TAB2", "TCOF1", "TIMELESS", "TIPIN", "TOP2A", "UBE2C", "USP1"),
    "het70" = c("AHCYL1", "AKT3", "ANO10", "ANTXR1", "ATP6V0E1", "ATXN1", "B4GALT2", "BASP1", "BHLHE40", "BLVRA", "CALU", "CAP1", "CAST", "CAV1", "CLIC4", "CTSL1", "CYB5R3", "ELOVL1", "EMP3", "FKBP14", "FN1", "FST", "GNA12", "GOLT1B", "HECTD3", "HEG1", "HOMER3", "IGFBP3", "IL6ST", "ITCH", "LEPRE1", "LEPREL1", "LEPROT", "LGALS1", "LIMA1", "LPP", "MED8", "MMP2", "MUL1", "MYO10", "NAGK", "NR1D2", "NRIP3", "P4HA2", "PKIG", "PLOD2", "PMP22", "POFUT2", "POMGNT1", "PRKAR2A", "RAGE", "RHOC", "RRAGC", "SEC22B", "SERPINB8", "SPAG9", "SQSTM1", "TIMP2", "TMEM111", "TRIM16", "TRIO", "TUBB2A", "VEGFC", "VIM", "WASL", "YIPF5", "YKT6", "ZBTB38", "ZCCHC24", "ZMPSTE24"),
    "bucc_up" = c("HPDL", "L1CAM", "SMPD4", "ELF4", "IRAK1", "ISG20L2", "ZNF512B", "TMEM164", "RCE1", "AMPD2", "SCXB", "HSF1", "PLXNA3", "RBP4", "MRPL11", "NSUN5", "PNMA5", "FADS3", "TRAIP", "IGSF9", "TFAP4", "BOP1", "FAM131A", "EIF1AD", "FAM122C", "SLC10A3", "CCNK", "ELK1", "DGAT1", "DHCR7", "GNAZ", "SCAMP3", "ZNF7", "CCDC86", "SLC35A2", "CPSF4", "USP39", "FTSJD2", "DBF4B", "MEPCE", "ZBTB2", "SEMA4D", "LRRC14", "YKT6", "METTL7B", "SCAMP5", "IP6K1", "MTCH1", "SLC6A9", "SMCR7L", "SAMD4B", "CD3EAP", "PPP1R16A", "TMCO6", "U2AF2", "POM121", "TSR2", "UNKL", "IFRD2", "CYHR1", "MAGEA9B", "TBC1D25", "ZNF275", "TJAP1", "CLDN4", "SPNS1", "TMEM120A", "CPSF1", "FABP3", "FAM58A"),
    "bucc_down" = c("CYTIP", "HLA-DMB", "FYCO1", "NUBP1", "AGPS", "CADPS2", "C3orf52", "FGL2", "PDIA3P", "IL18", "APOBEC3G", "CA8", "HNRNPA1L2", "CD44", "TMEM2", "ELL2", "FCGBP", "TAPT1", "AHR", "HLCS", "AEN", "AGR3", "STEAP1", "TLR4", "HNRNPA1", "HLA-DRB5", "SLITRK6", "LPXN", "HLA-DRB1", "BACE2", "SYTL1", "KYNU", "RPS27L", "CD96", "HLA-DPA1", "TRA2A", "FBXO25", "HFE", "DRAM1", "DYRK4", "FBXL5", "GLUL", "HLA-DRA", "KCNMA1", "IL18R1", "SIPA1L2", "STEAP2", "TP53I3", "LIMA1", "BMX", "HLA-DMA", "C17orf97", "RPH3AL", "SLC37A1", "PHLDA3", "TK2", "NHS", "TMEM229B", "CYFIP1", "B3GALT5", "FAM107B", "ADH6", "CCND1", "LCK", "RABL2B", "SPA17", "GADD45A", "FHL2", "RPL36AL", "CCDC60", "MUC5B", "PHLDA1", "ZG16B", "SIDT1", "BTD", "TTC9", "GSTZ1", "ITGA6", "C10orf32", "FAM46C", "SPDEF", "ASRGL1", "TCN1", "NUDT2", "CDKN1A", "POLR3GL", "C4BPA", "LAMC2", "SLC6A14")
  )

  seu <- AddModuleScore(seu, cin_scores, name = "cin")

  names(seu@meta.data)[which(names(seu@meta.data) %in% paste0("cin", seq(1, length(cin_scores))))] <- names(cin_scores)

  cin_score_fplot <- FeaturePlot(seu, names(cin_scores)) +
    plot_annotation(title = sample_id)

  dir_create("results/effect_of_regression/cin_scores/")

  plot_path <- glue("results/effect_of_regression/cin_scores/{sample_id}_cin_scores.pdf")

  ggsave(plot_path, height = 8, width = 8)

  return(plot_path)
}

score_stachelek <- function(seu_path, stachelek_scores_table) {

  sample_id <- str_extract(seu_path, "SRR[0-9]*")

  seu <- readRDS(seu_path)

  stachelek_scores <-
    stachelek_scores_table %>%
    list_flatten() %>%
    map(dplyr::distinct, symbol) %>%
    map(pull, symbol)


  seu <- AddModuleScore(seu, stachelek_scores, name = "stachelek")

  names(seu@meta.data)[which(names(seu@meta.data) %in% paste0("stachelek", seq(1, length(stachelek_scores))))] <- names(stachelek_scores)

  # stachelek_score_fplot <- FeaturePlot(seu, names(stachelek_scores)) +
  #   plot_annotation(title = sample_id)

  # stachelek_score_vlnplot <- VlnPlot(seu, names(stachelek_scores), group.by = "scna", ncol = 8) +
  #   plot_annotation(title = sample_id)

  stachelek_score_vlnplot <- VlnPlot(seu, names(stachelek_scores), group.by = "abbreviation", ncol = 8) +
    plot_annotation(title = sample_id)

  dir_create("results/effect_of_regression/stachelek/")

  plot_path <- glue("results/effect_of_regression/stachelek/{sample_id}_stachelek_scores.pdf")

  ggsave(plot_path, height = 8, width = 24)

  return(plot_path)
}


plot_plae_celltype_expression <- function(mygenes = c("RXRG", "NRL"), plot_type="box"){
  # browser()

  celltypes = c("Retinal Ganglion Cell", "Amacrine Cell", "Horizontal Cell", "RPC", "Early RPC", "Muller Glia", "Bipolar Cell", "Late RPC", "Neurogenic Cell", "B-Cell", "Rod", "Photoreceptor Precursor", "Cone", "RPE", "Microglia", "Red Blood Cell", "Astrocyte", "Rod Bipolar Cell")

  # pseudo_meta <- read_tsv("/dataVolume/storage/scEiad/human_pseudobulk/4000-counts-universe-study_accession-scANVIprojection-15-5-20-50-0.1-CellType-Homo_sapiens.meta.tsv.gz") %>%
  #   tidyr::unite(study_type, study_accession, CellType)

  sub_annotable <-
    annotables::grch38 %>%
    dplyr::filter(symbol %in% mygenes)

  pseudo_counts <-
    "data/plae_pseudobulk_counts.csv" %>%
    read_csv() %>%
    dplyr::inner_join(sub_annotable, by = c("Gene" = "ensgene"), relationship = "many-to-many") %>%
    dplyr::distinct(study, type, symbol, .keep_all = TRUE) %>%
    # dplyr::filter(type %in% celltypes) %>%
    identity()

  if(plot_type == "box"){
    exp_plot <- ggplot(pseudo_counts, aes(x = type,
                                                y = counts,
                                                color = study)) +
      geom_boxplot(color = 'black', outlier.shape = NA) +
      ggbeeswarm::geom_quasirandom(groupOnX = TRUE) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_radius(range=c(2, 6)) +
      scale_colour_manual(values = rep(c(pals::alphabet() %>% base::unname()), 20)) +
      theme(legend.position="bottom") +
      facet_wrap(ncol = 2, scales = 'free_y', ~symbol) +
      NULL
  } else if(plot_type == "hmap"){
    exp_plot <-
      pseudo_counts %>%
      group_by(symbol, type) %>%
      dplyr::summarize(median_counts = median(counts)) %>%
      ggplot(aes(x = symbol,
                 y = type,
                 fill = median_counts)) +
      geom_tile() +
      # scale_fill_gradient(name = "median_count", trans = "log") +
      NULL
  }

    return(exp_plot)
  }

reference_plae_celltypes <- function(seu_path, mycluster = "SCT_snn_res.0.4", mygenes){

  # con <- dbConnect(RSQLite::SQLite(), "/dataVolume/storage/scEiad/human_pseudobulk/diff_resultsCellType.sqlite")
  #
  # mytable <-
  #   tbl(con, "diffex") %>%
  #   dplyr::filter(Against =="All") %>%
  #   dplyr::filter(Base %in% celltypes) %>%
  #   dplyr::filter(padj < 0.05, abs(log2FoldChange) > 2) %>%
  #   dplyr::filter(Organism == "Homo sapiens") %>%
  #   collect()

  seu_path <- "output/seurat/SRR13884242_regressed_seu.rds"

  seu <- readRDS(seu_path)

  cluster_markers <-
    seu@misc$markers[[mycluster]][["presto"]] %>%
    dplyr::group_by(Cluster) %>%
    dplyr::slice_head(n=30) %>%
    dplyr::rename(symbol = `Gene.Name`)

  mytable <-
    "data/plae_top_diffex.csv" %>%
    read_csv() %>%
    dplyr::group_by(Base) %>%
    dplyr::arrange(desc(log2FoldChange)) %>%
    dplyr::mutate(symbol = str_remove(Gene, " \\(.*\\)")) %>%
    dplyr::filter(!is.na(symbol)) %>%
    dplyr::left_join(cluster_markers, by = "symbol", relationship = "many-to-many") %>%
    dplyr::filter(!is.na(Cluster)) %>%
    dplyr::filter(log2FoldChange > 0) %>%
    dplyr::arrange(Cluster, Base) %>%
    dplyr::select(Cluster, Base, everything()) %>%
    # dplyr::mutate(checked_gene = case_when(symbol %in% mygenes ~ 1)) %>%
    identity()

  mytable0 <-
    mytable %>%
    dplyr::group_by(Cluster, Base) %>%
    dplyr::summarise(mean_fc = mean(log2FoldChange)) %>%
    dplyr::arrange(Cluster, desc(mean_fc))

  # janitor::tabyl(mytable, Cluster, Base)

}

drop_bad_cells <- function(seu_path, bad_cell_types = c("RPCs", "Late RPCs", c("Red Blood Cells", "Microglia", "Muller Glia", "RPE", "Horizontal Cells", "Rod Bipolar Cells", "Pericytes", "Bipolar Cells", "Astrocytes", "Endothelial", "Schwann", "Fibroblasts"))) {
  browser()
  seu <- readRDS(seu_path)

  seu <- seu[,!seu$type %in% bad_cell_types]

  saveRDS(seu, str_replace(seu_path, "_seu.rds", "_dropped_cells_seu.rds"))

  retainedcells_dimplot <- DimPlot(
    seu,
    group.by = "type"
  ) +
    labs(title = fs::path_file(seu_path))

  return(seu_path)

}

generate_plae_ref <- function(plae_seu_path = "data/plae_human_fetal_seu.rds"){

  plae_human_fetal_seu <- readRDS(plae_seu_path)

  plae_ref <- AggregateExpression(plae_human_fetal_seu, group.by = "CellType_predict") %>%
    pluck("RNA")

  return(plae_ref)
}


plot_celltype_predictions <- function(seu_path, sample_id, plae_ref = NULL, group.by = "SCT_snn_res.0.4", query_genes = NULL) {
  # browser()

  seu <- readRDS(seu_path)

  celltypes = c("Amacrine Cells", "Bipolar Cells", "Cones",
                "Early RPCs", "Horizontal Cells", "Late RPCs",
                "Neurogenic Cells", "Photoreceptor Precursors",
                "Retinal Ganglion Cells", "Rods", "RPCs", "RPE"
  )

  if(is.null(query_genes)){
    query_genes <- VariableFeatures(seu)
  }

  seu_mat <- GetAssayData(seu, slot = "data")

  sub_annotable <-
    annotables::grch38 %>%
    dplyr::filter(symbol %in% query_genes)

  if(is.null(plae_ref)){
    plae_ref <-
      "data/plae_pseudobulk_counts.csv" %>%
      read_csv() %>%
      dplyr::inner_join(sub_annotable, by = c("Gene" = "ensgene"), relationship = "many-to-many") %>%
      dplyr::distinct(study, type, symbol, .keep_all = TRUE) %>%
      dplyr::filter(type %in% str_remove(celltypes, "s$")) %>%
      dplyr::group_by(symbol, type) %>%
      dplyr::summarize(total_counts = mean(counts)) %>%
      tidyr::pivot_wider(names_from = "type", values_from = "total_counts") %>%
      tibble::column_to_rownames("symbol") %>%
      as.matrix() %>%
      identity()
  } else {
    plae_ref = plae_ref[rownames(plae_ref) %in% query_genes,colnames(plae_ref) %in% celltypes]
  }

  res <- clustify(
    input = seu_mat,
    metadata = seu[[group.by]][[1]],
    ref_mat = plae_ref,
    query_genes = query_genes
  )

  cor_to_call(res)

  res2 <- cor_to_call(
    cor_mat = res,                  # matrix correlation coefficients
    cluster_col = group.by # name of column in meta.data containing cell clusters
  )

  res3 <- cor_to_call_rank(
    cor_mat = res,                  # matrix correlation coefficients
    cluster_col = group.by # name of column in meta.data containing cell clusters
  ) %>%
    dplyr::mutate(cluster = .data[[group.by]]) %>%
    dplyr::group_by(cluster) %>%
    dplyr::arrange(rank)

  dir_create("results/clustify")
  table_path = glue("results/clustify/{sample_id}_clustifyr.csv")
  write_csv(res3, table_path)


  # Insert into original metadata as "type" column
  seu@meta.data <- call_to_metadata(
    res = res2,                     # data.frame of called cell type for each cluster
    metadata = seu@meta.data,           # original meta.data table containing cell clusters
    cluster_col = group.by # name of column in meta.data containing cell clusters
  )

  neurogenic_table <- seu@meta.data %>%
    dplyr::mutate(cluster = .data[[group.by]]) %>%
    janitor::tabyl(cluster, type) %>%
    dplyr::mutate(sample = sample_id) %>%
    # dplyr::mutate(percent_neurogenic = `Neurogenic Cells`/(Cones+`Neurogenic Cells`)) %>%
    identity()

  allcells_dimplot <- DimPlot(
    seu,
    group.by = "type",
    split.by = group.by
  ) +
    plot_annotation(title = sample_id) +
    NULL

  dir_create("results/clustify")
  plot_path = glue("results/clustify/{sample_id}_clustifyr.pdf")
  ggsave(plot_path)

  return(list("plot" = plot_path, "table" = neurogenic_table))

}

read_zinovyev_genes <- function(zinovyev_file = "data/zinovyev_cc_genes.tsv"){
  zinovyev_cc_genes <-
    read_tsv(zinovyev_file) %>%
    dplyr::group_by(term) %>%
    tidyr::nest(data = symbol) %>%
    # split(.$term) %>%
    tibble::deframe() %>%
    map(tibble::deframe) %>%
    identity()
}

read_giotti_genes <- function(cc_file = "data/giotti_cc_genes.tsv"){
  giotti_cc_genes <-
    read_tsv(cc_file) %>%
    dplyr::filter(!(term %in% c("Function known but link to cell division not well established", "Uncharacterized"))) %>%
    dplyr::mutate(term = str_wrap(term, width = 10)) %>%
    identity()
}

giotti_genes <- read_giotti_genes()


plot_phase_distribution_of_all_samples_by_scna <- function(seu_paths, selected_samples = c("SRR13884243", "SRR13884249", "SRR14800534",
                                                                                           "SRR14800535", "SRR14800536", "SRR14800543",
                                                                                           "SRR17960481")){
  # browser()

  plot_phase_distribution_by_scna <- function(seu_path, seu_name){
    seu <- readRDS(seu_path)
    seu$Phase <- factor(seu$Phase, levels = c("G1", "S", "G2M"))
    plot_distribution_of_clones_across_clusters(seu, seu_name, var_x = "scna", var_y = "Phase", both_ways = FALSE)
  }

  seu_paths <-
    seu_paths %>%
    set_names(str_extract(., "SRR[0-9]*"))

  if(!is.null(selected_samples)){

    seu_paths <- seu_paths[selected_samples]
  }

  seu_plots <-
    seu_paths %>%
    imap(plot_phase_distribution_by_scna)

  test0 <- wrap_plots(seu_plots, ncol = 3) +
    plot_layout(guides = 'collect') +
    plot_annotation(title = 'Proliferation varies directly with scna abundance',
                    theme = theme(plot.title = element_text(size = 18)))

  plot_path = "results/relation_between_scna_and_phase.pdf"
  ggsave(plot_path, height= 8, width = 12)

  return(plot_path)

}


plot_filtering_timeline <- function(all_cells_meta, scna_meta, qc_meta, cell_type_meta, sample_id) {
  # browser()

  n_cells <- nrow(all_cells_meta)

  all_cells_meta$sample_id <- sample_id

  all_cells_bar <-
    all_cells_meta %>%
    tibble::rownames_to_column("cell") %>%
    ggplot(fill = "gray", aes(x = sample_id)) +
    geom_bar(position = "stack", width = 0.1) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin = margin(0, 0, 0, 0, "pt"),
          legend.position = "bottom") +
    ylim(0,n_cells) +
    labs(title = "all cells") +
    scale_x_discrete(expand = c(0,0)) +
    # scale_y_log10() +
    NULL

  # scna_meta <- dplyr::mutate(scna_meta, scna = ifelse(scna == "", "none", scna))
  #
  scna_levels <- levels(factor(scna_meta$scna))

  scna_pal <- scales::hue_pal()(length(scna_levels)) %>%
    set_names(scna_levels) %>%
    identity()

  scna_bar <-
    scna_meta %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(scna = wrap_scna_labels(scna)) %>%
    ggplot(aes(fill = scna, x = sample_id)) +
    geom_bar(position = "stack", width = 0.1) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.margin = margin(0, 0, 0, 0, "pt"),
          legend.position = "bottom") +
    ylim(0,n_cells) +
    # labs(title = "scna") +
    scale_x_discrete(expand = c(0,0)) +
    scale_fill_manual(values = scna_pal) +
    # scale_y_log10() +
    NULL

  qc_bar <-
    qc_meta %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(scna = wrap_scna_labels(scna)) %>%
    ggplot(aes(fill = scna, x = sample_id)) +
    geom_bar(position = "stack", width = 0.1) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          plot.margin = margin(0, 0, 0, 0, "pt"),
          legend.position = "bottom") +
    ylim(0,n_cells) +
    # labs(title = "qc") +
    scale_x_discrete(expand = c(0,0)) +
    scale_fill_manual(values = scna_pal) +
    # scale_y_log10() +
    NULL

  cell_type_bar <-
    cell_type_meta %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::mutate(scna = wrap_scna_labels(scna)) %>%
    ggplot(aes(fill = .data[["scna"]], x = sample_id)) +
    geom_bar(position = "stack", width = 0.2) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          plot.margin = margin(0, 0, 0, 0, "pt"),
          # panel.grid = element_blank(),
          # panel.border = element_blank(),
          legend.position = "bottom") +
    ylim(0,n_cells) +
    # labs(title = "bad cells") +
    scale_x_discrete(expand = c(0,0)) +
    scale_fill_manual(values = scna_pal) +
    theme(legend.position = "none") +
    # scale_y_log10() +
    NULL

  wrap_plots(scna_bar, qc_bar, cell_type_bar, nrow = 1) +
    plot_layout(guides = "collect") +
    plot_annotation(title = sample_id)

}

pull_common_markers <- function(filtered_seus, gene_lists){

  annotated_genes <-
    gene_lists %>%
    tibble::enframe("mp", "symbol") %>%
    tidyr::unnest(symbol) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    dplyr::arrange(mp) %>%
    identity()

  names(filtered_seus) <- str_extract(filtered_seus, "SRR[0-9]*")

  # load filtered_seus ------------------------------
  # find cluster markers of every seu, compare with zinovyev markers
  # find common markers for clusters
  pull_cluster_markers <- function(seu_path){
    seu <- readRDS(seu_path)

    table_cluster_markers(seu)

  }

  my_cluster_markers <- map(filtered_seus, pull_cluster_markers)

  common_genes <-
    map(my_cluster_markers, "SCT_snn_res.0.4") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    dplyr::filter(abs(Average.Log.Fold.Change) > 0.5) %>%
    dplyr::arrange(Gene.Name) %>%
    dplyr::group_by(Gene.Name) %>%
    dplyr::filter(dplyr::n() > 3) %>%
    dplyr::filter(!str_detect(Gene.Name, pattern = "^RP.*")) %>%
    dplyr::filter(!str_detect(Gene.Name, pattern = "^MT.*")) %>%
    dplyr::distinct(Gene.Name) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(annotated_genes, by = c("Gene.Name" = "symbol")) %>%
    # dplyr::slice_sample(n =50) %>%
    # dplyr::pull(Gene.Name) %>%
    # unique() %>%
    # sample(50) %>%
    identity()

  return(common_genes)

}

heatmap_marker_genes <- function(seu_path, common_genes, gene_lists, label="", sample_id = NULL, marker_col = "SCT_snn_res.0.4", group.by = c("SCT_snn_res.0.4", "scna", "subtype1", "subtype2"), col_arrangement = c("SCT_snn_res.0.4", "scna")){
  # browser()
  if(is.null(sample_id)){
    sample_id <- str_extract(seu_path, "SRR[0-9]*")
  }

  seu <- readRDS(seu_path)

  seu <- Seurat::AddModuleScore(seu, features = gene_lists, name = "subtype")

  module_names = paste0("subtype", seq(length(gene_lists)))

  seu@meta.data[names(gene_lists)] <- seu@meta.data[module_names]

  test0 <- seu@misc$markers[[marker_col]]$presto %>%
    dplyr::group_by(Cluster) %>%
    dplyr::slice_head(n=10) %>%
    dplyr::select(Gene.Name, Cluster) %>%
    dplyr::inner_join(common_genes, by = "Gene.Name") %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Cluster) %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE) %>%
    dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
    dplyr::mutate(mp = replace_na(mp, "")) %>%
    identity()

  mymarkers <-
    test0 %>%
    dplyr::pull(Gene.Name)

  seu$scna <- factor(seu$scna)
  levels(seu$scna)[1] <- "none"

  row_ha = ComplexHeatmap::rowAnnotation(mp = rev(test0$mp))

  ggplotify::as.ggplot(seu_complex_heatmap(seu, features = mymarkers, group.by = group.by, col_arrangement = col_arrangement, cluster_rows = FALSE, right_annotation = row_ha)) +
    labs(title = sample_id)


  # seu_heatmap <- Seurat::DoHeatmap(seu, features = mymarkers, group.by = "SCT_snn_res.0.4")

  heatmap_file <- glue("results/{sample_id}_{label}heatmap.pdf")

  ggsave(heatmap_file, height = 8, width = 8)

  return(heatmap_file)

}

wrap_scna_labels <- function(scna_labels){
  str_wrap(str_replace_all(scna_labels, ",", " "), width = 10, whitespace_only = FALSE)
}

make_pairwise_plots <- function(clone_set, clone_comparison, seu, sample_id){
	# browser()
	pair_seu <- seu[,seu$scna %in% clone_set]
	
	clone_ratio = janitor::tabyl(as.character(pair_seu$scna))$percent[[2]]
	
	comparison_scna <-
		janitor::tabyl(as.character(pair_seu$scna))[2,1]
	
	myplot <- plot_distribution_of_clones_across_clusters(
		pair_seu, seu_name = clone_comparison, var_x = "scna", var_y = "clusters", avg_line = clone_ratio, signif = TRUE, plot_type = "clone"
	)
	
	cluster_values <-
		pair_seu@meta.data %>%
		dplyr::group_by(clusters, scna) %>%
		dplyr::summarize(value = dplyr::n())
	
	all_plot_table <- 
		cluster_values %>% 
		dplyr::group_by(scna) %>% 
		dplyr::summarize(value = sum(value)) %>% 
		mutate(percent = value/sum(value)) %>%
		dplyr::select(-value) %>%
		dplyr::mutate(clusters = "all") %>% 
		identity()
	
	mytable <- 
		cluster_values %>% 
		mutate(percent = value/sum(value)) %>%
		dplyr::select(-value) %>%
		dplyr::bind_rows(all_plot_table) %>% 
		tidyr::pivot_wider(names_from = "scna", values_from = "percent") %>%
		dplyr::mutate(sample_id = sample_id) %>% 
		identity()
	
	mytable <-
		mytable %>%
		dplyr::mutate(up = ifelse(.data[[as.character(comparison_scna)]] > (clone_ratio+0.03), 1, 0)) %>%
		dplyr::mutate(down = ifelse(.data[[as.character(comparison_scna)]] < (clone_ratio-0.03), 1, 0)) %>%
		identity()
	
	return(
		list(
			"plot" = myplot,
			"table" = mytable)
	)
	
}

plot_clone_pearls <- function(seu_path, cluster_order, phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star")){
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	cluster_order = cluster_order[[sample_id]]
	
	seu <- readRDS(seu_path)
	seu$scna <- factor(seu$scna)
	levels(seu$scna)[levels(seu$scna) == ""] <- "diploid"
	
	if(!is.null(cluster_order)){
		
		group.by = unique(cluster_order$resolution)
		
		cluster_order <-
			cluster_order %>%
			dplyr::mutate(order = dplyr::row_number()) %>%
			dplyr::filter(!is.na(clusters)) %>%
			dplyr::mutate(clusters = as.character(clusters))
		
		seu@meta.data$clusters = seu@meta.data[[group.by]]
		
		seu_meta <- seu@meta.data %>%
			tibble::rownames_to_column("cell") %>%
			dplyr::left_join(cluster_order, by = "clusters") %>%
			dplyr::select(-clusters) %>%
			dplyr::rename(clusters = phase) %>%
			identity()
		
		phase_levels = phase_levels[phase_levels %in% unique(seu_meta$clusters)]
		
		seu_meta <-
			seu_meta %>%
			tidyr::unite("new_clusters", all_of(c("clusters", group.by)), remove = FALSE) %>%
			dplyr::arrange(clusters, new_clusters) %>%
			dplyr::mutate(clusters = factor(new_clusters, levels = unique(new_clusters))) %>%
			tibble::column_to_rownames("cell") %>%
			identity()
		
		seu@meta.data <- seu_meta[rownames(seu@meta.data),]
		
	}
	
	pearls_plot <- plot_distribution_of_clones_pearls(seu, seu_name = sample_id, var_x = "scna", var_y = "clusters") + 
		labs(title= sample_id)
	
	n_scnas = n_distinct(pearls_plot$data$scna)
	
	plot_width = 3*n_scnas
	
	plot_path = glue("results/{sample_id}_pearls.pdf")
	
	# browser()
	ggsave(plot_path, pearls_plot, width = plot_width, height = 6)
	
	return(
		plot_path
	)
	
}

calculate_clone_distribution <- function(seu_path, cluster_order, pairwise = FALSE, phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star")){
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")

  cluster_order = cluster_order[[sample_id]]

  seu <- readRDS(seu_path)
  seu$scna <- factor(seu$scna)
  levels(seu$scna)[levels(seu$scna) == ""] <- "diploid"

  if(!is.null(cluster_order)){
  	
  	group.by = unique(cluster_order$resolution)
  	
  	cluster_order <-
  		cluster_order %>%
  		dplyr::mutate(order = dplyr::row_number()) %>%
  		dplyr::filter(!is.na(clusters)) %>%
  		dplyr::mutate(clusters = as.character(clusters))

    seu@meta.data$clusters = seu@meta.data[[group.by]]

    seu_meta <- seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::left_join(cluster_order, by = "clusters") %>%
      dplyr::select(-clusters) %>%
      dplyr::rename(clusters = phase) %>%
      identity()

    phase_levels = phase_levels[phase_levels %in% unique(seu_meta$clusters)]

    seu_meta <-
      seu_meta %>%
      tidyr::unite("new_clusters", all_of(c("clusters", group.by)), remove = FALSE) %>%
      dplyr::arrange(clusters, new_clusters) %>%
      dplyr::mutate(clusters = factor(new_clusters, levels = unique(new_clusters))) %>%
      tibble::column_to_rownames("cell") %>%
      identity()

    seu@meta.data <- seu_meta[rownames(seu@meta.data),]

  }

  all_seu_plot <- plot_distribution_of_clones_across_clusters(seu, seu_name = sample_id, var_x = "scna", var_y = "clusters")

  if(pairwise){
    pairwise_seu_plots <- list()

    scna_clones <- unique(sort(as.factor(seu@meta.data$scna)))

    pairwise_clone_vectors <-
      bind_cols(scna_clones[-length(scna_clones)], scna_clones[-1]) %>%
      t() %>%
      as.data.frame() %>%
      as.list() %>%
      map(as.character) %>%
        identity()

    names(pairwise_clone_vectors) <- map(pairwise_clone_vectors, ~paste(., collapse = "_v_"))
    
    pairwise_seu_plots <- imap(pairwise_clone_vectors, make_pairwise_plots, seu, sample_id)
  }

  all_seu_plot <- plot_distribution_of_clones_across_clusters(seu, seu_name = sample_id, var_x = "scna", var_y = "clusters", plot_type = "clone")

  plot_path = glue("results/numbat_sridhar/{sample_id}_clone_distribution_filtered.pdf")

  pdf(plot_path, width = 8)
  print(all_seu_plot)
  if(pairwise){
    print(pairwise_seu_plots)
  }
  dev.off()
  
  if(pairwise){
  	plot_tables <-
  		pairwise_seu_plots %>%
  		map("table") %>%
  		set_names(str_replace_all(names(.), "\\n", "_")) %>%
  		set_names(str_replace_all(names(.), "\\s", "_")) %>%
  		map(dplyr::mutate, sample_id = sample_id) %>% 
  		identity()	
  	
  	table_path = write_xlsx(plot_tables, glue("results/{sample_id}_pairwise_clone_distribution.xlsx"))
  	
  } else {
  	cluster_values <-
  		seu@meta.data %>%
  		dplyr::group_by(clusters, scna) %>%
  		dplyr::summarize(value = dplyr::n())
  	
  	all_plot_table <- 
  		cluster_values %>% 
  		dplyr::group_by(scna) %>% 
  		dplyr::summarize(value = sum(value)) %>% 
  		mutate(percent = value/sum(value)) %>%
  		dplyr::select(-value) %>%
  		# dplyr::mutate(scna = na_if(scna, "")) %>%
  		# dplyr::mutate(scna = replace_na(scna, "diploid")) %>%
  		dplyr::mutate(clusters = "all") %>% 
  	identity()
  	
  	plot_tables <- 
  		cluster_values %>% 
  		mutate(percent = value/sum(value)) %>%
  		dplyr::select(-value) %>%
  		# dplyr::mutate(scna = na_if(scna, "")) %>%
  		# dplyr::mutate(scna = replace_na(scna, "diploid")) %>%
  		dplyr::bind_rows(all_plot_table) %>% 
  		tidyr::pivot_wider(names_from = "scna", values_from = "percent") %>%
  		dplyr::mutate(sample_id = sample_id) %>% 
  	identity()

  	table_path = write_xlsx(plot_tables, glue("results/{sample_id}_clone_distribution.xlsx"))
  }


  return(
    list(
      "plot" = plot_path,
      "table" = plot_tables
      )
    )

}

# speckle ------------------------------
run_speckle <- function(seu_path = "output/seurat/SRR14800534_filtered_seu.rds", group.by = "SCT_snn_res.0.6") {
  # browser()
  seu <- readRDS(seu_path)
  sample_id = str_extract(seu_path, "SRR[0-9]*")

  if(!"sample_id" %in% colnames(seu@meta.data)){
    seu <- AddMetaData(seu, sample(2, ncol(seu), replace = TRUE), "sample_id")
  }

  seu_meta <-
    seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    dplyr::select(clusters = .data[[group.by]], group = scna, sample_id, `orig.ident`) %>%
    # dplyr::mutate(group = as.numeric(factor(group))) %>%
    dplyr::mutate(sample_group = paste0(sample_id, "_", group)) %>%
    identity()

  # Run propeller testing for cell type proportion differences between the two
  # groups
  mytable <- propeller(clusters = seu_meta$clusters, sample = seu_meta$sample_group,
                       group = seu_meta$group) %>%
    tibble::rownames_to_column("cluster")

  seu_meta$clusters <- factor(seu_meta$clusters, levels = mytable$cluster)

  cell_prop_plot <- plotCellTypeProps(sample=seu_meta$clusters, clusters=seu_meta$group) +
    labs(title = sample_id)

  return(list("table" = mytable, "plot" = cell_prop_plot))
}

run_speckle_for_set <- function(seu_paths, group.by = "SCT_snn_res.0.6"){
  # browser()
  dir_create("results/speckle")
  test0 <- map(seu_paths, safe_run_speckle, group.by = group.by)

  plotlist <- compact(map(test0, c("result", "plot")))

  pdf(glue("results/speckle/{group.by}.pdf"))
  print(plotlist)
  dev.off()

  tables <- compact(map(test0, c("result", "table")))

  write_xlsx(tables, glue("results/speckle/{group.by}.xlsx"))

}

make_clustrees_for_sample <- function(seu_path, mylabel = "sample_id", assay = "SCT", resolutions = seq(0.2, 2.0, by = 0.2), fisher_p_val_threshold = 0.1){
	sample_id <- str_extract(seu_path, "SRR[0-9]*")
	sample_id <- mylabel
	
	seu <- readRDS(seu_path)
	seu@meta.data <- seu@meta.data %>%
		dplyr::mutate(scna = na_if(scna, "")) %>%
		dplyr::mutate(scna = replace_na(scna, ".diploid")) %>%
		identity()

	scna_counts <-
		seu@meta.data %>%
		dplyr::mutate(scna = factor(scna))
	
	scna_clones <- levels(scna_counts$scna)
	
	pairwise_clone_vectors <-
		bind_cols(scna_clones[-length(scna_clones)], scna_clones[-1]) %>%
		t() %>%
		as.data.frame() %>%
		as.list() %>%
		map(as.character) %>%
		identity()
	
	names(pairwise_clone_vectors) <- map(pairwise_clone_vectors, ~paste(., collapse = "_v_"))
	
	possible_make_clustree_for_clone_comparison <- possibly(make_clustree_for_clone_comparison)
	
	clustree_output <- map(pairwise_clone_vectors, ~possible_make_clustree_for_clone_comparison(seu, sample_id, .x, mylabel, assay, resolutions, fisher_p_val_threshold))
	
	return(clustree_output)
	
}

make_clustree_for_clone_comparison <- function(seu, sample_id, clone_set, mylabel = "sample_id", assay = "SCT", resolutions = seq(0.2, 2.0, by = 0.2), fisher_p_val_threshold = 0.1){
	
	seu <- seu[,seu$scna %in% clone_set]
	
	resolutions <-
		glue("{assay}_snn_res.{resolutions}") %>%
		set_names(.)
	
	seu_meta <- seu@meta.data %>%
		dplyr::mutate(scna = na_if(scna, "")) %>%
		dplyr::mutate(scna = replace_na(scna, "diploid")) %>%
		identity()
	
	speckle_proportions <- map(resolutions, ~seu_meta[,c(.x, "scna")]) %>%
		map(set_names, c("cluster", "scna")) %>%
		map(janitor::tabyl, cluster, scna) %>%
		map(janitor::adorn_percentages) %>%
		map(dplyr::rename, samples = cluster) %>%
		imap(~dplyr::mutate(.x, clusters = paste0(.y, "C", samples))) %>%
		dplyr::bind_rows() %>%
		janitor::clean_names() %>%
		na.omit() %>% 
		identity()
	
	clustree_plot <- clustree::clustree(seu, assay = assay, show_axis = TRUE)
	
	clustree_meta <- seu@meta.data[,str_subset(colnames(seu@meta.data), "SCT_snn.res.*")]
	
	clustree_graph <- clustree:::build_tree_graph(clusterings = clustree_meta,
												  prefix = "SCT_snn_res.",
												  metadata = clustree_meta,
												  node_aes_list = list(colour = list(value = "SCT_snn_res.", aggr = NULL), size = list(
												  	value = "size", aggr = NULL), alpha = list(value = 1, aggr = NULL)),
												  prop_filter = 0.1,
												  count_filter = 0)
	
	daughter_clusters <-
		clustree_graph %>%
		tidygraph::activate(edges) %>%
		data.frame() %>%
		dplyr::group_by(from_SCT_snn_res., from_clust) %>%
		dplyr::arrange(from_SCT_snn_res., from_clust, to_clust) %>%
		dplyr::filter(n_distinct(to_clust) > 1) %>%
		dplyr::mutate(from_SCT_snn_res. = as.character(from_SCT_snn_res.)) %>%
		dplyr::mutate(from_clust = as.character(from_clust)) %>%
		split(.$from_SCT_snn_res.) %>%
		map(~split(.x, .x$from_clust)) %>%
		identity()
	
	for(resolution in names(daughter_clusters)){
		for(from_clust in names(daughter_clusters[[resolution]])){
			daughter_clusters[[resolution]][[from_clust]] = chi_sq_daughter_clusters(seu, daughter_clusters, resolution = resolution, from_clust = from_clust)
		}
	}
	
	daughter_clusters <-
		daughter_clusters %>%
		map(dplyr::bind_rows) %>%
		dplyr::bind_rows() %>%
		dplyr::filter(p.value < fisher_p_val_threshold) %>%
		# dplyr::group_by(clone_comparison) %>%
		dplyr::ungroup() %>%
		dplyr::arrange(clone_comparison, p.value) %>%
		dplyr::slice_min(p.value, n = 20, by = clone_comparison) %>%
		dplyr::distinct(from_clust, to_clust, clone_comparison, from_SCT_snn_res., .keep_all = TRUE) %>%
		dplyr::select(all_of(c("to_clust", "to_SCT_snn_res.", "from_clust", "from_SCT_snn_res.",
							   "count", "in_prop", "clone_comparison", "p.value", "method"))) %>%
		dplyr::group_by(from_clust, from_SCT_snn_res., clone_comparison) %>%
		dplyr::mutate(to_clust = paste(to_clust, collapse = "_")) %>%
		dplyr::distinct(from_clust, from_SCT_snn_res., clone_comparison, .keep_all = TRUE) %>%
		tidyr::pivot_wider(names_from = "clone_comparison", values_from = "p.value") %>%
		identity()
	
	clone_comparisons <- str_subset(colnames(daughter_clusters), "_v_") %>%
		set_names(.)
	
	clustree_plot$data <-
		dplyr::left_join(clustree_plot$data, speckle_proportions, by = c("node" = "clusters")) %>%
		dplyr::left_join(daughter_clusters, by = c("SCT_snn_res." = "from_SCT_snn_res.", "cluster" = "from_clust")) %>%
		# dplyr::mutate(signif = ifelse(is.na(method), 0, 1)) %>%
		identity()
	
	clustree_plot$layers[[2]] <- NULL
	
	clustree_res <- map(clone_comparisons, plot_clustree_per_comparison, clustree_plot, speckle_proportions, sample_id)
	
	clustree_plots <- map(clustree_res, "plot")
	
	clustree_plot_path = glue("results/{mylabel}_{clone_comparisons}_clustree.pdf")
	
	pdf(clustree_plot_path, width = 8, height = 10)
	print(clustree_plots)
	dev.off()
	
	clustree_table_path = glue("results/{mylabel}_{clone_comparisons}_clustree.xlsx")
	
	clustree_tables <- map(clustree_res, "table") %>%
		purrr::flatten()
	
	writexl::write_xlsx(clustree_tables, clustree_table_path)
	
	return(clustree_plot_path)
	
}

plot_clustree_per_comparison <- function(mylabel, clustree_plot, speckle_proportions, sample_id) {
	# browser()
	
	comparison_clones <- str_split(mylabel, pattern = "_v_") %>%
		unlist()
	
	comparison_clones[comparison_clones == "x"] <- "diploid"
	
	clone_names <- colnames(speckle_proportions)[!colnames(speckle_proportions) %in% c("samples", "clusters")]
	
	brewer_palettes <- c("Reds", "Greens", "Blues", "Purples", "Oranges")
	
	clone_colors <- brewer_palettes[1:length(clone_names)] %>%
		set_names(clone_names)
	
	clone_colors <- clone_colors[names(clone_colors) %in% comparison_clones]
	
	clustree_res <- imap(clone_colors, ~color_clustree_by_clone(clustree_plot, .x, .y, mylabel = mylabel, sample_id))
	
	return(
		list(
			"plot" = map(clustree_res, "plot"),
			"table" = map(clustree_res, "table")
		))
	
}

color_clustree_by_clone <- function(clustree_plot, mycolor, myclone, mylabel = "asdf", sample_id){
	# browser()
	
	label_data <-
		clustree_plot$data %>%
		dplyr::filter(!!sym(mylabel) < 0.05) %>%
		# dplyr::filter(x12p_16q_1q_2p_v_x12p_16q_1q_2p_11p_8p < 0.05) %>%
		identity()
	
	myplot <- clustree_plot +
		ggraph::geom_node_point(aes(colour = .data[[myclone]], size = size)) +
		ggraph::geom_node_text(aes(label = cluster)) +
		ggraph::geom_node_label(data = label_data, aes(label = cluster)) +
		labs(title = mylabel, subtitle = myclone, caption = sample_id, colour = "clone %") +
		scale_color_distiller(palette = mycolor, direction = 1) +
		NULL
	
	return(list("plot" = myplot, "table" = clustree_plot$data))
}

chi_sq_daughter_clusters <- function(seu, daughter_clusters, resolution, from_clust) {
	# browser()
	from_resolution = glue("SCT_snn_res.{resolution}")
	
	to_clusts <- daughter_clusters[[as.character(resolution)]][[as.character(from_clust)]][["to_clust"]]
	to_resolution <- unique(daughter_clusters[[as.character(resolution)]][[as.character(from_clust)]][["to_SCT_snn_res."]])
	to_resolution = glue("SCT_snn_res.{to_resolution}")
	
	message(glue("resolution: {from_resolution} from_clust: {from_clust}"))
	test_seu <- seu[,seu[[]][[to_resolution]] %in% to_clusts]
	
	test_seu$scna <- janitor::make_clean_names(test_seu$scna, allow_dupes = TRUE)
	
	test_seu$clusters <- as.numeric(test_seu@meta.data[[to_resolution]])
	
	if(length(unique(test_seu$scna)) > 1){
		scna_counts <-
			test_seu@meta.data %>%
			dplyr::mutate(scna = factor(scna))
		
		scna_clones <- levels(scna_counts$scna)
		
		pairwise_clone_vectors <-
			bind_cols(scna_clones[-length(scna_clones)], scna_clones[-1]) %>%
			t() %>%
			as.data.frame() %>%
			as.list() %>%
			map(as.character) %>%
			identity()
		
		names(pairwise_clone_vectors) <- map(pairwise_clone_vectors, ~paste(., collapse = "_v_"))
		
		scna_counts <-
			scna_counts %>%
			janitor::tabyl(clusters, scna) %>%
			tibble::column_to_rownames("clusters") %>%
			as.matrix() %>%
			identity()
		
		fisher_results <- pairwise_clone_vectors %>%
			map(~scna_counts[,.x]) %>%
			map(fisher.test, simulate.p.value=TRUE, B=1e5) %>%
			map(broom::tidy) %>%
			dplyr::bind_rows(.id = "clone_comparison") %>%
			dplyr::mutate(from_clust = from_clust)
		
		test0 <- dplyr::left_join(daughter_clusters[[resolution]][[from_clust]], fisher_results, by = "from_clust")
		
	} else {
		return(daughter_clusters[[resolution]][[from_clust]])
	}
	
}

pull_cluster_orders <- function(raw_cluster_file = "data/raw_cluster_ids.csv"){
	
  cluster_ids <-
    raw_cluster_file %>%
    read_csv() %>%
    janitor::clean_names() %>%
    group_by(resolution, sample_id) %>%
    dplyr::transmute(across(-any_of(c("resolution", "sample_id")), as.character))
  
  cluster_id_list <- 
  	cluster_ids %>%
    tidyr::pivot_longer(-any_of(c("resolution", "sample_id", "tumor_id", "preferred_resolution")), names_to = "phase", values_to = "clusters") %>%
    dplyr::mutate(clusters = str_split(clusters, pattern = "_")) %>%
    tidyr::unnest(clusters) %>%
    dplyr::mutate(phase = factor(phase, levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star"))) %>%
    split(.$sample_id) %>%
  	map(~split(.x, .x$preferred_resolution)) %>%
    identity()
  
  return(cluster_id_list)

}

plot_seu_marker_heatmap <- function(seu_path = NULL, cluster_order = NULL, nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", label = "_filtered_", height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL) {

  kept_phases <- kept_phases %||% phase_levels

  # browser()
  
  tumor_id <- str_extract(seu_path, "SRR[0-9]*")
  
  sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
  
  message(sample_id)
  cluster_order = cluster_order[[sample_id]]

  seu <- readRDS(seu_path)

  nb_paths <- nb_paths %>% 
  	set_names(str_extract(., "SRR[0-9]*"))
  
  nb_path <- nb_paths[[tumor_id]]
  
  if(!is.null(cluster_order)){
    group.by = unique(cluster_order$resolution)
  }

  if(equalize_scna_clones){
    seu_meta <- seu@meta.data %>%
      tibble::rownames_to_column("cell")

    clones <- table(seu_meta$scna)

    min_clone_num <- clones[which.min(clones)]

    selected_cells <-
      seu_meta %>%
      dplyr::group_by(scna) %>%
      slice_sample(n = min_clone_num) %>%
      pull(cell)

    seu <- seu[,selected_cells]

  }

  if(!is.null(cluster_order)){
  	
  	group.by = unique(cluster_order$resolution)
  	
  	cluster_order <-
  		cluster_order %>%
  		dplyr::mutate(order = dplyr::row_number()) %>%
  		dplyr::filter(!is.na(clusters)) %>%
  		dplyr::mutate(clusters = as.character(clusters))

    seu@meta.data$clusters = seu@meta.data[[group.by]]

    seu_meta <- seu@meta.data %>%
    	dplyr::select(-any_of(c("phase_level", "order"))) %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::left_join(cluster_order, by = "clusters") %>%
      dplyr::select(-clusters) %>%
      dplyr::rename(phase_level = phase) %>%
      identity()

    phase_levels = phase_levels[phase_levels %in% unique(seu_meta$phase_level)]

    seu_meta <-
      seu_meta %>%
      tidyr::unite("clusters", all_of(c("phase_level", group.by)), remove = FALSE) %>%
      dplyr::arrange(phase_level, order) %>%
    	# dplyr::arrange(phase_level) %>%
      dplyr::mutate(clusters = factor(clusters, levels = unique(clusters))) %>%
      tibble::column_to_rownames("cell") %>%
      identity()

    seu@meta.data <- seu_meta[rownames(seu@meta.data),]

    seu <-
      seu[,seu$phase_level %in% kept_phases] %>%
      find_all_markers(metavar = "clusters", seurat_assay = "SCT") %>%
      identity()

    seu@meta.data$clusters <- forcats::fct_drop(seu@meta.data$clusters)

    heatmap_features  <-
      seu@misc$markers[["clusters"]][["presto"]] %>%
      dplyr::filter(Gene.Name %in% VariableFeatures(seu))

    tidy_eval_arrange <- function(.data, ...) {
      .data %>%
        arrange(...)
    }

    cluster_order_vec <-
      seu@meta.data %>%
      dplyr::select(clusters, !!group.by) %>%
      dplyr::arrange(clusters, !!sym(group.by)) %>%
      dplyr::pull(!!group.by) %>%
      unique() %>%
      as.character() %>%
      identity()

    heatmap_features[["Cluster"]] <-
      factor(heatmap_features[["Cluster"]], levels = levels(seu_meta$clusters))

    heatmap_features <-
      heatmap_features %>%
      dplyr::arrange(Cluster) %>%
      dplyr::group_by(Cluster) %>%
      slice_max(Average.Log.Fold.Change, n = 5) %>%
      identity()

  } else {

    heatmap_features  <-
      seu@misc$markers[[group.by]][["presto"]]

    cluster_order <- levels(seu@meta.data[[group.by]]) %>%
      set_names(.)

    seu@meta.data[[group.by]] <-
      factor(seu@meta.data[[group.by]], levels = cluster_order)

    group_by_clusters <- seu@meta.data[[group.by]]

    seu@meta.data$clusters <- names(cluster_order[group_by_clusters])

    seu@meta.data$clusters <- factor(seu@meta.data$clusters, levels = unique(setNames(names(cluster_order), cluster_order)[levels(seu@meta.data[[group.by]])]))

    heatmap_features <-
      heatmap_features %>%
      dplyr::arrange(Cluster) %>%
      group_by(Cluster) %>%
      slice_head(n=6) %>%
      dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
      dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
      identity()
  }

  seu$scna <- factor(seu$scna)
  levels(seu$scna)[1] <- "none"

  heatmap_features <-
    heatmap_features %>%
    dplyr::ungroup() %>%
    left_join(giotti_genes, by = c("Gene.Name" = "symbol")) %>%
    # select(Gene.Name, term) %>%
    dplyr::mutate(term = replace_na(term, "")) %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE)

  row_ha = ComplexHeatmap::rowAnnotation(term = rev(heatmap_features$term))

  seu_heatmap <- ggplotify::as.ggplot(
    seu_complex_heatmap(seu,
                        features = heatmap_features$Gene.Name,
                        group.by = c("G2M.Score", "S.Score", "scna", "clusters"),
                        col_arrangement = c("clusters", "scna"),
                        cluster_rows = FALSE,
                        column_split =  sort(seu@meta.data$clusters),
                        row_split = rev(heatmap_features$Cluster),
                        row_title_rot = 0,
                        # row_split = sort(seu@meta.data$clusters)
    )) +
    labs(title = sample_id) +
    theme()


  # browser()
  labels <- data.frame(clusters=unique(seu[[]][["clusters"]]), label =unique(seu[[]][["clusters"]])) %>%
    # dplyr::rename({{group.by}} := cluster) %>%
    identity()

  cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))

  centroid_data <-
    cc_data %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
    dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
    dplyr::mutate(centroid = "centroids") %>%
    identity()

  centroid_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["clusters"]])) +
    geom_point(size = 0.1) +
    theme_light() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
    NULL


  facet_cell_cycle_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
    geom_point(size = 0.1) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
    facet_wrap(~.data[["clusters"]], ncol = 2) +
    theme_light() +
    geom_label(data = labels,
               aes(label = label),
               # x = Inf,
               # y = -Inf,
               x = max(cc_data$S.Score)+0.05,
               y = max(cc_data$G2M.Score)-0.1,
               hjust=1,
               vjust=1,
               inherit.aes = FALSE) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    # guides(color = "none") +
    NULL

  appender <- function(string) str_wrap(string, width = 40)

  labels <- data.frame(scna=unique(seu$scna), label=str_replace(unique(seu$scna), "^$", "diploid"))

  clone_distribution_plot <-
    plot_distribution_of_clones_across_clusters(seu, tumor_id, var_x = "scna", var_y = "clusters")

  umap_plots <- DimPlot(seu, group.by = c("scna", "clusters"), combine = FALSE) %>% 
  	# map(~(.x + theme(legend.position = "bottom"))) %>% 
  	wrap_plots(ncol = 1)
  
  if(!is.null(nb_path)){
    clone_tree_plot <-
      plot_clone_tree(seu, tumor_id, nb_path, clone_simplifications, sample_id = sample_id, legend = FALSE, horizontal = TRUE)

    collage_plots <- list(seu_heatmap, facet_cell_cycle_plot, umap_plots, clone_distribution_plot, clone_tree_plot)

    layout <- "
            AAAAAAAAAAAEEEECCCC
            AAAAAAAAAAABBBBCCCC
            AAAAAAAAAAABBBBDDDD
            AAAAAAAAAAABBBBDDDD
            AAAAAAAAAAABBBBDDDD
            "

    wrap_plots(collage_plots) +
      # plot_layout(widths = c(16, 4)) +
      plot_layout(design = layout) +
    	plot_annotation(tag_levels = "A") +
      NULL

  } else {
    layout <- "
            AAAAAAAAAABBBBCCCC
            AAAAAAAAAABBBBCCCC
            AAAAAAAAAABBBBDDDD
            AAAAAAAAAABBBBDDDD
            AAAAAAAAAABBBBDDDD
    "

    collage_plots <- list(seu_heatmap, facet_cell_cycle_plot, centroid_plot, clone_distribution_plot)

    wrap_plots(collage_plots) +
      # plot_layout(widths = c(16, 4)) +
      plot_layout(design = layout) +
    	plot_annotation(tag_levels = "A") +
      NULL


  }
  
  file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")

  ggsave(glue("results/{file_slug}_{label}heatmap_phase_scatter_patchwork.pdf"), height = height, width = width)

}

make_faded_umap_plots <- function(full_seu, retained_clones, group_by = "clusters"){
	scna_plot <- DimPlot(full_seu, group.by = c("scna")) + 
		aes(alpha = alpha_var) + 
		NULL
	scna_plot[[1]]$layers[[1]]$aes_params$alpha =  ifelse ( full_seu@meta.data$clone_opt %in% retained_clones, 1, .05 )
	
	cluster_plot <- DimPlot(full_seu, group.by = group_by) + 
		aes(alpha = alpha_var) + 
		NULL
	cluster_plot[[1]]$layers[[1]]$aes_params$alpha =  ifelse ( full_seu@meta.data$clone_opt %in% retained_clones, 1, .05 )
	
	umap_plots <- scna_plot + 
		cluster_plot + 
		plot_layout(ncol = 1)
	
	return(umap_plots)
	
}

plot_seu_marker_heatmap_by_scna_every_res <- function(resolution, ...){
	plot_seu_marker_heatmap_by_scna(resolution = resolution, ...)
}

plot_seu_marker_heatmap_by_scna <- function(seu_path = NULL, cluster_order = NULL, nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", resolution = 1, height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, rb_scna_samples, large_clone_comparisons, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, split_columns = "clusters") {
	
	kept_phases <- kept_phases %||% phase_levels
	
	# browser()
	
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	message(sample_id)
	cluster_order = cluster_order[[sample_id]]
	
	full_seu <- readRDS(seu_path)
	
	# subset by retained clones ------------------------------
	clone_comparisons = names(large_clone_comparisons[[sample_id]])
	clone_comparison = clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
	retained_clones <- clone_comparison %>% 
		str_extract("[0-9]_v_[0-9]") %>%
		str_split("_v_", simplify = TRUE)
	
	nb_paths <- nb_paths %>% 
		set_names(str_extract(., "SRR[0-9]*"))
	
	nb_path <- nb_paths[[tumor_id]]
	
	plot_paths <- vector(mode = "list", length = length(cluster_order))
	names(plot_paths) <- names(cluster_order)
	
	file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	plot_path <- glue("results/{file_slug}_{scna_of_interest}_heatmap_phase_scatter_patchwork.pdf")
	
	pdf(plot_path, height = height, width = width)
	
	for(resolution in names(cluster_order)){
  	# start loop ------------------------------
  	
		seu <- full_seu[,full_seu$clone_opt %in% retained_clones]
		
  	if(!is.null(cluster_order)){
  		
  		single_cluster_order <- cluster_order[[resolution]]
  		
  		group.by = unique(single_cluster_order$resolution)
  	}
  	
  	if(equalize_scna_clones){
  		seu_meta <- seu@meta.data %>%
  			tibble::rownames_to_column("cell")
  		
  		clones <- table(seu_meta$scna)
  		
  		min_clone_num <- clones[which.min(clones)]
  		
  		selected_cells <-
  			seu_meta %>%
  			dplyr::group_by(scna) %>%
  			slice_sample(n = min_clone_num) %>%
  			pull(cell)
  		
  		seu <- seu[,selected_cells]
  		
  	}
  	
  	
  	if(!is.null(single_cluster_order)){
  		
  		single_cluster_order <- 
  			single_cluster_order |> 
  			dplyr::mutate(order = dplyr::row_number()) %>%
  			dplyr::filter(!is.na(clusters)) %>%
  			dplyr::mutate(clusters = as.character(clusters))	
  		
  		group.by = unique(single_cluster_order$resolution)
  		
  		seu@meta.data$clusters = seu@meta.data[[group.by]]
  		
  		seu_meta <- seu@meta.data %>%
  			tibble::rownames_to_column("cell") %>%
  			dplyr::select(-any_of(c("phase_level", "order"))) %>%
  			dplyr::left_join(single_cluster_order, by = "clusters") %>%
  			dplyr::select(-clusters) %>%
  			dplyr::rename(phase_level = phase) %>%
  			identity()
  		
  		phase_levels = phase_levels[phase_levels %in% unique(seu_meta$phase_level)]
  		
  		seu_meta <-
  			seu_meta %>%
  			tidyr::unite("clusters", all_of(c("phase_level", group.by)), remove = FALSE) %>%
  			dplyr::arrange(phase_level, order) %>%
  			dplyr::mutate(clusters = factor(clusters, levels = unique(clusters))) %>%
  			tibble::column_to_rownames("cell") %>%
  			identity()
  		
  		seu@meta.data <- seu_meta[rownames(seu@meta.data),]
  		
  		seu <-
  			seu[,seu$phase_level %in% kept_phases] %>%
  			find_all_markers(metavar = "clusters", seurat_assay = "SCT") %>%
  			identity()
  		
  		seu@meta.data$clusters <- forcats::fct_drop(seu@meta.data$clusters)
  		
  		# mysec ------------------------------
  		
  		heatmap_features  <-
  			seu@misc$markers[["clusters"]][["presto"]] %>%
  			dplyr::filter(Gene.Name %in% VariableFeatures(seu))
  		
  		tidy_eval_arrange <- function(.data, ...) {
  			.data %>%
  				arrange(...)
  		}
  		# browser()
  		single_cluster_order_vec <-
  			seu@meta.data %>%
  			dplyr::select(clusters, !!group.by) %>%
  			dplyr::arrange(clusters, !!sym(group.by)) %>%
  			dplyr::select(clusters, !!group.by) |> 
  			dplyr::distinct(.data[[group.by]], .keep_all = TRUE) |> 
  			dplyr::mutate(!!group.by := as.character(.data[[group.by]])) |> 
  			tibble::deframe() |> 
  			identity()
  		
  		heatmap_features[["Cluster"]] <-
  			factor(heatmap_features[["Cluster"]], levels = levels(seu_meta$clusters))
  		
  		heatmap_features <-
  			heatmap_features %>%
  			dplyr::arrange(Cluster) %>%
  			dplyr::group_by(Cluster) %>%
  			slice_max(Average.Log.Fold.Change, n = 5) %>%
  			identity()
  		
  	} else {
  		
  		heatmap_features  <-
  			seu@misc$markers[[group.by]][["presto"]]
  		
  		single_cluster_order <- levels(seu@meta.data[[group.by]]) %>%
  			set_names(.)
  		
  		seu@meta.data[[group.by]] <-
  			factor(seu@meta.data[[group.by]], levels = single_cluster_order)
  		
  		group_by_clusters <- seu@meta.data[[group.by]]
  		
  		seu@meta.data$clusters <- names(single_cluster_order[group_by_clusters])
  		
  		seu@meta.data$clusters <- factor(seu@meta.data$clusters, levels = unique(setNames(names(single_cluster_order), single_cluster_order)[levels(seu@meta.data[[group.by]])]))
  		
  		heatmap_features <-
  			heatmap_features %>%
  			dplyr::arrange(Cluster) %>%
  			group_by(Cluster) %>%
  			slice_head(n=6) %>%
  			dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
  			dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
  			identity()
  	}
  	
  	large_enough_clusters <- 
  		seu@meta.data %>% 
  		dplyr::group_by(clusters) %>% 
  		dplyr::count()
  	
  	large_enough_clusters <- 
  		large_enough_clusters	%>% 
  		dplyr::filter(n >= min_cells_per_cluster) %>% 
  		dplyr::pull(clusters)
  	
  	seu <- seu[,seu$clusters %in% large_enough_clusters]
  	
  	seu$scna[seu$scna == ""] <- ".diploid"
  	seu$scna <- factor(seu$scna)
  	# levels(seu$scna)[1] <- "none"
  	
  	heatmap_features <-
  		heatmap_features %>%
  		dplyr::ungroup() %>%
  		left_join(giotti_genes, by = c("Gene.Name" = "symbol")) %>%
  		# select(Gene.Name, term) %>%
  		dplyr::mutate(term = replace_na(term, "")) %>%
  		dplyr::distinct(Gene.Name, .keep_all = TRUE)
  	
  	row_ha = ComplexHeatmap::rowAnnotation(term = rev(heatmap_features$term))
  	
  	if(!is.null(split_columns)){
  		column_split = sort(seu@meta.data[[split_columns]])
  		column_title = unique(column_split)
  	} else {
  		column_split = split_columns
  		column_title = NULL
  	}
  	
  	seu_heatmap <- ggplotify::as.ggplot(
  		seu_complex_heatmap(seu,
  												features = heatmap_features$Gene.Name,
  												group.by = c("G2M.Score", "S.Score", "scna", "clusters"),
  												col_arrangement = c("clusters", "scna"),
  												cluster_rows = FALSE,
  												column_split =  column_split,
  												row_split = rev(heatmap_features$Cluster),
  												row_title_rot = 0,
  												column_title = column_title,
  												column_title_rot = 90
  		)) +
  		labs(title = sample_id) +
  		theme()
  	
  	labels <- data.frame(clusters=unique(seu[[]][["clusters"]]), label =unique(seu[[]][["clusters"]])) %>%
  		# dplyr::rename({{group.by}} := cluster) %>%
  		identity()
  	
  	cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))
  	
  	centroid_data <-
  		cc_data %>%
  		dplyr::group_by(clusters) %>%
  		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
  		dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
  		dplyr::mutate(centroid = "centroids") %>%
  		identity()
  	
  	centroid_plot <-
  		cc_data %>%
  		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["clusters"]])) +
  		geom_point(size = 0.1) +
  		theme_light() +
  		theme(
  			strip.background = element_blank(),
  			strip.text.x = element_blank()
  		) +
  		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
  		guides(fill = "none", color = "none") +
  		NULL
  	
  	
  	facet_cell_cycle_plot <-
  		cc_data %>%
  		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
  		geom_point(size = 0.1) +
  		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
  		facet_wrap(~.data[["clusters"]], ncol = 2) +
  		theme_light() +
  		geom_label(data = labels,
  							 aes(label = label),
  							 # x = Inf,
  							 # y = -Inf,
  							 x = max(cc_data$S.Score)+0.05,
  							 y = max(cc_data$G2M.Score)-0.1,
  							 hjust=1,
  							 vjust=1,
  							 inherit.aes = FALSE) +
  		theme(
  			strip.background = element_blank(),
  			strip.text.x = element_blank()
  		) +
  		# guides(color = "none") +
  		NULL
  	
  	appender <- function(string) str_wrap(string, width = 40)
  	
  	labels <- data.frame(scna=unique(seu$scna), label=str_replace(unique(seu$scna), "^$", "diploid"))
  	
  	# browser()
  	
  	clone_ratio = janitor::tabyl(as.character(seu$scna))$percent[[2]]
  	
  	comparison_scna <-
  		janitor::tabyl(as.character(seu$scna))[2,1]
  	
  	clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
  		seu, seu_name = glue("{tumor_id} {comparison_scna}"), var_x = "scna", var_y = "clusters", signif = TRUE, plot_type = "clone"
  	)
  	
  	# umap_plots <- DimPlot(full_seu, group.by = c("scna", "clusters"), combine = FALSE) %>% 
  	# 	# map(~(.x + theme(legend.position = "bottom"))) %>% 
  	# 	wrap_plots(ncol = 1)
  	# full_seu$clusters
  	# full_seu[[group.by]] <- 
  	full_seu@meta.data[[group.by]] <- factor(full_seu@meta.data[[group.by]], levels = single_cluster_order_vec)
  	levels(full_seu@meta.data[[group.by]]) <- names(single_cluster_order_vec)
  	umap_plots <- make_faded_umap_plots(full_seu, retained_clones, group_by = group.by)
  	
  	if(!is.null(nb_path)){
  		clone_tree_plot <-
  			plot_clone_tree(seu, tumor_id, nb_path, clone_simplifications, sample_id = sample_id, legend = FALSE, horizontal = FALSE)
  		
  		collage_plots <- list(
  			"seu_heatmap" = seu_heatmap, 
  			"facet_cell_cycle_plot" = facet_cell_cycle_plot, 
  			"umap_plots" = umap_plots, 
  			"clone_distribution_plot" = clone_distribution_plot, 
  			"clone_tree_plot" = clone_tree_plot, 
  			"centroid_plot" = centroid_plot)
  		
  		layout <- "
              AAAAAAAAAAAAEEFFCCCC
              AAAAAAAAAAAABBBBCCCC
              AAAAAAAAAAAABBBBDDDD
              AAAAAAAAAAAABBBBDDDD
              AAAAAAAAAAAABBBBDDDD
              "
  		
  		plot_collage <- wrap_plots(collage_plots) +
  			# plot_layout(widths = c(16, 4)) +
  			plot_layout(design = layout) +
  			plot_annotation(tag_levels = "A") +
  			NULL
  		
  	} else {
  		layout <- "
              AAAAAAAAAABBBBCCCC
              AAAAAAAAAABBBBCCCC
              AAAAAAAAAABBBBDDDD
              AAAAAAAAAABBBBDDDD
              AAAAAAAAAABBBBDDDD
      "
  		
  		collage_plots <- list("seu_heatmap" = seu_heatmap, 
  													"facet_cell_cycle_plot" = facet_cell_cycle_plot, 
  													"centroid_plot" = centroid_plot, 
  													"clone_distribution_plot" = clone_distribution_plot)
  		
  		plot_collage <- wrap_plots(collage_plots) +
  			# plot_layout(widths = c(16, 4)) +
  			plot_layout(design = layout) +
  			plot_annotation(tag_levels = "A") +
  			NULL
  		
  		
  	}
  	
  	print(plot_collage)
  	# end loop------------------------------
	}
	
	dev.off()
	
	return(plot_path)
}

plot_seu_marker_heatmap_by_scna_ara <- function(seu_path = NULL, cluster_order = NULL, nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", label = "_filtered_", height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL, rb_scna_samples, large_clone_comparisons, scna_of_interest = "1q", min_cells_per_cluster = 50, return_plots = FALSE, column_split = "clusters") {
	
	kept_phases <- kept_phases %||% phase_levels
	
	# browser()
	
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	message(sample_id)
	cluster_order = cluster_order[[sample_id]]
	
	full_seu <- readRDS(seu_path)
	
	# subset by retained clones ------------------------------
	clone_comparisons = names(large_clone_comparisons[[sample_id]])
	clone_comparison = clone_comparisons[str_detect(clone_comparisons, scna_of_interest)]
	retained_clones <- clone_comparison %>% 
		str_extract("[0-9]_v_[0-9]") %>%
		str_split("_v_", simplify = TRUE)
	seu <- full_seu[,full_seu$clone_opt %in% retained_clones]
	
	
	if(!is.null(cluster_order)){
		
		heatmap_features  <-
			seu@misc$markers[["clusters"]][["presto"]] %>%
			dplyr::filter(Gene.Name %in% VariableFeatures(seu))
		
		heatmap_features[["Cluster"]] <-
			factor(heatmap_features[["Cluster"]], levels = levels(seu$clusters))
		
		heatmap_features <-
			heatmap_features %>%
			dplyr::arrange(Cluster) %>%
			dplyr::group_by(Cluster) %>%
			slice_max(Average.Log.Fold.Change, n = 5) %>%
			identity()
		
	} else {
		
		heatmap_features <-
			heatmap_features %>%
			dplyr::arrange(Cluster) %>%
			group_by(Cluster) %>%
			slice_head(n=6) %>%
			dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
			dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
			identity()
	}
	
	large_enough_clusters <- 
		seu@meta.data %>% 
		dplyr::group_by(clusters) %>% 
		dplyr::count()
	
	large_enough_clusters <- 
		large_enough_clusters	%>% 
		dplyr::filter(n >= min_cells_per_cluster) %>% 
		dplyr::pull(clusters)
	
	seu <- seu[,seu$clusters %in% large_enough_clusters]
	
	seu$scna[seu$scna == ""] <- ".diploid"
	seu$scna <- factor(seu$scna)
	# levels(seu$scna)[1] <- "none"
	
	heatmap_features <-
		heatmap_features %>%
		dplyr::ungroup() %>%
		dplyr::distinct(Gene.Name, .keep_all = TRUE)
	
	if(!is.null(column_split)){
		column_split = sort(seu@meta.data[[column_split]])
		column_title = unique(column_split)
	} else {
		column_title = NULL
	}
	
	# browser()
	labels <- data.frame(clusters=unique(seu[[]][["clusters"]]), label =unique(seu[[]][["clusters"]])) %>%
		# dplyr::rename({{group.by}} := cluster) %>%
		identity()
	
	cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))
	
	centroid_data <-
		cc_data %>%
		dplyr::group_by(clusters) %>%
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
		dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
		dplyr::mutate(centroid = "centroids") %>%
		identity()
	
	centroid_plot_by_cluster <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["clusters"]])) +
		geom_point(size = 0.1) +
		theme_light() +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank(),
			legend.text = element_text(size = 14), # increase legend text size
			legend.title = element_text(size = 16)
		) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
		NULL
	
	centroid_plot_by_scna <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
		geom_point(size = 0.1) +
		theme_light() +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank(),
			legend.text = element_text(size = 14), # increase legend text size
			legend.title = element_text(size = 16)
		) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
		guides(colour = guide_legend(override.aes = list(size=10)),
					 fill = FALSE) +
		NULL
	
	
	facet_cell_cycle_plot <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
		geom_point(size = 0.1) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
		facet_wrap(~.data[["clusters"]], ncol = 2) +
		theme_light() +
		geom_label(data = labels,
							 aes(label = label),
							 # x = Inf,
							 # y = -Inf,
							 x = max(cc_data$S.Score)+0.05,
							 y = max(cc_data$G2M.Score)-0.1,
							 hjust=1,
							 vjust=1,
							 inherit.aes = FALSE) +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank(),
			legend.text = element_text(size = 14), # increase legend text size
			legend.title = element_text(size = 16),
			legend.position = "none"
		) +
		guides(colour = guide_legend(override.aes = list(size=10))) +
		# guides(color = "none") +
		NULL
	
	comparison_scna <-
		janitor::tabyl(as.character(seu$scna))[2,1]
	
	clone_distribution_plot <- plot_distribution_of_clones_across_clusters(
		seu, seu_name = glue("{tumor_id} {comparison_scna}"), var_x = "scna", var_y = "clusters", signif = TRUE, plot_type = "clone"
	)
	
	collage_plots <- list(
		"centroid_plot_by_cluster" = centroid_plot_by_cluster,
		"centroid_plot_by_scna" = centroid_plot_by_scna,
		"facet_cell_cycle_plot" = facet_cell_cycle_plot,
		plot_spacer(),
		"clone_distribution_plot" = clone_distribution_plot,
		plot_spacer())
	
	layout <- "
		EECCCDDDDDD
		EECCCDDDDDD
		AACCCDDDDDD
		AACCCDDDDDD
		BBCCCDDDDDD
		BBCCCDDDDDD
		"
	
	plot_collage <- wrap_plots(collage_plots) +
		# plot_layout(widths = c(16, 4)) +
		plot_layout(design = layout) +
		plot_annotation(tag_levels = "A") +
		NULL
	
	file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	if(is.character(return_plots)){
		plot_path <- ggsave(glue("results/{file_slug}_{scna_of_interest}_{return_plots}_ara.pdf"), collage_plots[[return_plots]], height = height, width = width)
	} else {
		plot_path <- ggsave(glue("results/{file_slug}_{scna_of_interest}_heatmap_phase_scatter_patchwork_ara.pdf"), plot_collage, height = height, width = width)
	}
	
}

plot_seu_marker_heatmap_integrated <- function(seu_path = NULL, cluster_order = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", label = "_filtered_", height = 10, width = 18, equalize_scna_clones = FALSE, phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star"), kept_phases = NULL) {
	
	kept_phases <- kept_phases %||% phase_levels
	
	# browser()
	
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	message(sample_id)
	
	seu <- readRDS(seu_path)
	
	
	if(!is.null(cluster_order)){
		group.by = unique(cluster_order$resolution)
	}
	
	if(equalize_scna_clones){
		seu_meta <- seu@meta.data %>%
			tibble::rownames_to_column("cell")
		
		clones <- table(seu_meta$scna)
		
		min_clone_num <- clones[which.min(clones)]
		
		selected_cells <-
			seu_meta %>%
			dplyr::group_by(scna) %>%
			slice_sample(n = min_clone_num) %>%
			pull(cell)
		
		seu <- seu[,selected_cells]
		
	}
	
	
	if(!is.null(cluster_order)){
		
		group.by = unique(cluster_order$resolution)
		
		cluster_order <-
			cluster_order %>%
			dplyr::mutate(order = dplyr::row_number()) %>%
			dplyr::filter(!is.na(clusters)) %>%
			dplyr::mutate(clusters = as.character(clusters))
		
		seu@meta.data$clusters = seu@meta.data[[group.by]]
		
		seu_meta <- seu@meta.data %>%
			tibble::rownames_to_column("cell") %>%
			dplyr::left_join(cluster_order, by = "clusters") %>%
			dplyr::select(-clusters) %>%
			dplyr::select(-any_of(c("phase_level"))) %>%
			dplyr::rename(phase_level = phase) %>%
			identity()
		
		phase_levels = phase_levels[phase_levels %in% unique(seu_meta$phase_level)]
		
		seu_meta <-
			seu_meta %>%
			tidyr::unite("clusters", all_of(c("phase_level", group.by)), remove = FALSE) %>%
			dplyr::arrange(phase_level, order) %>%
			dplyr::mutate(clusters = factor(clusters, levels = unique(clusters))) %>%
			tibble::column_to_rownames("cell") %>%
			identity()
		
		seu@meta.data <- seu_meta[rownames(seu@meta.data),]
		
		seu <-
			seu[,seu$phase_level %in% kept_phases] %>%
			find_all_markers(metavar = "clusters", seurat_assay = "SCT") %>%
			identity()
		
		seu@meta.data$clusters <- forcats::fct_drop(seu@meta.data$clusters)
		
		heatmap_features  <-
			seu@misc$markers[["clusters"]][["presto"]] %>%
			dplyr::filter(Gene.Name %in% VariableFeatures(seu))
		
		tidy_eval_arrange <- function(.data, ...) {
			.data %>%
				arrange(...)
		}
		
		cluster_order_vec <-
			seu@meta.data %>%
			dplyr::select(clusters, !!group.by) %>%
			dplyr::arrange(clusters, !!sym(group.by)) %>%
			dplyr::pull(!!group.by) %>%
			unique() %>%
			as.character() %>%
			identity()
		
		heatmap_features[["Cluster"]] <-
			factor(heatmap_features[["Cluster"]], levels = levels(seu_meta$clusters))
		
		heatmap_features <-
			heatmap_features %>%
			dplyr::arrange(Cluster) %>%
			dplyr::group_by(Cluster) %>%
			slice_max(Average.Log.Fold.Change, n = 5) %>%
			identity()
		
	} else {
		
		heatmap_features  <-
			seu@misc$markers[[group.by]][["presto"]]
		
		if(!is.ordered(seu@meta.data[[group.by]])){ 
			seu@meta.data[[group.by]] <- factor(as.numeric(seu@meta.data[[group.by]]))
			}
		
		cluster_order <- levels(seu@meta.data[[group.by]]) %>%
			set_names(.)
		
		seu@meta.data[[group.by]] <-
			factor(seu@meta.data[[group.by]], levels = cluster_order)
		
		group_by_clusters <- seu@meta.data[[group.by]]
		
		seu@meta.data$clusters <- names(cluster_order[group_by_clusters])
		
		seu@meta.data$clusters <- factor(seu@meta.data$clusters, levels = unique(setNames(names(cluster_order), cluster_order)[levels(seu@meta.data[[group.by]])]))
		
		heatmap_features <-
			heatmap_features %>%
			dplyr::arrange(Cluster) %>%
			group_by(Cluster) %>%
			slice_head(n=6) %>%
			dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
			dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
			identity()
	}
	
	seu$scna <- factor(seu$scna)
	levels(seu$scna)[1] <- "none"
	
	heatmap_features <-
		heatmap_features %>%
		dplyr::ungroup() %>%
		left_join(giotti_genes, by = c("Gene.Name" = "symbol")) %>%
		# select(Gene.Name, term) %>%
		dplyr::mutate(term = replace_na(term, "")) %>%
		dplyr::distinct(Gene.Name, .keep_all = TRUE)
	
	row_ha = ComplexHeatmap::rowAnnotation(term = rev(heatmap_features$term))
	
	seu_heatmap <- ggplotify::as.ggplot(
		seu_complex_heatmap(seu,
												features = heatmap_features$Gene.Name,
												group.by = c("G2M.Score", "S.Score", "scna", "clusters"),
												col_arrangement = c("clusters", "scna"),
												cluster_rows = FALSE,
												column_split =  sort(seu@meta.data$clusters),
												row_split = rev(heatmap_features$Cluster),
												row_title_rot = 0,
												# row_split = sort(seu@meta.data$clusters)
		)) +
		labs(title = sample_id) +
		theme()
	
	
	# browser()
	labels <- data.frame(clusters=unique(seu[[]][["clusters"]]), label =unique(seu[[]][["clusters"]])) %>%
		# dplyr::rename({{group.by}} := cluster) %>%
		identity()
	
	cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))
	
	centroid_data <-
		cc_data %>%
		dplyr::group_by(clusters) %>%
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
		dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
		dplyr::mutate(centroid = "centroids") %>%
		identity()
	
	centroid_plot <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
		geom_point(size = 0.1) +
		theme_light() +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank()
		) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
		NULL
	
	
	facet_cell_cycle_plot <-
		cc_data %>%
		ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[["scna"]])) +
		geom_point(size = 0.1) +
		geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
		facet_wrap(~.data[["clusters"]], ncol = 2) +
		theme_light() +
		geom_label(data = labels,
							 aes(label = label),
							 # x = Inf,
							 # y = -Inf,
							 x = max(cc_data$S.Score)+0.05,
							 y = max(cc_data$G2M.Score)-0.1,
							 hjust=1,
							 vjust=1,
							 inherit.aes = FALSE) +
		theme(
			strip.background = element_blank(),
			strip.text.x = element_blank()
		) +
		# guides(color = "none") +
		NULL
	
	appender <- function(string) str_wrap(string, width = 40)
	
	labels <- data.frame(scna=unique(seu$scna), label=str_replace(unique(seu$scna), "^$", "diploid"))
	
	clone_distribution_plot <-
		plot_distribution_of_clones_across_clusters(seu, tumor_id, var_x = "scna", var_y = "clusters")
	
	umap_plots <- DimPlot(seu, group.by = c("scna", "clusters"), combine = FALSE) %>% 
		# map(~(.x + theme(legend.position = "bottom"))) %>% 
		wrap_plots(ncol = 1)
	
		layout <- "
            AAAAAAAAAABBBBCCCC
            AAAAAAAAAABBBBCCCC
            AAAAAAAAAABBBBDDDD
            AAAAAAAAAABBBBDDDD
            AAAAAAAAAABBBBDDDD
    "
		
		collage_plots <- list(seu_heatmap, facet_cell_cycle_plot, centroid_plot, clone_distribution_plot)
		
		wrap_plots(collage_plots) +
			# plot_layout(widths = c(16, 4)) +
			plot_layout(design = layout) +
			plot_annotation(tag_levels = "A") +
			NULL

	
	file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	ggsave(glue("results/{file_slug}_{label}heatmap_phase_scatter_patchwork.pdf"), height = height, width = width)
	
}

plot_seu_marker_heatmap_all_resolutions <- function(seu_path = NULL, nb_paths = NULL, clone_simplifications = NULL, cluster_orders = NULL){
	
  sample_id <- str_extract(seu_path, "SRR[0-9]*")
  
  nb_paths <- 
  	nb_paths %>% 
  	set_names(str_extract(., "SRR[0-9]*"))
  
  nb_path <- nb_paths[[sample_id]]
  
  if(!is.null(cluster_orders)){
  	cluster_orders = cluster_orders[[sample_id]]
  	
  	resolutions <- unique(cluster_orders$resolution)
  	
  	cluster_orders <- split(cluster_orders, cluster_orders$resolution)
  	
  	cluster_orders <- map(cluster_orders, ~list(.x)) |> 
  		map(set_names, sample_id)
  	
  	collages <- map2(resolutions, cluster_orders, ~plot_seu_marker_heatmap(seu = seu_path, cluster_order = .y, nb_path = nb_path, clone_simplifications = clone_simplifications, group.by = .x, label = .x))
  	
  } else {
  	resolutions = glue("SCT_snn_res.{seq(0.2, 2.0, by = 0.2)}") %>%
  		set_names(.)
  	
  	collages <- map(resolutions, ~plot_seu_marker_heatmap(seu = seu_path, nb_path = nb_path, clone_simplifications = clone_simplifications, group.by = .x, label = .x))
  	
  }

  file_slug <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")

  qpdf::pdf_combine(collages, glue("results/{file_slug}_collage_all_resolutions.pdf"))
}


compare_markers <- function(seu_path, cluster_order = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", mygene =  "EZH2", label = "_filtered_", height = 14, width = 22, phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star")) {
  # browser()
  sample_id <- str_extract(seu_path, "SRR[0-9]*")
  message(sample_id)
  cluster_order = cluster_order[[sample_id]]

  seu <- readRDS(seu_path)

    group.by = unique(cluster_order$resolution)

    group.by = unique(cluster_order$resolution)

    cluster_order <-
      cluster_order %>%
      dplyr::filter(!is.na(clusters)) %>%
      dplyr::mutate(clusters = as.character(clusters))

    seu@meta.data$clusters = seu@meta.data[[group.by]]

    seu_meta <- seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::left_join(cluster_order, by = "clusters") %>%
      dplyr::select(-clusters) %>%
      dplyr::rename(clusters = phase) %>%
      identity()

    phase_levels = phase_levels[phase_levels %in% unique(seu_meta$clusters)]

    seu_meta <-
      seu_meta %>%
      tidyr::unite("new_clusters", all_of(c("clusters", group.by)), remove = FALSE) %>%
      dplyr::arrange(clusters, new_clusters) %>%
      dplyr::mutate(clusters = factor(new_clusters, levels = unique(new_clusters))) %>%
      tibble::column_to_rownames("cell") %>%
      identity()

    seu@meta.data <- seu_meta[rownames(seu@meta.data),]

    g1_seu <- seu[,str_detect(seu$clusters, "g1_[0-9]")] %>%
      # find_all_markers(metavar = "clusters", seurat_assay = "SCT") %>%
      identity()

  # test0 <-
  #   g1_seu@misc$markers$clusters$presto %>%
  #   dplyr::filter(str_detect(Cluster, "g1_[0-9]")) %>%
  #   group_by(Cluster) %>%
  #   slice_head(n=30) %>%
  #   dplyr::filter(Gene.Name %in% rownames(seu@assays$SCT@scale.data)) %>%
  #   dplyr::filter(Gene.Name %in% VariableFeatures(seu)) %>%
  #   identity()

  return(g1_seu)

}

compare_cluster_continuous_var <- function(seu_path, cluster_order = NULL, continuous_var = "TFF1", group.by = "SCT_snn_res.0.6", assay = "SCT", mygene =  "EZH2", label = "_filtered_", height = 14, width = 22, phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star"), gene_lists = NULL) {
  # browser()
  sample_id <- str_extract(seu_path, "SRR[0-9]*")
  message(sample_id)
  cluster_order = cluster_order[[sample_id]]

  seu <- readRDS(seu_path)

  heatmap_features  <-
    table_cluster_markers(seu, assay = assay)

    group.by = unique(cluster_order$resolution)

    cluster_order <-
      cluster_order %>%
      dplyr::filter(!is.na(clusters)) %>%
      dplyr::mutate(clusters = as.character(clusters))

    seu@meta.data$clusters = seu@meta.data[[group.by]]

    seu_meta <- seu@meta.data %>%
      tibble::rownames_to_column("cell") %>%
      dplyr::left_join(cluster_order, by = "clusters") %>%
      dplyr::select(-clusters) %>%
      dplyr::rename(clusters = phase) %>%
      identity()

    phase_levels = phase_levels[phase_levels %in% unique(seu_meta$clusters)]

    seu_meta <-
      seu_meta %>%
      tidyr::unite("new_clusters", all_of(c("clusters", group.by)), remove = FALSE) %>%
      dplyr::arrange(clusters, new_clusters) %>%
      dplyr::mutate(clusters = factor(new_clusters, levels = unique(new_clusters))) %>%
      tibble::column_to_rownames("cell") %>%
      identity()

    seu@meta.data <- seu_meta[rownames(seu@meta.data),]

    gene_lists <- map(gene_lists, ~(.x[.x %in% VariableFeatures(seu)]))

    seu <- Seurat::AddModuleScore(seu, features = gene_lists, name = "subtype")

    module_names = paste0("subtype", seq(length(gene_lists)))

    seu@meta.data[names(gene_lists)] <- seu@meta.data[module_names]

    cont_tibble <- FetchData(seu, vars = c(continuous_var, "clusters"))

}

collate_clone_distribution_tables <- function(clone_distribution_plots, interesting_samples){
  names(clone_distribution_plots) <- interesting_samples

  cluster_by_clone_tables <- clone_distribution_plots %>%
    map("table") %>%
    map("cluster_by_clone") %>%
    map(dplyr::select, -value) %>%
    map(dplyr::mutate, scna = na_if(scna, "")) %>%
    map(dplyr::mutate, scna = replace_na(scna, "diploid")) %>%
    map(tidyr::pivot_wider, names_from = "scna", values_from = "percent") %>%
    writexl::write_xlsx("results/cluster_by_clone_tables.xlsx") %>%
    identity()

  clone_by_cluster_tables <- clone_distribution_plots %>%
    map("table") %>%
    map("clone_by_cluster") %>%
    map(dplyr::select, -value) %>%
    map(dplyr::mutate, scna = na_if(scna, "")) %>%
    map(dplyr::mutate, scna = replace_na(scna, "diploid")) %>%
    map(tidyr::pivot_wider, names_from = "scna", values_from = "percent") %>%
    writexl::write_xlsx("results/clone_by_cluster_tables.xlsx")

  return(list(
    "cluster_by_clone" = "results/cluster_by_clone_tables.xlsx",
    "clone_by_cluster" = "results/clone_by_cluster_tables.xlsx",
         ))

}


score_clusters_up_down <- function(clone_distribution_plots, interesting_samples, csv_file = "results/cluster_up_down_scorecard.csv"){
  test0 <-
    clone_distribution_plots %>%
    map("table") %>%
    set_names(interesting_samples) %>%
    identity()

  summarize_clusters <- function(df){
    df %>%
      dplyr::select(clusters, up, down) %>%
      tidyr::pivot_longer(-clusters, names_to = "direction", values_to = "whether") %>%
      dplyr::filter(whether == 1) %>%
      dplyr::select(-whether) %>%
      identity()
  }

  test1 <- map(test0,
               ~map(.x, summarize_clusters)
  ) %>%
    map(dplyr::bind_rows, .id = "comparison") %>%
    dplyr::bind_rows(.id = "sample_id") %>%
    write_csv(csv_file) %>%
    identity()
}

plot_clone_tree_from_path <- function(seu_path, nb_paths, clone_simplifications, label = "_clone_tree", ...){
	
	seu <- readRDS(seu_path)
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	nb_paths <- nb_paths %>% 
		set_names(str_extract(., "SRR[0-9]*"))
	
	nb_path <- nb_paths[[tumor_id]]
	
	clone_tree <- plot_clone_tree(seu, tumor_id = tumor_id, nb_path, clone_simplifications, sample_id = sample_id, ...)
	
	return(clone_tree)
}

save_clone_tree_from_path <- function(seu_path, nb_paths, clone_simplifications, label = "_clone_tree", ...){
	
	seu <- readRDS(seu_path)
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	nb_paths <- nb_paths %>% 
		set_names(str_extract(., "SRR[0-9]*"))
	
	nb_path <- nb_paths[[tumor_id]]
	
	plot_clone_tree(seu, tumor_id = tumor_id, nb_path, clone_simplifications, sample_id = sample_id, ...)
	
	plot_path <- ggsave(glue("results/{sample_id}{label}.pdf"), width = 4, height = 4)
	return(plot_path)
}

save_cc_space_plot_from_path <- function(seu_path, clone_simplifications, label = "_clone_tree", ...){
	
	seu <- readRDS(seu_path)
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	plot_cc_space_plot(seu, tumor_id = tumor_id, sample_id = sample_id, ...)
	
	plot_path <- ggsave(glue("results/{sample_id}{label}.pdf"), width = 4, height = 4)
	return(plot_path)
}


plot_clone_tree <- function(seu, tumor_id, nb_path, clone_simplifications, sample_id = NULL, ...){
  
  mynb <- readRDS(nb_path)
  
  mynb$mut_graph <- tidygraph::as_tbl_graph(mynb$mut_graph) %>% 
  	tidygraph::activate(nodes) %>% 
  	dplyr::filter(clone %in% unique(seu$clone_opt)) %>%
  	as.igraph() %>% 
  	identity()
  
  mynb$clone_post <- dplyr::filter(mynb$clone_post, cell %in% colnames(seu))
  
  # mynb$mut_graph <- induced_subgraph(mynb$mut_graph, retained_vertices)
  
  ## clone tree ------------------------------

  nclones <- length(unique(seu$clone_opt))

  mypal <- scales::hue_pal()(nclones) %>%
    set_names(1:nclones)

  rb_scnas = clone_simplifications[[tumor_id]]

  mynb <- simplify_gt(mynb, rb_scnas)
  
  plot_title = ifelse(is.null(sample_id), tumor_id, sample_id)

  mynb$plot_mut_history(pal = mypal, ...) +
    labs(title = plot_title) + 
  	theme(plot.title = element_text(hjust = 0.5))
}


find_diffex_from_clustree <- function(to_SCT_snn_res. = 1, to_clust = "1_10", sample_id, tumor_id, seu, mynb, ...){
	to_clust = str_split(to_clust, pattern = "_") %>% 
		unlist()
	
	to_SCT_snn_res. = glue("SCT_snn_res.{to_SCT_snn_res.}")
	
	divergent_diffex = find_diffex_bw_divergent_clusters(sample_id, tumor_id, seu, mynb, to_SCT_snn_res., to_clust, ...) 
	
	# browser()
	
	divergent_diffex <- 
		divergent_diffex %>% 
		# compact() %>%
		map(dplyr::bind_rows, .id = "location") %>% 
		bind_rows(.id = "clone_comparison") %>%
		dplyr::mutate(sample_id := {{sample_id}}) %>%
		# dplyr::arrange(cluster, p_val_adj) %>%
		identity() %>% 
		group_by(to_clust, clone_comparison, location)
	
	group_names <- 
		divergent_diffex %>% 
		group_keys() %>% 
		dplyr::mutate(plot_label = glue("clusters: {to_clust}; clones: {clone_comparison}; location: {location}")) %>% 
		dplyr::select(plot_label) %>% 
		tibble::deframe()
	
	volcano_plots <- 
		divergent_diffex %>% 
		group_split() %>% 
		set_names(group_names) %>% 
		map(tibble::column_to_rownames, "symbol") %>%
		map(dplyr::mutate, diffex_comparison = to_clust) %>% 
		imap(make_volcano_plots, sample_id = sample_id) %>%
		identity()
	
	pdf_path <- glue("results/divergent_cluster_diffex_{sample_id}_{to_SCT_snn_res.}_{paste(to_clust, collapse = '_')}.pdf")
	pdf(pdf_path)
	print(volcano_plots)
	dev.off()
	
	# enrichment_table <- 
	# 	diffex %>% 
	# 	dplyr::distinct(symbol, .keep_all = TRUE) %>% 
	# 	tibble::column_to_rownames("symbol") %>% 
	# 	dplyr::select(-any_of(colnames(annotables::grch38))) %>% 
	# 	enrichment_analysis() %>% 
	# 	setReadable(org.Hs.eg.db::org.Hs.eg.db, keyType = "auto")
	# 
	# enrichment_plot <- ggplotify::as.ggplot(
	# 	plot_enrichment(enrichment_table)
	# ) + 
	# 	labs(title = glue("{sample_id}_{unique(diffex$to_clust)}_{unique(diffex$to_SCT_snn_res.)}"))
	
	
	# return(list("diffex" = divergent_diffex, "enrichment_table" = enrichment_table, "enrichment_plot" = enrichment_plot))
	
	return(list("diffex" = divergent_diffex, "plot" = pdf_path))
	
}

find_all_diffex_from_clustree <- function(table_set, debranched_seus, ...) {
	
	tumor_id <- str_extract(names(table_set), "SRR[0-9]*")
	message(tumor_id)
	
	sample_id <- names(table_set)
	message(sample_id)
	
	table_set <- table_set[[1]] %>% 
		dplyr::distinct(to_SCT_snn_res., to_clust, .keep_all = TRUE)
	
	debranched_seus <- 
		debranched_seus %>% 
		unlist() %>% 
		set_names(str_remove(fs::path_file(.), "_filtered_seu.rds"))
	
	seu <- readRDS(debranched_seus[[sample_id]])
	
	mynb <- readRDS(glue("output/numbat_sridhar/{tumor_id}_numbat.rds"))
	
	check_table_set <- function(to_clust, to_SCT_snn_res., seu){
		# browser()
		idents <- 
			to_clust %>% 
			str_split(pattern = "_") %>%
			unlist()
		
		to_SCT_snn_res. <- glue("SCT_snn_res.{to_SCT_snn_res.}")
		
		all(idents %in% seu@meta.data[[to_SCT_snn_res.]])
	}
	
	table_set <- table_set %>% 
		dplyr::rowwise() %>% 
		dplyr::mutate(good_set = check_table_set(to_clust, to_SCT_snn_res., seu)) %>% 
		identity()
	
	message("running comparison")
	test0 <- purrr::map2(table_set$to_SCT_snn_res., table_set$to_clust, find_diffex_from_clustree, sample_id, tumor_id, seu, mynb, ...)
	
	return(test0)
}

pull_clustree_tables <- function(clustrees, divergent_cluster_file = "data/clustree_divergent_clusters.csv") {
	
	clustrees <- map(clustrees, purrr::compact)
	
	divergent_clusters <- divergent_cluster_file %>% 
		read_csv() %>% 
		identity()
	
	clustree_tables <- 
		clustrees %>% 
		unlist() %>% 
		str_replace(".pdf", ".xlsx") %>%
		set_names(str_remove(fs::path_file(.), "_x.*")) %>%
		set_names(str_remove(names(.), "_diploid.*")) %>%
		map(myreadxl) %>%
		map(bind_rows, .id = "clone") %>% 
		bind_rows(.id = "sample_id") %>% 
		dplyr::mutate(from_SCT_snn_res. = as.double(SCT_snn_res.), from_clust = cluster) %>% 
		dplyr::mutate(from_clust = as.integer(str_extract(node, "[0-9]$"))) %>%
		dplyr::select(-to_clust) %>%
		dplyr::inner_join(divergent_clusters, by = c("sample_id", "to_SCT_snn_res.", "from_SCT_snn_res.", "from_clust")) %>%
		dplyr::filter(!is.na(from_clust)) %>%
		dplyr::select(to_clust, everything()) %>%
		dplyr::arrange(sample_id, to_clust) %>%
		dplyr::select(-clone_comparison) %>%
		identity()
	
	
	clustree_tables <- 
		clustree_tables %>% 
		tidyr::pivot_longer(contains('_v_'), names_to = "clone_comparison", values_to = "p.value") %>%
		dplyr::filter(!is.na(p.value)) %>%
		dplyr::filter(!is.na(from_SCT_snn_res.)) %>%
		dplyr::mutate(clone = str_remove(clone, "_[0-9]$")) %>% 
		dplyr::select(-clone) %>% 
		dplyr::distinct(.keep_all = TRUE) %>%
		dplyr::select(-starts_with("x")) %>%
		dplyr::select(sample_id, clone_comparison, from_SCT_snn_res., from_clust, to_SCT_snn_res., to_clust, p.value, everything()) %>% 
		dplyr::arrange(sample_id, clone_comparison) %>% 
		split(.$sample_id) %>% 
		identity()
	
	return(clustree_tables)
}

pull_branches <- function(branch_dictionary_file = "data/branch_dictionary.csv"){
	branch_dictionary <- read_csv(branch_dictionary_file) %>% 
		dplyr::mutate(branch_members = str_split(branch_members, "_")) %>% 
		group_by(sample_id) %>% 
		tidyr::nest(branch_termination, branch_members) %>%
		# group_split() %>% 
		# dplyr::select(sample_id, branch_members) %>% 
		tibble::deframe() %>%
		map(deframe) %>% 
		identity()
	
	# branch_dictionary <- branch_dictionary[map(branch_dictionary, length) > 1]
	
	return(branch_dictionary)
}


debranch_seus <- function(filtered_seus, branch_dictionary, ...){
	
	filtered_seus <- filtered_seus %>% 
		purrr::set_names(str_extract(., "SRR[0-9]*"))
	
	debranched_seus <- map(filtered_seus, split_seu_by_branch, branch_dictionary, ...)
	
	debranched_seus <- unlist(debranched_seus)
	
	return(debranched_seus)
}

prep_seu_branch <- function(debranched_seu, ...){
	DefaultAssay(debranched_seu) <- "gene"
	
	debranched_seu <- 
		debranched_seu %>% 
		seurat_preprocess() %>% 
		seurat_reduce_dimensions()
	
	debranched_seu <- FindNeighbors(debranched_seu, dims = 1:30, verbose = FALSE)
	
	debranched_seu <- seurat_cluster(seu = debranched_seu, resolution = seq(0.2, 2.0, by = 0.2),
																	 reduction = "pca")
	
	debranched_seu <- find_all_markers(debranched_seu, seurat_assay = "gene")
	
	DefaultAssay(debranched_seu) <- "SCT"
	
	debranched_seu <- seurat_cluster(seu = debranched_seu, resolution = seq(0.2, 2.0, by = 0.2),
																	 reduction = "pca")
	
	debranched_seu <- find_all_markers(debranched_seu, seurat_assay = "SCT")
	
	debranched_seu <- assign_phase_clusters(debranched_seu, ...)
	
	return(debranched_seu)
	
}

split_seu_by_branch <- function(seu_path, branches, ...){
	
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	branches <-  branches[[tumor_id]]
	
	if(length(branches) == 1){
		seu_paths <- seu_path
	} else {
		seu <- readRDS(seu_path)
		seu_paths <- list()
		for(terminal_clone in names(branches)){
			debranched_seu = seu[,seu$clone_opt %in% branches[[terminal_clone]]]
			
			sample_id = glue("{tumor_id}_branch_{terminal_clone}")
			
			prep_seu_branch(debranched_seu, tumor_id = tumor_id, sample_id = sample_id, ...)
			
			seu_paths[[terminal_clone]] <- str_replace(seu_path, "_filtered_seu.rds", glue("_branch_{terminal_clone}_filtered_seu.rds"))
			saveRDS(debranched_seu, seu_paths[[terminal_clone]])
		}	
	}
	
	return(seu_paths)
}

find_candidate_cis_in_clustree_diffexes <- function(clustree_diffexes){
	
	sample_id_list <- clustree_diffexes %>% 
		map("diffex") %>% 
		map("sample_id") %>% 
		map_chr(unique)
	
	to_clust_list <- clustree_diffexes %>% 
		map("diffex") %>% 
		map("to_clust") %>% 
		map_chr(unique) %>% 
		identity()
	
	table_names <- glue("{sample_id_list}_{to_clust_list}")
	
	test1 <- 
		clustree_diffexes %>% 
		map("diffex") %>% 
	map(function(mytable){
		test0 <- 
			mytable %>% 
			dplyr::arrange(to_clust, avg_log2FC) %>% 
			dplyr::filter(p_val_adj < 0.1, location == "cis") %>% 
			dplyr::slice_max(abs(avg_log2FC), n = 5) %>% 
			identity()
		
	}) %>% 
		set_names(table_names)
	
	table_path <- writexl::write_xlsx(test1, "results/cis_divergent_cluster_diffex.xlsx")
	
	clustree_plots <- clustree_diffexes %>% 
		map("plot")
	
	fs::dir_create("results/cis_clustree_diffex")
	
	cis_clustree_plots <- clustree_plots %>% 
		str_replace(".pdf", "_cis.pdf") %>% 
		str_replace("results", "results/cis_clustree_diffex")
	
	cis_pages_to_extract <- 
	clustree_plots %>% 
		map(qpdf::pdf_length) %>% 
		map(~seq(2, .x, 3))
	
	pmap(list(clustree_plots, cis_pages_to_extract, cis_clustree_plots), qpdf::pdf_subset)
	
	plot_path <- glue("results/divergent_cluster_diffex_cis.pdf")
	
	qpdf::pdf_combine(cis_clustree_plots, plot_path)
	
	return(list("table" = table_path, "plot" = plot_path))
	
}

find_candidate_trans_in_clustree_diffexes <- function(clustree_diffexes, gene_location = "trans"){
	
	sample_id_list <- clustree_diffexes %>% 
		map("diffex") %>% 
		map("sample_id") %>% 
		map_chr(unique)
	
	to_clust_list <- clustree_diffexes %>% 
		map("diffex") %>% 
		map("to_clust") %>% 
		map_chr(unique) %>% 
		identity()
	
	table_names <- glue("{sample_id_list}_{to_clust_list}")
	
	test1 <- 
		clustree_diffexes %>% 
		map("diffex") %>% 
		map(function(mytable){
			test0 <- 
				mytable %>% 
				dplyr::arrange(to_clust, avg_log2FC) %>% 
				dplyr::filter(p_val_adj < 0.1, location == gene_location) %>% 
				dplyr::slice_max(abs(avg_log2FC), n = 5) %>% 
				identity()
			
		}) %>% 
		set_names(table_names)
	
	table_path <- writexl::write_xlsx(test1, "results/trans_divergent_cluster_diffex.xlsx")
	
	clustree_plots <- clustree_diffexes %>% 
		map("plot")
	
	fs::dir_create("results/trans_clustree_diffex")
	
	trans_clustree_plots <- clustree_plots %>% 
		str_replace(".pdf", "_trans.pdf") %>% 
		str_replace("results", "results/trans_clustree_diffex")
	
	trans_pages_to_extract <- 
		clustree_plots %>% 
		map(qpdf::pdf_length) %>% 
		map(~seq(3, .x, 3))
	
	pmap(list(clustree_plots, trans_pages_to_extract, trans_clustree_plots), qpdf::pdf_subset)
	
	plot_path <- glue("results/divergent_cluster_diffex_trans.pdf")
	
	qpdf::pdf_combine(trans_clustree_plots, plot_path)
	
	return(list("table" = table_path, "plot" = plot_path))
	
}

find_candidate_all_in_clustree_diffexes <- function(clustree_diffexes, gene_location = "all"){
	
	sample_id_list <- clustree_diffexes %>% 
		map("diffex") %>% 
		map("sample_id") %>% 
		map_chr(unique)
	
	to_clust_list <- clustree_diffexes %>% 
		map("diffex") %>% 
		map("to_clust") %>% 
		map_chr(unique) %>% 
		identity()
	
	table_names <- glue("{sample_id_list}_{to_clust_list}")
	
	test1 <- 
		clustree_diffexes %>% 
		map("diffex") %>% 
		map(function(mytable){
			test0 <- 
				mytable %>% 
				dplyr::arrange(to_clust, avg_log2FC) %>% 
				dplyr::filter(p_val_adj < 0.1, location == gene_location) %>% 
				dplyr::slice_max(abs(avg_log2FC), n = 5) %>% 
				identity()
			
		}) %>% 
		set_names(table_names)
	
	table_path <- writexl::write_xlsx(test1, "results/all_divergent_cluster_diffex.xlsx")
	
	clustree_plots <- clustree_diffexes %>% 
		map("plot")
	
	fs::dir_create("results/all_clustree_diffex")
	
	all_clustree_plots <- clustree_plots %>% 
		str_replace(".pdf", "_all.pdf") %>% 
		str_replace("results", "results/all_clustree_diffex")
	
	all_pages_to_extract <- 
		clustree_plots %>% 
		map(qpdf::pdf_length) %>% 
		map(~seq(1, .x, 3))
	
	pmap(list(clustree_plots, all_pages_to_extract, all_clustree_plots), qpdf::pdf_subset)
	
	plot_path <- glue("results/divergent_cluster_diffex_all.pdf")
	
	qpdf::pdf_combine(all_clustree_plots, plot_path)
	
	return(list("table" = table_path, "plot" = plot_path))
	
}


plot_seu_clusters_and_markers <- function(seu_path, cluster_order, phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star")){
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	cluster_order = cluster_order[[tumor_id]]
	
	seu <- readRDS(seu_path)
	
	if(!is.null(cluster_order)){
		group.by = unique(cluster_order$resolution)
		
		cluster_order <-
			cluster_order %>%
			dplyr::mutate(order = dplyr::row_number()) %>%
			dplyr::filter(!is.na(clusters)) %>%
			dplyr::mutate(clusters = as.character(clusters))
		
		seu@meta.data$clusters = seu@meta.data[[group.by]]
		
		seu_meta <- seu@meta.data %>%
			tibble::rownames_to_column("cell") %>%
			dplyr::left_join(cluster_order, by = "clusters") %>%
			dplyr::select(-clusters) %>%
			dplyr::rename(clusters = phase) %>%
			identity()
		
		phase_levels = phase_levels[phase_levels %in% unique(seu_meta$clusters)]
		
		seu_meta <-
			seu_meta %>%
			tidyr::unite("new_clusters", all_of(c("clusters", group.by)), remove = FALSE) %>%
			dplyr::arrange(clusters, new_clusters) %>%
			dplyr::mutate(clusters = factor(new_clusters, levels = unique(new_clusters))) %>%
			tibble::column_to_rownames("cell") %>%
			identity()
		
		seu@meta.data <- seu_meta[rownames(seu@meta.data),]
		
	}
	
	dplot <- Seurat::DimPlot(seu, group.by = "clusters")
	
	
	mplot <- plot_markers(seu, metavar = "clusters", marker_method = "presto", return_plotly = FALSE, hide_technical = "all") +
		ggplot2::scale_y_discrete(position = "left") +
		labs(title = sample_id)
	
	mypatch = dplot + mplot
	
	plot_path = glue("results/{sample_id}_clusters_and_markers.pdf")
	
	ggsave(plot_path, mypatch, height = 8, width = 12)
	
	return(
		plot_path
	)
	
}

select_genes_from_arbitrary_diffex <- function(oncoprint_input_by_scna_unfiltered){
	
	samples_1q_without_preceding_16q <- c("SRR14800534", "SRR14800535", "SRR14800536")
	
	arrange_by_recurrence <- function(df){
		df %>% 
			dplyr::group_by(symbol) %>% 
			dplyr::mutate(abs_mean_FC = abs(mean(avg_log2FC))) %>% 
			dplyr::mutate(recurrence = dplyr::n()) %>% 
			dplyr::arrange(desc(recurrence), desc(abs_mean_FC))
	}
	
	
	test0 <- 
		oncoprint_input_by_scna_unfiltered$cis$`1q+` %>% 
		dplyr::filter(sample_id %in% samples_1q_without_preceding_16q) %>% 
		arrange_by_recurrence() %>% 
		identity()
	
}

read_liu_lu_supp_tables <- function(){
	myreadxl <- function(excel_file, skip = 0){
		sheets <- readxl::excel_sheets(excel_file) %>% set_names(.)
		
		map(sheets, ~readxl::read_excel(excel_file, .x, skip = skip))
	}
	
	liu_supp_tables <- "data/liu_lu_supp_data/supp_table_4.xlsx" %>% 
		myreadxl(skip = 1) %>% 
		map(janitor::clean_names) %>% 
		identity()
}


plot_cc <- function(object = NULL, group.by = "gene_snn_res.0.6", assay = "gene", 
					label = "_filtered_", color.by = "batch", mytitle = "Title", 
					faceted = TRUE, pt.size = 1, phase_levels = c("g1", "g1_s", 
																												"s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", 
																												"other", "s_star"), kept_phases = NULL) 
{
	kept_phases <- kept_phases %||% phase_levels
	heatmap_features <- object@misc$markers[[group.by]][["presto"]]
	
	if(is.factor(object@meta.data[[group.by]])){
		cluster_order <- levels(object@meta.data[[group.by]]) %>% 
			set_names(.)
	} else {
		cluster_order <- unique(object@meta.data[[group.by]]) %>% 
			set_names(.)
	}

	
	object@meta.data[[group.by]] <- factor(object@meta.data[[group.by]], 
																				 levels = cluster_order)
	group_by_clusters <- object@meta.data[[group.by]]
	object@meta.data$clusters <- names(cluster_order[group_by_clusters])
	object@meta.data$clusters <- factor(object@meta.data$clusters, 
																			levels = unique(setNames(names(cluster_order), cluster_order)[levels(object@meta.data[[group.by]])]))
	heatmap_features <- heatmap_features %>% dplyr::arrange(Cluster) %>% 
		group_by(Cluster) %>% slice_head(n = 6) %>% dplyr::filter(Gene.Name %in% 
																																rownames(GetAssayData(object, layer = "gene", slot = "scale.data"))) %>% 
		dplyr::filter(Gene.Name %in% VariableFeatures(object)) %>% 
		identity()
	heatmap_features <- heatmap_features %>% dplyr::ungroup() %>% 
		dplyr::distinct(Gene.Name, .keep_all = TRUE)
	object_heatmap <- ggplotify::as.ggplot(seu_complex_heatmap(object, 
																														 features = heatmap_features$Gene.Name, group.by = c("G2M.Score", 
																														 																										"S.Score", "clusters"), col_arrangement = c("clusters"), 
																														 cluster_rows = FALSE, column_split = sort(object@meta.data$clusters), 
																														 row_split = rev(heatmap_features$Cluster), row_title_rot = 0, 
	)) + labs(title = mytitle) + theme()
	labels <- data.frame(clusters = unique(object[[]][["clusters"]]), 
											 label = unique(object[[]][["clusters"]])) %>% identity()
	cc_data <- FetchData(object, c("clusters", "G2M.Score", "S.Score", 
																 "Phase", color.by))
	centroid_data <- cc_data %>% dplyr::group_by(clusters) %>% 
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>% 
		dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>% 
		dplyr::mutate(centroid = "centroids") %>% identity()
	centroid_plot <- cc_data %>% ggplot(aes(x = S.Score, y = G2M.Score, 
																					group = .data[["clusters"]], color = .data[[color.by]])) + 
		geom_point(size = pt.size) + theme_light() + theme(strip.background = element_blank(), 
																											 strip.text.x = element_blank()) + geom_point(data = centroid_data, 
																											 																						 aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), 
																											 																						 size = 6, alpha = 0.7, shape = 23, colour = "black") + 
		NULL
	facet_cell_cycle_plot <- cc_data %>% ggplot(aes(x = S.Score, 
																									y = G2M.Score, group = .data[["clusters"]], color = .data[[color.by]])) + 
		geom_point(size = pt.size) + geom_point(data = centroid_data, 
																						aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), 
																						size = 6, alpha = 0.7, shape = 23, colour = "black") + 
		facet_wrap(~.data[["clusters"]], ncol = 2) + theme_light() + 
		geom_label(data = labels, aes(label = label), x = max(cc_data$S.Score) + 
							 	0.05, y = max(cc_data$G2M.Score) - 0.1, hjust = 1, 
							 vjust = 1, inherit.aes = FALSE) + theme(strip.background = element_blank(), 
							 																				strip.text.x = element_blank()) + NULL
	appender <- function(string) str_wrap(string, width = 40)
	if (faceted) {
		return(facet_cell_cycle_plot)
	}
	else {
		return(centroid_plot)
	}
}

set_final_seus <- function(interesting_samples){
	final_seus <- fs::dir_ls("output/seurat/", regexp = ".*SRR[0-9]*_filtered_seu.rds") %>% 
		set_names(str_extract(., "SRR[0-9]*"))
	
	final_seus[names(final_seus) %in% interesting_samples]
}

run_velocity <- function(debranched_seu){
	adata = chevreul::run_scvelo(debranched_seu)
	chevreul::plot_scvelo(adata, mode = "dynamical")
}


assign_phase_clusters <- function(seu = NULL, tumor_id, sample_id, cluster_orders = NULL, nb_paths = NULL, clone_simplifications = NULL, group.by = "SCT_snn_res.0.6", assay = "SCT", label = "_filtered_", height = 10, width = 18) {
	
	phase_levels <- c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "pm", "hsp", "hypoxia", "other", "s_star")
	
	# browser()
	
	message(sample_id)
	
	cluster_order <-
		cluster_orders[[sample_id]]
	
	nb_paths <- nb_paths %>% 
		set_names(str_extract(., "SRR[0-9]*"))
	
	nb_path <- nb_paths[[tumor_id]]
	
	if(!is.null(cluster_order)){
		group.by = unique(cluster_order$resolution)
	}
	
		group.by = unique(cluster_order$resolution)
		
		
		cluster_order <-
			cluster_order %>%
			dplyr::mutate(order = dplyr::row_number()) %>%
			dplyr::filter(!is.na(clusters)) %>%
			dplyr::mutate(clusters = as.character(clusters))
		
		seu@meta.data$clusters = seu@meta.data[[group.by]]
		
		seu_meta <- seu@meta.data %>%
			tibble::rownames_to_column("cell") %>%
			dplyr::left_join(cluster_order, by = "clusters") %>%
			dplyr::select(-clusters) %>%
			dplyr::select(-any_of(c("phase_level"))) %>%
			dplyr::rename(phase_level = phase) %>%
			identity()
		
		phase_levels = phase_levels[phase_levels %in% unique(seu_meta$phase_level)]
		
		seu_meta <-
			seu_meta %>%
			tidyr::unite("clusters", all_of(c("phase_level", group.by)), remove = FALSE) %>%
			dplyr::arrange(phase_level, order) %>%
			dplyr::mutate(clusters = factor(clusters, levels = unique(clusters))) %>%
			tibble::column_to_rownames("cell") %>%
			identity()
		
		seu@meta.data <- seu_meta[rownames(seu@meta.data),]
		
		seu <- find_all_markers(seu, metavar = "clusters", seurat_assay = "SCT")
		
		seu@meta.data$clusters <- forcats::fct_drop(seu@meta.data$clusters)
		
		return(seu)
		
}

merge_orders <- function(...){
	
	overall_orders <- c(...)
	
	overall_orders <- overall_orders[!duplicated(names(overall_orders))]
	
	overall_orders <- overall_orders[order(names(overall_orders))]
	
}

clean_diffex <- function(diffex, p_val_limit = 0.05, log2fc_limit = 1){
	diffex %>% 
		tibble::rownames_to_column("symbol") %>% 
		dplyr::left_join(annotables::grch38, by = "symbol") %>% 
		dplyr::select(description, everything()) %>% 
		dplyr::distinct(symbol, .keep_all = TRUE) %>% 
		dplyr::filter(abs(avg_log2FC) > log2fc_limit) %>% 
		dplyr::filter(p_val_adj < p_val_limit) %>% 
		dplyr::arrange(desc(avg_log2FC), p_val_adj) %>% 
		identity()
	
}

table_and_plot_enrichment <- function(ident.1, ident.2, scna_of_interest, seu, sample_id) {
	# browser()
	my_diffex <- Seurat::FindMarkers(seu, group.by = "clusters", ident.1 = ident.1, ident.2 = ident.2)
	
	my_enrichment <- 
		my_diffex %>% 
		# clean_diffex() %>% 
		enrichment_analysis(gene_set = "hallmark") %>% 
		DOSE::setReadable(OrgDb = org.Hs.eg.db, keyType="ENTREZID")
	
	core_enriched_genes <- my_enrichment %>% 
		as_tibble() %>% 
		dplyr::mutate(core_enrichment = str_split(core_enrichment, "\\/")) %>% 
		dplyr::rowwise() %>% 
		dplyr::mutate(core_enrichment = list(core_enrichment[core_enrichment %in% rownames(my_diffex)])) %>% 
		dplyr::mutate(core_enrichment = list(paste(core_enrichment, collapse = "/"))) %>% 
		dplyr::pull(core_enrichment) %>% 
		identity()
	
	simplified_enrichment <- my_enrichment
	
	simplified_enrichment@result$core_enrichment <- core_enriched_genes
	
	my_enrichment_plot <- 
		my_enrichment %>% 
		plot_enrichment() + 
		labs(title = glue("{scna_of_interest}: {ident.1} v. {ident.2}"))
	
	netplot <- simplified_enrichment %>% 
		enrichplot::cnetplot(node_label="all", 
											showCategory = 10, 
											# cex.params = list(gene_node = 2), 
											cex_category= 0.5,
											cex_gene= 0.5,
											cex_label_category= 0.6,
											cex_label_gene= 0.3
	) + 
		labs(title = glue("{scna_of_interest}: {ident.1} v. {ident.2}"))
	
	return(
		list("diffex" = clean_diffex(my_diffex),
				 "enrichment_table" = as_tibble(my_enrichment), 
				 "dotplot" = my_enrichment_plot,
				 "netplot" = netplot
		)
	)
}

pull_cluster_comparisons <- function(cluster_comparisons_file){
	cluster_comparisons <- read_csv(cluster_comparisons_file) %>% 
		split(.$sample_id)
}

make_cluster_comparisons_by_phase_for_disctinct_clones <- function(cluster_comparison, debranched_seus){
	
	sample_id <- unique(cluster_comparison[[1]][["sample_id"]])
	
	tumor_id <- str_extract(sample_id, "SRR[0-9]*")
	
	debranched_seus <- debranched_seus %>% 
		set_names(str_extract(., "SRR[0-9]*.*(?=_filtered_seu.rds)"))
	
	seu_path <- debranched_seus[[sample_id]]
	
	test0 <- 
		cluster_comparison[[1]] %>% 
		split(.$scna_of_interest) %>% 
		map(dplyr::select, all_of(c("ident.1", "ident.2", "scna_of_interest"))) %>% 
		as.list() %>% 
		identity()
	
	seu <- readRDS(seu_path)
	
	safe_table_and_plot_enrichment <- safely(table_and_plot_enrichment)
	
	plot_outcome <- map(test0, ~pmap(.x, safe_table_and_plot_enrichment, seu, sample_id))
	
	outfile = glue("results/cluster_comparisons_by_phase_for_disctinct_clones_{sample_id}.pdf")
	
	slugs <- 
		test0 %>% 
		map(dplyr::select, ident.1, ident.2) %>% 
		map(tidyr::unite, "slug", everything(), sep = "-v-") %>% 
		map(dplyr::pull, slug) %>% 
		identity()
		
	plot_outcome <- map2(plot_outcome, slugs, ~set_names(.x, .y))
	
	results <- 
		plot_outcome %>% 
		list_flatten() %>% 
		map("result") %>% 
		purrr::compact() %>%
		map(~.x[c("diffex", "enrichment_table")]) %>%
		imap(~write_xlsx(.x, glue("results/{sample_id}_{.y}.xlsx"))) %>%
		identity()
	
	pdf(outfile, height = 8, width = 8)
	print(plot_outcome)
	dev.off()

	return(outfile)
	
	
}

print_table_tally <- function(myvec){
	table_tally <- table(myvec)
	message(glue("{names(table_tally)}: {table_tally}; "))
}


poster_plot_markers <- function(seu, metavar = "batch", num_markers = 5, selected_values = NULL, 
						return_plotly = FALSE, marker_method = "presto", seurat_assay = "gene", 
						hide_technical = NULL, unique_markers = FALSE, p_val_cutoff = 1, 
						...) 
	{
		Idents(seu) <- seu[[]][[metavar]]
		seu <- find_all_markers(seu, metavar, seurat_assay = seurat_assay, 
														p_val_cutoff = p_val_cutoff)
		marker_table <- seu@misc$markers[[metavar]][[marker_method]]
		markers <- marker_table %>% seuratTools::enframe_markers() %>% dplyr::mutate(dplyr::across(.fns = as.character))
		if (!is.null(hide_technical)) {
			markers <- purrr::map(markers, c)
			if (hide_technical == "pseudo") {
				markers <- purrr::map(markers, ~.x[!.x %in% pseudogenes[[seurat_assay]]])
			}
			else if (hide_technical == "mito_ribo") {
				markers <- purrr::map(markers, ~.x[!str_detect(.x, 
																											 "^MT-")])
				markers <- purrr::map(markers, ~.x[!str_detect(.x, 
																											 "^RPS")])
				markers <- purrr::map(markers, ~.x[!str_detect(.x, 
																											 "^RPL")])
			}
			else if (hide_technical == "all") {
				markers <- purrr::map(markers, ~.x[!.x %in% pseudogenes[[seurat_assay]]])
				markers <- purrr::map(markers, ~.x[!str_detect(.x, 
																											 "^MT-")])
				markers <- purrr::map(markers, ~.x[!str_detect(.x, 
																											 "^RPS")])
				markers <- purrr::map(markers, ~.x[!str_detect(.x, 
																											 "^RPL")])
			}
			min_length <- min(purrr::map_int(markers, length))
			markers <- purrr::map(markers, head, min_length) %>% 
				dplyr::bind_cols()
		}
		if (unique_markers) {
			markers <- markers %>% dplyr::mutate(precedence = row_number()) %>% 
				pivot_longer(-precedence, names_to = "group", values_to = "markers") %>% 
				dplyr::arrange(markers, precedence) %>% dplyr::group_by(markers) %>% 
				dplyr::filter(row_number() == 1) %>% dplyr::arrange(group, 
																														precedence) %>% tidyr::drop_na() %>% dplyr::group_by(group) %>% 
				dplyr::mutate(precedence = row_number()) %>% tidyr::pivot_wider(names_from = "group", 
																																				values_from = "markers") %>% dplyr::select(-precedence)
		}
		sliced_markers <- markers %>% dplyr::slice_head(n = num_markers) %>% 
			tidyr::pivot_longer(everything(), names_to = "group", 
													values_to = "feature") %>% dplyr::arrange(group) %>% 
			dplyr::distinct(feature, .keep_all = TRUE) %>% identity()
		if (!is.null(selected_values)) {
			seu <- seu[, Idents(seu) %in% selected_values]
			sliced_markers <- sliced_markers %>% dplyr::filter(group %in% 
																												 	selected_values) %>% dplyr::distinct(feature, .keep_all = TRUE)
		}
		vline_coords <- head(cumsum(table(sliced_markers$group)) + 
												 	0.5, -1)
		sliced_markers <- dplyr::pull(sliced_markers, feature)
		seu[[metavar]][is.na(seu[[metavar]])] <- "NA"
		Idents(seu) <- metavar
		markerplot <- DotPlot(seu, assay = seurat_assay, features = sliced_markers, 
													group.by = metavar, dot.scale = 3) + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, 
																																																									angle = 45, vjust = 1, hjust = 1), axis.text.y = ggplot2::element_text(size = 10)) + 
			ggplot2::scale_y_discrete(position = "left") + ggplot2::scale_x_discrete(limits = sliced_markers) + 
			ggplot2::geom_vline(xintercept = vline_coords, linetype = 2) + 
			NULL
		if (return_plotly == FALSE) {
			return(markerplot)
		}
		plot_height <- (150 * num_markers)
		plot_width <- (100 * length(levels(Idents(seu))))
		markerplot <- plotly::ggplotly(markerplot, height = plot_height, 
																	 width = plot_width) %>% plotly_settings() %>% plotly::toWebGL() %>% 
			identity()
		return(list(plot = markerplot, markers = marker_table))
	}


calc_silhouette <- function(seu_path, assay = "SCT", reduction = "pca", dims = 1:30, ...) {
	# browser()
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	plot_path <- glue("results/{sample_id}_silhouette.pdf")
	
	seu <- readRDS(seu_path)
	
	resolutions <- glue("{assay}_snn_res.{seq(0.2, 2.0, by = 0.2)}") %>%
		set_names(.)
	
	for(resolution in resolutions){
		dist.matrix <- dist(x = Embeddings(object = seu[[reduction]])[, dims])
		clusters <- seu@meta.data[[resolution]]
		sil <- cluster::silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
		
		sil_res <- str_replace(resolution, "SCT_snn_res.", "sil_")
		
		seu[[sil_res]] <- sil[, 3]	
	}
	plot_list <- vector(mode = "list", length = length(resolutions)) |> 
		set_names(resolutions)
	
	metrics <- vector(mode = "list", length = length(resolutions)) |> 
		set_names(resolutions)
	
	for(resolution in resolutions){
		# browser()
		
		mod_out <- utils::capture.output(Seurat::FindClusters(seu, resolution = as.numeric(str_replace(resolution, "SCT_snn_res.", ""))), type =  "output")
		modularity <- as.numeric(unlist(strsplit(x = mod_out[7], split = ":"))[2])
		
		sil_res <- str_replace(resolution, "SCT_snn_res.", "sil_")
		
		seu_meta <- 
			seu@meta.data |> 
			tibble::rownames_to_column("cell") |> 
			dplyr::arrange(.data[[resolution]], desc(.data[[sil_res]])) |> 
			dplyr::mutate(cell = factor(cell, levels = unique(cell)))
		
		mean_sil <- mean(seu_meta[[sil_res]])
		
		plot_list[[resolution]] = ggplot(data = seu_meta, aes(x = .data[["cell"]], y = .data[[sil_res]], fill = .data[[resolution]])) + 
			geom_col(outlier.size = 0.1)  +
			geom_hline(aes(yintercept = mean_sil)) + 
			xlab("Method") + ylab("Silhoutte Metric") +
			labs(title = glue("{sample_id} {resolution}"), subtitle = glue("silhouette mean: {mean_sil}; modularity: {modularity}")) + 
			NULL
		
		metrics[[resolution]] <- c(
			"resolution" = as.numeric(str_replace(resolution, "SCT_snn_res.", "")),
			"silhouette_mean" = mean_sil
		)
	}
	
	metrics_plot <- dplyr::bind_rows(metrics) |> 
		ggplot(aes(x = resolution, y = silhouette_mean)) +
		geom_point() + 
		geom_line() + 
		labs(title = sample_id)
	
	pdf(plot_path, ...)
	print(metrics_plot)
	print(plot_list)
	dev.off()
	
	return(plot_path)
	
}


find_cc_genes_by_arm <- function(){
	
	cc_genes <- Seurat::cc.genes.updated.2019 %>%
		tibble::enframe("phase_of_gene", "symbol") %>%
		tidyr::unnest(symbol) |> 
		dplyr::left_join(annotables::grch38, by = "symbol") |> 
		dplyr::mutate(seqnames = chr) |> 
		dplyr::filter(!is.na(start)) |> 
		as_granges() |>
		identity()
	
	arms_df <- read_tsv("http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz", 
											col_names = c("chrom","chromStart","chromEnd","name","gieStain")) |> 
		mutate(arm = substring(name, 1, 1)) |> 
		group_by(chrom, arm) |> 
		summarise(start = min(chromStart),
							end = max(chromEnd),
							length = end - start) |> 
		dplyr::mutate(chrom = str_remove(chrom, "chr")) |> 
		dplyr::mutate(seqnames = chrom) |> 
		as_granges() |> 
		join_overlap_intersect(cc_genes) |> 
		as_tibble() |> 
		dplyr::mutate(seqnames = str_pad(seqnames, 2, pad= "0")) |> 
		dplyr::arrange(seqnames, arm) |> 
		dplyr::select(seqnames, arm, everything())
	
}

#' Convert seurat object to seurat V5 format
#'
#' Convert seurat object from v5 to v3 format
#'
#' @param seu_v5 a version 5 seurat object
#'
#' @return
#' @export
#'
#' @examples
#' convert_v5_to_v3(human_gene_transcript_seu)
convert_v5_to_v3 <- function(seu_v5) {
	
		meta <- seu_v5@meta.data
		
		seu_v3 <- CreateSeuratObject(counts = seu_v5$gene@counts, data = seu_v5$gene@data, assay = "gene", meta.data = meta)
		
		# transcript_assay.v5 <- CreateAssay5Object(counts = seu_v3$transcript@counts, data = seu_v3$transcript@data)
		# seu_v5$transcript <- transcript_assay.v5
		
		seu_v3$gene <- seurat_preprocess(seu_v5$gene, normalize = FALSE)
		
		# seu_v5 <- clustering_workflow(seu_v5)
		seu_v3@reductions <- seu_v5@reductions
		seu_v3@graphs <- seu_v5@graphs
		seu_v3@neighbors <- seu_v5@neighbors
		
		seu_v3@misc <- seu_v5@misc
		
		Idents(seu_v3) <- Idents(seu_v5)
	
	return(seu_v3)
}


seu_factor_heatmap <- function(seu, features = NULL, group.by = "ident", cells = NULL, 
															 layer = "scale.data", assay = NULL, group.bar.height = 0.01, 
															 column_split = NULL, col_arrangement = "ward.D2", mm_col_dend = 30, 
															 embedding = "pca", factor_cols = NULL, ...) 
{
	# browser()
	if(!is.null(factor_cols)){
		data <- seu@meta.data[,factor_cols]	
	} else {
		data <- seu@meta.data[,str_detect(colnames(seu@meta.data), pattern = "^k[0-9]*_[0-9]*$")]	
	}
	
	if (any(col_arrangement %in% c("ward.D", "single", "complete", 
																 "average", "mcquitty", "median", "centroid", "ward.D2"))) {
		
		cluster_columns <- function(m) as.dendrogram(cluster::agnes(m), 
																								 method = col_arrangement)
	} else {
		cells <- seu %>% Seurat::FetchData(vars = col_arrangement) %>% 
			dplyr::arrange(across(all_of(col_arrangement))) %>% 
			rownames()
		data <- data[cells, ]
		group.by <- base::union(group.by, col_arrangement)
		cluster_columns <- FALSE
	}
	group.by <- group.by %||% "ident"
	groups.use <- seu[[group.by]][cells, , drop = FALSE]
	# groups.use <- seu[[group.by]]
	groups.use <- groups.use %>% tibble::rownames_to_column("sample_id") %>% 
		dplyr::mutate(across(where(is.character), ~str_wrap(str_replace_all(.x, ",", " "), 10))) %>% 
		dplyr::mutate(across(where(is.character), as.factor)) %>% 
		data.frame(row.names = 1) %>% 
		identity()
	
	groups.use.factor <- groups.use[sapply(groups.use, is.factor)]
	ha_cols.factor <- NULL
	if (length(groups.use.factor) > 0) {
		ha_col_names.factor <- lapply(groups.use.factor, levels)
		ha_cols.factor <- purrr::map(ha_col_names.factor, ~(scales::hue_pal())(length(.x))) %>% 
			purrr::map2(ha_col_names.factor, purrr::set_names)
	}
	groups.use.numeric <- groups.use[sapply(groups.use, is.numeric)]
	ha_cols.numeric <- NULL
	if (length(groups.use.numeric) > 0) {
		numeric_col_fun <- function(myvec, color) {
			circlize::colorRamp2(range(myvec), c("white", color))
		}
		ha_col_names.numeric <- names(groups.use.numeric)
		ha_col_hues.numeric <- (scales::hue_pal())(length(ha_col_names.numeric))
		ha_cols.numeric <- purrr::map2(groups.use[ha_col_names.numeric], 
																	 ha_col_hues.numeric, numeric_col_fun)
	}
	ha_cols <- c(ha_cols.factor, ha_cols.numeric)
	column_ha <- ComplexHeatmap::HeatmapAnnotation(df = groups.use, 
																								 height = grid::unit(group.bar.height, "points"), col = ha_cols)
	
	hm <- ComplexHeatmap::Heatmap(t(data), name = "normalized usage", 
																top_annotation = column_ha, cluster_columns = cluster_columns,
																cluster_rows = FALSE,
																show_column_names = FALSE, column_dend_height = grid::unit(mm_col_dend, 
																																													 "mm"), column_split =  sort(seu@meta.data$clusters), column_title = NULL, 
																...)
	return(hm)
}


annotate_cell_cycle_without_1q <- function(seu, organism = "human", ...){
	
	cc_genes_by_arm <- find_cc_genes_by_arm() |> 
		dplyr::distinct(symbol, .keep_all = TRUE) |> 
		group_by(phase_of_gene) |>
			identity()
	
	cc_wo_1q <- 
		cc_genes_by_arm |> 
		dplyr::filter(!(seqnames == "01" & arm == "q")) |>
		identity()
			
	cc_wo_1q <- split(cc_wo_1q, cc_wo_1q$phase_of_gene) |> 
		map(pull, symbol)
	
	s_genes <- cc_wo_1q$s.genes
	g2m_genes <- cc_wo_1q$g2m.genes
	if (organism == "mouse") {
		s_genes <- dplyr::filter(human_to_mouse_homologs, HGNC.symbol %in% 
														 	s_genes) %>% dplyr::pull(MGI.symbol)
		g2m_genes <- dplyr::filter(human_to_mouse_homologs, HGNC.symbol %in% 
															 	g2m_genes) %>% dplyr::pull(MGI.symbol)
	}
	seu <- CellCycleScoring(seu, s.features = s_genes, g2m.features = g2m_genes, 
													set.ident = FALSE)
}

plot_phase_wo_arm <- function(seu_path, pdf_path = NULL) {
	
	pdf_path <- pdf_path %||% str_replace("seu_path", ".rds", "_cc_wo_1q.pdf")
	
	seu <- readRDS(seu_path)
	
	ccplot1 <- DimPlot(seu, group.by = "Phase") + 
		labs(subtitle = paste(glue("{names(table(seu$Phase))}: {table(seu$Phase)}"), collapse = "; "))
	
	seu <- annotate_cell_cycle_without_1q(seu)
	
	ccplot2 <- DimPlot(seu, group.by = "Phase") +
		labs(subtitle = paste(glue("{names(table(seu$Phase))}: {table(seu$Phase)}"), collapse = "; "))
	
	ccplot1 + ccplot2
	
	ggsave(pdf_path)
}

plot_cc_space_plot <- function(seu_path = "output/seurat/SRR14800534_filtered_seu.rds", facet = FALSE, group_by = "clusters", color_by = "clusters"){
	
	tumor_id <- str_extract(seu_path, "SRR[0-9]*")
	
	sample_id <- str_remove(fs::path_file(seu_path), "_filtered_seu.rds")
	
	seu <- readRDS(seu_path)
	
	cc_data <- FetchData(seu, c("clusters", "G2M.Score", "S.Score", "Phase", "scna"))
	
	centroid_data <-
		cc_data %>%
		dplyr::group_by(clusters) %>%
		dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
		dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
		dplyr::mutate(centroid = "centroids") %>%
		identity()
	
	if(!facet){
		centroid_plot <-
			cc_data %>%
			ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[group_by]], color = .data[[color_by]])) +
			geom_point(size = 0.1) +
			theme_light() +
			theme(
				strip.background = element_blank(),
				strip.text.x = element_blank()
			) +
			geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
			labs(title = sample_id)
		
	} else{
	
		facet_cell_cycle_plot <-
			cc_data %>%
			ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[[group_by]], color = .data[[color_by]])) +
			geom_point(size = 0.1) +
			geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[[group_by]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
			facet_wrap(~.data[["clusters"]], ncol = 2) +
			theme_light() +
			geom_label(data = labels,
								 aes(label = label),
								 # x = Inf,
								 # y = -Inf,
								 x = max(cc_data$S.Score)+0.05,
								 y = max(cc_data$G2M.Score)-0.1,
								 hjust=1,
								 vjust=1,
								 inherit.aes = FALSE) +
			theme(
				strip.background = element_blank(),
				strip.text.x = element_blank()
			) +
			labs(title = sample_id) + 
			# guides(color = "none") +
			NULL
	}
		
		plot_path <- ggsave(glue("results/{sample_id}_cc_space_plot.pdf"))
	
		return(plot_path)
	
}

seu_complex_heatmap <- function(seu, features = NULL, group.by = "ident", cells = NULL, 
																layer = "scale.data", assay = NULL, group.bar.height = 0.01, 
																column_split = NULL, col_arrangement = "ward.D2", mm_col_dend = 30, 
																embedding = "pca", ...) 
{
	if (length(GetAssayData(seu, layer = "scale.data")) == 0) {
		message("seurat object has not been scaled. Please run `Seurat::ScaleData` to view a scaled heatmap; showing unscaled expression data")
		layer <- "data"
	}
	cells <- cells %||% colnames(x = seu)
	if (is.numeric(x = cells)) {
		cells <- colnames(x = seu)[cells]
	}
	assay <- assay %||% Seurat::DefaultAssay(object = seu)
	Seurat::DefaultAssay(object = seu) <- assay
	features <- features %||% VariableFeatures(object = seu)
	features <- rev(x = unique(x = features))
	possible.features <- rownames(x = GetAssayData(object = seu, 
																								 layer = layer))
	if (any(!features %in% possible.features)) {
		bad.features <- features[!features %in% possible.features]
		features <- features[features %in% possible.features]
		if (length(x = features) == 0) {
			stop("No requested features found in the ", layer, 
					 " layer for the ", assay, " assay.")
		}
		warning("The following features were omitted as they were not found in the ", 
						layer, " layer for the ", assay, " assay: ", paste(bad.features, 
																															 collapse = ", "))
	}
	data <- as.data.frame(x = t(x = as.matrix(x = GetAssayData(object = seu, 
																														 layer = layer)[features, cells, drop = FALSE])))
	seu <- suppressMessages(expr = StashIdent(object = seu, 
																						save.name = "ident"))
	if (any(col_arrangement %in% c("ward.D", "single", "complete", 
																 "average", "mcquitty", "median", "centroid", "ward.D2"))) {
		if ("pca" %in% Seurat::Reductions(seu)) {
			cluster_columns <- Seurat::Embeddings(seu, embedding) %>% 
				dist() %>% hclust(col_arrangement)
		}
		else {
			message(glue("{embedding} not computed for this dataset; cells will be clustered by displayed features"))
			cluster_columns <- function(m) as.dendrogram(cluster::agnes(m), 
																									 method = col_arrangement)
		}
	}
	else {
		cells <- seu %>% Seurat::FetchData(vars = col_arrangement) %>% 
			dplyr::arrange(across(all_of(col_arrangement))) %>% 
			rownames()
		data <- data[cells, ]
		group.by <- base::union(group.by, col_arrangement)
		cluster_columns <- FALSE
	}
	group.by <- group.by %||% "ident"
	groups.use <- seu[[group.by]][cells, , drop = FALSE]
	groups.use <- groups.use %>% tibble::rownames_to_column("sample_id") %>% 
		dplyr::mutate(across(where(is.character), ~str_wrap(str_replace_all(.x, 
																																				",", " "), 10))) %>% dplyr::mutate(across(where(is.character), 
																																																									as.factor)) %>% data.frame(row.names = 1) %>% identity()
	groups.use.factor <- groups.use[sapply(groups.use, is.factor)]
	ha_cols.factor <- NULL
	if (length(groups.use.factor) > 0) {
		ha_col_names.factor <- lapply(groups.use.factor, levels)
		ha_cols.factor <- purrr::map(ha_col_names.factor, ~(scales::hue_pal())(length(.x))) %>% 
			purrr::map2(ha_col_names.factor, purrr::set_names)
	}
	groups.use.numeric <- groups.use[sapply(groups.use, is.numeric)]
	ha_cols.numeric <- NULL
	if (length(groups.use.numeric) > 0) {
		numeric_col_fun <- function(myvec, color) {
			circlize::colorRamp2(range(myvec), c("white", color))
		}
		ha_col_names.numeric <- names(groups.use.numeric)
		ha_col_hues.numeric <- (scales::hue_pal())(length(ha_col_names.numeric))
		ha_cols.numeric <- purrr::map2(groups.use[ha_col_names.numeric], 
																	 ha_col_hues.numeric, numeric_col_fun)
	}
	ha_cols <- c(ha_cols.factor, ha_cols.numeric)
	column_ha <- ComplexHeatmap::HeatmapAnnotation(df = groups.use, 
																								 height = grid::unit(group.bar.height, "points"), col = ha_cols)
	hm <- ComplexHeatmap::Heatmap(t(data), name = "log expression", 
																top_annotation = column_ha, cluster_columns = cluster_columns, 
																show_column_names = FALSE, column_dend_height = grid::unit(mm_col_dend, 
																																													 "mm"), column_split = column_split,
																...)
	return(hm)
}

all_same_sign <- function(x) {
	OR  <- `||`
	AND <- `&&`
	
	OR(length(x) <= 1L,
		 {
		 	if (anyNA(x1 <- x[1L])) {
		 		return(NA)
		 	}
		 	if (x1 == 0) {
		 		AND(min(x) == 0,
		 				max(x) == 0)
		 	} else if (x1 > 0) {
		 		min(x) > 0
		 	} else {
		 		max(x) < 0
		 	}
		 })
}