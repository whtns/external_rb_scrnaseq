#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(targets)
})

usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript export_targets_rulegraph.R [output_html] [store_path] [targets_only] [arrangement] [view]",
      "",
      "Arguments:",
      "  output_html  Output HTML path (default: doc/targets_rulegraph.html)",
      "  store_path   targets store path (default: _targets_r431)",
      "  targets_only TRUE/FALSE to hide helper nodes (default: TRUE)",
      "  arrangement  one of: hierarchical-lr, hierarchical-ud, force (default: hierarchical-lr)",
      "  view         one of: full, compact (default: full)",
      sep = "\n"
    ),
    "\n"
  )
}

parse_bool <- function(x, default = TRUE) {
  if (is.null(x) || !nzchar(x)) {
    return(default)
  }
  x <- tolower(x)
  if (x %in% c("true", "t", "1", "yes", "y")) {
    return(TRUE)
  }
  if (x %in% c("false", "f", "0", "no", "n")) {
    return(FALSE)
  }
  stop("Invalid boolean value for targets_only: ", x)
}

apply_arrangement <- function(graph, arrangement) {
  arrangement <- tolower(arrangement)

  if (arrangement == "hierarchical-lr") {
    return(
      graph |>
        visNetwork::visHierarchicalLayout(
          enabled = TRUE,
          direction = "LR",
          sortMethod = "directed",
          levelSeparation = 260,
          nodeSpacing = 220,
          treeSpacing = 320
        ) |>
        visNetwork::visPhysics(enabled = FALSE)
    )
  }

  if (arrangement == "hierarchical-ud") {
    return(
      graph |>
        visNetwork::visHierarchicalLayout(
          enabled = TRUE,
          direction = "UD",
          sortMethod = "directed",
          levelSeparation = 210,
          nodeSpacing = 180,
          treeSpacing = 240
        ) |>
        visNetwork::visPhysics(enabled = FALSE)
    )
  }

  if (arrangement == "force") {
    return(
      graph |>
        visNetwork::visPhysics(
          enabled = TRUE,
          solver = "forceAtlas2Based",
          stabilization = list(enabled = TRUE, iterations = 250)
        )
    )
  }

  stop(
    "Invalid arrangement: ", arrangement,
    ". Use one of hierarchical-lr, hierarchical-ud, force"
  )
}

apply_view <- function(graph, view) {
  view <- tolower(view)
  if (view == "full") {
    return(graph)
  }

  if (view == "compact") {
    # In this pipeline DAG, "pattern" nodes are square branching helpers.
    # Keep stem targets and drop pattern helpers for a cleaner high-level graph.
    nodes <- graph$x$nodes
    if (!is.null(nodes$type)) {
      keep <- is.na(nodes$type) | nodes$type != "pattern"
      nodes_kept <- nodes[keep, , drop = FALSE]

      if (nrow(nodes_kept) == 0) {
        stop("Compact view removed all nodes; refusing to write empty graph")
      }

      keep_ids <- nodes_kept$id
      graph$x$nodes <- nodes_kept

      if (!is.null(graph$x$edges) && nrow(graph$x$edges) > 0) {
        graph$x$edges <- graph$x$edges[
          graph$x$edges$from %in% keep_ids & graph$x$edges$to %in% keep_ids,
          ,
          drop = FALSE
        ]
      }
    }

    if (!is.null(graph$x$legend) && !is.null(graph$x$legend$nodes)) {
      status_labels <- c("Outdated", "Up to date", "Errored", "Started")
      legend_nodes <- graph$x$legend$nodes
      graph$x$legend$nodes <- legend_nodes[legend_nodes$label %in% status_labels, , drop = FALSE]
    }

    return(graph)
  }

  stop("Invalid view: ", view, ". Use one of full, compact")
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0 && args[[1]] %in% c("-h", "--help")) {
  usage()
  quit(status = 0)
}

output_html <- if (length(args) >= 1) args[[1]] else "doc/targets_rulegraph.html"
store_path <- if (length(args) >= 2) args[[2]] else "_targets_r431"
targets_only <- parse_bool(if (length(args) >= 3) args[[3]] else NULL, default = TRUE)
arrangement <- if (length(args) >= 4) args[[4]] else "hierarchical-lr"
view <- if (length(args) >= 5) args[[5]] else "full"

if (!dir.exists(store_path)) {
  stop("Store path does not exist: ", store_path)
}

if (!requireNamespace("visNetwork", quietly = TRUE)) {
  stop(
    "Package 'visNetwork' is required to save the graph. ",
    "Install with: install.packages('visNetwork')"
  )
}

out_dir <- dirname(output_html)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

tar_config_set(store = store_path)
graph <- tar_visnetwork(targets_only = targets_only)
graph <- apply_view(graph, view)
graph <- apply_arrangement(graph, arrangement)

font_size <- if (tolower(view) == "compact") 14 else 22
graph <- graph |>
  visNetwork::visNodes(font = list(size = font_size)) |>
  visNetwork::visEdges(arrows = "to", smooth = FALSE)
visNetwork::visSave(graph, file = output_html)

cat("Saved targets rulegraph to:", output_html, "\n")
cat("Using store:", store_path, "\n")
cat("targets_only:", targets_only, "\n")
cat("arrangement:", arrangement, "\n")
cat("view:", view, "\n")