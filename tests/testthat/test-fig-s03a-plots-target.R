test_that("fig_s03a_plots target is configured for unfiltered numbat heatmaps", {
  path_candidates <- c(
    "_targets_BB_edit.R",
    "../_targets_BB_edit.R",
    "../../_targets_BB_edit.R"
  )
  targets_file <- path_candidates[file.exists(path_candidates)][1]
  expect_true(file.exists(targets_file))

  lines <- readLines(targets_file, warn = FALSE)
  text <- paste(lines, collapse = "\n")

  # Verify target name and core command contract.
  expect_match(text, "tar_target\\(fig_s03a_plots,", perl = TRUE)
  expect_match(text, "make_numbat_heatmaps\\(original_seus, numbat_rds_files, p_min = 0\\.9, line_width = 0\\.1, extension = \"_unfiltered\"\\)", perl = TRUE)

  # Verify dynamic branching shape for this target.
  expect_match(text, "pattern = map\\(original_seus\\)", perl = TRUE)
  expect_match(text, "iteration = \"list\"", perl = TRUE)
})
