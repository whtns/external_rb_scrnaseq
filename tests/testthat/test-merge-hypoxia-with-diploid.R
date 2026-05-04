make_minimal_seu <- function(n_genes = 20, n_cells = 10, cell_prefix = "cell") {
  counts <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes, ncol = n_cells,
    dimnames = list(
      paste0("gene", seq_len(n_genes)),
      paste0(cell_prefix, seq_len(n_cells))
    )
  )
  SeuratObject::CreateSeuratObject(counts)
}

test_that("merge_hypoxia_with_diploid returns path matching expected pattern", {
  seu_hyp <- make_minimal_seu(cell_prefix = "hyp")
  seu_dip <- make_minimal_seu(cell_prefix = "dip")

  withr::with_tempdir({
    dir.create("output/seurat", recursive = TRUE)
    hyp_path <- "SRR12345678_low_hypoxia.rds"
    dip_path <- "diploid_seu.rds"
    saveRDS(seu_hyp, hyp_path)
    saveRDS(seu_dip, dip_path)

    result <- merge_hypoxia_with_diploid(hyp_path, dip_path, slug = "low_hypoxia")

    expect_match(result, "SRR12345678_low_hypoxia_diploid_merged\\.rds$")
    expect_match(result, "^output/seurat/")
  })
})

test_that("merge_hypoxia_with_diploid writes a readable merged Seurat object", {
  seu_hyp <- make_minimal_seu(cell_prefix = "hyp")
  seu_dip <- make_minimal_seu(cell_prefix = "dip")

  withr::with_tempdir({
    dir.create("output/seurat", recursive = TRUE)
    saveRDS(seu_hyp, "SRR12345678_low_hypoxia.rds")
    saveRDS(seu_dip, "diploid_seu.rds")

    result <- merge_hypoxia_with_diploid(
      "SRR12345678_low_hypoxia.rds", "diploid_seu.rds", slug = "low_hypoxia"
    )

    expect_true(file.exists(result))
    merged <- readRDS(result)
    expect_s4_class(merged, "Seurat")
    expect_equal(ncol(merged), ncol(seu_hyp) + ncol(seu_dip))
    expect_equal(nrow(merged), nrow(seu_hyp))
  })
})

test_that("merge_hypoxia_with_diploid incorporates slug in output filename", {
  seu_hyp <- make_minimal_seu(cell_prefix = "hyp")
  seu_dip <- make_minimal_seu(cell_prefix = "dip")

  withr::with_tempdir({
    dir.create("output/seurat", recursive = TRUE)
    saveRDS(seu_hyp, "SRR12345678_high_hypoxia.rds")
    saveRDS(seu_dip, "diploid_seu.rds")

    result <- merge_hypoxia_with_diploid(
      "SRR12345678_high_hypoxia.rds", "diploid_seu.rds", slug = "high_hypoxia"
    )

    expect_match(result, "high_hypoxia_diploid_merged")
  })
})

test_that("merge_hypoxia_with_diploid joins Assay5 layers before merging", {
  skip_if_not(
    utils::packageVersion("SeuratObject") >= "5.0.0",
    "Assay5 requires SeuratObject >= 5.0.0"
  )

  counts <- matrix(
    rpois(200, lambda = 5), nrow = 20, ncol = 10,
    dimnames = list(paste0("gene", 1:20), paste0("hyp", 1:10))
  )
  seu_hyp <- SeuratObject::CreateSeuratObject(counts)
  # Split into layers to simulate an un-joined Assay5 object
  seu_hyp[["RNA"]] <- SeuratObject::split(seu_hyp[["RNA"]], f = rep(c("a", "b"), each = 5))

  seu_dip <- make_minimal_seu(cell_prefix = "dip")

  withr::with_tempdir({
    dir.create("output/seurat", recursive = TRUE)
    saveRDS(seu_hyp, "SRR99999999_high_hypoxia.rds")
    saveRDS(seu_dip, "diploid_seu.rds")

    result <- merge_hypoxia_with_diploid(
      "SRR99999999_high_hypoxia.rds", "diploid_seu.rds", slug = "high_hypoxia"
    )

    merged <- readRDS(result)
    expect_s4_class(merged, "Seurat")
    expect_equal(ncol(merged), 20)
  })
})
