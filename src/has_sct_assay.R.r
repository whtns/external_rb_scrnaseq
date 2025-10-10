has_sct_assay <- function(rds_path) {
  obj <- readRDS(rds_path)
  sample_id <- str_extract(rds_path, "SRR[0-9]*")
  print(sample_id)
  "SCT" %in% names(obj@assays)
}

map(seus, has_sct_assay)
