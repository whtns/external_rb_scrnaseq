setwd('/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj')
suppressPackageStartupMessages(source('packages.R'))
lapply(list.files('R', full.names = TRUE, pattern = '[.]R$'), source) |> invisible()
targets::tar_load(cluster_dictionary, store = '_targets_r431')

for (srr in c('SRR27187900', 'SRR27187901')) {
  cat('\n=== Testing sample:', srr, '===\n')
  path <- paste0('output/seurat/', srr, '_unfiltered_seu.rds')
  filtered_seus_vals <- targets::tar_read(filtered_seus, store = '_targets_r431')
  filtered_path <- unlist(filtered_seus_vals)[grepl(srr, unlist(filtered_seus_vals))]
  if (length(filtered_path) == 0) filtered_path <- NULL
  cat('unfiltered_path:', path, '\n')
  cat('filtered_path:', filtered_path, '\n')
  cat('Memory before:', pryr::mem_used(), '\n')
  result <- tryCatch(
    plot_effect_of_filtering(path, filtered_path, cluster_dictionary = cluster_dictionary),
    error = function(e) { cat('ERROR:', conditionMessage(e), '\n'); NULL }
  )
  cat('Result:', class(result), '\n')
  rm(result)
  gc()
  cat('Memory after gc:', pryr::mem_used(), '\n')
}
