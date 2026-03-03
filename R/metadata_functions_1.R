# Parameter and metadata functions (11)

#' Add batch hash metadata column
#'
#' @param seu Seurat object with a 'batch' metadata column
#' @return Seurat object with added 'hash' metadata column containing digest of each batch
#' @export
add_batch_hash_metadata <- function(filepath = NULL, seu = NULL, sqlite_path = "batch_hashes.sqlite") {
  if (is.null(seu)) {
    seu <- readRDS(filepath)
  }
  # Split by batch
  seu_list <- Seurat::SplitObject(seu, split.by = "batch")
  
  # Compute hash for each batch
  batch_hashes <- seu_list  |> 
  map(colnames)  |> 
  map(str_remove, "_.*")  |> 
  purrr::map_chr(digest::digest)
  
  # Create a lookup data frame
  hash_lookup <- data.frame(
    batch = names(batch_hashes),
    batch_hash = base::unname(batch_hashes),
    stringsAsFactors = FALSE
  )
  
  # Add hash to metadata by matching batch
  seu$batch_hash <- hash_lookup$batch_hash[match(seu$batch, hash_lookup$batch)]
  # add integration hash
  hash <- digest::digest(colnames(seu))
  seu$hash <- hash

  if (is.null(filepath)){
    filepath <- glue::glue("output/seurat/{hash}_seu.rds")
  }
  
  # Write to sqlite
  # DBI::dbExecute <- get("dbExecute", asNamespace("DBI")) # for safe use in targets
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  # Create table if not exists
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS hashes (filepath TEXT PRIMARY KEY, hash TEXT)")
  # Insert or replace
  DBI::dbExecute(con, "INSERT OR REPLACE INTO hashes (filepath, hash) VALUES (?, ?)", params = list(filepath, hash))
  DBI::dbDisconnect(con)
  saveRDS(seu, filepath)
  return(filepath)

  return(seu)

}

add_hash_metadata <- function(filepath = NULL, seu = NULL, sqlite_path = "batch_hashes.sqlite"){
  if (is.null(seu)) {
    seu <- readRDS(filepath)
  }
  hash <- digest::digest(colnames(seu))
  seu$hash <- hash

  if (is.null(filepath)){
    filepath <- glue::glue("output/seurat/{hash}_seu.rds")
  }

  # Write to sqlite
  # DBI::dbExecute <- get("dbExecute", asNamespace("DBI")) # for safe use in targets
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  # Create table if not exists
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS hashes (filepath TEXT PRIMARY KEY, hash TEXT)")
  # Insert or replace
  DBI::dbExecute(con, "INSERT OR REPLACE INTO hashes (filepath, hash) VALUES (?, ?)", params = list(filepath, hash))
  DBI::dbDisconnect(con)
  saveRDS(seu, filepath)
  return(filepath)

  return(seu)
}

read_seu_hash <- function(sqlite_path = "batch_hashes.sqlite", filepath = NULL) {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  
  if (!is.null(filepath)) {
    # Query for specific filepath
    hashes_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes WHERE filepath = ?", params = list(filepath))
  } else {
    # Query all records
    hashes_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes")
  }
  
  DBI::dbDisconnect(con)
  return(hashes_df)
}

make_filepaths_unique_in_hashes_table <- function(){
  con <- DBI::dbConnect(RSQLite::SQLite(), "batch_hashes.sqlite")

  hashes_df <-
  tbl(con, "hashes")  |> 
  dplyr::distinct(hash, .keep_all = TRUE)

  DBI::dbWriteTable(
      con,
      "hashes", # The name of the table to replace
      hashes_df,
      overwrite = TRUE,
      temporary = FALSE # Set to FALSE for a persistent table
    )
  
  DBI::dbDisconnect(con)
  message("Made filepaths unique in hashes table")
}

read_seu_path <- function(hash = NULL, sqlite_path = "batch_hashes.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)

  if (!is.null(hash)) {
    # Query for specific hash
    filepath_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes WHERE hash = ?", params = list(hash))
  } else {
    # Query all records
    filepath_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes")
  }
  
  DBI::dbDisconnect(con)

  filepath_df <- dplyr::slice_tail(filepath_df, n = 1)

  return(filepath_df$filepath)
}

read_batch_hashes <- function(filepath = NULL, sqlite_path = "batch_hashes.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  
  if (!is.null(filepath)) {
    # Query for specific filepath
    hashes_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes WHERE filepath = ?", params = list(filepath))
    seu <- readRDS(hashes_df$filepath)
    return(unique(seu$batch_hash))

    batches_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes WHERE filepath = ?", params = list(filepath))

  } 

  DBI::dbDisconnect(con)

}

encode_cluster_order_to_hash_table <- function(cluster_orders, sqlite_path = "batch_hashes.sqlite"){
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  # Create table if not exists
  DBI::dbExecute(con, "CREATE TABLE IF NOT EXISTS cluster_orders (file_id TEXT PRIMARY KEY, cluster_order TEXT)")
  
  for (file_id in names(cluster_orders)) {
    # Convert to JSON
    cluster_order_json <- jsonlite::toJSON(cluster_orders[[file_id]], auto_unbox = TRUE)
    # Insert or replace
    DBI::dbExecute(con, "INSERT OR REPLACE INTO cluster_orders (file_id, cluster_order) VALUES (?, ?)", params = list(file_id, cluster_order_json))
  }
  
  DBI::dbDisconnect(con)
}

read_cluster_orders_table <- function(sqlite_path = "batch_hashes.sqlite", file_id = NULL, hash = NULL, as_list = TRUE) {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  
  if (!is.null(file_id)) {
    # Query for specific file_id
    cluster_orders_df <- DBI::dbGetQuery(con, "SELECT * FROM cluster_orders WHERE file_id = ?", params = list(file_id))
  } else if (!is.null(hash)) {
    # Query for specific hash
    filepaths_df <- DBI::dbGetQuery(con, "SELECT * FROM hashes WHERE hash = ?", params = list(hash))
    cluster_orders_df <- DBI::dbGetQuery(con, "SELECT * FROM cluster_orders WHERE file_id = ?", params = list(fs::path_file(filepaths_df$filepath)))
  } else {
    # Query all records
    cluster_orders_df <- DBI::dbGetQuery(con, "SELECT * FROM cluster_orders")
  }
  
  DBI::dbDisconnect(con)
  
  # Parse JSON if requested
  if (as_list && nrow(cluster_orders_df) > 0) {
    cluster_orders_list <- lapply(cluster_orders_df$cluster_order, function(x) {
      tryCatch({
        jsonlite::fromJSON(x)
      }, error = function(e) {
        # If JSON parsing fails, try to parse as R expression (legacy format)
        warning("Non-JSON format detected for entry, attempting to parse as R expression")
        eval(parse(text = x))
      })
    })
    names(cluster_orders_list) <- cluster_orders_df$file_id
    return(cluster_orders_list)
  }
  
  return(cluster_orders_df)
}

#' Drop and recreate cluster_orders table
#' Use this if you need to clear old non-JSON data
drop_cluster_orders_table <- function(sqlite_path = "batch_hashes.sqlite") {
  con <- DBI::dbConnect(RSQLite::SQLite(), sqlite_path)
  DBI::dbExecute(con, "DROP TABLE IF EXISTS cluster_orders")
  DBI::dbDisconnect(con)
  message("Dropped cluster_orders table from ", sqlite_path)
}

#' Perform retrieve snakemake params operation
#'
#' @param numbat_rds_file File path
#' @return List object
#' @export
retrieve_snakemake_params <- function(numbat_rds_file) {
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

#' Perform retrieve current param operation
#'
#' @param current_params Parameter for current params
#' @param myparam Parameter for myparam
#' @return Function result
#' @export
retrieve_current_param <- function(current_params, myparam) {
  sample_ids <- map(current_params, 1)
  param_values <- map(current_params, 2) %>%
    map(myparam) %>%
    set_names(sample_ids)
}