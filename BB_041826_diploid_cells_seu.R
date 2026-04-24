library(Seurat)

# Directory containing Seurat objects
seurat_dir <- "/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/seurat"

# Get all .rds files in the directory
rds_files <- list.files(seurat_dir, pattern = "^SRR[0-9]+_filtered_seu\\.rds$", full.names = TRUE, ignore.case = TRUE)
# Output directory for diploid-subsetted objects
output_dir <- "/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/seurat/diploid_subsets"
diploid_paths <- c()

for (file_path in rds_files) {
    message("Checking: ", basename(file_path))
    
    tryCatch({
        seurat_obj <- readRDS(file_path)
        
        # Check if SCNA column exists in metadata
        if ("scna" %in% colnames(seurat_obj@meta.data)) {
            
            # Check if any cell has SCNA == "diploid"
            if (any(seurat_obj@meta.data$scna == "", na.rm = TRUE)) {
                diploid_paths <- c(diploid_paths, file_path)
                message("  ✔ Contains diploid cells: ", file_path)
                
                # Subset to diploid cells only
                diploid_obj <- subset(seurat_obj, subset = scna == "")
                message("  → Diploid cells: ", ncol(diploid_obj), " / ", ncol(seurat_obj))
                
                # Save subsetted object to output directory
                out_filename <- sub("\\.rds$", "_diploid.rds", basename(file_path))
                out_path <- file.path(output_dir, out_filename)
                saveRDS(diploid_obj, file = out_path)
                message("  → Saved to: ", out_path)
                
                rm(diploid_obj)
            } else {
                message("  ✘ No diploid cells found")
            }
            
        } else {
            message("  ! 'SCNA' column not found in metadata")
        }
        
        # Remove object from memory after each iteration
        rm(seurat_obj)
        gc()
        
    }, error = function(e) {
        message("  ERROR reading file: ", file_path, "\n  ", e$message)
    })
}

message("\nDone! Found ", length(diploid_paths), " Seurat object(s) with diploid cells.")
print(diploid_paths)


### Integrate these seurat objects


path <- "/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/seurat/diploid_subsets"

if (!file.exists(path)) {
    stop("Path does not exist: ", path)
} else {
    files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
    seurat_list <- setNames(lapply(files, readRDS), 
                            tools::file_path_sans_ext(basename(files)))
}
integration_workflow()

