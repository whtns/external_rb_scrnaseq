library(Seurat)
library(seuratTools)

# Directory containing Seurat objects
seurat_dir <- "/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/seurat"

# Get all .rds files in the directory
rds_files <- list.files(seurat_dir, pattern = "^SRR[0-9]+_filtered_seu\\.rds$", full.names = TRUE, ignore.case = TRUE)
# Output directory for diploid objects
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
    
    # Filter out Seurat objects with fewer than 100 cells
    seurat_list <- Filter(function(obj) ncol(obj) >= 100, seurat_list)
}


library(purrr)
# run integration
test <- integration_workflow(seurat_list)


debug(find_all_markers)

# seurat_culster ran; it successfully clustered cells at lower resolutions
# but the higher resolution clustering was incomplete.
# therefore the problem is that the find_all_markers's stash_marker_features is failing
# new_markers <- purrr::map(metavar, stash_marker_features, 
# seu, seurat_assay = seurat_assay, ...)


sapply(clusters, function(x) sum(is.na(x)))
# gene_snn_res.0.2  gene_snn_res.0.4  gene_snn_res.0.6 
# 0                 0                 0 
# gene_snn_res.0.8    gene_snn_res.1  gene_snn_res.1.2 
# 3544              3544              3544 
# gene_snn_res.1.4  gene_snn_res.1.6  gene_snn_res.1.8 
# 3544              3544              3544 
# gene_snn_res.2 gene_snn_res.0.15 
# 3544              6346 

# during debug; before the stash_marker_features step
# I removed all the resolutions with any NAs using
metavar <- metavar[sapply(metavar, function(col) {
    sum(is.na(clusters[[col]])) == 0
})]

# with this you still have all cells, 
# you're just finding markers only at resolutions where every cell has a valid cluster assignment.

saveRDS(test, "/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/seurat/diploid_subsets/diploid_seu.rds")

# $diploid_seu
# An object of class Seurat 
# 84133 features across 6465 samples within 4 assays 
# Active assay: integrated (2000 features, 2000 variable features)
# 2 layers present: data, scale.data
# 3 other assays present: gene, SCT, RNA
# 3 dimensional reductions calculated: pca, tsne, umap

