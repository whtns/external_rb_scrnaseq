# R Script for Seurat Layer to BPCells Conversion
# This function requires Seurat (v5) and BPCells packages.

# Function: convert_seurat_layers_to_bpcells
# Description:
#   Reads a Seurat object, converts specified (or all) in-memory layers of the 
#   active assay to on-disk BPCells matrices, and saves the modified object. 
#   This is ideal for reducing memory usage when working with large datasets.
#
# Arguments:
#   input_path (character): Path to the existing Seurat object (.rds file).
#   output_path (character): Path where the new Seurat object with BPCells 
#                            layers will be saved (.rds file).
#   bpcells_dir (character): Directory where the on-disk BPCells matrix files 
#                            will be stored. This directory must be kept with 
#                            the saved Seurat object.
#   assay_name (character/NULL): The name of the assay to process (e.g., "RNA").
#                                If NULL (default), the function processes the
#                                active assay of the Seurat object.
#
# Returns:
#   Invisibly returns the modified Seurat object.

library(tictoc)

convert_seurat_layers_to_bpcells <- function(
  input_path, 
  output_path, 
  bpcells_dir, 
  assay_name = NULL
) {
  
  # --- 1. Load and Verify Packages ---
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required but not installed.")
  }
  if (!requireNamespace("BPCells", quietly = TRUE)) {
    stop("Package 'BPCells' is required but not installed.")
  }
  
  message("--- Starting BPCells Conversion Process ---")
  message(paste("1. Reading Seurat object from:", input_path))
  
  # --- 2. Read Object and Setup Directories ---
  tryCatch({
    seurat_obj <- readRDS(input_path)
  }, error = function(e) {
    stop(paste("Failed to read Seurat object:", e$message))
  })
  
  # Set up the assay to process
  if (is.null(assay_name)) {
    assay_name <- Seurat::DefaultAssay(seurat_obj)
    message(paste("   Using Default Assay:", assay_name))
  } else if (!assay_name %in% Seurat::Assays(seurat_obj)) {
    stop(paste("Assay '", assay_name, "' not found in the Seurat object."))
  }
  
  # Create the BPCells storage directory
  if (!dir.exists(bpcells_dir)) {
    message(paste("   Creating BPCells storage directory:", bpcells_dir))
    dir.create(bpcells_dir, recursive = TRUE)
  }
  
  # --- 3. Identify and Convert Layers ---
  
  # Get all layers in the target assay
  all_layers <- SeuratObject::Layers(seurat_obj[[assay_name]])
  if (length(all_layers) == 0) {
    message(paste("   Assay '", assay_name, "' has no layers to process. Exiting."))
    return(invisible(seurat_obj))
  }
  
  message(paste("2. Converting layers for assay '", assay_name, "': ", paste(all_layers, collapse = ", ")))
  
  for (layer in all_layers) {
    # Extract the layer data
    layer_data <- SeuratObject::LayerData(seurat_obj, assay = assay_name, layer = layer)
    
    # Check if the layer is already an on-disk BPCells matrix
    if (inherits(layer_data, "BPCellsMatrix")) {
      message(paste("   [SKIP]", layer, "is already an on-disk BPCells matrix."))
      next
    }
    
    # Convert the in-memory matrix to BPCells format on disk
    matrix_storage_path <- file.path(bpcells_dir, paste0(assay_name, "_", layer))
    message(paste("   [WRITE]", layer, "-> Saving to:", matrix_storage_path))
    
    tryCatch({
      # Use BPCells function to write the matrix to disk
      bpcells_matrix <- BPCells::write_matrix_dir(
        mat = layer_data, 
        dir = matrix_storage_path
      )
      
      # Replace the in-memory data in the Seurat object with the BPCells matrix pointer
      SeuratObject::LayerData(seurat_obj, assay = assay_name, layer = layer) <- bpcells_matrix
      message(paste("   [SUCCESS]", layer, "conversion complete."))
      
    }, error = function(e) {
      warning(paste("   [FAILED] Conversion for layer", layer, "failed:", e$message))
    })
  }
  
  # --- 4. Save the Modified Object ---
  message(paste("3. Saving modified Seurat object to:", output_path))
  tryCatch({
    # Use saveRDS, which stores the Seurat object metadata and the BPCells file paths
    saveRDS(object = seurat_obj, file = output_path)
    message("--- Conversion and Saving Complete ---")
  }, error = function(e) {
    stop(paste("Failed to save Seurat object:", e$message))
  })
  
  return(invisible(seurat_obj))
}

# Example Usage (Demonstration - requires a sample Seurat object and the packages)
# NOTE: This example is commented out and requires the user to replace paths and have
# a large Seurat object available for a real test.

# # -------------------------------------------------------------------------
# # Example Setup (Run these lines in R console before running the function)
# # -------------------------------------------------------------------------
# # 1. Install necessary packages if not present:
# # install.packages("devtools")
# # devtools::install_github("satijalab/seurat", ref = "seurat5")
# # install.packages("BPCells")
# 
2. Create a dummy Seurat object for demonstration (optional)
library(Seurat)
set.seed(42)
data <- matrix(rnbinom(50000, mu = 5, size = 10), ncol = 10000)
rownames(data) <- paste0("Gene", 1:nrow(data))
colnames(data) <- paste0("Cell", 1:ncol(data))
obj <- Seurat::CreateSeuratObject(counts = data, assay = "RNA")
file.remove("dummy_input.rds")
saveRDS(obj, file = "dummy_input.rds")

tic()
obj <- readRDS("dummy_input.rds")
toc()

# Define paths for the demonstration
INPUT_FILE <- "dummy_input.rds"
OUTPUT_FILE <- "modified_bpcells_output.rds"
BPCELLS_DIR <- "bpcells_storage_files"

fs::dir_delete(BPCELLS_DIR)
file.remove(OUTPUT_FILE)

# -------------------------------------------------------------------------
# Example Function Call:
# -------------------------------------------------------------------------
converted_obj <- convert_seurat_layers_to_bpcells(
  input_path = INPUT_FILE,
  output_path = OUTPUT_FILE,
  bpcells_dir = BPCELLS_DIR,
  assay_name = "RNA" # Optional: specify the assay
)

# Verify the conversion (check the layer class)
Layers(converted_obj$RNA)
class(SeuratObject::LayerData(converted_obj, assay = "RNA", layer = "counts"))

tic()
obj <- readRDS("dummy_input.rds")
toc()
