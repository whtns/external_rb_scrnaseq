module load r/4.4.1
export LD_LIBRARY_PATH="$HOME/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"

Rscript -e "
library(targets)
library(stringr)
library(numbatHelpers)

# Load targets data
tar_config_set(store = '_targets_r431')
cat('Loading unfiltered_seus...\n')
tar_load(unfiltered_seus, store = '_targets_r431')
tar_load(numbat_rds_files, store = '_targets_r431')

# Get first sample
seu_path <- unfiltered_seus[[1]]
cat('Testing with seu_path:', seu_path, '\n')

# Extract sample ID
sample_id <- str_extract(seu_path, 'SRX[0-9]+')
cat('Sample ID:', sample_id, '\n')

# Run make_numbat_heatmaps
tryCatch({
  cat('\nCalling make_numbat_heatmaps...\n')
  result <- make_numbat_heatmaps(
    seu_path,
    numbat_rds_files,
    p_min = 0.9,
    line_width = 0.1,
    extension = '_unfiltered_test',
    show_segment_names_on_x = TRUE
  )
  cat('SUCCESS! Result:\n')
  print(result)
  
  # Check if files were created
  for (file in result) {
    if (!is.na(file)) {
      cat('Checking file:', file, '\n')
      if (file.exists(file)) {
        cat('  ✓ File exists (', file.size(file), 'bytes)\n')
      } else {
        cat('  ✗ File NOT found\n')
      }
    }
  }
}, error = function(e) {
  cat('ERROR:', conditionMessage(e), '\n')
  traceback()
})
"
