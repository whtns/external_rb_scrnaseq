# This function identifies function definitions present in one R script but not another.
# It sources each file into a temporary, isolated environment to avoid conflicts
# with your current R session.

compare_sourced_functions <- function(file1, file2) {
  # Step 1: Input validation
  # Check if both files exist.
  if (!file.exists(file1)) {
    stop(paste("Error: file not found at path:", file1))
  }
  if (!file.exists(file2)) {
    stop(paste("Error: file not found at path:", file2))
  }

  # Step 2: Source each file into a clean, temporary environment.
  # This prevents functions from the sourced files from polluting the global environment.
  env1 <- new.env()
  env2 <- new.env()

  # Safely source the files.
  tryCatch(
    lapply(list.files("./R", full.names = TRUE), source, local = env1),
    error = function(e) {
      stop(paste("Error sourcing", file1, ":", e$message))
    }
  )

  tryCatch(
    source("functions.R", local = env2),
    error = function(e) {
      stop(paste("Error sourcing", file2, ":", e$message))
    }
  )

  # Step 3: Identify all function objects in each environment.
  # We use `sapply` to check if each object is a function.
  functions1 <- names(which(sapply(ls(env1), function(x) is.function(get(x, envir = env1)))))
  functions2 <- names(which(sapply(ls(env2), function(x) is.function(get(x, envir = env2)))))

  # Step 4: Compare the two lists of function names.
  # `setdiff` finds elements in the first vector that are not in the second.
  unique_to_file1 <- setdiff(functions1, functions2)
  unique_to_file2 <- setdiff(functions2, functions1)

  # Step 5: Return a clear, named list with the results.
  return(list(
    unique_to_file1 = unique_to_file1,
    unique_to_file2 = unique_to_file2
  ))
}

# --- Example Usage ---

# Step 1: Create two temporary R files for demonstration.
# File 1 will have functions A, B, and C.
# File 2 will have functions A, B, and D.
cat("
# File 1 content
func_A <- function() {
  print('This is function A.')
}
func_B <- function() {
  print('This is function B.')
}
func_C <- function() {
  print('This is function C, unique to File 1.')
}
", file = "temp_file1.R")

cat("
# File 2 content
func_A <- function() {
  print('This is function A.')
}
func_B <- function() {
  print('This is function B.')
}
func_D <- function() {
  print('This is function D, unique to File 2.')
}
", file = "temp_file2.R")

# Step 2: Call the function to compare the two temporary files.
comparison_results <- compare_sourced_functions("temp_file1.R", "temp_file2.R")

# Step 3: Print the results to the console.
print("Comparison Results:")
print("Functions unique to temp_file1.R:")
print(comparison_results$unique_to_file1)
print("---------------------------------")
print("Functions unique to temp_file2.R:")
print(comparison_results$unique_to_file2)

# Step 4: Clean up the temporary files.
unlink(c("temp_file1.R", "temp_file2.R"))
