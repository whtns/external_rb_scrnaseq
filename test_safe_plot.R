module <- "r/4.4.1"
system(paste("module load", module))
Sys.setenv(LD_LIBRARY_PATH = paste(Sys.getenv("HOME"), "/lib", sep = "", Sys.getenv("LD_LIBRARY_PATH")))

library(numbatHelpers)

# Test that safe_plot_numbat returns a list with otherwise
result <- safe_plot_numbat(NULL, NULL, NULL, NULL)
cat("Result class:", class(result), "\n")
cat("Result structure:\n")
str(result)

# Test the extraction
extracted <- result[["result"]]
cat("Extracted result:", extracted, "\n")
