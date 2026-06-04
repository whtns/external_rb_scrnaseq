# Generate roxygen2 documentation
# This script will process all the roxygen comments and generate proper documentation

# Install required packages if not already installed
if (!require("roxygen2")) {
  install.packages("roxygen2")
}

if (!require("devtools")) {
  install.packages("devtools")
}

library(roxygen2)
library(devtools)

# Generate documentation
cat("Generating roxygen2 documentation...\n")

# Set the working directory to the package root (current directory)
setwd(".")

# Generate documentation using roxygen2
roxygen2::roxygenise(package.dir = ".", roclets = c("rd", "collate", "namespace"))

cat("Documentation generation completed!\n")
cat("Check the 'man/' directory for generated .Rd files\n")

# If you want to build a proper package, uncomment the lines below:
# devtools::document()
# devtools::check()