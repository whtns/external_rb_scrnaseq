library(jsonlite)
library(readr)
library(purrr)
library(stringr)
library(dplyr)

parse_value <- function(x) {
  if (is.na(x) || x == "") {
    return(NULL)
  }
  # If all parts are digits and separated by underscores
  if (str_detect(x, "_") && all(str_detect(str_split(x, "_")[[1]], "^\\d+$"))) {
    return(as.integer(str_split(x, "_")[[1]]))
  }
  # Try to convert to integer
  if (str_detect(x, "^\\d+$")) {
    return(as.integer(x))
  }
  return(x)
}

input_csv <- "data/scna_cluster_order.csv"
output_json <- "data/scna_cluster_order.json"

df <- read_csv(input_csv, show_col_types = FALSE)
parsed <- df %>%
  mutate(across(everything(), ~map(.x, parse_value)))

json_data <- toJSON(parsed, pretty = TRUE, na = "null")
write(json_data, file = output_json)

cat("Converted", input_csv, "to", output_json, "\n")