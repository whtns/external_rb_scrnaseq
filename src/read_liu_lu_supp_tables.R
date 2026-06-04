#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

myreadxl <- function(excel_file, skip = 0){
	sheets <- readxl::excel_sheets(excel_file) %>% set_names(.)
	
	map(sheets, ~readxl::read_excel(excel_file, .x, skip = skip))
}

liu_supp_tables <- "data/liu_lu_supp_data/supp_table_4.xlsx" %>% 
	myreadxl(skip = 1) %>% 
	map(janitor::clean_names) %>% 
	identity()