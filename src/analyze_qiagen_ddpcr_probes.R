#!/usr/bin/env Rscript

library('tidyverse')
library('fs')
library('readxl')

qiagen_probes <- readxl::read_xlsx("data/qiagen_ddpcr_probes.xlsx") %>% 
	janitor::clean_names() %>% 
	dplyr::mutate(symbol = str_remove(name, ".* for ")) %>% 
	dplyr::left_join(annotables::grch38, by = "symbol")
