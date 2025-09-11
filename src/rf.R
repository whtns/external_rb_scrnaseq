#!/usr/bin/env Rscript   

library(targets)
suppressPackageStartupMessages(source("packages.R"))
source("functions.R")
library(plotgardener)
library(tidymodels)

seu <- readRDS("/dataVolume/storage/single_cell_projects/resources/external_rb_scrnaseq_proj/output/seurat/integrated_1q_16q/integrated_seu_1q_afterall.rds")

data <- seu@assays$SCT@scale.data %>%
    t() |>
    as.data.frame() |>
    identity()

# mt_genes <- str_subset(rownames(seu), "MT-.*")
mt_genes <- c("MT-CO3", "MT-ND3", "MT-CYB", "MT-ATP6", "MT-CO2")

# data <- data[,colnames(data) %in% unique(seu@misc$markers$clusters$presto$Gene.Name)]
# data <- data[,colnames(data) %in% VariableFeatures(seu)]

data$clusters <- seu$clusters

data$cell <- colnames(seu)

set.seed(123)

data_split <- initial_split(data, strata = "clusters")
data_train <- training(data_split)
data_test <- testing(data_split)

# 2 fold cross validation
data_fold <- vfold_cv(data_train, v = 2)

rf_recipe <- 
    recipe(formula = clusters ~ ., data = data_train) %>%
    update_role(cell, new_role = "ID") %>%
    step_zv(all_predictors())

## feature importance sore to TRUE
rf_spec <- rand_forest() %>%
    set_engine("randomForest", importance = TRUE) %>%
    set_mode("classification")

rf_workflow <- workflow() %>% 
    add_recipe(rf_recipe) %>% 
    add_model(rf_spec)

rf_fit <- fit(rf_workflow, data = data_train)

saveRDS(rf_fit, "tmp.rds")

rf_fit <- readRDS("tmp.rds")

## confusion matrix, perfect classification! 
predict(rf_fit, new_data = data_test) %>%
    bind_cols(data_test %>% select(clusters)) %>%
    conf_mat(truth = clusters, estimate = .pred_class)

# use vip https://koalaverse.github.io/vip/articles/vip.html
# also read https://github.com/tidymodels/parsnip/issues/311
rf_fit %>%
    extract_fit_parsnip() %>%
    vip::vip(geom = "col", num_features = 25) + 
    theme_bw(base_size = 14)+
    labs(title = "Random forest variable importance") 

ggsave("tmp.pdf")