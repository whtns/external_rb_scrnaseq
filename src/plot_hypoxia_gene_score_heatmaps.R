source("functions.R")
source("src/heatmap_functions.R")

rename_hyp_columns <- function(seu){
    seu$hypoxia <- seu$hypoxia1
    seu$MT <- seu$hypoxia2
    
    return(seu)
}

all_mt_genes <- c("MT-CO1", "MT-CO2", "MT-CYB", "MT-CO3", "MT-ND3", "MT-ATP6", 
                  "MT-ND4", "MT-ND2", "MT-ND1", "MT-ND6", "MT-ND5", "MT-ND4L")

five_mt_genes <- c("MT-CO3", "MT-ND3", "MT-CYB", "MT-ATP6", "MT-CO2")

mt_genes <- c(all_mt_genes[!all_mt_genes %in% five_mt_genes], five_mt_genes)

hypoxia_genes <- c("BNIP3", "GAS5", "EPB41L4A-AS1")

seu_1q <- subset_to_1q(integrated_seu, file_id = "integrated_seu_1q_afterall6.rds" , cluster_orders)

# debug(heatmap_marker_genes_debug)

test12 <- seu_1q |> 
    rename_hyp_columns() |> 
    heatmap_marker_genes_debug(list("hypoxia" = hypoxia_genes, "MT" = mt_genes), "filtered_", marker_col = "clusters", group.by = c("clusters", "scna", "hypoxia_score", "MT", "hypoxia"), col_arrangement = c("clusters"), height = 8, width = 12, title = "No hypoxia threshold")

browseURL(test12)

test13 <- 
    readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_less_0_5.rds") |> 
    rename_hyp_columns() |> 
    heatmap_marker_genes_debug(list("hypoxia" = hypoxia_genes, "MT" = mt_genes), "filtered_", marker_col = "clusters", group.by = c("clusters", "scna", "hypoxia_score", "MT", "hypoxia"), col_arrangement = c("clusters"), height = 8, width = 12, title = "Hypoxia score less than 0.5")

test14 <- 
    readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_more_0_5.rds") |> 
    rename_hyp_columns() |> 
    heatmap_marker_genes_debug(list("hypoxia" = hypoxia_genes, "MT" = mt_genes), "filtered_", marker_col = "clusters", group.by = c("clusters", "scna", "hypoxia_score", "MT", "hypoxia"), col_arrangement = c("clusters"), height = 8, width = 12, title = "Hypoxia score greater than 0.5")

test15 <- 
    readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_less_0_4.rds") |> 
    rename_hyp_columns() |> 
    heatmap_marker_genes_debug(list("hypoxia" = hypoxia_genes, "MT" = mt_genes), "filtered_", marker_col = "clusters", group.by = c("clusters", "scna", "hypoxia_score", "MT", "hypoxia"), col_arrangement = c("clusters"), height = 8, width = 12, title = "Hypoxia score less than 0.4")

test16 <- 
    readRDS("output/seurat/integrated_1q_16q/integrated_seu_1q_afterall6_more_0_4.rds") |> 
    rename_hyp_columns() |> 
    heatmap_marker_genes_debug(list("hypoxia" = hypoxia_genes, "MT" = mt_genes), "filtered_", marker_col = "clusters", group.by = c("clusters", "scna", "hypoxia_score", "MT", "hypoxia"), col_arrangement = c("clusters"), height = 8, width = 12, title = "Hypoxia score greater than 0.4")

qpdf::pdf_combine(list(test12, test13, test14, test15, test16), output = "results/expression_of_hypoxia_genes.pdf") |> browseURL()
