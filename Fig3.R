library(dplyr)
library(data.table)
library(Seurat)
library(ggplot2)
library(viridis)
library(Nebulosa)
options(Seurat.object.assay.version = "v3")

#set the working directory to the desired folder path
setwd("~/cross_species_analysis")

#load integrated plaque object ------------------------------------------------
plaques.integrated <- readRDS("plaques.integrated.rds")

# Figure 3a & b----------------------------------------------------------------
md <- plaques.integrated@meta.data %>% as.data.table
cluster_abundance_sample <- md[, .N, by = c("orig.ident", "lineages")] %>% 
  dcast(., orig.ident ~ lineages, value.var = "N")

#Figure 3c & d-----------------------------------------------------------------
Idents(plaques.integrated) <- "orig.ident"
levels(plaques.integrated)
plaques.integrated <- RenameIdents(object = plaques.integrated, 
                                   'one'='human', 
                                   'two'='human', 
                                   'three'='human', 
                                   'four'='human', 
                                   'five'='human', 
                                   'six'='human',
                                   'seven'='human',
                                   'eight'='human',
                                   'apoe_hfd'='mouse', 
                                   "ldlr_chow"='mouse',
                                   "ldlr_hfd"='mouse',
                                   'apoe_chow'='mouse',
                                   'ad_wt'='mouse',
                                   'ad_apoe'='mouse',
                                   'base'='mouse')
plaques.integrated$species <- Idents(plaques.integrated)
levels(plaques.integrated)

Idents(plaques.integrated) <- "species"
DefaultAssay(plaques.integrated) <- "RNA"
plaques.integrated <- NormalizeData(plaques.integrated, verbose = FALSE, assay= "RNA")
Idents(plaques.integrated) <- "lineages"
levels(plaques.integrated)
folder_names <- c("CD4 T cell", "CD8 T cell", "NK cell", "NK T cell", "B cell",
                  "Neutrophil", "Monocyte", "Myeloid DC", "Macrophage", "Plasmocytoid DC",
                  "Other T cell", "Proliferating", "ILC", "Mast cell")

# Specify the parent directory
parent_directory <- "~/cross_species_analysis/"

# Create subfolders based on the vector
for (folder_name in folder_names) {
  subfolder_path <- file.path(parent_directory, folder_name)
  dir.create(subfolder_path, recursive = TRUE, showWarnings = T)
  if (file.exists(subfolder_path)) {
    cat("Subfolder", folder_name, "created successfully.\n")
  } else {
    warning("Failed to create subfolder", folder_name)
  }
}

# Loop over cell types
for (folder_name in folder_names) {
  
  subfolder_path <- file.path(parent_directory, folder_name)
  
  if (!dir.exists(subfolder_path)) {
    dir.create(subfolder_path, recursive = TRUE, showWarnings = TRUE)
  }
  
  cell_type_subset <- subset(plaques.integrated, idents = folder_name)
  
  de_results <- FindMarkers(cell_type_subset, 
                            ident.1 = "human",
                            ident.2 = "mouse", 
                            group.by = "species")
  
  de_results$gene <- rownames(de_results)
  de_results <- subset(de_results, !grepl("^HLA-", gene))
  de_results <- subset(de_results, !grepl("^H2", gene))
  de_results <- subset(de_results, !grepl("^RPL", gene))
  de_results <- subset(de_results, !grepl("^RPS", gene))
  de_results <- subset(de_results, !grepl("^MT-", gene))
  de_results <- subset(de_results, !grepl("^ATP", gene))
  de_results <- subset(de_results, !(de_results$pct.1=="0"))
  de_results <- subset(de_results, !(de_results$pct.2=="0"))
  
  result_file <- file.path(subfolder_path, paste0(folder_name, "_DEGs.csv"))
  write.csv(de_results, file = result_file, row.names = T)
  
  Up_in_human <- subset(de_results, p_val_adj <= 0.05 & avg_log2FC > 0)
  Up_in_human$Up_in_human <- rownames(Up_in_human)
  result_file <- file.path(subfolder_path, "Up_in_human.csv")
  write.csv(Up_in_human, file = result_file, row.names = T)
  
  Up_in_mouse <- subset(de_results, p_val_adj <= 0.05 & avg_log2FC < 0)
  Up_in_mouse$Up_in_mouse <- rownames(Up_in_mouse)
  result_file <- file.path(subfolder_path, "Up_in_mouse.csv")
  write.csv(Up_in_mouse, file = result_file, row.names = T)
  
  cat("Analysis for", folder_name, "completed.\n")
}

#Figure 3e & f-----------------------------------------------------------------
Idents(plaques.integrated) <- "orig.ident"
levels(plaques.integrated)
plaques.integrated <- RenameIdents(plaques.integrated, 
                                   `one` = "Human", 
                                   `two` = "Human",
                                   `three` = "Human",
                                   `four` = "Human",
                                   `five` = "Human",
                                   `six` = "Human",
                                   `seven` = "Human",
                                   `eight` = "Human",
                                   `apoe_hfd` = "Apoe WD",
                                   `ldlr_chow` = "rest",
                                   `ldlr_hfd` = "LDLr HFD",
                                   `apoe_chow` = "Apoe CHOW",
                                   `ad_wt` = "rest",
                                   `ad_apoe` = "rest",
                                   `base` = "rest"
)
levels(plaques.integrated)
plaques.integrated[["models"]] <- Idents(object = plaques.integrated)
levels(plaques.integrated)
plaques.subset <- subset(plaques.integrated, idents = "rest", invert = T)
levels(plaques.subset)
plaques.subset$models <- factor(plaques.subset$models, 
                                         levels = c("Apoe WD",
                                                    "Apoe CHOW",
                                                    "LDLr HFD",
                                                    "Human"))

DefaultAssay(plaques.subset) <- "RNA"
plaques.subset <- NormalizeData(plaques.subset, verbose = T, assay="RNA")

species_markers <- FindAllMarkers(plaques.subset, only.pos = TRUE, min.pct = 0)
species_markers <- species_markers[- grep("HLA-", species_markers$gene),]
species_markers <- species_markers[- grep("H2", species_markers$gene),]
species_markers <- species_markers[- grep("RPL", species_markers$gene),]
species_markers <- species_markers[- grep("RPS", species_markers$gene),]
species_markers <- species_markers[- grep("MT-", species_markers$gene),]
species_markers <- species_markers[- grep("ATP", species_markers$gene),]
species_markers <- species_markers[!(species_markers$pct.1=="0"),]
species_markers <- species_markers[!(species_markers$pct.2=="0"),]

#Figure 3g-i-------------------------------------------------------------------

#All CD4
Idents(plaques.subset) <- "lineages"
levels(plaques.subset)
CD4 <- subset(plaques.subset, idents = "CD4 T cell") 
Idents(CD4) <- "models"
levels(CD4)
DefaultAssay(CD4) <- "RNA"
CD4 <- NormalizeData(CD4)
CD4_markers <- FindAllMarkers(CD4, only.pos = T, logfc.threshold = 0)
CD4_markers <- as_tibble(CD4_markers)
CD4_markers <- na.omit(CD4_markers)
CD4_markers <- subset(CD4_markers,
                             !grepl("^MT-", gene))
CD4_markers <- subset(CD4_markers,
                             !grepl("^RPS", gene))
CD4_markers <- subset(CD4_markers,
                             !grepl("^RPL", gene))
CD4_markers <- subset(CD4_markers,
                             !grepl("^HLA", gene))
CD4_markers <- subset(CD4_markers,
                             !grepl("^TCR", gene))
CD4_markers <- subset(CD4_markers,
                             !grepl("^FC", gene))
CD4_markers <- subset(CD4_markers,
                             !grepl("^GM", gene))
CD4_markers <- subset(CD4_markers,
                             !grepl("^HLA-", gene))
CD4_markers <- subset(CD4_markers,
                             !grepl("^H2", gene))
CD4_markers <- subset(CD4_markers,
                             !grepl("^ATP", gene))
head(CD4_markers)

#export Up in Aopoe WD
CD4_Up_in_ApoeWD <- CD4_markers[CD4_markers$cluster == "Apoe WD", ]
CD4_Up_in_ApoeWD <- subset(CD4_Up_in_ApoeWD, p_val_adj <= 0.05)

#export Up in Aopoe CHOW
CD4_Up_in_ApoeCHOW <- CD4_markers[CD4_markers$cluster == "Apoe CHOW", ]
CD4_Up_in_ApoeCHOW <- subset(CD4_Up_in_ApoeCHOW, p_val_adj <= 0.05)

#export Up in LDLr HFD
CD4_Up_in_LDLrHFD <- CD4_markers[CD4_markers$cluster == "LDLr HFD", ]
CD4_Up_in_LDLrHFD <- subset(CD4_Up_in_LDLrHFD, p_val_adj <= 0.05)

#export Up in Human
CD4_Up_in_Human <- CD4_markers[CD4_markers$cluster == "Human", ]
CD4_Up_in_Human <- subset(CD4_Up_in_Human, p_val_adj <= 0.05)

#All Macrophage
Idents(plaques.subset) <- "lineages"
levels(plaques.subset)
Macrophage <- subset(plaques.subset, idents = "Macrophage") 
Idents(Macrophage) <- "models"
levels(Macrophage)
DefaultAssay(Macrophage) <- "RNA"
Macrophage <- NormalizeData(Macrophage)
Macrophage_markers <- FindAllMarkers(Macrophage, only.pos = T, logfc.threshold = 0)
Macrophage_markers <- as_tibble(Macrophage_markers)
Macrophage_markers <- na.omit(Macrophage_markers)
Macrophage_markers <- subset(Macrophage_markers,
                             !grepl("^MT-", gene))
Macrophage_markers <- subset(Macrophage_markers,
                             !grepl("^RPS", gene))
Macrophage_markers <- subset(Macrophage_markers,
                             !grepl("^RPL", gene))
Macrophage_markers <- subset(Macrophage_markers,
                             !grepl("^HLA", gene))
Macrophage_markers <- subset(Macrophage_markers,
                             !grepl("^TCR", gene))
Macrophage_markers <- subset(Macrophage_markers,
                             !grepl("^FC", gene))
Macrophage_markers <- subset(Macrophage_markers,
                             !grepl("^GM", gene))
Macrophage_markers <- subset(Macrophage_markers,
                             !grepl("^HLA-", gene))
Macrophage_markers <- subset(Macrophage_markers,
                             !grepl("^H2", gene))
Macrophage_markers <- subset(Macrophage_markers,
                             !grepl("^ATP", gene))
head(Macrophage_markers)

#export Up in Aopoe WD
Macrophage_Up_in_ApoeWD <- Macrophage_markers[Macrophage_markers$cluster == "Apoe WD", ]
Macrophage_Up_in_ApoeWD <- subset(Macrophage_Up_in_ApoeWD, p_val_adj <= 0.05)

#export Up in Aopoe CHOW
Macrophage_Up_in_ApoeCHOW <- Macrophage_markers[Macrophage_markers$cluster == "Apoe CHOW", ]
Macrophage_Up_in_ApoeCHOW <- subset(Macrophage_Up_in_ApoeCHOW, p_val_adj <= 0.05)

#export Up in LDLr HFD
Macrophage_Up_in_LDLrHFD <- Macrophage_markers[Macrophage_markers$cluster == "LDLr HFD", ]
Macrophage_Up_in_LDLrHFD <- subset(Macrophage_Up_in_LDLrHFD, p_val_adj <= 0.05)

#export Up in Human
Macrophage_Up_in_Human <- Macrophage_markers[Macrophage_markers$cluster == "Human", ]
Macrophage_Up_in_Human <- subset(Macrophage_Up_in_Human, p_val_adj <= 0.05)

#All B_cell
Idents(plaques.subset) <- "lineages"
levels(plaques.subset)
B_cell <- subset(plaques.subset, idents = "B cell") 
Idents(B_cell) <- "models"
levels(B_cell)
DefaultAssay(B_cell) <- "RNA"
B_cell <- NormalizeData(B_cell)
B_cell_markers <- FindAllMarkers(B_cell, only.pos = T, logfc.threshold = 0)
B_cell_markers <- as_tibble(B_cell_markers)
B_cell_markers <- na.omit(B_cell_markers)
B_cell_markers <- subset(B_cell_markers,
                         !grepl("^MT-", gene))
B_cell_markers <- subset(B_cell_markers,
                         !grepl("^RPS", gene))
B_cell_markers <- subset(B_cell_markers,
                         !grepl("^RPL", gene))
B_cell_markers <- subset(B_cell_markers,
                         !grepl("^HLA", gene))
B_cell_markers <- subset(B_cell_markers,
                         !grepl("^TCR", gene))
B_cell_markers <- subset(B_cell_markers,
                         !grepl("^FC", gene))
B_cell_markers <- subset(B_cell_markers,
                         !grepl("^GM", gene))
B_cell_markers <- subset(B_cell_markers,
                         !grepl("^HLA-", gene))
B_cell_markers <- subset(B_cell_markers,
                         !grepl("^H2", gene))
B_cell_markers <- subset(B_cell_markers,
                         !grepl("^ATP", gene))
head(B_cell_markers)

#export Up in Aopoe WD
B_cell_Up_in_ApoeWD <- B_cell_markers[B_cell_markers$cluster == "Apoe WD", ]
B_cell_Up_in_ApoeWD <- subset(B_cell_Up_in_ApoeWD, p_val_adj <= 0.05)

#export Up in Aopoe CHOW
B_cell_Up_in_ApoeCHOW <- B_cell_markers[B_cell_markers$cluster == "Apoe CHOW", ]
B_cell_Up_in_ApoeCHOW <- subset(B_cell_Up_in_ApoeCHOW, p_val_adj <= 0.05)

#export Up in LDLr HFD
B_cell_Up_in_LDLrHFD <- B_cell_markers[B_cell_markers$cluster == "LDLr HFD", ]
B_cell_Up_in_LDLrHFD <- subset(B_cell_Up_in_LDLrHFD, p_val_adj <= 0.05)

#export Up in Human
B_cell_Up_in_Human <- B_cell_markers[B_cell_markers$cluster == "Human", ]
B_cell_Up_in_Human <- subset(B_cell_Up_in_Human, p_val_adj <= 0.05)

#All CD8
Idents(plaques.subset) <- "lineages"
levels(plaques.subset)
CD8 <- subset(plaques.subset, idents = "CD8 T cell") 
Idents(CD8) <- "models"
levels(CD8)
DefaultAssay(CD8) <- "RNA"
CD8 <- NormalizeData(CD8)
CD8_markers <- FindAllMarkers(CD8, only.pos = T, logfc.threshold = 2.5)
CD8_markers <- as_tibble(CD8_markers)
CD8_markers <- na.omit(CD8_markers)
CD8_markers <- subset(CD8_markers,
                      !grepl("^MT-", gene))
CD8_markers <- subset(CD8_markers,
                      !grepl("^RPS", gene))
CD8_markers <- subset(CD8_markers,
                      !grepl("^RPL", gene))
CD8_markers <- subset(CD8_markers,
                      !grepl("^HLA", gene))
CD8_markers <- subset(CD8_markers,
                      !grepl("^TCR", gene))
CD8_markers <- subset(CD8_markers,
                      !grepl("^FC", gene))
CD8_markers <- subset(CD8_markers,
                      !grepl("^GM", gene))
CD8_markers <- subset(CD8_markers,
                      !grepl("^HLA-", gene))
CD8_markers <- subset(CD8_markers,
                      !grepl("^H2", gene))
CD8_markers <- subset(CD8_markers,
                      !grepl("^ATP", gene))
head(CD8_markers)

#export Up in Aopoe WD
CD8_Up_in_ApoeWD <- CD8_markers[CD8_markers$cluster == "Apoe WD", ]
CD8_Up_in_ApoeWD <- subset(CD8_Up_in_ApoeWD, p_val_adj <= 0.05)

#export Up in Aopoe CHOW
CD8_Up_in_ApoeCHOW <- CD8_markers[CD8_markers$cluster == "Apoe CHOW", ]
CD8_Up_in_ApoeCHOW <- subset(CD8_Up_in_ApoeCHOW, p_val_adj <= 0.05)

#export Up in LDLr HFD
CD8_Up_in_LDLrHFD <- CD8_markers[CD8_markers$cluster == "LDLr HFD", ]
CD8_Up_in_LDLrHFD <- subset(CD8_Up_in_LDLrHFD, p_val_adj <= 0.05)

#export Up in Human
CD8_Up_in_Human <- CD8_markers[CD8_markers$cluster == "Human", ]
CD8_Up_in_Human <- subset(CD8_Up_in_Human, p_val_adj <= 0.05)


