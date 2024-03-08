library(dplyr)
library(data.table)
library(Seurat)
library(ggplot2)
library(viridis)
library(Nebulosa)
options(Seurat.object.assay.version = "v3")

#set the working directory to the desired folder path
setwd("~/cross_species_analysis")

# Run CCA integration----------------------------------------------------------

# load data files
one.data <- read.table(file =  "_T_one.csv", header = T, row.names = 1, sep=",", as.is = T)
two.data <- read.table(file =  "_T_two.csv", header = T, row.names = 1, sep=",", as.is = T)
three.data <- read.table(file =  "_T_three.csv", header = T, row.names = 1, sep=",", as.is = T)
four.data <- read.table(file =  "_T_four.csv", header = T, row.names = 1, sep=",", as.is = T)
five.data <- read.table(file =  "_T_five.csv", header = T, row.names = 1, sep=",", as.is = T)
six.data <- read.table(file =  "_T_six.csv", header = T, row.names = 1, sep=",", as.is = T)
seven.data <- read.table(file =  "_T_seven.csv", header = T, row.names = 1, sep=",", as.is = T)
eight.data <- read.table(file =  "_T_eight.csv", header = T, row.names = 1, sep=",", as.is = T)
hfd.data <- read.table(file =  "_T_apoe_hfd.csv", header = T, row.names = 1, sep=",", as.is = T)
chow.data <- read.table(file =  "_T_apoe_chow.csv", header = T, row.names = 1, sep=",", as.is = T)
chow_LDLR.data <- read.table(file =  "_T_ldlr_chow.csv", header = T, row.names = 1, sep=",", as.is = T)
HFD_LDLR.data <- read.table(file =  "_T_ldlr_hfd.csv", header = T, row.names = 1, sep=",", as.is = T)
ad_wt.data <- read.table(file =  "_T_ad_wt.csv", header = T, row.names = 1, sep=",", as.is = T)
ad_apoe.data <- read.table(file =  "_T_ad_apoe.csv", header = T, row.names = 1, sep=",", as.is = T)
baseline.data <- read.table(file =  "_T_base.csv", header = T, row.names = 1, sep=",", as.is = T)

# create seurat objects
base <- CreateSeuratObject(counts = baseline.data, project = "base", min.cells = 3, min.features = 200)
one <- CreateSeuratObject(counts = one.data, project = "one", min.cells = 3, min.features = 200)
two <- CreateSeuratObject(counts = two.data, project = "two", min.cells = 3, min.features = 200)
three <- CreateSeuratObject(counts = three.data, project = "three", min.cells = 3, min.features = 200)
four <- CreateSeuratObject(counts = four.data, project = "four", min.cells = 3, min.features = 200)
five <- CreateSeuratObject(counts = five.data, project = "five", min.cells = 3, min.features = 200)
six <- CreateSeuratObject(counts = six.data, project = "six", min.cells = 3, min.features = 200)
seven <- CreateSeuratObject(counts = seven.data, project = "seven", min.cells = 3, min.features = 200)
eight <- CreateSeuratObject(counts = eight.data, project = "eight", min.cells = 3, min.features = 200)
apoehfd <- CreateSeuratObject(counts = hfd.data, project = "apoehfd", min.cells = 3, min.features = 200)
ldlrchow <- CreateSeuratObject(counts = chow_LDLR.data, project = "ldlrchow", min.cells = 3, min.features = 200)
ldlrhfd <- CreateSeuratObject(counts = HFD_LDLR.data, project = "ldlrhfd", min.cells = 3, min.features = 200)
apoechow <- CreateSeuratObject(counts = chow.data, project = "apoechow", min.cells = 3, min.features = 200)
adapoe <- CreateSeuratObject(counts = ad_apoe.data, project = "adapoe", min.cells = 3, min.features = 200)
adwt <- CreateSeuratObject(counts = ad_wt.data, project = "adwt", min.cells = 3, min.features = 200)

#merge
all <- merge(one, y = c(two, three, four, five, six, seven, eight, apoehfd, ldlrchow, ldlrhfd, apoechow, adwt, adapoe, base),project = "PLAQUE")

#split in list for normalization and integration 
T_cell.list <- SplitObject(all, split.by = "orig.ident")

#Normalize using vst 
for (i in 1:length(T_cell.list)) {
  T_cell.list[[i]] <- NormalizeData(T_cell.list[[i]], verbose = T)
  T_cell.list[[i]] <- FindVariableFeatures(T_cell.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = T)}

#perform integration 
T_cell.anchors <- FindIntegrationAnchors(object.list = T_cell.list, dims = 1:30, k.filter = 70)
T_cell_reclustering <- IntegrateData(anchorset = T_cell.anchors, dims = 1:30, k.weight = 70)

#Run standard downstream analysis workflow 
T_cell_reclustering <- ScaleData(T_cell_reclustering, verbose = T)
T_cell_reclustering <- RunPCA(T_cell_reclustering, npcs = 30, verbose = T)
T_cell_reclustering <- RunUMAP(T_cell_reclustering, reduction = "pca", dims = 1:30)
T_cell_reclustering <- FindNeighbors(T_cell_reclustering , dims = 1:20)
T_cell_reclustering <- FindClusters(T_cell_reclustering , resolution = 1.8)

#Normalize and Scale RNA assay
DefaultAssay(T_cell_reclustering) <- "RNA"
T_cell_reclustering <- NormalizeData(T_cell_reclustering, verbose = T)
T_cell_reclustering <- FindVariableFeatures(object = T_cell_reclustering)
T_cell_reclustering <- ScaleData(T_cell_reclustering, verbose = T)

#Figure 4a---------------------------------------------------------------------
Idents(T_cell_reclustering) <- "seurat_clusters"
levels(T_cell_reclustering)
T_cell_reclustering <- RenameIdents(object = T_cell_reclustering, 
                                    '0'= "TCL_1", 
                                    '1'= "TCL_2", 
                                    '2'= "TCL_3", 
                                    '3'= "TCL_4", 
                                    '4'= "TCL_5",
                                    '5'= "TCL_6", 
                                    '6'= "TCL_7", 
                                    '7'= "TCL_8", 
                                    '8'= "TCL_9", 
                                    '9'= "TCL_10",
                                    '10'= "TCL_11",
                                    '11'= "TCL_12",
                                    '12'= "TCL_13",
                                    '13'= "TCL_14",
                                    '14'= "TCL_15"
)

T_cell_reclustering[["TCLs"]] <- Idents(object = T_cell_reclustering)

DimPlot(T_cell_reclustering, group.by = "TCLs", label = F)

#Figure 4b
DefaultAssay(T_cell_reclustering) <- "RNA"
FeaturePlot(T_cell_reclustering, features = c("CD3E","CD8A","CD4", "IFIT1","RAG1"),
            ncol = 3, cols = c("lightblue","orange","red"), order = T)

Idents(T_cell_reclustering) <- "seurat_clusters"
T_cell_reclustering <- RenameIdents(object = T_cell_reclustering, 
                                   '0'='CD4', 
                                   '1'='CD4', 
                                   '2'='CD8', 
                                   '3'='CD8', 
                                   '4'='CD8',
                                   '5'='CD8', 
                                   '6'='Immature', 
                                   '7'='CD4', 
                                   '8'='CD4', 
                                   '9'='CD4',
                                   '10'='CD4',
                                   '11'='IFN', 
                                   '12'='CD8', 
                                   '13'='Immature', 
                                   '14'='Immature'
)

T_cell_reclustering$major <- Idents(T_cell_reclustering)
levels(T_cell_reclustering)
DimPlot(T_cell_reclustering, group.by = "major", reduction = "umap", raster = F, 
        cols = c("Immature" = "#fc8d62", 
                 "CD8" = "#66c2a5", 
                 "CD4" = "#8da0cb", 
                 "IFN" = "#e78ac3"), pt.size = 3)
                                    

#Figure 4c---------------------------------------------------------------------
Idents(T_cell_reclustering) <- "major"
DefaultAssay(T_cell_reclustering) <- "RNA"
markers <- FindAllMarkers(T_cell_reclustering, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- markers %>% group_by(cluster) %>% top_n(20,avg_log2FC)
cluster.averages <- AverageExpression(object = T_cell_reclustering,assay="RNA", return.seurat = T)
DoHeatmap(cluster.averages,  features = top10$gene, size=2, angle=0, raster=F, assay="RNA", 
          slot="scale.data", draw.lines=FALSE, label=FALSE) + scale_fill_viridis_c(option="magma")

#Figure 4e---------------------------------------------------------------------
md <- T_cell_reclustering@meta.data %>% as.data.table
cluster_abundance_sample <- md[, .N, by = c("orig.ident", "major")] %>% 
  dcast(., orig.ident ~ major, value.var = "N")

#Figure 4f---------------------------------------------------------------------
DefaultAssay(T_cell_reclustering) <- "RNA"
plot_density(T_cell_reclustering, features = "SELL")
plot_density(T_cell_reclustering, features = "CD44")
plot_density(T_cell_reclustering, features = "CCR7")
plot_density(T_cell_reclustering, features = "IL7R")
plot_density(T_cell_reclustering, features = "KLRB1")
plot_density(T_cell_reclustering, features = "FOXP3")

#Figure 4g & h-----------------------------------------------------------------
Idents(T_cell_reclustering) <- "seurat_clusters"
T_cell_reclustering <- RenameIdents(object = T_cell_reclustering, 
                                   '0'='TEM', 
                                   '1'='naive', 
                                   '2'='TEM', 
                                   '3'='naive', 
                                   '4'='TEM',
                                   '5'='TEM', 
                                   '6'='n/a', 
                                   '7'='TEM', 
                                   '8'='TCM', 
                                   '9'='TCM',
                                   '10'='Treg',
                                   '11'='n/a', 
                                   '12'='TEM', 
                                   '13'='n/a', 
                                   '14'='n/a'
)

T_cell_reclustering$activation <- Idents(T_cell_reclustering)

md <- T_cell_reclustering@meta.data %>% as.data.table
cluster_abundance_sample <- md[, .N, by = c("orig.ident", "activation")] %>% 
  dcast(., orig.ident ~ activation, value.var = "N")

#Figure 4i---------------------------------------------------------------------
Idents(T_cell_reclustering) <- "major"
CD4 <- subset(x = T_cell_reclustering, idents= c('CD4'))
Idents(CD4) <- "TCLs"
Treg_score <- c("FOXP3", "CTLA4", "IRF4", "BATF", "TNFRSF18", "TOX2")
CD4 <- AddModuleScore(object = CD4, features = Treg_score, ctrl = 5, name = 'Treg_score')
Th17_score <- c("IRF4", "IL17A", "IL17F", "RORA", "PTGFRN", "IL21", "MAF")
CD4 <- AddModuleScore(object = CD4, features = Th17_score, ctrl = 5, name = 'Th17_score')
meta.data <- as.data.frame(CD4@meta.data)

#Figure 4j---------------------------------------------------------------------
Idents(T_cell_reclustering) <- "major"
CD4 <- subset(x = T_cell_reclustering, idents= c('CD4'))

Idents(CD4) <- "orig.ident"
levels(CD4)
CD4 <- RenameIdents(object = CD4, 
                    'one'='human', 
                    'two'='human', 
                    'three'='human', 
                    'four'='human', 
                    'five'='human', 
                    'six'='human',
                    'seven'='human',
                    'eight'='human',
                    'apoehfd'='mouse', 
                    "ldlrchow"='mouse',
                    "ldlrhfd"='mouse',
                    'apoechow'='mouse',
                    'adwt'='mouse',
                    'adapoe'='mouse',
                    'base'='mouse')
CD4$species <- Idents(CD4)
levels(CD4)

Idents(CD4) <- "TCLs"
DefaultAssay(CD4) <- "RNA"
features.TH <- c("CXCR3", "IFNG", "CCR5", "HAVCR2", "TBX21", "CXCR6", "ID2", "IL12RB1", "IFNGR1", "IFNGR2", "BCL6", "ASAP1", "CXCR5",
                 "IRF4", "IL17A", "IL17F", "RORA", "PTGFRN", "IL21", "MAF", "GATA3", "IL5", "BATF", "CCR3", "IL13", "IL4", "RUNX3",
                 "IL2RA", "ICOS", "CTLA4", "FOXP3", "IL10", "TIGIT", "TGFA", "IL2RB", "LILRB4A", "NRP1", "TNFRSF4", "TNFRSF18") 
Idents(CD4)<- "species"
DotPlot(CD4, features = features.TH, cols = "RdYlBu", split.by = "species", group.by = "TCLs") & RotatedAxis()







