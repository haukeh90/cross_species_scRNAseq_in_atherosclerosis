library(dplyr)
library(data.table)
library(Seurat)
library(ggplot2)
library(viridis)
options(Seurat.object.assay.version = "v3")

#set the working directory to the desired folder path
setwd("~/cross_species_analysis")

# Run CCA integration----------------------------------------------------------

# load data files
one.data <- read.table(file = "two.csv", header = T, row.names = 1, sep=",", as.is = T)
two.data <- read.table(file = "three.csv", header = T, row.names = 1, sep=",", as.is = T)
three.data <- read.table(file = "four.csv", header = T, row.names = 1, sep=",", as.is = T)
four.data <- read.table(file = "five.csv", header = T, row.names = 1, sep=",", as.is = T)
five.data <- read.table(file = "ten.csv", header = T, row.names = 1, sep=",", as.is = T)
six.data <- read.table(file = "eleven.csv", header = T, row.names = 1, sep=",", as.is = T)
seven.data <- read.table(file = "twelve.csv", header = T, row.names = 1, sep=",", as.is = T)
eight.data <- read.table(file = "thirteen.csv", header = T, row.names = 1, sep=",", as.is = T)

# create seurat objects
one <- CreateSeuratObject(counts = one.data, project = "one", min.cells = 3, min.features = 200)
two <- CreateSeuratObject(counts = two.data, project = "two", min.cells = 3, min.features = 200)
three <- CreateSeuratObject(counts = three.data, project = "three", min.cells = 3, min.features = 200)
four <- CreateSeuratObject(counts = four.data, project = "four", min.cells = 3, min.features = 200)
five <- CreateSeuratObject(counts = five.data, project = "five", min.cells = 3, min.features = 200)
six <- CreateSeuratObject(counts = six.data, project = "six", min.cells = 3, min.features = 200)
seven <- CreateSeuratObject(counts = seven.data, project = "seven", min.cells = 3, min.features = 200)
eight <- CreateSeuratObject(counts = eight.data, project = "eight", min.cells = 3, min.features = 200)

#merge
all <- merge(one, y = c(two, three, four, five, six, seven, eight), 
             add.cell.ids = c("one", "two", "three", "four", "five", "six", "seven", "eight"),
             project = "PLAQUE")

#split in list for normalization and integration 
plaques.list <- SplitObject(all, split.by = "ident")

#Normalize using vst 
for (i in 1:length(plaques.list)) {
  plaques.list[[i]] <- NormalizeData(plaques.list[[i]], verbose = T)
  plaques.list[[i]] <- FindVariableFeatures(plaques.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = T)}

#perform integration 
plaques.anchors <- FindIntegrationAnchors(object.list = plaques.list, dims = 1:30, verbose = T)
plaques.integrated <- IntegrateData(anchorset = plaques.anchors, dims = 1:30, verbose = T)

#Run standard downstream analysis workflow 
DefaultAssay(plaques.integrated) <- "integrated"
plaques.integrated <- ScaleData(plaques.integrated, verbose = T)
plaques.integrated <- RunPCA(plaques.integrated, npcs = 30, verbose = T)
plaques.integrated <- RunUMAP(plaques.integrated, reduction = "pca", dims = 1:30, verbose = T)
plaques.integrated <- FindNeighbors(plaques.integrated , dims = 1:20)
plaques.integrated <- FindClusters(plaques.integrated , resolution = 0.8)

# curated annotations
Idents(plaques.integrated) <- "seurat_clusters"
plaques.integrated <- RenameIdents(object = plaques.integrated, 
                                   '0'='1', 
                                   '1'='2', 
                                   '2'='3', 
                                   '3'='4', 
                                   '4'='5',
                                   '5'='6', 
                                   '6'='7', 
                                   '7'='8', 
                                   '8'='9', 
                                   '9'='10',
                                   '10'='11',
                                   '11'='12',
                                   '12'='13',  
                                   '13'='14', 
                                   '14'='15', 
                                   '15'='16', 
                                   '16'='17', 
                                   '17'='18', 
                                   '18'='19', 
                                   '19'='20', 
                                   '20'='21', 
                                   '21'='22', 
                                   '22'='23'
)
                                
plaques.integrated$seurat_clusters <- Idents(plaques.integrated)
levels(plaques.integrated)

#Normalize and Scale RNA assay
DefaultAssay(plaques.integrated) <- "RNA"
plaques.integrated <- NormalizeData(plaques.integrated, verbose = T)
plaques.integrated <- FindVariableFeatures(object = plaques.integrated)
plaques.integrated <- ScaleData(plaques.integrated, verbose = T)

#Figure 1c---------------------------------------------------------------------
DimPlot(plaques.integrated, reduction = "umap", label= T) 
ggsave("umap_1c.pdf", units = "cm", height = 12, width = 16)

#Figure 1d---------------------------------------------------------------------
DefaultAssay(plaques.integrated) <- "RNA"
FeaturePlot(plaques.integrated, features = c("CD3E", "CD8A", "CD4", "LYZ", "CD68", "CD19"),
            ncol = 3, cols=c("lightblue","orange","red"), order = TRUE)

#Figure 1e---------------------------------------------------------------------
Idents(plaques.integrated) <- "seurat_clusters"
levels(plaques.integrated)
plaques.integrated <- RenameIdents(object = plaques.integrated, 
                                   '1'='CD4+ TCM', 
                                   '2'='CD4+ naïve I', 
                                   '3'='CD8+ TEM/TCM', 
                                   '4'='CD8+ naïve', 
                                   '5'='CD16+ NK cell',
                                   '6'='CD16- NK cell', 
                                   '7'='Classical Monocyte', 
                                   '8'='Treg', 
                                   '9'='Natural killer T cell', 
                                   '10'='B cell',
                                   '11'='CD4+ TEM I',
                                   '12'='pDC',
                                   '13'='cDC 1',
                                   '14'='CD4+ naïve II',
                                   '15'='C1Q+ Macrophage', 
                                   '16'='Proliferating', 
                                   '17'='TREM2+ Macrophage', 
                                   '18'='CD4+ TEM II',
                                   '19'='Non-classical Monocyte',
                                   '20'='Plasma cell',
                                   '21'='cDC1',
                                   '22'='ILC',
                                   '23'='Non-leukocyte'
)

plaques.integrated$labels <- Idents(plaques.integrated)

VlnPlot(plaques.integrated, features = c("PTPRC", "CD3E","CD8A",
                                         "CD4", "NKG7","FLT3",
                                         "CD79A","CD68", "CD14", "C1QA"), 
        ncol= 1, pt.size = 0, group.by = "labels")
                                          
#Figure 1f---------------------------------------------------------------------
DefaultAssay(plaques.integrated) <- "RNA"
Idents(plaques.integrated) <- "seurat_clusters"
plaques.integrated <- NormalizeData(plaques.integrated, verbose = T)
plaques.integrated <- FindVariableFeatures(object = plaques.integrated)
plaques.integrated <- ScaleData(plaques.integrated, verbose = T)
markers <- FindAllMarkers(plaques.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(plaques.integrated, features = top10$gene, size=1, angle=0, raster=FALSE) & NoLegend() &
  scale_fill_viridis_c(option="inferno",na.value = "white")

#Figure 1g & h-----------------------------------------------------------------
md <- plaques.integrated@meta.data %>% as.data.table
cluster_abundance_sample <- md[, .N, by = c("orig.ident", "labels")] %>% 
  dcast(., orig.ident ~ labels, value.var = "N")






