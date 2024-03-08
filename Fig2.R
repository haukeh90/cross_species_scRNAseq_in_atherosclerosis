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
one.data <- read.table(file =  "two.csv", header = T, row.names = 1, sep=",", as.is = T)
two.data <- read.table(file =  "three.csv", header = T, row.names = 1, sep=",", as.is = T)
three.data <- read.table(file =  "four.csv", header = T, row.names = 1, sep=",", as.is = T)
four.data <- read.table(file =  "five.csv", header = T, row.names = 1, sep=",", as.is = T)
five.data <- read.table(file =  "ten.csv", header = T, row.names = 1, sep=",", as.is = T)
six.data <- read.table(file =  "eleven.csv", header = T, row.names = 1, sep=",", as.is = T)
seven.data <- read.table(file =  "twelve.csv", header = T, row.names = 1, sep=",", as.is = T)
eight.data <- read.table(file =  "thirteen.csv", header = T, row.names = 1, sep=",", as.is = T)
hfd.data <- read.table(file =  "apoe_hfd.csv", header = T, row.names = 1, sep=",", as.is = T)
chow.data <- read.table(file =  "apoe_chow.csv", header = T, row.names = 1, sep=",", as.is = T)
chow_LDLR.data <- read.table(file =  "ldlr_chow.csv", header = T, row.names = 1, sep=",", as.is = T)
HFD_LDLR.data <- read.table(file =  "ldlr_hfd.csv", header = T, row.names = 1, sep=",", as.is = T)
ad_wt.data <- read.table(file =  "ad_wt_new.csv", header = T, row.names = 1, sep=",", as.is = T)
ad_apoe.data <- read.table(file =  "ad_apoe_new.csv", header = T, row.names = 1, sep=",", as.is = T)
baseline.data <- read.table(file =  "base.csv", header = T, row.names = 1, sep=",", as.is = T)

# create seurat objects
one <- CreateSeuratObject(counts = one.data, project = "one", min.cells = 3, min.features = 200)
two <- CreateSeuratObject(counts = two.data, project = "two", min.cells = 3, min.features = 200)
three <- CreateSeuratObject(counts = three.data, project = "three", min.cells = 3, min.features = 200)
four <- CreateSeuratObject(counts = four.data, project = "four", min.cells = 3, min.features = 200)
five <- CreateSeuratObject(counts = five.data, project = "five", min.cells = 3, min.features = 200)
six <- CreateSeuratObject(counts = six.data, project = "six", min.cells = 3, min.features = 200)
seven <- CreateSeuratObject(counts = seven.data, project = "seven", min.cells = 3, min.features = 200)
eight <- CreateSeuratObject(counts = eight.data, project = "eight", min.cells = 3, min.features = 200)
base <- CreateSeuratObject(counts = baseline.data, project = "base", min.cells = 3, min.features = 200)
apoe_hfd <- CreateSeuratObject(counts = hfd.data, project = "apoe_hfd", min.cells = 3, min.features = 200)
ldlr_chow <- CreateSeuratObject(counts = chow_LDLR.data, project = "ldlr_chow", min.cells = 3, min.features = 200)
ldlr_hfd <- CreateSeuratObject(counts = HFD_LDLR.data, project = "ldlr_hfd", min.cells = 3, min.features = 200)
apoe_chow <- CreateSeuratObject(counts = chow.data, project = "apoe_chow", min.cells = 3, min.features = 200)
ad_apoe <- CreateSeuratObject(counts = ad_apoe.data, project = "ad_apoe", min.cells = 3, min.features = 200)
ad_wt <- CreateSeuratObject(counts = ad_wt.data, project = "ad_wt", min.cells = 3, min.features = 200)

#merge
all <- merge(one, y = c(two, three, four, five, six, seven, eight,
                        apoe_hfd, ldlr_chow, ldlr_hfd, apoe_chow, ad_wt, ad_apoe, base), 
             add.cell.ids = c("one", "two", "three", "four", "five", "six", "seven", "eight",
                              "apoe_hfd", "ldlr_chow", "ldlr_hfd", "apoe_chow", "ad_wt", "ad_apoe", "base"), project = "PLAQUE")

#split in list for normalization and integration 
plaques.list <- SplitObject(all, split.by = "orig.ident")

#Normalize using vst 
for (i in 1:length(plaques.list)) {
  plaques.list[[i]] <- NormalizeData(plaques.list[[i]], verbose = T)
  plaques.list[[i]] <- FindVariableFeatures(plaques.list[[i]], selection.method = "vst", 
                                            nfeatures = 2000, verbose = T)}

#perform integration 
plaques.anchors <- FindIntegrationAnchors(object.list = plaques.list, dims = 1:30)
plaques.integrated <- IntegrateData(anchorset = plaques.anchors, dims = 1:30)

#Run standard downstream analysis workflow 
DefaultAssay(plaques.integrated) <- "integrated"
plaques.integrated <- ScaleData(plaques.integrated, verbose = T)
plaques.integrated <- RunPCA(plaques.integrated, npcs = 30, verbose = T)
plaques.integrated <- RunUMAP(plaques.integrated, reduction = "pca", dims = 1:30, verbose = T)
plaques.integrated <- FindNeighbors(plaques.integrated , dims = 1:30)
plaques.integrated <- FindClusters(plaques.integrated , resolution = 3)

#Normalize and Scale RNA assay
DefaultAssay(plaques.integrated) <- "RNA"
plaques.integrated <- NormalizeData(plaques.integrated, verbose = T)
plaques.integrated <- FindVariableFeatures(object = plaques.integrated)
plaques.integrated <- ScaleData(plaques.integrated, verbose = T)

#Figure 2b---------------------------------------------------------------------
Idents(plaques.integrated) <- "seurat_clusters"
plaques.integrated <- RenameIdents(object = plaques.integrated, 
                                   '0'='CD4 T cell', 
                                   '1'='CD8 T cell', 
                                   '2'='CD4 T cell', 
                                   '3'='CD8 T cell', 
                                   '4'='CD4 T cell',
                                   '5'='NK cell', 
                                   '6'='CD4 T cell', 
                                   '7'='NK T cell', 
                                   '8'='NK cell', 
                                   '9'='B cell',
                                   '10'='Neutrophil',
                                   '11'='CD8 T cell',
                                   '12'='CD4 T cell', 
                                   '13'='CD4 T cell', 
                                   '14'='B cell', 
                                   '15'='NK T cell', 
                                   '16'='CD4 T cell', 
                                   '17'='Monocyte', 
                                   '18'='CD8 T cell', 
                                   '19'='NK cell', 
                                   '20'='Monocyte', 
                                   '21'='CD4 T cell', 
                                   '22'='Myeloid DC', 
                                   '23'='Macrophage', 
                                   '24'='Plasmocytoid DC', 
                                   '25'='Macrophage', 
                                   '26'='Other T cell', 
                                   '27'='B cell', 
                                   '28'='CD8 T cell', 
                                   '29'='Plasmocytoid DC', 
                                   '30'='Proliferating',
                                   '31'='Other T cell',
                                   '32'='Other T cell',
                                   '33'='CD4 T cell',
                                   '34'='B cell',
                                   '35'='Monocyte',
                                   '36'='CD8 T cell',
                                   '37'='CD4 T cell',
                                   '38'='CD4 T cell',
                                   '39'='Proliferating',
                                   '40'='NK cell',
                                   '41'='Myeloid DC',
                                   '42'='NK cell',
                                   '43'='ILC',
                                   '44'='CD8 T cell',
                                   '45'='Monocyte',
                                   '46'='Mast cell',
                                   '47'='Monocyte',
                                   '48'='Macrophage',
                                   '49'='Myeloid DC',
                                   '50'='ILC'
)
plaques.integrated$lineages <- Idents(plaques.integrated)

cols <- c('#16989F','#1C70A7','#B883FE', '#0C2666','#37C628', '#942192' ,'#FC67A2','#E58709','#FA0008', '#108002','#0433FF','#15AFF5','#F78386','#D6D6D6')

DimPlot(plaques.integrated, group.by = "lineages", reduction = "umap", label = F, cols = cols, pt.size= 2) 

#Figure 2c---------------------------------------------------------------------
Idents(plaques.integrated) <- "orig.ident"
plaques.integrated <- RenameIdents(object = plaques.integrated, 'one'='human', 'two'='human', 'three'='human', 'four'='human',
                                   'five'='human', 'six'='human','seven'='human','eight'='human', 'apoe_hfd'='apoe', "ldlr_chow"='ldlr',
                                   "ldlr_hfd"='ldlr', 'apoe_chow'='apoe', 'ad_wt'='adv', 'ad_apoe'='adv', 'base'='apoe')
plaques.integrated$grouping <- Idents(plaques.integrated)
DimPlot(plaques.integrated, reduction = "umap", split.by='grouping', ncol=2, pt.size = 1, cols=c('#24318A', '#CE313F', '#FF7F00', '#43942F'))

#Figure 2d---------------------------------------------------------------------
DefaultAssay(plaques.integrated) <- "RNA"
FeaturePlot(plaques.integrated, features = c("CD14", "CD3E", "CD68", "CD79A"),
            ncol = 2, cols = c("lightblue","orange","red"), order = T)

#Figure 2e---------------------------------------------------------------------
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(plaques.integrated) <- "RNA"
plaques.integrated <- NormalizeData(plaques.integrated, verbose = T, assay="RNA")
plaques.integrated <- CellCycleScoring(plaques.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)

DimPlot(plaques.integrated, reduction = "umap", label=F)
plaques.integrated[["percent.mt"]] <- PercentageFeatureSet(plaques.integrated, pattern = "^MT-")
FeaturePlot(plaques.integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
            ncol = 3, cols = c("lightblue","orange","red"), order = T)

#Figure 2f---------------------------------------------------------------------
VlnPlot(plaques.integrated, group.by = "lineages", features =  c("nFeature_RNA", "percent.mt", "S.Score", "G2M.Score"),
        ncol = 1, pt.size = 0, cols = cols)

#Figure 2g---------------------------------------------------------------------
features <- c("PTPRC", "CD3E","CD8A","CD4", "IL7R","NKG7", "ITGAX", "FLT3", "AREG", "GATA3", "CD19", "CD79A", "ITGAM",
                  "ADGRE1", "CSF1R", "CD68", "CD14", "C1QC", "CD9", "LY6G")
VlnPlot(plaques.integrated, group.by = "lineages", features =  features,
        ncol = 1, pt.size = 0, cols = cols)

#Figure 2h---------------------------------------------------------------------
DefaultAssay(plaques.integrated) <- "RNA"
Idents(plaques.integrated) <- "lineages"
plaques.integrated <- NormalizeData(plaques.integrated, verbose = T)
plaques.integrated <- FindVariableFeatures(object = plaques.integrated)
plaques.integrated <- ScaleData(plaques.integrated, verbose = T)
markers <- FindAllMarkers(plaques.integrated, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(plaques.integrated, features = top10$gene, size = 1, angle = 0, raster = F) & NoLegend() &
  scale_fill_viridis_c(option="inferno",na.value = "white")
