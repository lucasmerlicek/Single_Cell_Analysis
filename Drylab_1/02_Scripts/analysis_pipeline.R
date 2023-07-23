# Load Libraries for the Analysis 
library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Matrix)
library(dplyr)
library(svglite)
#library(sctransform)

###################################### Step 0 General information from the Paper 
# Eventually we were given the data after correction -> 32464 cells that were doublet correction
# and cells with > 30% MT were removed
# Matrices were merged for the: 10,327 cells for telencephalic sample,
#                               5,347 cells for ChP sample 1,
#                               8,603 cells for ChP sample 2,
#                               8,573 cells for ChP sample 3

# Eventually just start of analysis should be the normalization Step 3
# svg for vector graphics works for everything except heat maps 
# use ggsave for heatmaps with the png device 

###################################### Step 0 Define Variables

Analysis_Number <- 10
PCA <- 10 # Define the number of PCAs to take into account for the analysis
min_feature <- 1000 # min number of different genes a cell must have to be included 
max_feature <- 5000 # max number of genes a cell can have to be included 
mt_percent <- 10 # max percent of genes a cell can have, that are form the mitochondrial transcriptome
variable_Features <- 3000 # Number of variable genes to include into differential expression analysis
Resolution <- 0.3 # Set resolution for the clustering

dir.create(file.path("03_PLots",paste0("Analysis_N",Analysis_Number)))

df <- data.frame(var1 = c(PCA, min_feature, max_feature,mt_percent, variable_Features, Resolution),
                 row.names = c("PCA", "min_feature", "max_feature",
                               "mt_percent", "variable_Feature", "Resolution"))
write.table(df, file = paste0("../03_Plots/Analysis_N",Analysis_Number,
                              "/Analysis_Parameters_N",Analysis_Number), row.names = TRUE)

###################################### Step 1 create a Seurat object 
counts <- readMM("../01_Raw_Data/counts.mtx.gz")
features <- read.csv("../01_Raw_Data/features.tsv.gz", stringsAsFactors=F, sep="\t", header=F)
#barcodes <- read.table("01_Raw_Data/metadata.tsv.gz", stringsAsFactors=F)[,1]
rownames(counts) <- make.unique(features[,1])

# Create barcodes from metaData
metaData <- read.csv("../01_Raw_Data/metadata.tsv.gz", stringsAsFactors=F, sep="\t", header=T)
metaData %>% row.names() -> barcodes
colnames(counts) <- barcodes

# Create Seurat Object with the metaData for the different samples (Tel,ChP1,ChP2,ChP3)
seurat <- CreateSeuratObject(counts, project="Data1", meta.data = metaData)


###################################### Step 2 Quality Control 
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT[-\\.]")
plot1 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot2 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0)

# save violin plots as pdf graphics
svg(paste0("../03_Plots/Analysis_N",Analysis_Number,"/QC_Violin_Plot_Points.svg"))
plot1
dev.off() 

svg(paste0("../03_Plots/Analysis_N",Analysis_Number,"/QC_Violin_Plot.svg"))
plot2
dev.off() 

# create correlation plots for mt % vs RNA count and Genes(Feature) vs RNA count
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1$layers[[1]]$aes_params$alpha <- 0.5

plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2$layers[[1]]$aes_params$alpha <- 0.5

# save correlation plots as pdf graphics 
svg(paste0("../03_Plots/Analysis_N",Analysis_Number,"/Correlation_Plot_PCA.svg"))
plot1 + plot2
dev.off()

# Eventually not necessary if we want to copy the analysis in the paper, just in case use the same threshold as used in the paper 30%
seurat <- subset(seurat, subset = nFeature_RNA > min_feature & nFeature_RNA < max_feature & percent.mt < mt_percent)

# How the data looks after removing potentially apoptotic cells 
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot1$layers[[1]]$aes_params$alpha <- 0.5
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2$layers[[1]]$aes_params$alpha <- 0.5

svg(paste0("../03_Plots/Analysis_N",Analysis_Number,"/Correlation_Plot_After_Correction_PCA.svg"))
plot1 + plot2
dev.off()

# Optional could do DoubletFinder but not sure it is necessary since the RNA counts aren't that high but there are certain populations visible
# if True paper already performed doublet removal


###################################### Step 3 Normalization using TPM afterwards log transformed with (X+1)

#seurat <- NormalizeData(seurat)
#Integration
seurat_tel <- subset(seurat, subset = orig.ident=="T")
seurat_1 <- subset(seurat, subset = orig.ident=="1")
seurat_2 <- subset(seurat, subset = orig.ident=="2")
seurat_3 <- subset(seurat, subset = orig.ident=="3")

#Seurat
seurat_tel <- NormalizeData(seurat_tel) %>% FindVariableFeatures(nfeatures = 3000)
seurat_1 <- NormalizeData(seurat_1) %>% FindVariableFeatures(nfeatures = 3000)
seurat_2 <- NormalizeData(seurat_2) %>% FindVariableFeatures(nfeatures = 3000)
seurat_3 <- NormalizeData(seurat_3) %>% FindVariableFeatures(nfeatures = 3000)

seurat_objs <- list(tel = seurat_tel, chp1 = seurat_1, chp2=seurat_2, chp3=seurat_3)
anchors <- FindIntegrationAnchors(object.list = seurat_objs, dims = 1:30)
seurat <- IntegrateData(anchors, dims = 1:30)


seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat, npcs = 50)
seurat <- RunUMAP(seurat, dims = 1:20)
seurat <- FindNeighbors(seurat, dims = 1:20) %>% FindClusters(resolution = 0.3)
saveRDS(seurat, file="integrated_seurat.rds")

DefaultAssay(seurat) <- "RNA"
plot1 <- UMAPPlot(seurat, group.by="orig.ident")
plot2 <- UMAPPlot(seurat, label = T)
plot3 <- FeaturePlot(seurat, c("FOXG1","EMX1","DLX2","LHX9"), ncol=2, pt.size = 0.1)
((plot1 / plot2) | plot3) + plot_layout(width = c(1,2))


###################################### Step 4 Feature selection for following heterogeneity analysis 

seurat <- FindVariableFeatures(seurat, nfeatures = variable_Features)

top_features <- head(VariableFeatures(seurat), 20)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top_features, repel = TRUE)

svg(paste0("../03_Plots/Analysis_N",Analysis_Number,"/Variable_Genes_Plot_PCA.svg"))
plot1 + plot2
dev.off()


###################################### Step 5 Data scaling 

# Data scaling with regressing out the number of features per cell and the percentage of mitochondrial RNA
seurat <- ScaleData(seurat, vars.to.regress = c("nFeature_RNA", "percent.mt"))

###################################### Step 6 Linear dimensionality reduction using principal component analysis (PCA)
# Dimension reduction using PCA -> gives the first 50 PCAs
seurat <- RunPCA(seurat, npcs = 30)
# Shows the PCs and the variation they explain 
plot1 <- ElbowPlot(seurat, ndims = ncol(Embeddings(seurat, "pca")))

svg(paste0("../03_Plots/Analysis_N",Analysis_Number,"/Elbowplot_PCA",PCA,".svg"))
plot1
dev.off()
# Heatmaps for the top 20 pCAs
plot1 <- PCHeatmap(seurat, dims = 1:PCA, cells = 500, balanced = TRUE, ncol = 4, fast = FALSE)

plot1 
ggsave(filename = paste0("Heatmap_PCA_new_",PCA,".png"), plot = plot1 ,device = "png", width = 15, height = 10, path = paste0("03_Plots/Analysis_N",Analysis_Number,"/"))


###################################### Step 7 Non-linear dimension reduction for visualization 

seurat <- RunTSNE(seurat, dims = 1:PCA)
seurat <- RunUMAP(seurat, dims = 1:PCA)

plot1 <- TSNEPlot(seurat)
plot2 <- UMAPPlot(seurat)

svg(paste0("../03_Plots/Analysis_N",Analysis_Number,"/2D_Reduction_Plots_PCA",PCA ,".svg"))
plot1 + plot2
dev.off()

# Data is from Pellegrini_2020_Science

# Markers for different cell types
# Cell cylce ("MKI67") -> There seems to be a cluster for this marker, maybe need to regress cell cycle genes out 
# Neurons ("DCX", )
# Choroid plexus ChP ("CLIC6", "HTR2C", "TTR")
# ChP development ("MSX1", "PAX6") stromal ("COL1A1")
# ChP immature/hem ("OTX2", "RSPO3", "PAX6")
# ChP mature ("TTR", "KRT18", "NME5")
# ChP stromal ("LUM", "DCN", "DLK1")
# Neural prog /ImmN ("NES", "SOX2")

# plot_Neurons <- FeaturePlot(seurat, c("MKI67","DCX"),
#                             ncol=2, reduction = "umap")
# 
# plot_Neural_Prog <- FeaturePlot(seurat, c("MKI67","NES", "SOX2"),
#                                 ncol=2, reduction = "umap")
# 
# plot_ChP <- FeaturePlot(seurat, c("MKI67","CLIC6", "HTR2C", "TTR"),
#                         ncol=2, reduction = "umap")
# 
# plot_ChP_Imm <- FeaturePlot(seurat, c("MKI67","OTX2", "RSPO3", "PAX6"),
#                             ncol=2, reduction = "umap")
# 
# plot_ChP_Mature <- FeaturePlot(seurat, c("MKI67","TTR", "KRT18", "NME5"),
#                                ncol=2, reduction = "umap")
# 
# plot_ChP_Stromal <- FeaturePlot(seurat, c("MKI67","LUM", "DCN", "DLK1"),
#                                 ncol=2, reduction = "umap")
# 
# plot_Neurons
# 
# plot_Neural_Prog
# 
# plot_ChP_Imm
# 
# plot_ChP_Mature
# 
# plot_ChP_Stromal


###################################### Step 8 Cluster the cells

seurat <- FindNeighbors(seurat, dims = 1:PCA)
seurat <- FindClusters(seurat, resolution = Resolution)

# Visualize the clustering 
plot1 <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot2 <- DimPlot(seurat, reduction = "umap", label = TRUE)

svg(paste0("../03_Plots/Analysis_N",Analysis_Number,"/Clustering_Plot_PCA",PCA ,".svg"))
plot1 + plot2
dev.off()

###################################### Step 9 Annotate cell clusters

Cluster_Markers <- c("DCX", # Neurons
                     "CLIC6", "HTR2C", "TTR", # ChP
                     "MSX1", # ChP developmental stage
                     "OTX2", "RSPO3", "PAX6", # ChP immature/hem
                     "TTR", "KRT18", "NME5", # ChP mature 
                     "LUM", "DCN", "DLK1","COL1A1", # ChP stromal 
                     "NES", "SOX2") # NPC

plot1 <- DoHeatmap(seurat, features = Cluster_Markers) + NoLegend()

ggsave(filename = paste0("Heatmap_Markers_Paper_PCA",PCA,".png"), plot = plot1,device = "png", width = 15, height = 10, path = paste0("../03_Plots/Analysis_N",Analysis_Number,"/"))

# Identification of marker genes from clustering
cl_markers <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.2))
cl_markers %>% group_by(cluster) %>% top_n( n = 10, wt = avg_log2FC) %>% print(n = 100)

top10_cl_markers <- cl_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
plot1 <- DoHeatmap(seurat, features = top10_cl_markers$gene) + NoLegend()

ggsave(filename = paste0("Heatmap_Markers_Automatic_PCA",PCA,".png"), plot = plot1,device = "png", width = 18, height = 12, path = paste0("../03_Plots/Analysis_N",Analysis_Number,"/"))

#comparison of simmilar clusters
plot1 <- FeaturePlot(seurat, c("STMN2","DCX"), ncol = 2)
plot2 <- VlnPlot(seurat, features = c("STMN2","DCX"), pt.size = 0, ncol = 1)
plot_tsne <- DimPlot(seurat, reduction = "tsne", label = TRUE)
plot3 <- plot1 + plot2 + plot_tsne + plot_layout(widths = c(1, 1))

ggsave(filename = paste0("Cell_Types_STMN2_DCX",PCA,".png"), plot = plot3,device = "png", width = 15, height = 10, path = paste0("../03_Plots/Analysis_N",Analysis_Number,"/"))

saveRDS(seurat, file="integrated_seurat_clusters.rds")
seurat <- readRDS("integrated_seurat_clusters.rds")

#Mannualy determining marker genes

"DCX", # Neurons
"CLIC6", "HTR2C", "TTR", # ChP
"MSX1", # ChP developmental stage
"OTX2", "RSPO3", "PAX6", # ChP immature/hem
"TTR", "KRT18", "NME5", # ChP mature 
"LUM", "DCN", "DLK1","COL1A1", # ChP stromal 
"NES", "SOX2") # NPC

#Paper--
#"DCX", # Neurons: Cluster 2
#"CLIC6", "HTR2C", "TTR", # ChP: 0, 4 ,11
#"MSX1", # ChP developmental stage: 0, 4+, 8+, 9
#"OTX2", "RSPO3", "PAX6", # ChP immature/hem: 6+, 8, 12
#"TTR", "KRT18", "NME5", # ChP mature: 0+, 4, 11+ 
#"LUM", "DCN", "DLK1","COL1A1", # ChP stromal: 6-, 7+, 13+
#"NES", "SOX2") # NPC: 9, 11

#VoxHunt mapping
# Load VoxHunt
library(voxhunt)


# Point VoxHunt to ABA expression data
load_aba_data('../01_Raw_Data/voxhunt_data/')

# Find 300 most variable genes from the E13.5 mouse brain
genes_use <- variable_genes('E13', 300)$gene

# Calculate the similarity map of a seurat object to the E13.5 mouse brain
vox_map <- voxel_map(seurat, genes_use=genes_use)

# Plot the result
vox_map_plot <- plot_map(vox_map)
ggsave(filename = paste0("Vox_Map_E13.png"), plot = vox_map_plot ,device = "png", width = 15, height = 10, path = paste0("../03_Plots/Analysis_N",Analysis_Number,"/"))

#Rough brain regions
p1 <- voxhunt::plot_annotation('E13', show_legend=T)
ggsave(filename = paste0("Brain_Map_E13.png"), plot = p1 ,device = "png", width = 15, height = 10, path = paste0("../03_Plots/Analysis_N",Analysis_Number,"/"))


#Finding regional markers in the ABA dataset
regional_markers <- structure_markers('E13') %>%
  group_by(group) %>%
  top_n(10, auc) %>% 
  {unique(.$gene)}
head(regional_markers)

vox_map_2 <- voxel_map(
  example_seurat, 
  stage = 'E13', 
  group_name = 'cluster', 
  genes_use = regional_markers
)

ggsave(filename = paste0("Vox_cell_types.png"), plot = plot_structure_similarity(vox_map_2, cluster=F) ,device = "png", width = 6, height = 4, path = paste0("../03_Plots/Analysis_N",Analysis_Number,"/"))



ggsave(filename = paste0("Vox_cell_types_map_2.png"), plot = plot_map(vox_map_2, nrow=2) & no_legend() ,device = "png", width = 6, height = 4, path = paste0("../03_Plots/Analysis_N",Analysis_Number,"/"))


plot3d <- plot_map_3d(
  vox_map_2, 
  show_group = 'ctx_cerebral', 
  width = 800, 
  height = 600
)

#save 3d plot and open in browser
temp_file <- tempfile(fileext = ".html")
htmlwidgets::saveWidget(plot3d, file = temp_file, selfcontained = TRUE)

cell_assignment <- assign_cells(vox_map)
head(cell_assignment)

# Test for markers 
# plot1 <- FeaturePlot(seurat,c("NEUROD2","NEUROD6"),ncol = 1)
# plot2 <- VlnPlot(seurat,features = c("NEUROD2", "NEUROD6"),pt.size = 0)
# plot1 + plot2 + plot_layout(widths = c(1,2))

# new_ident <- setNames(c("Dorsal telen. NPC",
#                         "Midbrain-hindbrain boundary neuron",
#                         "Dorsal telen. neuron",
#                         "Dien. and midbrain excitatory neuron",
#                         "MGE-like neuron","G2M dorsal telen. NPC",
#                         "Dorsal telen. IP","Dien. and midbrain NPC",
#                         "Dien. and midbrain IP and excitatory early neuron",
#                         "G2M Dien. and midbrain NPC",
#                         "G2M dorsal telen. NPC",
#                         "Dien. and midbrain inhibitory neuron",
#                         "Dien. and midbrain IP and early inhibitory neuron",
#                         "Ventral telen. neuron",
#                         "Unknown 1",
#                         "Unknown 2"),
#                       levels(seurat))
# seurat <- RenameIdents(seurat, new_ident)
# DimPlot(seurat, reduction = "umap", label = TRUE) + NoLegend()

################################ Cell-cycling
# seurat <- CellCycleScoring(seurat,
#                                   s.features = cc.genes$s.genes,
#                                   g2m.features = cc.genes$g2m.genes,
#                                   set.ident = TRUE)
# seurat <- ScaleData(seurat, vars.to.regress = c("S.Score", "G2M.Score"))



