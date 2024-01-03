#(Reference) Single Cell Genomics, Transcriptomics & Proteomics
#(Datasets)
# https://www.10xgenomics.com/resources/datasets/mouse-brain-section-coronal-1-standard
# - 10x Genomics, spatial gene expression
# - Mouse Brain coronal section (v1)
# 1. Feature / barcode matrix HDF5(filtered)
# 2. Spatial imaging data

library(Seurat)
library(tidyverse)
library(patchwork)

#load data (https://satijalab.org/seurat/reference/load10x_spatial)

brain <- Load10X_Spatial('Data/', 
                         filename = "Visium_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5",
                         assay = 'spatial',
                         slice = 'slice1',
                         filter.matrix = TRUE,
                         to.upper = FALSE,
                         image = NULL)

View(brain@meta.data)

# nCount_spatial: total number of transcripts detected in each spatial location
# nFeature_spatial: total number of unique genes detected in each spatial location

# VlnPlot & SpatialFeaturePlot
plot1 <- VlnPlot(brain, features = "nCount_spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features =
                              "nCount_spatial") + theme(legend.position = "right")

# wrap_plots
wrap_plots(plot1, plot2)

plot3 <- VlnPlot(brain, features = 'nFeature_spatial', pt.size = 0.1) + NoLegend()
plot4 <- SpatialFeaturePlot(brain, features = "nFeature_spatial") + theme(legend.position = 'right')

wrap_plots(plot3, plot4)

# RNA-seq data pre-processing
range(brain@meta.data$nCount_spatial)
range(brain@meta.data$nFeature_spatial)

# scRNA-seq
brain <- NormalizeData(brain)
brain <- FindVariableFeatures(brain)
brain <- ScaleData(brain)

# Spatial RNA-seq Data pre-processing
## Used for normalization and variance stabilization

brain <- SCTransform(brain, assay = 'spatial', verbose = FALSE)

# Dimensionality reduction, clustering, and visualization
brain <- RunPCA(brain, asay = 'SCT', verbose = FALSE)
brain <- FindNeighbors(brain, reduction = 'pca', dims =1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = 'pca', dims =1:30)

p1 <-DimPlot(brain, reduction ='umap', label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1+p2

