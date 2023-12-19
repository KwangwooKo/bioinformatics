# Reference: GSE132771 (NML1, NML2, NML3)

library(Seurat)
library(tidyverse)
library(Matrix)
library(patchwork)
library(DT)


NML1 <- Read10X("Data/NML1/")
NML2 <- Read10X("Data/NML2/")
NML3 <- Read10X("Data/NML3/")



NML1_seurat <- CreateSeuratObject(counts=NML1,
                                  min.cells = 3,
                                  project = 'NML1_seurat')

NML2_seurat <- CreateSeuratObject(counts=NML2,
                                  min.cells = 3,
                                  project = 'NML2_seurat')

NML3_seurat <- CreateSeuratObject(counts=NML3,
                                  min.cells = 3,
                                  project = 'NML3_seurat')

# Merge all

MergedNML <- merge(NML1_seurat, y=c(NML2_seurat, NML3_seurat),
                   add.cell.ids =c('NML1', 'NML2', 'NML3'),
                   project='MergedNML.seurat')

MergedNML[['percentage.mt']] <- PercentageFeatureSet(MergedNML, pattern='^MT-')
MergedNML <- subset(MergedNML,
                      subset = nFeature_RNA >200 &
                        percentage.mt <10)

MergedNML <- NormalizeData(MergedNML)
MergedNML <- FindVariableFeatures(MergedNML,
                                    selection.method = 'vst',
                                    nfeatures = 2000)
MergedNML <- ScaleData(MergedNML)

VlnPlot(MergedNML,
        features = c('nFeature_RNA', 'nCount_RNA', 'percentage.mt'),
        pt.size = 0.25,
        ncol=3)

# Dimensionality reduction

MergedNML <- RunPCA(MergedNML)
MergedNML <- FindNeighbors(MergedNML, reduction = 'pca', dims = 1:50)

# Higher resolution => more clusters
## I think that it would be better to start from 0.05
MergedNML <- FindClusters(MergedNML, resolution = 0.01)
MergedNML <- RunUMAP(MergedNML, reduction = 'pca', dims = 1:50)

DimPlot(MergedNML, reduction = 'umap', label = TRUE)
DimPlot(MergedNML, reduction = 'umap', label = TRUE, group.by ='orig.ident')

View(MergedNML@meta.data)

# Find markers and differentially expressed genes
# Marker genes for Epithelial, endothelial, mesenthymal, and Immune cells
FeaturePlot(MergedNML, features = c("EPCAM", "CLDN5", "COL1A2", 'PTPRC'),
            cols = c('lightgrey', 'blue'))

# DefaultAssay(MergedNML) <- 'integrated'
DefaultAssay(MergedNML) <- 'RNA'
MergedNML <- JoinLayers(MergedNML)

# Find markers for cluster 2
cluster2.marker <- FindMarkers(MergedNML, ident.1 = 2, min.pct = 0.5,
                               only.pos = TRUE, logfc.threshold = 0.25)

# Find top 10 genes of cluster 2 
cluster2.top10 <- cluster2.marker %>% arrange(desc(avg_log2FC))
cluster2.top10 <- dplyr::slice(cluster2.top10, 1:10)

# Make data table of cluster 2
datatable(cluster2.top10,
          extensions =c('KeyTable', 'FixedHeader'),
          caption = 'Table 1: Cluster 2 genes',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10,
                         lengthMenu = c("10", "15", "50", "70"))) %>%
  formatRound(columns = c(2:6), digits=2)


# Find markers of all clusters
All.markers <- FindAllMarkers(MergedNML, min.pct = 0.5, only.pos = TRUE,
                              logfc.threshold = 0.25)

All.top10 <- All.markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))
All.top10 <- dplyr::slice(All.top10, 1:10)

# After running All.top10, you will notice that gene column is the last column
## so, to make gene column to be the first,
All.top10 <- All.top10 %>% select(gene, everything())

DoHeatmap(MergedNML, features = All.top10$gene)

datatable(All.top10,
          extensions =c('KeyTable', 'FixedHeader'),
          caption = 'Table 1: All Clusters : genes',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10,
                         lengthMenu = c("10", "15", "50", "70"))) %>%
  formatRound(columns = c(2:6), digits=2)

# Find markers for cluster 0
cluster0.marker <- FindMarkers(MergedNML, ident.1 = 0, min.pct = 0.5,
                               only.pos = TRUE, logfc.threshold = 0.25)
# Find top 10 genes of cluster 0 
cluster0.top10 <- cluster0.marker %>% arrange(desc(avg_log2FC))
cluster0.top10 <- dplyr::slice(cluster0.top10, 1:10)

# Make data table of cluster 0
datatable(cluster0.top10,
          extensions =c('KeyTable', 'FixedHeader'),
          caption = 'Table 1: Cluster 0 genes',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10,
                         lengthMenu = c("10", "15", "50", "70"))) %>%
  formatRound(columns = c(2:6), digits=2)

# Visualize marker gene expression for cluster 2
FeaturePlot(MergedNML, features = c("SFTPC", "SFTPB", "NAPSA", 'SCGB3A2'),
            cols = c('lightgrey', 'blue'))
DotPlot(MergedNML, features = c("SFTPC", "SFTPB", "NAPSA", 'SCGB3A2','AGER'))
VlnPlot(MergedNML, features = c("SFTPC", "SFTPB", "NAPSA", 'SCGB3A2','AGER'))

# Differential gene expression
## Find markers distinguishing cluster 0 from cluster 3
cluster0_3.markers <- FindMarkers(MergedNML,
                                  ident.1 = 0, ident.2 = 3, min.pct = 0.25)
# (note) positive of avg.log2FC means that a gene is more expressed in cluster 0
## than cluster 3

# Find all markers distinguishing cluster 0 from clusters 3 and 6
cluster0_3and6.markers <- FindMarkers(MergedNML,
                                      ident.1 = 0, ident.2 = c(3,6), min.pct = 0.25)

# Find markers of cluster 5
cluster5.markers <- FindMarkers(MergedNML, ident.1 = 5, min.pct = 0.5,
                                only.pos = TRUE, logfc.threshold = 0.25)

FeaturePlot(MergedNML, features=c("MS4A2", "CPA3", "TPSAB1", "TPSB2"),
            cols = c('grey', 'red')) # Mast cells, # cluster 5 = Immune cells

# Find markers of cluster 6
cluster6.markers <- FindMarkers(MergedNML, ident.1 = 6, min.pct = 0.5,
                                only.pos = TRUE, logfc.threshold = 0.25)

FeaturePlot(MergedNML, features=c("CCL21", "NNMT", "TFF3", "MMRN1"),
            cols = c('grey', 'red')) # cluster 6 = Lymphatic Endothelial Cells

MergedNML <- RenameIdents(MergedNML, 
                              `0` = "Immune", `3` = "Immune", `5` = "Immune", 
                              `1` = "Endothelial", `6` = "Lymphatic_Endo",  
                              `2` = "Epithelial",     
                              `4` = "Mensenchymal")

DimPlot(MergedNML, reduction = "umap", label = TRUE)

# Chose your favorite color for the cell clusters, 
# https://htmlcolorcodes.com/color-names/

P1 <- DimPlot(MergedNML, reduction = "umap", label = TRUE, label.size = 3, cols =
                c("Red", "Green", "DeepPink", "Blue"," Orange", "DeepSkyBlue"))

P2 <- DimPlot(MergedNML, reduction = "umap", label = TRUE, label.size = 3)

P1+P2
