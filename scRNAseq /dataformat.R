# (Reference) 
# https://github.com/kpatel427/YouTubeTutorials/blob/main/singleCell_standard_workflow.R
# https://www.youtube.com/watch?v=3xcTpqQzUwQ&t=154s

setwd("~/Desktop/Coding_practice/scRNAseq")

# Install BiocManager
install.packages("BiocManager")

# Install the seurat package
BiocManager::install("Seurat")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

#load libraries
library(dplyr)
library(tidyverse)
library(Seurat)
library(SeuratDisk)

# script to demonstrate reading single cell matrices in various format 
# and converting to seurat object

# .RDS format
## https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/QB5CC8

rds_obj <- readRDS('Data/ependymal_cells.rds')

# 10X CellRanger .HDF5 format 
hdf5_obj <- Read10X_h5(filename = "Data/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)
seurat_hdf5 <- CreateSeuratObject(counts = hdf5_obj)

# .mtx file (Unfortunately, I could not find the mtx files)
mtx_obj <- ReadMtx(mtx = "raw_feature_bc_matrix/matrix.mtx.gz",
                   features = "raw_feature_bc_matrix/features.tsv.gz",
                   cells = "raw_feature_bc_matrix/barcodes.tsv.gz")
seurat_mtx <- CreateSeuratObject(counts = mtx_obj)

# .loom files (Unfortunately, I could not find the mtx files)
loom_oj <- Connect(filename = "adult-hem-organs-10X-bone-marrow.loom", mode = 'r')
seurat_loom <- as.Seurat(loom_oj)

# .h5ad format 
# step 1: convert AnnData object to an h5Seurat file
Convert("Data/adata_SS2_for_download.h5ad", dest = "h5seurat", overwrite = TRUE)

# step 2: Load h5Seurat file into a Seurat object
seurat_anndata <- LoadH5Seurat("Data/adata_SS2_for_download.h5seurat")
