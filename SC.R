library(dplyr)
library(Seurat)
library(patchwork)

#### Load the PBMC dataset ####
# Load the PBMC dataset from 10X genomics (Cell Ranger Matrix) data
pbmc.data <- Read10X(data.dir = "/Users/Zen/Downloads/filtered_gene_bc_matrices/hg19")
# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k",
                           min.cells = 3, min.features = 200)
pbmc #check object class

#### Preprocessing ####

# QC

# Calculate the percentage of mitochondrial genes per cell
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(pbmc@meta.data, 5)

# Visualize QC metrics as a violin plot and a scatter plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter cells based on number of genes expressed and mitochondrial content
# Example here uses simple cutoffs on each, but you can use any stats you like
# Cells that have unique feature counts over 2,500 or less than 200 +
# >5% mitochondrial counts are filtered
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Normalize data
# By default it is using log-normalization with a scale factor of 10,000
# They are stored in "pbmc[["RNA"]]$data"
# ??Normalize the data to 10,000 unique molecular identifiers (UMIs) per cell, so that counts become comparable among cells.
pbmc <- NormalizeData(pbmc)

## Feature Selection
# Identify variable features (feature selection)
# Seurat directly models the mean-variance relationship inherent in SC data
# By default, the 2,000 most highly variable features are identified
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

## Scaling
# Scale the data
# Scale the data to remove unwanted sources of variation, such as differences in sequencing depth
# ScaleData normalizes the gene expression measurements for each cell by the total expression, multiplies by a scale factor (default is 10,000), and log-transforms the result.
# Also prevents domination of highly expressed genes in downstream analyses
# The scaled data is stored in the "RNA" assay slot of the Seurat object, "pbmc[["RNA"]]$scale.data"
# SCTransform() normalization is better in removing unwanted sources of variation
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#### Principal Component Analysis (PCA) ####
# Run PCA  note: ALWAYS SCALE DATA BEFORE RUNNING PCA
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
# No 1. prints out the first 5 principal components with 5 features
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
# No 2. plot the first 2 principal components
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
# No 3. plot a scatterplot for first 2 principal components (the "PCA plot")
DimPlot(pbmc, reduction = "pca") + NoLegend()
# No 4. plot a heatmap from 500 SC that is:
# A. contributing to the first principal component
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
# B. contributing to all principal component from 1 to 15 (plots 15 heatmaps)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

# Determine the "dimentiality"
# The elbow plot is a plot of variance% explained by each principal component
ElbowPlot(pbmc)
# When the SD flattens (elbow), those are not selected for use (cutoff)
# In this example PC7-12 is the cutoff point
# However, more PC used will not dramatically affect results after cutoff
# BUT using too few PC will adversely affects results
# Highest SD is the most important and can be used in GSEA

#### Clustering ####

# Find neighbours
# This use graph based clustering algorithm utilizing KNN
# Then it is refined by Jaccard similarity
pbmc <- FindNeighbors(pbmc, dims = 1:10)

# Find clusters
# This use Louvain algorithm to find clusters
# Alternative option is SLM (smart local moving) algorithm
# Resolution parameter sets the ‘granularity’ of the downstream clustering
# Increased resolution values leads to a greater number of clusters
# For 3k cells: 0.4-1.2 is a good range, increase for larger datasets
pbmc <- FindClusters(pbmc, resolution = 0.5)

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

#### Dimensionality Reduction ####

# Two types: tSNE and UMAP
# Good for visualization but don't draw conclusions just from them
# for UMAP run:
pbmc <- RunUMAP(pbmc, dims = 1:10)
# for tSNE run:
pbmc <- RunTSNE(pbmc, dims = 1:10)

# Plot UMAP
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(pbmc, reduction = "tsne")
# Note: for UMAP, use `reduction = "umap"`, and tsne for tSNE
# Note: you can also use the `patchwork` package to combine plots[?]

# Run code below to save the output file if needed
#saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")

#### Biomarker Identification (DE) ####
# Finds +/- markers for each cluster compared to all other cells
# specify ident.2 when comparing clusters
# FindAllMarkers automatically tests different clusters
# Means, you can compares cluster groups or against cells

# find all markers of cluster 2
cluster2.markers <- FindMarkers(pbmc, ident.1 = 2)
head(cluster2.markers, n = 5)

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3))
head(cluster5.markers, n = 5)

# find markers for every cluster compared to all remaining cells,
# report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# Lets talk about DE testing!
# The default test is the Wilcoxon Rank Sum test (wilcox), but you can also use:
# 1. `wilcox_limma` 2. `bimod` 3. `roc` 4. `t` 5. `LR` 6. `MAST` 7. `DESeq2`
# Use only for UMI-based datasets: 1. `negbinom` 2. `poisson`
# P.S. you have Presto! It will run wilcox using presto as default
# Not: `roc` is the most sensitive test for small datasets
cluster0.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# Visualize the expression of a marker gene
# Methods: FeaturePlot, VlnPlot, RidgePlot, DotPlot, CellScatter, and DoHeatmap

# VlnPlot shows expression probability distributions across clusters
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

# you can plot raw counts as well
VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

# FeaturePlot visualizes feature expression on a tSNE or PCA plot
FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
                               "CD8A"))

# DotPlot is a good way to visualize the expression of a gene across clusters
# Here, map generated from given cells and features for top 20 markers
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#### Cell Identity Assigning ####
# Assign cell identities to clusters

# For this is simple and just need to match markets to known cell types!
# Cluster ID	Markers	        Cell Type
# 0         	IL7R, CCR7	    Naive CD4+ T
# 1	          CD14, LYZ	      CD14+ Mono
# 2	          IL7R, S100A4	  Memory CD4+
# 3	          MS4A1	          B
# 4	          CD8A	          CD8+ T
# 5	          FCGR3A, MS4A7 	FCGR3A+ Mono
# 6	          GNLY, NKG7	    NK
# 7	          FCER1A, CST3	  DC
# 8	          PPBP	          Platelet

new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(pbmc, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
plot

# Save the output plots
#ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
#saveRDS(pbmc, file = "../output/pbmc3k_final.rds")