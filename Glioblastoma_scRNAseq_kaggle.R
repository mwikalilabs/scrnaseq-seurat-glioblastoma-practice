setwd("C:\\Users\\PC\\OneDrive\\Attachments\\Desktop\\Bioinformatics\\practice_projects\\")

# load required libararies
library(Seurat) # scRNAseq analysis toolset
suppressPackageStartupMessages(library(dplyr)) # data manipulation
library(patchwork) # data visualization in grids

# Load the dataset
sc3_data <- Read10X(data.dir = "../practice_projects/human-glioblastoma-dataset/")

#dimension(genes and cells)
cat("Original dimension of the count matrix:", dim(sc3_data))

dim(sc3_data)

#to createseuratobject,  
#filter for genes in at least 3 cells and a sinle cell expresses at least 200genes
sc3 <- CreateSeuratObject(counts = sc3_data, project = "SC3v3", min.cells = 3, min.features = 200)
sc3

#Result:: filtered  1186−1138=48genes and  5697−4433=1264cells(Original dimension of the count matrix was: 1186 5697)

# Look at first three genes in the first thirty cells
sc3_data[1:3, 1:30]

# calculate mitochondrial percentage
sc3[["percent.mt"]] <- PercentageFeatureSet(sc3, pattern = "^MT-")

# check the dataframe with the new column, cells with the most mitochondrial contamination first
# bc this is a filtered dataset, it should be free of contamination
sc3@meta.data %>%
  arrange(desc(percent.mt)) %>%
  head(n = 10)
#visualize the quality of our dataset with Vlnplot()
# custom plot size
# visualize QC metrics
VlnPlot(sc3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

cat("The dataset contains cells with", min(sc3@meta.data$nFeature_RNA), "up to", max(sc3@meta.data$nFeature_RNA), 
    "features (genes) and", min(sc3@meta.data$nCount_RNA), "up to", max(sc3@meta.data$nCount_RNA), "neuro-specific target gene copies.")

#Result::The top 10 most variable genes among glioblastoma cells are: COL4A1 APOD COL4A2 MT1H PLA2G2A CLDN5 MMP9 EGFL7 CXCL10 MARCO

#log-normalization
sc3 <- NormalizeData(sc3, normalization.method = "LogNormalize", scale.factor = 10000)
head(sc3)

# return the 500 most variable features
sc3 <- FindVariableFeatures(sc3, selection.method = "vst", nfeatures = 500)

# Identify the 10 most variable genes
top10_genes <- head(VariableFeatures(sc3), 10)
cat("The top 10 most variable genes among glioblastoma cells are:\n", top10_genes)

#visualite these genes
#VariableFeaturePlot--creates the plot
plot1 = VariableFeaturePlot(sc3)
plot1

#variable features with labels of the top 10 genes
#LabelPoints--labels the variable genes
plot2 = LabelPoints(plot1, points = top10_genes, repel = TRUE, xnudge = 0, ynudge = 0)
plot2

head(sc3)
rownames(sc3)
colnames(sc3)
sc3

# scaling the data for PCA
genes_names <- rownames(sc3)
sc3 <- ScaleData(sc3, features = genes_names, verbose = FALSE)
head(sc3)

#PCA
sc3 <- RunPCA(sc3, features = VariableFeatures(object = sc3),npcs = 50 , verbose = FALSE)

# 5 most important positive and negative expression features in the first 5 principle components:
print(sc3[["pca"]], dims = 1:5, nfeatures = 5)

# plot PCA, first 5 PCs
VizDimLoadings(sc3, dims = 1:6, reduction = "pca", ncol = 3)

# % variance explained by each PC/individual principles components
mat = GetAssayData(sc3, assay = "RNA", layer = "scale.data")
pca = sc3[["pca"]]

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues <- (pca@stdev)^2  ## EigenValues
varExplained <- eigValues / total_variance

# calculate the % variance explained by the first two PCs
cat("The first two principle components explain", 
    scales::label_percent()(round(sum(varExplained[c(1,2)]), 4)), "of the total variance in the dataset.")
# store varEplained in a dataframe alongside the PC number
var_explained_df <- tibble(pc_number = 1:length(varExplained), var_explained = cumsum(varExplained))

# load ggplot2 package
library(ggplot2)

#plot
ggplot(var_explained_df, aes(x = pc_number, y = var_explained)) +
  geom_point(size = 2) +
  geom_line(size = 0.8) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Cumulative Variance Explained by the first 50 PCs",
       x = "PC",
       y = "Cumulative Variance explained") +
  theme_bw() +
  theme(plot.title = element_text(size = 22),
        plot.subtitle = element_text(size = 16),
        axis.text.x= element_text(size = 15),
        axis.text.y= element_text(size = 15), 
        axis.title=element_text(size = 18))

# dimplot, PC1 & PC2
DimPlot(sc3, dims = c(1, 2), reduction = "pca")   

#Results:: suggests that there might be at least three distinct clusters of glioblastoma cells with few cells in between

# heatmap with 500 cells
DimHeatmap(sc3, dims = 1:15, cells = 500, balanced = TRUE)   

# jackstraw method
sc3 <- JackStraw(sc3, num.replicate = 100, prop.freq = 0.05)
sc3 <- ScoreJackStraw(sc3, dims = 1:20)
JackStrawPlot(sc3, dims = 1:20)

# elbow plot method
ElbowPlot(sc3, ndims = 30)

# cell clustering
sc3 <- FindNeighbors(sc3, dims = 1:18)
sc3 <- FindClusters(sc3, algorithm = 1, resolution = 0.2)

# cluster identity of the first 10 cells
head(Idents(sc3), 10)

# PCA plot again with clusters encoded as colors
DimPlot(sc3, reduction = "pca", dims = c(1,2))

#Results::found 6 clusters of glioblastoma cells within our dataset

# find all markers of cluster 3
cluster3_markers <- FindMarkers(sc3, ident.1 = 3, min.pct = 0.3)
head(cluster3_markers, n = 5)

# find all markers distinguishing cluster 0 from cluster 4
cluster0_markers <- FindMarkers(sc3, ident.1 = 0, ident.2 = 4, min.pct = 0.3)
head(cluster0_markers, n = 5)

# find all markers distinguishing cluster 2 from clusters 0 & 4
cluster2_markers <- FindMarkers(sc3, ident.1 = 2, ident.2 = c(0, 4), min.pct = 0.3)
head(cluster2_markers, n = 5)

# tSNE
sc3 <- RunTSNE(sc3, reduction = "pca", dims = 1:18, tsne.method = "Rtsne")
# plot tSNE 
DimPlot(sc3, reduction = "tsne", dims = c(1,2))
