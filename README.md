## scrnaseq-seurat-glioblastoma-practice
Hands-on practice: single-cell RNA-seq workflow in R (Seurat) on a human glioblastoma dataset, following a public kaggle tutorial.

## Credit
This project follows the tutorial **"Seurat single-cell RNA sequencing tutorial: Exploring the Human Glioblastoma scRNAseq dataset"** by Reminho on Kaggle.The workflow and dataset come from that tutorial. I ran it myself to gain practical experience, so this is a learning exercise, not original research.

## Dataset
Human glioblastoma scRNA-seq data from 10X Genomics, as provided in the tutorial: https://www.kaggle.com/code/reminho/basic-scrnaseq-tutorial/notebook#Seurat-single-cell-RNA-sequencing-Tutorial:-Exploring-the-Human-Glioblastoma-scRNAseq-dataset

## What I practiced
1. Loaded the data and created a Seurat object
2. Quality control and cell filtering
3. Normalisation
4. Finding highly variable genes
5. Scaling and PCA
6. Choosing the number of PCs (DimHeatmap, JackStraw, elbow plot)
7. Clustering cells
8. Finding marker genes (differential expression) per cluster
9. Visualisation with tSNE and UMAP

### Repository contents
R code and plots
