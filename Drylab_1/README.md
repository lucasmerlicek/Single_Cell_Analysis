# scRNAseq Pipeline

This repository contains a detailed scRNAseq pipeline implemented in R, leveraging the power of the `Seurat` package for single-cell RNA sequencing data analysis. The pipeline is carefully designed to perform various steps of scRNAseq data analysis, including normalization, quality control, feature selection, data scaling, linear and non-linear dimensionality reduction, clustering, and cell type annotation.

## Features

1. Load and pre-process scRNAseq data
2. Perform quality control and remove outliers
3. Normalize data and identify highly variable features
4. Scale data and regress out unwanted sources of variation
5. Perform linear (PCA) and non-linear (t-SNE, UMAP) dimensionality reduction
6. Cluster cells based on their gene expression profiles
7. Annotate cell types based on cluster-specific marker genes
8. Generate high-quality visualizations such as violin plots, scatterplots, heatmaps, PCA, t-SNE and UMAP plots, etc.

## Requirements

To successfully run this pipeline, you need to have the following R packages installed:

- Seurat
- patchwork
- tidyverse
- ggplot2
- Matrix
- dplyr
- svglite

## Usage

To run the pipeline, just execute the script in R or Rstudio. Ensure that the input data paths are correctly set and the necessary R libraries are installed.

## Data 

The pipeline uses an example scRNAseq dataset derived from telencephalic sample and choroid plexus (ChP) samples. Note that users will need to replace the data loading steps with paths to their own scRNAseq datasets. Raw Data is not included

## Contributions 

This pipeline was writen with the help of a tutorial by Zhisong He of the Barbara Treutlein lab. 
GitHub: https://github.com/quadbio/scRNAseq_analysis_vignette