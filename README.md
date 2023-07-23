# Single-Cell Analysis Repo - README

This repository contains two comprehensive R scripts that perform single-cell analysis, including both single-cell RNA sequencing (scRNAseq) and integrative analysis of the single-cell multiomic (RNA-ATAC scMultiome) data.

## Project 1: scRNAseq Pipeline

A detailed pipeline to perform single-cell RNA sequencing data analysis using the `Seurat` package. This pipeline is equipped to handle data loading and pre-processing, quality control, normalization, feature selection, data scaling, linear and non-linear dimensionality reduction, clustering, and cell type annotation. It also generates high-quality visualizations for better understanding of the data.

**Key Features**
- Load and pre-process scRNAseq data.
- Perform quality control and remove outliers.
- Normalize data and identify highly variable features.
- Scale data and regress out unwanted sources of variation.
- Perform linear (PCA) and non-linear (t-SNE, UMAP) dimensionality reduction.
- Cluster cells based on their gene expression profiles.
- Annotate cell types based on cluster-specific marker genes.
- Generate high-quality visualizations.

## Project 2: sc Multiomic Analysis

This script provides a comprehensive approach to analyze both mono- and bi-modal data from a single-cell multiomic analysis of a 35-week-old human retinal organoid. It employs the Seurat package and uses RNA-ATAC scMultiome data. The pipeline identifies cell type heterogeneity and putative cis- and trans-regulatory elements.

**Key Features**
- Mono-modal integrative analysis of the RNA-ATAC scMultiome data.
- Bi-modal integrative analysis of the RNA-ATAC scMultiome data.
- Gene regulatory network reconstruction.

## Requirements

To run these pipelines, make sure you have the following R packages installed:

- Seurat
- patchwork
- tidyverse
- ggplot2
- Matrix
- dplyr
- svglite

Ensure that the input data paths are correctly set and the necessary R libraries are installed to successfully run the scripts.

## Usage

To use the pipelines, execute the scripts in R or RStudio.

## Data 

While specific datasets are used in the scripts (e.g., scRNAseq dataset derived from telencephalic sample and choroid plexus (ChP) samples, and a 35-week-old human retinal organoid for sc Multiomic analysis), these are not included in the repo. Users need to replace the data loading steps with paths to their own data. 

## Contributions 

The scRNAseq pipeline was written with the help of a tutorial by Zhisong He of the Barbara Treutlein lab.
