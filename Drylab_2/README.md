# Single-Cell Multiomic Analysis - README

This R script provides an integrative approach to analyze mono- and bi-modal data from the single-cell multiomic analysis of a 35-week-old human retinal organoid. It employs the Seurat package and uses RNA-ATAC scMultiome data. The study aims to identify cell type heterogeneity and putative cis- and trans-regulatory elements critical for cell type/state identity establishment, maintenance, and transition.

## Section 1: Mono-modal Integrative Analysis

- **Step 1:** Data loading and creation of the Seurat object.
- **Step 2:** Quality control of the RNA and ATAC assays.
- **Step 3:** Analysis of the RNA assay which includes RNA analysis pipeline, data normalization, identification of highly variable genes, clustering, and plotting of clusters.
- **Step 4:** Analysis of the ATAC assay including feature selection, normalization, linear and non-linear dimension reduction, and UMAP embedding for visualization.

## Section 2: Bi-modal Integrative Analysis

- **Step 1:** Weighted nearest neighbor analysis using combined RNA and ATAC assays. UMAP with weighted nearest neighbor is used for visualization. 
- **Step 2:** Cell type gene/peak marker identification and visualization of the chromatin. Both RNA and ATAC differential expression and accessibility analysis are conducted.
- **Step 3:** TF binding motif enrichment analysis using JASPAR db to identify TF binding motifs.
- **Step 4:** ChromVAR analysis to estimate motif activity score for each motif on the cell level and differential activity analysis.

## Section 3: Gene Regulatory Network Reconstruction

- **Step 1:** Initialization and fusion of the Pando object.
- **Step 2:** Gene regulatory network (GRN) generation.
- **Step 3:** Visualization and analysis of regulons.
- **Step 4:** Visualization of the GRN.

Note: Some steps have alternative methods provided in the script. Refer to the in-script comments for more information. The script also uses marker genes from Wahle et al. 2023 and takes care of removing irrelevant cell types due to possible cell death in the organoid sample.