# RNA-Seq Analysis Shiny App  

This repository contains an interactive R Shiny application for RNA-Seq analysis. The app provides functionality for exploring sample metadata, analyzing gene expression counts matrices, performing differential expression analysis (DEA), and conducting Gene Set Enrichment Analysis (GSEA).  

---

## **Features**
1. **Sample Information Exploration**:  
   - Upload and explore sample metadata.  
   - Summarize dataset structure and visualize distributions by groups.  

2. **Counts Matrix Exploration**:  
   - Upload normalized counts matrices.  
   - Filter genes by variance and non-zero thresholds.  
   - Visualize diagnostic plots, clustered heatmaps, and PCA results.  

3. **Differential Expression Analysis (DEA)**:  
   - Load precomputed DESeq2 results or upload your own.  
   - Explore results interactively with thresholds for significance.  
   - Visualize p-value and log2 fold change distributions, volcano plots, and top gene counts.  

4. **Gene Set Enrichment Analysis (GSEA)**:  
   - Upload fgsea results and explore enriched pathways.  
   - Visualize top pathways with bar plots and scatter plots.  
   - Download filtered enrichment results.  

---

## **Dataset Used**
The application uses RNA-Seq data from post-mortem human brain tissue of Huntington’s Disease patients and neurologically healthy controls.  

- **Metadata**: Prepared from the supplementary files of the dataset.  
- **Normalized Counts**: Provided with the dataset.  
- **DESeq2 Results**: Precomputed differential expression results.  
- **GSEA Input**: Generated using fgsea with DESeq2 results.  

---

## **Installation and Requirements**
### **1. Install R and RStudio**  
Ensure you have the latest versions of [R](https://cran.r-project.org/) and [RStudio](https://rstudio.com/).  


### **2. Install Required Packages**  
Run the following commands to install all necessary libraries:  
```R
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("shiny", "shinydashboard", "DT", "ggplot2", 
                       "fgsea", "msigdbr", "dplyr", "data.table", 
                       "pheatmap", "viridis", "ggbeeswarm"))
