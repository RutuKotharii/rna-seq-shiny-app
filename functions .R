# Installing and loading required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Installing Bioconductor packages if not already installed
BiocManager::install(c("fgsea", "msigdbr", "dplyr", "ggplot2"))

# Loading necessary libraries
library(DT)
library(ggplot2)
library(pheatmap)
library(stats)
library(reshape2)
library(dplyr)
library(tibble)
library(viridis)
library(scales)
library(ggbeeswarm)

### PART 1: Handling Sample Information and Data Summarization ###

# Loading sample data from a CSV file
# Parameters:
# - file_path: Specifying the path to the CSV file
# Returning: Loaded data as a data frame
load_sample_data <- function(file_path) {
  data <- read.csv(file_path, stringsAsFactors = TRUE)
  
  # Ensuring Vonsattel_Grade is treated as a factor
  if ("Vonsattel_Grade" %in% colnames(data)) {
    data$Vonsattel_Grade <- as.factor(data$Vonsattel_Grade)
  }
  
  return(data)
}

# Generating summary statistics for each column in the dataset
# Parameters:
# - data: Specifying the data frame to summarize
# Returning: A data frame summarizing column types and statistics
generate_summary <- function(data) {
  data.frame(
    Column_Name = colnames(data),
    Type = sapply(data, class),
    Summary = sapply(data, function(col) {
      if (is.numeric(col)) {
        # Calculating mean and SD grouped by Diagnosis
        group_stats <- tapply(col, data$Diagnosis, function(x) {
          mean_val <- round(mean(x, na.rm = TRUE), 2)
          sd_val <- round(sd(x, na.rm = TRUE), 2)
          paste("Mean:", mean_val, "SD:", sd_val)
        })
        
        # Creating a formatted string that lists mean and SD for each group
        summary_text <- sapply(names(group_stats), function(group) {
          # Extracting and formatting the mean and standard deviation
          paste(group, ":", group_stats[[group]])
        })
        
        # Combining all formatted statistics into one string
        return(paste(summary_text, collapse = "; "))
      } else if (is.factor(col) || is.character(col)) {
        # For categorical columns, listing unique values
        unique_values <- unique(col[!is.na(col)])
        return(paste(unique_values, collapse = ", "))
      } else {
        # For other column types, counting non-missing values
        return(paste0("Non-missing values: ", sum(!is.na(col))))
      }
    }),
    row.names = NULL
  )
}

# Creating a density plot for a specified column
# Parameters:
# - data: Specifying the data frame
# - x_col: Selecting the column to plot on the x-axis
# - group_col: Selecting the column used to group data (fill color)
# Returning: A ggplot2 density plot
create_density_plot <- function(data, x_col, group_col) {
  ggplot(data, aes_string(x = x_col, fill = group_col)) +
    geom_density(alpha = 0.7) +
    labs(title = paste("Distribution of", x_col, "by", group_col),
         x = x_col,
         y = "Density",
         fill = group_col) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, face = "bold", color = "black"),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.grid.major = element_blank(),   # Remove grid lines
      panel.grid.minor = element_blank(),   # Remove minor grid lines
      panel.border = element_rect(color = "black", fill = NA, size = 1), # Add border around the plot
      axis.line = element_line(color = "black") # Add black axis lines
    )
}


### PART 2: Processing Counts Matrix ###

# Loading a counts matrix from a CSV file
# Parameters:
# - file_path: Specifying the path to the CSV file
# Returning: The counts matrix as a data frame
load_counts_matrix <- function(file_path) {
  data <- read.csv(file_path, stringsAsFactors = FALSE, row.names = 1)
  return(as.data.frame(data))
}

# Filtering the counts matrix based on variance and non-zero counts
# Parameters:
# - counts: Specifying the counts matrix
# - var_percentile: Setting the percentile threshold for variance filtering
# - min_nonzero: Setting the minimum number of non-zero values required per gene
# Returning: A list with filtered counts matrix and summary statistics
filter_counts <- function(counts, var_percentile, min_nonzero) {
  # Calculating variance, median, and nonzero counts for each gene
  gene_variance <- apply(counts, 1, var, na.rm = TRUE)
  gene_nonzero <- apply(counts, 1, function(x) sum(x > 0))
  
  # Determining the variance threshold
  var_threshold <- quantile(gene_variance, var_percentile)
  
  # Identifying filtered genes based on variance and nonzero counts
  filtered_genes <- (gene_variance >= var_threshold) & (gene_nonzero >= min_nonzero)
  
  # Filtering genes based on criteria
  filtered_counts <- counts[filtered_genes, ]
  
  # Number of samples (columns) in the counts matrix
  num_samples <- ncol(counts)
  
  # Creating a summary statistics data frame
  stats <- data.frame(
    Metric = c("Total Genes", "Filtered Genes", "Filtered Percentage",
               "Unfiltered Genes", "Unfiltered Percentage", "Number of Samples"),
    Value = c(
      nrow(counts),                              # Total number of genes
      sum(filtered_genes),                       # Number of filtered genes
      round(sum(filtered_genes) / nrow(counts) * 100, 2), # Filtered percentage
      nrow(counts) - sum(filtered_genes),        # Number of unfiltered genes
      round((nrow(counts) - sum(filtered_genes)) / nrow(counts) * 100, 2), # Unfiltered percentage
      num_samples                                # Number of samples
    )
  )
  
  # Returning filtered counts and summary statistics
  return(list(
    filtered_counts = filtered_counts,
    stats = stats
  ))
}

# Function to create diagnostic plots dynamically
create_diagnostic_plot <- function(counts, x_col, y_col, filtered_indices, log_transform, num_zeros) {
  # Computing metrics
  gene_median <- apply(counts, 1, median, na.rm = TRUE)
  gene_variance <- apply(counts, 1, var, na.rm = TRUE)
  
  # Preparing data frame for plotting
  data <- data.frame(
    Median = gene_median,
    Variance = gene_variance,
    Zeros = num_zeros,  
    Filtered = filtered_indices
  )
  
  # Generating plot dynamically based on selected metrics
  plot <- ggplot(data, aes_string(x = x_col, y = y_col, color = "Filtered")) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "lightgrey")) +
    labs(
      title = paste(y_col, "vs", x_col),
      x = if (log_transform) paste0(x_col, " (log scale)") else x_col,
      y = if (log_transform) paste0(y_col, " (log scale)") else y_col,
      color = "Passed"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, size = 1), # Add border around the plot
      axis.line = element_line(color = "black") # Add black axis lines
    )
  
  # Applying log-transform to axes 
  if (log_transform) {
    plot <- plot +
      scale_x_continuous(trans = 'log10') +
      scale_y_continuous(trans = 'log10')
  }
  
  return(plot)
}

# Function to generate a heatmap with customization options
generate_heatmap <- function(filtered_counts, top_genes, log_transform, color_palette) {
  # Calculating variance for top genes selection
  gene_variance <- apply(filtered_counts, 1, var, na.rm = TRUE)
  top_genes_indices <- order(gene_variance, decreasing = TRUE)[1:top_genes]
  selected_counts <- filtered_counts[top_genes_indices, ]
  
  # Log transform if enabled
  if (log_transform) {
    selected_counts <- log10(selected_counts + 1)
  }
  
  # Choosing color palette
  palette <- switch(color_palette,
                    "Viridis" = viridis::viridis(100),
                    "Blue-Red" = colorRampPalette(c("blue", "white", "red"))(100),
                    "Heat" = heat.colors(100),
                    "Coolwarm" = colorRampPalette(c("blue", "white", "orange", "red"))(100))
  
  # Generating heatmap
  pheatmap(
    selected_counts,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = palette,
    main = "Clustered Heatmap of Filtered Counts",
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = TRUE
  )
}

# Function to perform PCA on the counts matrix
perform_pca <- function(counts, top_genes = NULL) {
  # Optionally filter to the top most variable genes
  if (!is.null(top_genes)) {
    gene_variance <- apply(counts, 1, var, na.rm = TRUE)
    top_genes_indices <- order(gene_variance, decreasing = TRUE)[1:top_genes]
    counts <- counts[top_genes_indices, ]
  }
  
  # Performing PCA on the transposed counts matrix (samples as rows, genes as columns)
  pca <- prcomp(t(counts), scale. = TRUE)
  
  # Calculating the variance explained by each principal component
  explained_variance <- (pca$sdev^2 / sum(pca$sdev^2)) * 100
  
  # Returning PCA object and explained variance
  return(list(pca = pca, explained_variance = explained_variance))
}

# Function to create PCA scatter plots
create_pca_scatter <- function(pca_result, pc_x, pc_y, explained_variance) {
  # Extract PCA components and format into a data frame
  pca_data <- data.frame(pca_result$x)
  pca_data$Sample <- rownames(pca_data)
  
  # Generating scatter plot
  plot <- ggplot(pca_data, aes_string(x = pc_x, y = pc_y, label = "Sample")) +
    geom_point(size = 3, color = "darkblue", alpha = 0.7) +
    labs(
      title = paste("PCA Plot:", pc_x, "vs", pc_y),
      x = paste0(pc_x, " (", round(explained_variance[as.numeric(gsub("PC", "", pc_x))], 2), "% variance)"),
      y = paste0(pc_y, " (", round(explained_variance[as.numeric(gsub("PC", "", pc_y))], 2), "% variance)")
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, size = 1), # Add border around the plot
      axis.line = element_line(color = "black") # Add black axis lines
    )
  
  return(plot)
}

# Function to create PCA beeswarm plot for multiple components
create_pca_beeswarm_multiple <- function(pca_result, components, explained_variance) {
  # Extracting PCA components and format into a data frame
  pca_data <- data.frame(pca_result$x)
  pca_data$Sample <- rownames(pca_data)
  
  # Reshaping data for beeswarm plotting
  pca_long <- tidyr::pivot_longer(
    pca_data,
    cols = components,
    names_to = "Component",
    values_to = "Value"
  )
  
  # Generating beeswarm plot
  plot <- ggplot(pca_long, aes(x = Component, y = Value, color = Component)) +
    geom_beeswarm(size = 3, alpha = 0.7) +
    labs(
      title = "Beeswarm Plot for PCA Components",
      x = "Principal Components",
      y = "PCA Values"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 1), # Add border around the plot
      axis.line = element_line(color = "black") # Add black axis lines
    ) +
    scale_color_brewer(palette = "Set1")
  # Adding explained variance as annotation
  for (comp in components) {
    variance <- explained_variance[as.numeric(gsub("PC", "", comp))]
    plot <- plot + annotate(
      "text", x = comp, y = max(pca_long$Value),
      label = paste0(comp, " (", round(variance, 2), "%)"),
      hjust = 0.5, vjust = -0.5, size = 3
    )
  }
  
  return(plot)
}

#3
# Function to preprocess DESeq2 results and add volcano plot status
process_deseq2_results <- function(deseq2_results) {
  deseq2_results <- deseq2_results %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(volc_plot_status = dplyr::case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "UP",
      padj < 0.05 & log2FoldChange < 0 ~ "DOWN",
      TRUE ~ "NS"
    ))
  return(deseq2_results)
}

# Function to plot histogram of unadjusted p-values
plot_pvalue_histogram <- function(deseq2_results) {
  ggplot(deseq2_results, aes(x = pvalue)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Histogram of P-Values",
      x = "P-Value",
      y = "Frequency"
    ) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
}

# Function to plot histogram of log2FoldChange for significant genes
plot_logfc_histogram <- function(deseq2_results, padj_threshold = 0.05) {
  significant_genes <- deseq2_results %>% filter(padj < padj_threshold)
  ggplot(significant_genes, aes(x = log2FoldChange)) +
    geom_histogram(bins = 100, fill = "skyblue", color = "black") +
    theme_minimal(base_size = 14) +
    labs(
      title = "Histogram of Log2 Fold Change (Significant Genes)",
      x = "Log2 Fold Change",
      y = "Frequency"
    ) +
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
}

# Function to plot normalized counts for top genes by adjusted p-value
plot_top_genes_counts <- function(counts_matrix, deseq2_results, top_n = 10) {
  top_genes <- deseq2_results %>%
    arrange(padj) %>%
    head(top_n) %>%
    pull(symbol)
  
  top_gene_counts <- counts_matrix %>%
    rownames_to_column("gene") %>%
    filter(gene %in% top_genes) %>%
    tidyr::gather(sample, count, -gene)
  
  ggplot(top_gene_counts, aes(x = sample, y = count, color = gene)) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste("Normalized Counts for Top", top_n, "Genes"),
      x = "Sample",
      y = "Normalized Count"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, size = 1),
      legend.title = element_blank()
    )
}

# Creating a volcano plot for visualizing differential expression results
# Parameters:
# - de_results: Specifying the data frame with differential expression results
# - pval_threshold: Setting the p-value threshold for significance
# - logfc_threshold: Setting the log2 fold-change threshold for significance
# Returning: A ggplot2 volcano plot
plot_volcano <- function(deseq2_results) {
  ggplot(deseq2_results, aes(x = log2FoldChange, y = -log10(padj), color = volc_plot_status)) +
    geom_point(alpha = 0.5, size = 2) +
    theme_minimal() +  
    labs(x = "Log2 Fold Change", y = "-Log10 (padj)", 
         title = "Volcano Plot of DESeq2 Differential Expression Results", 
         color = "Significance") +  # Label axes and title
    scale_color_manual(values = c("red", "green", "blue")) +  # Custom colors for status categories
    theme(legend.position = "right", panel.border = element_rect(color = "black", fill = NA, size = 1) 
                                  )
}

### PART 4: Gene Set Enrichment Analysis (GSEA) ###

# Function to create a barplot of top pathways by NES
plot_top_pathways_by_nes <- function(fgsea_results, top_n) {
  fgsea_results <- fgsea_results %>% 
    arrange(desc(abs(NES))) %>%  # Sorting by absolute NES
    head(top_n)
  
  ggplot(fgsea_results, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
    geom_bar(stat = "identity", color = "black", width = 0.8) +
    coord_flip() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, name = "NES") +
    theme_minimal(base_size = 14, base_family = "Arial") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      legend.position = "right"
    ) +
    labs(
      x = "Pathway",
      y = "Normalized Enrichment Score (NES)",
      title = paste("Top", top_n, "Pathways by NES")
    )
}

# Function to create a scatter plot of pathways
plot_scatter <- function(fgsea_results, padj_threshold) {
  fgsea_results <- fgsea_results %>% 
    mutate(
      color_group = case_when(
        padj > padj_threshold ~ "Below Threshold",
        padj <= padj_threshold ~ "Above Threshold"
      )
    )
  
  ggplot(fgsea_results, aes(x = NES, y = -log10(padj), color = color_group)) +
    geom_point(size = 3, alpha = 0.7) +
    scale_color_manual(values = c("Below Threshold" = "grey", "Above Threshold" = "red")) +
    theme_minimal(base_size = 14, base_family = "Arial") +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      legend.title = element_blank(),
      legend.position = "top"
    ) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = "-log10(Adjusted p-value)",
      title = "Scatter Plot of Pathways"
    ) +
    geom_hline(yintercept = -log10(padj_threshold), linetype = "dashed", color = "blue", linewidth = 0.8) +
    annotate("text", x = max(fgsea_results$NES) * 0.8, y = -log10(padj_threshold) + 0.2,
             label = paste("Threshold: p-value <", padj_threshold), color = "blue", size = 4)
}

# Function to filter FGSEA results by p-value and NES type
filter_fgsea_results <- function(fgsea_results, padj_filter, nes_filter) {
  filtered <- fgsea_results %>%
    filter(padj <= padj_filter)
  
  if (nes_filter == "pos") {
    filtered <- filtered %>% filter(NES > 0)
  } else if (nes_filter == "neg") {
    filtered <- filtered %>% filter(NES < 0)
  }
  
  return(filtered)
}

# Saving filtered FGSEA results to a file
save_filtered_results <- function(filtered_results, file) {
  write.csv(filtered_results, file, row.names = FALSE)
}