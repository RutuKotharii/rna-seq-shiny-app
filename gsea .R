# Load required libraries
library(dplyr)
library(tibble)
library(fgsea)
library(msigdbr)
library(readr)

# Function to prepare ranked list from differential expression results
prepare_ranked_list <- function(de_results, ranking_metric = "log2FoldChange") {
  ranked_list <- de_results %>%
    arrange(desc(!!sym(ranking_metric))) %>%
    select(symbol, !!sym(ranking_metric)) %>%
    distinct() %>%
    deframe()
  return(ranked_list)
}

# Function to load Hallmark gene sets using msigdbr
load_hallmark_gene_sets <- function() {
  gene_set_data <- msigdbr(species = "Homo sapiens", category = "H")
  gene_set_db <- split(gene_set_data$gene_symbol, gene_set_data$gs_name)
  return(gene_set_db)
}

# Function to run FGSEA
run_fgsea_analysis <- function(gene_set_db, ranked_list) {
  fgsea_results <- fgsea(
    pathways = gene_set_db,
    stats = ranked_list,
    minSize = 15,
    maxSize = 500,
    nperm = 1000
  ) %>%
    as_tibble() %>%
    arrange(padj)
  return(fgsea_results)
}

# Script Execution
main <- function() {
  # Load differential expression results
  # Replace with the actual file path to your DE results
  de_results_file <- "DESeq2_results.csv"
  de_results <- read_csv(de_results_file)
  
  # Prepare the ranked list
  ranked_list <- prepare_ranked_list(de_results, ranking_metric = "log2FoldChange")
  
  # Load Hallmark gene sets
  gene_set_db <- load_hallmark_gene_sets()
  
  # Run FGSEA
  fgsea_results <- run_fgsea_analysis(gene_set_db, ranked_list)
  fgsea_results <- as.data.frame(fgsea_results)
  # Save FGSEA results to file
  output_file <- "fgsea_results.csv"
  write_csv(fgsea_results, output_file)
  
  message("FGSEA results saved to: ", output_file)
}

# Run the script
main()
