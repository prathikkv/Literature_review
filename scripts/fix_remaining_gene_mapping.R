#!/usr/bin/env Rscript
#' Fix Gene Mapping for Remaining Datasets
#' 
#' Apply proper probe-to-gene mapping for GSE115574, GSE41177, GSE79768

# Load functions
source("functions/data_processing.R")
source("functions/camk_definitions.R")

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("=== Fixing Gene Mapping for Remaining Datasets ===\n\n")

# Datasets to fix
datasets_to_fix <- c("GSE115574", "GSE41177", "GSE79768")
camk_genes <- get_camk_gene_categories()$core

for (dataset_id in datasets_to_fix) {
  cat("--- Fixing", dataset_id, "---\n")
  
  # Find dataset file
  dataset_file <- list.files("cache", pattern = paste0(dataset_id, "_processed.rds"), 
                           recursive = TRUE, full.names = TRUE)[1]
  
  if (is.null(dataset_file) || !file.exists(dataset_file)) {
    cat("ERROR: Dataset file not found\n\n")
    next
  }
  
  # Load dataset
  dataset <- readRDS(dataset_file)
  
  if (!dataset$success || is.null(dataset$expression_matrix)) {
    cat("ERROR: Dataset failed or no expression matrix\n\n")
    next
  }
  
  expr_matrix <- dataset$expression_matrix
  cat("Original:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
  
  # Check if gene mapping is needed
  gene_names <- rownames(expr_matrix)
  if (is.null(gene_names)) {
    cat("ERROR: No rownames in expression matrix\n\n")
    next
  }
  
  # Check for CAMK genes before mapping
  camk_before <- intersect(gene_names, camk_genes)
  cat("CAMK genes before mapping:", length(camk_before), "/", length(camk_genes), "\n")
  
  if (length(camk_before) == length(camk_genes)) {
    cat("All CAMK genes already present, no mapping needed\n\n")
    next
  }
  
  # Apply gene symbol mapping
  cat("Applying gene symbol mapping...\n")
  
  tryCatch({
    # Use the existing function to apply gene mapping
    mapped_matrix <- apply_gene_symbol_mapping(expr_matrix, platform = "GPL570")
    
    cat("After mapping:", nrow(mapped_matrix), "genes x", ncol(mapped_matrix), "samples\n")
    
    # Check CAMK genes after mapping
    camk_after <- intersect(rownames(mapped_matrix), camk_genes)
    cat("CAMK genes after mapping:", length(camk_after), "/", length(camk_genes), "\n")
    cat("CAMK genes found:", paste(camk_after, collapse = ", "), "\n")
    
    if (length(camk_after) > length(camk_before)) {
      # Update the dataset with mapped matrix
      dataset$expression_matrix <- mapped_matrix
      dataset$n_genes <- nrow(mapped_matrix)
      
      # Save the updated dataset
      output_file <- paste0(dirname(dataset_file), "/", dataset_id, "_processed_FIXED.rds")
      saveRDS(dataset, output_file)
      
      cat("SUCCESS: Fixed", dataset_id, "saved to:", output_file, "\n")
      cat("- Genes mapped from", nrow(expr_matrix), "to", nrow(mapped_matrix), "\n")
      cat("- CAMK genes improved from", length(camk_before), "to", length(camk_after), "\n")
      
      # Replace original file
      file.rename(output_file, dataset_file)
      cat("- Original file updated\n")
      
    } else {
      cat("WARNING: No improvement in CAMK gene detection after mapping\n")
    }
    
  }, error = function(e) {
    cat("ERROR in gene mapping:", e$message, "\n")
  })
  
  cat("\n")
}

cat("=== Gene Mapping Fix Complete for All Datasets ===\n")