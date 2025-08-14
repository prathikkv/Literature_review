#!/usr/bin/env Rscript
#' Complete CAMK Analysis for ALL Available Datasets
#' 
#' This script runs CAMK analysis on all cached datasets regardless of location

# Load required functions
source("functions/data_processing.R")
source("functions/analysis.R")
source("functions/camk_definitions.R")
source("functions/common_utils.R")
source("scripts/enhanced_group_detection.R")

# Load required libraries
suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
})

cat("=== Complete CAMK Analysis for ALL Datasets ===\n\n")

# Get CAMK gene definitions
camk_gene_categories <- get_camk_gene_categories()
camk_core_genes <- camk_gene_categories$core
cat("Core CAMK genes:", paste(camk_core_genes, collapse = ", "), "\n\n")

# Find ALL processed datasets across all cache subdirectories
all_cache_files <- list.files("cache", pattern = "_processed\\.rds$", 
                             recursive = TRUE, full.names = TRUE)
cat("Found", length(all_cache_files), "processed datasets:\n")
for (file in all_cache_files) {
  dataset_id <- gsub(".*/(GSE[0-9]+)_processed\\.rds", "\\1", file)
  cat(" -", dataset_id, ":", file, "\n")
}
cat("\n")

# Function to analyze a single dataset
analyze_single_dataset <- function(file_path) {
  dataset_id <- gsub(".*/(GSE[0-9]+)_processed\\.rds", "\\1", file_path)
  
  cat("=== Analyzing", dataset_id, "===\n")
  
  tryCatch({
    # Load dataset
    dataset <- readRDS(file_path)
    
    if (!dataset$success || is.null(dataset$expression_matrix)) {
      cat("ERROR: Dataset failed or no expression matrix\n\n")
      return(NULL)
    }
    
    expr_matrix <- dataset$expression_matrix
    cat("Dataset:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
    
    # Check for gene names in rownames
    gene_names <- rownames(expr_matrix)
    if (is.null(gene_names) || length(gene_names) == 0) {
      cat("WARNING: No gene names in rownames - dataset may use probe IDs\n")
      
      # Try to use feature_data for gene mapping
      if (!is.null(dataset$feature_data) && nrow(dataset$feature_data) > 0) {
        cat("Attempting gene mapping from feature data...\n")
        
        # Look for gene symbol columns
        gene_cols <- c("Gene.Symbol", "gene_symbol", "GENE_SYMBOL", 
                      "Gene_Symbol", "Symbol", "symbol")
        gene_col <- NULL
        
        for (col in gene_cols) {
          if (col %in% names(dataset$feature_data)) {
            gene_col <- col
            break
          }
        }
        
        if (!is.null(gene_col)) {
          symbols <- dataset$feature_data[[gene_col]]
          if (length(symbols) == nrow(expr_matrix)) {
            rownames(expr_matrix) <- symbols
            gene_names <- symbols
            cat("SUCCESS: Mapped", length(unique(symbols)), "unique gene symbols\n")
          }
        }
      }
      
      # If still no gene names, skip this dataset
      if (is.null(gene_names) || length(gene_names) == 0) {
        cat("ERROR: Cannot proceed without gene names\n\n")
        return(NULL)
      }
    }
    
    # Check for CAMK genes
    camk_present <- intersect(gene_names, camk_core_genes)
    cat("CAMK genes found:", length(camk_present), "/", length(camk_core_genes), "\n")
    
    if (length(camk_present) == 0) {
      cat("WARNING: No CAMK genes detected - skipping\n\n")
      return(NULL)
    }
    
    cat("CAMK genes present:", paste(camk_present, collapse = ", "), "\n")
    
    # Try to detect groups using enhanced detection
    detected_groups <- enhanced_auto_detect_groups(dataset)
    
    if (is.null(detected_groups)) {
      cat("WARNING: No suitable groups detected\n\n")
      return(NULL)
    }
    
    cat("Groups detected:", detected_groups$pattern_type, "\n")
    cat("  Healthy:", sum(detected_groups$groups == "Healthy"), "samples\n")
    cat("  Disease:", sum(detected_groups$groups == "Disease"), "samples\n")
    
    # Prepare data for DGE
    if (!is.null(detected_groups$sample_indices)) {
      expr_filtered <- expr_matrix[, detected_groups$sample_indices, drop = FALSE]
      groups_vector <- detected_groups$groups
    } else {
      expr_filtered <- expr_matrix
      groups_vector <- detected_groups$groups
    }
    
    cat("Analysis matrix:", nrow(expr_filtered), "genes x", ncol(expr_filtered), "samples\n")
    
    # Perform DGE analysis
    cat("Running genome-wide DGE analysis...\n")
    
    # Create design matrix
    design <- model.matrix(~ groups_vector)
    colnames(design) <- c("Intercept", "Disease_vs_Healthy")
    
    # Fit linear model
    fit <- lmFit(expr_filtered, design)
    fit <- eBayes(fit)
    
    # Get all results
    all_dge_results <- topTable(fit, coef = "Disease_vs_Healthy", 
                               number = Inf, adjust.method = "BH")
    
    cat("Genome-wide DGE completed:", nrow(all_dge_results), "genes analyzed\n")
    
    # Filter for CAMK genes
    camk_dge_results <- all_dge_results[rownames(all_dge_results) %in% camk_core_genes, , drop = FALSE]
    
    if (nrow(camk_dge_results) > 0) {
      # Format results
      camk_dge_results$Gene_Symbol <- rownames(camk_dge_results)
      camk_dge_results$Dataset <- dataset_id
      camk_dge_results$Significant <- camk_dge_results$adj.P.Val < 0.05
      camk_dge_results$Regulation <- ifelse(camk_dge_results$logFC > 0, 
                                           "UP in Disease", "DOWN in Disease")
      
      cat("CAMK genes analyzed:", nrow(camk_dge_results), "\n")
      cat("CAMK significant (FDR < 0.05):", sum(camk_dge_results$Significant), "\n")
      
      # Display results
      for (i in 1:nrow(camk_dge_results)) {
        gene <- camk_dge_results$Gene_Symbol[i]
        logfc <- round(camk_dge_results$logFC[i], 3)
        adj_pval <- camk_dge_results$adj.P.Val[i]
        sig_marker <- if (adj_pval < 0.05) "*" else ""
        direction <- if (logfc > 0) "UP" else "DOWN"
        
        cat(sprintf("  %-8s: %4s logFC=%6.3f, FDR=%8.2e %s\n", 
                   gene, direction, logfc, adj_pval, sig_marker))
      }
      
      cat("SUCCESS: Analysis completed for", dataset_id, "\n\n")
      return(camk_dge_results)
    } else {
      cat("WARNING: No CAMK genes found in DGE results\n\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("ERROR in", dataset_id, ":", e$message, "\n\n")
    return(NULL)
  })
}

# Analyze all datasets
all_results <- list()
successful_datasets <- character(0)

for (file_path in all_cache_files) {
  result <- analyze_single_dataset(file_path)
  if (!is.null(result)) {
    dataset_id <- gsub(".*/(GSE[0-9]+)_processed\\.rds", "\\1", file_path)
    all_results[[dataset_id]] <- result
    successful_datasets <- c(successful_datasets, dataset_id)
  }
}

# Combine all results
if (length(all_results) > 0) {
  cat("=== FINAL RESULTS SUMMARY ===\n")
  cat("Successfully analyzed datasets:", length(all_results), "\n")
  cat("Datasets:", paste(names(all_results), collapse = ", "), "\n")
  
  # Combine into single dataframe
  combined_results <- do.call(rbind, all_results)
  rownames(combined_results) <- NULL
  
  cat("Total gene-dataset combinations:", nrow(combined_results), "\n")
  
  # Save results
  output_file <- "output/CAMK_focused_DGE_all_datasets_FIXED.csv"
  write.csv(combined_results, output_file, row.names = FALSE)
  cat("Results saved to:", output_file, "\n")
  
  # CAMK2D summary
  camk2d_results <- combined_results[combined_results$Gene_Symbol == "CAMK2D", ]
  cat("\nCAMK2D results across", nrow(camk2d_results), "datasets:\n")
  if (nrow(camk2d_results) > 0) {
    for (i in 1:nrow(camk2d_results)) {
      ds <- camk2d_results$Dataset[i]
      logfc <- round(camk2d_results$logFC[i], 3)
      pval <- camk2d_results$adj.P.Val[i]
      sig <- camk2d_results$Significant[i]
      reg <- camk2d_results$Regulation[i]
      
      cat(sprintf("  %-10s: %4s logFC=%6.3f, FDR=%8.2e [%s]\n", 
                 ds, ifelse(logfc > 0, "UP", "DOWN"), logfc, pval, 
                 ifelse(sig, "SIG", "NS")))
    }
  }
  
} else {
  cat("ERROR: No datasets could be successfully analyzed\n")
}

cat("\n=== Analysis Complete ===\n")