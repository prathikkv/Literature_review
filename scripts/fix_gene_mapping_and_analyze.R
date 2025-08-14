#!/usr/bin/env Rscript
#' Fix Gene Mapping and Analyze All Datasets
#' 
#' This script applies proper gene symbol mapping to all cached datasets
#' and performs comprehensive CAMK analysis

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

cat("=== Fix Gene Mapping and Analyze All Datasets ===\n\n")

# Get CAMK gene definitions
camk_gene_categories <- get_camk_gene_categories()
camk_core_genes <- camk_gene_categories$core
cat("Core CAMK genes:", paste(camk_core_genes, collapse = ", "), "\n\n")

# Find all unique datasets (avoid duplicates)
all_cache_files <- list.files("cache", pattern = "_processed\\.rds$", 
                             recursive = TRUE, full.names = TRUE)

# Get unique datasets (some datasets appear in multiple cache locations)
unique_datasets <- unique(gsub(".*/(GSE[0-9]+)_processed\\.rds", "\\1", all_cache_files))
cat("Unique datasets found:", length(unique_datasets), "\n")
cat("Datasets:", paste(unique_datasets, collapse = ", "), "\n\n")

# For each unique dataset, find the best cache file (prefer comprehensive, then microarray)
dataset_files <- character(length(unique_datasets))
names(dataset_files) <- unique_datasets

for (dataset_id in unique_datasets) {
  # Find all files for this dataset
  matching_files <- all_cache_files[grepl(paste0(dataset_id, "_processed\\.rds"), all_cache_files)]
  
  # Prefer comprehensive > microarray > others
  if (any(grepl("/comprehensive/", matching_files))) {
    dataset_files[dataset_id] <- matching_files[grepl("/comprehensive/", matching_files)][1]
  } else if (any(grepl("/microarray/", matching_files))) {
    dataset_files[dataset_id] <- matching_files[grepl("/microarray/", matching_files)][1]
  } else {
    dataset_files[dataset_id] <- matching_files[1]
  }
}

# Function to analyze a single dataset with proper gene mapping
analyze_dataset_with_mapping <- function(file_path, dataset_id) {
  cat("=== Analyzing", dataset_id, "===\n")
  
  tryCatch({
    # Load dataset
    dataset <- readRDS(file_path)
    
    if (!dataset$success || is.null(dataset$expression_matrix)) {
      cat("ERROR: Dataset failed or no expression matrix\n\n")
      return(NULL)
    }
    
    expr_matrix <- dataset$expression_matrix
    cat("Original dataset:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
    
    # Determine platform for gene mapping
    platform <- "GPL570"  # Default
    if (!is.null(dataset$dataset_info$platform)) {
      platform <- dataset$dataset_info$platform
    }
    cat("Platform:", platform, "\n")
    
    # Apply gene symbol mapping if needed
    if (!is.null(rownames(expr_matrix))) {
      # Check if we need mapping (look for probe ID patterns)
      probe_pattern <- "_at$|_s_at$|_x_at$|_a_at$"
      needs_mapping <- any(grepl(probe_pattern, rownames(expr_matrix)[1:min(10, nrow(expr_matrix))]))
      
      if (needs_mapping) {
        cat("Applying gene symbol mapping...\n")
        expr_matrix <- apply_gene_symbol_mapping(expr_matrix, platform)
        cat("After mapping:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
      } else {
        cat("Gene symbols already present\n")
      }
    } else {
      cat("ERROR: No rownames in expression matrix\n\n")
      return(NULL)
    }
    
    # Check for CAMK genes after mapping
    gene_names <- rownames(expr_matrix)
    camk_present <- intersect(gene_names, camk_core_genes)
    cat("CAMK genes found:", length(camk_present), "/", length(camk_core_genes), "\n")
    
    if (length(camk_present) == 0) {
      cat("WARNING: No CAMK genes detected after mapping\n\n")
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
      cat("Results:\n")
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

# Analyze all unique datasets
all_results <- list()
successful_datasets <- character(0)

for (dataset_id in unique_datasets) {
  file_path <- dataset_files[dataset_id]
  cat("Processing", dataset_id, "from:", file_path, "\n")
  
  result <- analyze_dataset_with_mapping(file_path, dataset_id)
  if (!is.null(result)) {
    all_results[[dataset_id]] <- result
    successful_datasets <- c(successful_datasets, dataset_id)
  }
}

# Combine all results
if (length(all_results) > 0) {
  cat("=== COMPREHENSIVE RESULTS SUMMARY ===\n")
  cat("Successfully analyzed datasets:", length(all_results), "\n")
  cat("Datasets:", paste(names(all_results), collapse = ", "), "\n")
  
  # Combine into single dataframe
  combined_results <- do.call(rbind, all_results)
  rownames(combined_results) <- NULL
  
  cat("Total gene-dataset combinations:", nrow(combined_results), "\n")
  
  # Save comprehensive results
  output_file <- "output/CAMK_focused_DGE_all_datasets_COMPLETE.csv"
  write.csv(combined_results, output_file, row.names = FALSE)
  cat("Complete results saved to:", output_file, "\n")
  
  # CAMK2D summary across all datasets
  camk2d_results <- combined_results[combined_results$Gene_Symbol == "CAMK2D", ]
  cat("\n=== CAMK2D Results Summary ===\n")
  cat("CAMK2D found in", nrow(camk2d_results), "datasets:\n")
  
  if (nrow(camk2d_results) > 0) {
    for (i in 1:nrow(camk2d_results)) {
      ds <- camk2d_results$Dataset[i]
      logfc <- round(camk2d_results$logFC[i], 4)
      pval <- camk2d_results$P.Value[i]
      adj_pval <- camk2d_results$adj.P.Val[i]
      sig <- camk2d_results$Significant[i]
      reg <- camk2d_results$Regulation[i]
      
      sig_marker <- if (sig) " [SIGNIFICANT]" else " [NS]"
      cat(sprintf("  %-10s: %4s logFC=%7.4f, p=%8.2e, FDR=%8.2e%s\n", 
                 ds, ifelse(logfc > 0, "UP", "DOWN"), logfc, pval, adj_pval, sig_marker))
    }
    
    # Overall CAMK2D statistics
    sig_datasets <- sum(camk2d_results$Significant)
    down_datasets <- sum(camk2d_results$logFC < 0)
    up_datasets <- sum(camk2d_results$logFC > 0)
    
    cat("\nCAMK2D Summary Statistics:\n")
    cat("  Total datasets with CAMK2D:", nrow(camk2d_results), "\n")
    cat("  Significant results:", sig_datasets, "/", nrow(camk2d_results), "\n")
    cat("  Downregulated:", down_datasets, "/", nrow(camk2d_results), "\n")
    cat("  Upregulated:", up_datasets, "/", nrow(camk2d_results), "\n")
  }
  
} else {
  cat("ERROR: No datasets could be successfully analyzed\n")
}

cat("\n=== Gene Mapping and Analysis Complete ===\n")