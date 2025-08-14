#!/usr/bin/env Rscript
#' Analyze Key Missing Datasets (GSE57338 and GSE115574)
#' 
#' Focused analysis on the two most important missing datasets

# Load required functions
source("functions/data_processing.R")
source("functions/analysis.R") 
source("functions/camk_definitions.R")
source("functions/common_utils.R")

# Load required libraries
suppressPackageStartupMessages({
  library(limma)
  library(tidyverse)
})

cat("=== Analyzing Key Missing Datasets ===\n\n")

# Get CAMK gene definitions
camk_core_genes <- get_camk_gene_categories()$core
cat("Core CAMK genes:", paste(camk_core_genes, collapse = ", "), "\n\n")

# Function to analyze GSE57338 specifically
analyze_gse57338 <- function() {
  cat("=== Analyzing GSE57338 (Primary Dataset) ===\n")
  
  # Load from comprehensive cache
  dataset <- readRDS("cache/comprehensive/GSE57338_processed.rds")
  
  if (!dataset$success || is.null(dataset$expression_matrix)) {
    cat("ERROR: GSE57338 dataset failed\n")
    return(NULL)
  }
  
  # The existing GSE57338 DGE file exists but uses probe IDs
  # Let's work with the processed dataset and map genes properly
  
  expr_matrix <- dataset$expression_matrix
  cat("Dataset:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
  
  # Simple manual mapping for the few genes we have
  # Based on the GSE57338 dataset probe mapping
  probe_to_camk <- list(
    # These are example mappings - we'd need to check actual probe IDs in the data
    "CAMK2D" = "CAMK2D",  # If gene symbols are already present
    "CAMK2A" = "CAMK2A",
    "CAMK2B" = "CAMK2B"
  )
  
  # Check actual gene names in the dataset
  actual_genes <- rownames(expr_matrix)
  cat("First 10 gene names:", paste(head(actual_genes, 10), collapse = ", "), "\n")
  
  # For now, let's use the existing GSE57338 DGE file but try to extract CAMK genes
  # by using the comprehensive mapping approach
  
  # Use existing DGE file - it has probe-level results
  existing_dge <- read.csv("output/GSE57338_healthy_vs_disease_DGE.csv", stringsAsFactors = FALSE)
  cat("Existing DGE results:", nrow(existing_dge), "probes\n")
  
  # The probe IDs are in the first column (row names)
  probe_ids <- existing_dge[,1]
  
  # Apply gene mapping to the probe IDs
  mapped_genes <- map_probe_ids_to_genes(probe_ids, platform = "GPL570")
  
  # Find CAMK genes in the mapped results
  camk_indices <- which(mapped_genes$gene_symbol %in% camk_core_genes)
  
  if (length(camk_indices) > 0) {
    cat("Found", length(camk_indices), "CAMK genes in GSE57338\n")
    
    # Extract CAMK results
    camk_results <- existing_dge[camk_indices, ]
    camk_results$Gene_Symbol <- mapped_genes$gene_symbol[camk_indices]
    camk_results$Dataset <- "GSE57338"
    camk_results$Significant <- camk_results$adj.P.Val < 0.05
    camk_results$Regulation <- ifelse(camk_results$logFC > 0, "UP in Disease", "DOWN in Disease")
    
    # Remove the row names column
    if ("X" %in% names(camk_results)) {
      camk_results$X <- NULL
    } else if (names(camk_results)[1] == "") {
      camk_results <- camk_results[, -1]
    }
    
    cat("CAMK genes found in GSE57338:\n")
    for (i in 1:nrow(camk_results)) {
      gene <- camk_results$Gene_Symbol[i]
      logfc <- round(camk_results$logFC[i], 4)
      pval <- camk_results$P.Value[i]
      adj_pval <- camk_results$adj.P.Val[i]
      sig <- camk_results$Significant[i]
      
      sig_marker <- if (sig) " [SIG]" else " [NS]"
      cat(sprintf("  %-8s: %4s logFC=%7.4f, p=%8.2e, FDR=%8.2e%s\n", 
                 gene, ifelse(logfc > 0, "UP", "DOWN"), logfc, pval, adj_pval, sig_marker))
    }
    
    return(camk_results)
  } else {
    cat("No CAMK genes found in GSE57338 mapping\n")
    return(NULL)
  }
}

# Function to analyze GSE115574 with proper group detection
analyze_gse115574 <- function() {
  cat("\n=== Analyzing GSE115574 ===\n")
  
  # Load dataset
  dataset <- readRDS("cache/microarray/GSE115574_processed.rds")
  
  if (!dataset$success || is.null(dataset$expression_matrix)) {
    cat("ERROR: GSE115574 dataset failed\n")
    return(NULL)
  }
  
  expr_matrix <- dataset$expression_matrix
  cat("Dataset:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
  
  # Apply gene symbol mapping
  cat("Applying gene symbol mapping...\n")
  expr_matrix_mapped <- apply_gene_symbol_mapping(expr_matrix, platform = "GPL570")
  cat("After mapping:", nrow(expr_matrix_mapped), "genes\n")
  
  # Check for CAMK genes
  camk_present <- intersect(rownames(expr_matrix_mapped), camk_core_genes)
  cat("CAMK genes found:", length(camk_present), "/", length(camk_core_genes), "\n")
  cat("CAMK genes:", paste(camk_present, collapse = ", "), "\n")
  
  if (length(camk_present) == 0) {
    return(NULL)
  }
  
  # Manual group detection for GSE115574 (AF vs SR)
  if (!is.null(dataset$phenotype_data)) {
    pheno <- dataset$phenotype_data
    cat("Phenotype columns:", paste(names(pheno), collapse = ", "), "\n")
    
    # Look for AF vs SR pattern
    group_vector <- NULL
    
    # Check common AF/SR columns
    if ("source_name_ch1" %in% names(pheno)) {
      source_names <- pheno$source_name_ch1
      cat("Source names:", paste(unique(source_names), collapse = ", "), "\n")
      
      # Create AF vs SR groups
      af_samples <- grepl("AFib|AF", source_names, ignore.case = TRUE)
      sr_samples <- grepl("SR", source_names, ignore.case = TRUE)
      
      if (any(af_samples) && any(sr_samples)) {
        group_vector <- ifelse(af_samples, "Disease", "Healthy")  # AF as Disease, SR as Healthy
        names(group_vector) <- colnames(expr_matrix_mapped)
        
        cat("Groups assigned:\n")
        cat("  AF (Disease):", sum(group_vector == "Disease"), "samples\n")
        cat("  SR (Healthy):", sum(group_vector == "Healthy"), "samples\n")
      }
    }
    
    if (is.null(group_vector)) {
      cat("Could not detect suitable AF vs SR groups\n")
      return(NULL)
    }
    
    # Perform DGE analysis
    cat("Running DGE analysis (AF vs SR)...\n")
    
    # Create design matrix
    design <- model.matrix(~ group_vector)
    colnames(design) <- c("Intercept", "AF_vs_SR")
    
    # Fit linear model
    fit <- lmFit(expr_matrix_mapped, design)
    fit <- eBayes(fit)
    
    # Get all results
    all_dge_results <- topTable(fit, coef = "AF_vs_SR", number = Inf, adjust.method = "BH")
    
    # Filter for CAMK genes
    camk_dge_results <- all_dge_results[rownames(all_dge_results) %in% camk_present, , drop = FALSE]
    
    if (nrow(camk_dge_results) > 0) {
      # Format results
      camk_dge_results$Gene_Symbol <- rownames(camk_dge_results)
      camk_dge_results$Dataset <- "GSE115574"
      camk_dge_results$Significant <- camk_dge_results$adj.P.Val < 0.05
      camk_dge_results$Regulation <- ifelse(camk_dge_results$logFC > 0, 
                                           "UP in Disease", "DOWN in Disease")
      
      cat("CAMK genes analyzed:", nrow(camk_dge_results), "\n")
      
      # Display results
      for (i in 1:nrow(camk_dge_results)) {
        gene <- camk_dge_results$Gene_Symbol[i]
        logfc <- round(camk_dge_results$logFC[i], 4)
        pval <- camk_dge_results$P.Value[i]
        adj_pval <- camk_dge_results$adj.P.Val[i]
        sig <- camk_dge_results$Significant[i]
        
        sig_marker <- if (sig) " [SIG]" else " [NS]"
        cat(sprintf("  %-8s: %4s logFC=%7.4f, p=%8.2e, FDR=%8.2e%s\n", 
                   gene, ifelse(logfc > 0, "UP", "DOWN"), logfc, pval, adj_pval, sig_marker))
      }
      
      return(camk_dge_results)
    }
  }
  
  return(NULL)
}

# Analyze both key datasets
results_list <- list()

# Try GSE57338
gse57338_results <- analyze_gse57338()
if (!is.null(gse57338_results)) {
  results_list[["GSE57338"]] <- gse57338_results
}

# Try GSE115574
gse115574_results <- analyze_gse115574()
if (!is.null(gse115574_results)) {
  results_list[["GSE115574"]] <- gse115574_results
}

# Combine with existing results
if (length(results_list) > 0) {
  cat("\n=== Integration with Existing Results ===\n")
  
  # Load existing results
  existing_results <- read.csv("output/CAMK_focused_DGE_all_datasets.csv", stringsAsFactors = FALSE)
  cat("Existing results:", nrow(existing_results), "rows from", length(unique(existing_results$Dataset)), "datasets\n")
  
  # Combine new results
  new_results <- do.call(rbind, results_list)
  rownames(new_results) <- NULL
  
  cat("New results:", nrow(new_results), "rows from", length(unique(new_results$Dataset)), "datasets\n")
  
  # Remove any existing entries for these datasets
  datasets_to_replace <- unique(new_results$Dataset)
  existing_filtered <- existing_results[!existing_results$Dataset %in% datasets_to_replace, ]
  
  # Combine
  complete_results <- rbind(existing_filtered, new_results)
  
  cat("Complete combined results:", nrow(complete_results), "rows from", length(unique(complete_results$Dataset)), "datasets\n")
  
  # Save
  output_file <- "output/CAMK_focused_DGE_all_datasets_UPDATED.csv"
  write.csv(complete_results, output_file, row.names = FALSE)
  cat("Updated results saved to:", output_file, "\n")
  
  # CAMK2D summary
  camk2d_all <- complete_results[complete_results$Gene_Symbol == "CAMK2D", ]
  cat("\nCAMK2D across all datasets:\n")
  print(camk2d_all[, c("Dataset", "logFC", "P.Value", "adj.P.Val", "Significant")])
  
} else {
  cat("ERROR: No new datasets could be analyzed\n")
}

cat("\n=== Key Dataset Analysis Complete ===\n")