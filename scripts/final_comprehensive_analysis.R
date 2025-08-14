#!/usr/bin/env Rscript
#' Final Comprehensive CAMK Analysis
#' 
#' Complete DGE and Meta-Analysis on 4 High-Quality Finalized Datasets:
#' - GSE57338 (313 samples, 20K+ genes) - PRIMARY DATASET
#' - GSE115574 (59 samples, 31K+ genes) 
#' - GSE41177 (38 samples, 31K+ genes)
#' - GSE79768 (26 samples, 31K+ genes)
#' Total: 436 samples

# Load required functions and libraries
source("functions/data_processing.R")
source("functions/analysis.R")
source("functions/camk_definitions.R")
source("functions/common_utils.R")
source("scripts/enhanced_group_detection.R")

suppressPackageStartupMessages({
  library(limma)
  library(metafor)
  library(tidyverse)
})

cat("=== FINAL COMPREHENSIVE CAMK ANALYSIS ===\n")
cat("Target: 4 high-quality datasets, 436 total samples\n\n")

# Get CAMK gene definitions
camk_core_genes <- get_camk_gene_categories()$core
cat("Core CAMK genes:", paste(camk_core_genes, collapse = ", "), "\n\n")

# Define finalized datasets
finalized_datasets <- c("GSE57338", "GSE115574", "GSE41177", "GSE79768")

# Function to analyze a single dataset with comprehensive logging
analyze_dataset_final <- function(dataset_id) {
  cat("=== ANALYZING", dataset_id, "===\n")
  
  # Find dataset file
  dataset_files <- list.files("cache", pattern = paste0(dataset_id, "_processed.rds"), 
                             recursive = TRUE, full.names = TRUE)
  
  if (length(dataset_files) == 0) {
    cat("ERROR: No file found for", dataset_id, "\n\n")
    return(NULL)
  }
  
  # Use first matching file
  dataset_file <- dataset_files[1]
  cat("Loading from:", dataset_file, "\n")
  
  tryCatch({
    # Load dataset
    dataset <- readRDS(dataset_file)
    
    if (!dataset$success || is.null(dataset$expression_matrix)) {
      cat("ERROR: Invalid dataset\n\n")
      return(NULL)
    }
    
    expr_matrix <- dataset$expression_matrix
    cat("Dataset dimensions:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
    
    # Verify CAMK genes
    camk_present <- intersect(rownames(expr_matrix), camk_core_genes)
    cat("CAMK genes detected:", length(camk_present), "/", length(camk_core_genes), "\n")
    cat("CAMK genes:", paste(camk_present, collapse = ", "), "\n")
    
    if (length(camk_present) < 5) {
      cat("ERROR: Insufficient CAMK genes detected\n\n")
      return(NULL)
    }
    
    # Enhanced group detection
    detected_groups <- enhanced_auto_detect_groups(dataset)
    
    if (is.null(detected_groups)) {
      cat("ERROR: No groups detected\n\n")
      return(NULL)
    }
    
    cat("Groups detected:", detected_groups$pattern_type, "\n")
    cat("Sample distribution:\n")
    print(table(detected_groups$groups))
    
    # Prepare expression data
    if (!is.null(detected_groups$sample_indices)) {
      expr_filtered <- expr_matrix[, detected_groups$sample_indices, drop = FALSE]
      groups_vector <- detected_groups$groups
    } else {
      expr_filtered <- expr_matrix
      groups_vector <- detected_groups$groups
    }
    
    cat("Analysis matrix:", nrow(expr_filtered), "genes x", ncol(expr_filtered), "samples\n")
    
    # Ensure we have proper contrasts
    if (length(unique(groups_vector)) != 2) {
      cat("ERROR: Need exactly 2 groups for analysis\n\n")
      return(NULL)
    }
    
    # Create design matrix
    design <- model.matrix(~ groups_vector)
    colnames(design) <- c("Intercept", "Disease_vs_Control")
    
    cat("Running genome-wide differential expression analysis...\n")
    
    # Fit linear model
    fit <- lmFit(expr_filtered, design)
    fit <- eBayes(fit)
    
    # Get all results
    all_dge_results <- topTable(fit, coef = "Disease_vs_Control", 
                               number = Inf, adjust.method = "BH")
    
    cat("Genome-wide analysis completed:", nrow(all_dge_results), "genes\n")
    
    # Filter for CAMK genes
    camk_dge_results <- all_dge_results[rownames(all_dge_results) %in% camk_present, , drop = FALSE]
    
    if (nrow(camk_dge_results) > 0) {
      # Format results
      camk_dge_results$Gene_Symbol <- rownames(camk_dge_results)
      camk_dge_results$Dataset <- dataset_id
      camk_dge_results$Significant <- camk_dge_results$adj.P.Val < 0.05
      camk_dge_results$Regulation <- ifelse(camk_dge_results$logFC > 0, 
                                           "UP in Disease", "DOWN in Disease")
      
      cat("CAMK genes in final results:", nrow(camk_dge_results), "\n")
      cat("Significant CAMK genes:", sum(camk_dge_results$Significant), "\n")
      
      # Display key results
      cat("Key CAMK results:\n")
      for (i in 1:nrow(camk_dge_results)) {
        gene <- camk_dge_results$Gene_Symbol[i]
        logfc <- round(camk_dge_results$logFC[i], 4)
        adj_pval <- camk_dge_results$adj.P.Val[i]
        sig_status <- if (adj_pval < 0.05) "SIG" else "NS"
        direction <- if (logfc > 0) "UP" else "DOWN"
        
        cat(sprintf("  %-8s: %4s logFC=%7.4f, FDR=%8.2e [%s]\n", 
                   gene, direction, logfc, adj_pval, sig_status))
      }
      
      cat("SUCCESS: Analysis completed for", dataset_id, "\n\n")
      return(camk_dge_results)
    } else {
      cat("ERROR: No CAMK results generated\n\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("ERROR analyzing", dataset_id, ":", e$message, "\n\n")
    return(NULL)
  })
}

# Analyze all finalized datasets
cat("Starting analysis of", length(finalized_datasets), "finalized datasets...\n\n")

all_results <- list()
successful_analyses <- 0
total_samples_analyzed <- 0

for (dataset_id in finalized_datasets) {
  result <- analyze_dataset_final(dataset_id)
  
  if (!is.null(result)) {
    all_results[[dataset_id]] <- result
    successful_analyses <- successful_analyses + 1
    
    # Count samples for this dataset
    dataset_file <- list.files("cache", pattern = paste0(dataset_id, "_processed.rds"), 
                              recursive = TRUE, full.names = TRUE)[1]
    dataset_info <- readRDS(dataset_file)
    total_samples_analyzed <- total_samples_analyzed + dataset_info$n_samples
  }
}

# Combine all DGE results
if (length(all_results) > 0) {
  cat("=== FINAL DGE RESULTS SUMMARY ===\n")
  cat("Successfully analyzed datasets:", length(all_results), "/", length(finalized_datasets), "\n")
  cat("Datasets included:", paste(names(all_results), collapse = ", "), "\n")
  cat("Total samples analyzed:", total_samples_analyzed, "\n")
  
  # Combine results
  final_dge_results <- do.call(rbind, all_results)
  rownames(final_dge_results) <- NULL
  
  cat("Total gene-dataset combinations:", nrow(final_dge_results), "\n")
  
  # Save final DGE results
  dge_output_file <- "output/CAMK_focused_DGE_all_datasets_FINAL.csv"
  write.csv(final_dge_results, dge_output_file, row.names = FALSE)
  cat("Final DGE results saved to:", dge_output_file, "\n")
  
  # Dataset-wise summary
  cat("\nDataset-wise CAMK gene coverage:\n")
  for (dataset in names(all_results)) {
    dataset_results <- all_results[[dataset]]
    sig_genes <- sum(dataset_results$Significant)
    cat(sprintf("  %-10s: %2d CAMK genes, %d significant\n", 
               dataset, nrow(dataset_results), sig_genes))
  }
  
  # CAMK2D across all datasets
  camk2d_results <- final_dge_results[final_dge_results$Gene_Symbol == "CAMK2D", ]
  cat("\nCAMK2D results across all datasets:\n")
  if (nrow(camk2d_results) > 0) {
    for (i in 1:nrow(camk2d_results)) {
      ds <- camk2d_results$Dataset[i]
      logfc <- round(camk2d_results$logFC[i], 4)
      pval <- camk2d_results$P.Value[i]
      adj_pval <- camk2d_results$adj.P.Val[i]
      sig <- camk2d_results$Significant[i]
      
      sig_marker <- if (sig) " [SIGNIFICANT]" else " [NS]"
      direction <- if (logfc > 0) "UP" else "DOWN"
      
      cat(sprintf("  %-10s: %4s logFC=%7.4f, p=%8.2e, FDR=%8.2e%s\n", 
                 ds, direction, logfc, pval, adj_pval, sig_marker))
    }
  }
  
} else {
  cat("ERROR: No datasets could be analyzed successfully\n")
  quit(status = 1)
}

cat("\n=== FINAL COMPREHENSIVE DGE ANALYSIS COMPLETE ===\n")
cat("Ready for meta-analysis with", total_samples_analyzed, "total samples\n")