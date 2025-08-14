#!/usr/bin/env Rscript
#' Analyze Missing Datasets for Complete CAMK Analysis
#' 
#' This script performs DGE analysis on all available but unanalyzed datasets

library(limma)
library(edgeR)
library(DESeq2)

# Load functions
source("functions/camk_definitions.R")
source("functions/data_processing.R")
source("functions/analysis.R")
source("functions/common_utils.R")

# Get CAMK genes
camk_genes <- get_camk_gene_categories()$core
cat("=== Analyzing Missing Datasets for CAMK Genes ===\n")
cat("Core CAMK genes:", paste(camk_genes, collapse = ", "), "\n\n")

# Function to perform DGE analysis on a dataset
analyze_dataset <- function(dataset_path, dataset_id) {
  cat("Analyzing", dataset_id, "\n")
  
  tryCatch({
    # Load dataset
    dataset <- readRDS(dataset_path)
    
    if (!dataset$success || is.null(dataset$expression_matrix)) {
      cat("  ERROR: Invalid dataset\n")
      return(NULL)
    }
    
    expr_matrix <- dataset$expression_matrix
    cat("  Dataset:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
    
    # Get gene names (try different methods)
    gene_names <- rownames(expr_matrix)
    if (is.null(gene_names)) {
      cat("  WARNING: No rownames found\n")
      return(NULL)
    }
    
    # Check for CAMK genes
    camk_present <- intersect(gene_names, camk_genes)
    cat("  CAMK genes found:", length(camk_present), "/", length(camk_genes), "\n")
    
    if (length(camk_present) == 0) {
      cat("  WARNING: No CAMK genes detected\n")
      return(NULL)
    }
    
    cat("  CAMK genes present:", paste(camk_present, collapse = ", "), "\n")
    
    # Try to detect groups
    if (is.null(dataset$phenotype_data)) {
      cat("  WARNING: No phenotype data available\n")
      return(NULL)
    }
    
    # Simple group detection - look for healthy vs disease patterns
    pheno <- dataset$phenotype_data
    group_col <- NULL
    
    # Common group column names
    group_candidates <- c("Group", "group", "condition", "Condition", 
                         "disease_state", "status", "Status", "type", "Type")
    
    for (col in group_candidates) {
      if (col %in% names(pheno)) {
        group_col <- col
        break
      }
    }
    
    if (is.null(group_col)) {
      cat("  WARNING: No suitable group column found\n")
      return(NULL)
    }
    
    groups <- pheno[[group_col]]
    unique_groups <- unique(groups)
    cat("  Groups found:", paste(unique_groups, collapse = ", "), "\n")
    
    if (length(unique_groups) != 2) {
      cat("  WARNING: Need exactly 2 groups for comparison\n")
      return(NULL)
    }
    
    # Create design matrix
    design <- model.matrix(~ groups)
    
    # Perform DGE analysis using limma
    fit <- lmFit(expr_matrix, design)
    fit2 <- eBayes(fit)
    
    # Extract results for all genes
    all_results <- topTable(fit2, coef = 2, number = Inf, sort.by = "none")
    
    # Filter for CAMK genes
    camk_results <- all_results[rownames(all_results) %in% camk_present, ]
    
    if (nrow(camk_results) > 0) {
      # Format results
      formatted_results <- data.frame(
        logFC = camk_results$logFC,
        AveExpr = camk_results$AveExpr,
        t = camk_results$t,
        P.Value = camk_results$P.Value,
        adj.P.Val = camk_results$adj.P.Val,
        B = camk_results$B,
        Gene_Symbol = rownames(camk_results),
        Dataset = dataset_id,
        Significant = camk_results$adj.P.Val < 0.05,
        Regulation = ifelse(camk_results$logFC > 0, "UP in Disease", "DOWN in Disease"),
        stringsAsFactors = FALSE
      )
      
      cat("  SUCCESS: DGE analysis completed,", nrow(formatted_results), "CAMK genes analyzed\n")
      return(formatted_results)
    }
    
    return(NULL)
    
  }, error = function(e) {
    cat("  ERROR:", e$message, "\n")
    return(NULL)
  })
}

# Find all datasets to analyze
missing_datasets <- c("GSE57338", "GSE115574", "GSE161472", "GSE224997", "GSE226282")
all_new_results <- list()

for (dataset_id in missing_datasets) {
  # Find cached file
  cache_files <- list.files("cache", 
                           pattern = paste0(dataset_id, "_processed.rds"),
                           recursive = TRUE, full.names = TRUE)
  
  if (length(cache_files) > 0) {
    cat("\n=== Processing", dataset_id, "===\n")
    results <- analyze_dataset(cache_files[1], dataset_id)
    if (!is.null(results)) {
      all_new_results[[dataset_id]] <- results
    }
  } else {
    cat("WARNING: No cached file found for", dataset_id, "\n")
  }
}

# Combine all new results
if (length(all_new_results) > 0) {
  new_dge_results <- do.call(rbind, all_new_results)
  rownames(new_dge_results) <- NULL
  
  cat("\n=== New DGE Results Summary ===\n")
  cat("New datasets analyzed:", length(all_new_results), "\n")
  cat("Total gene-dataset combinations:", nrow(new_dge_results), "\n")
  
  # Load existing results
  existing_file <- "output/CAMK_focused_DGE_all_datasets.csv"
  existing_results <- read.csv(existing_file, stringsAsFactors = FALSE)
  
  # Combine with existing
  complete_results <- rbind(existing_results, new_dge_results)
  
  # Save complete results
  output_file <- "output/CAMK_focused_DGE_all_datasets_COMPLETE.csv"
  write.csv(complete_results, output_file, row.names = FALSE)
  
  cat("Complete DGE results saved to:", output_file, "\n")
  cat("Total datasets now:", length(unique(complete_results$Dataset)), "\n")
  cat("Total combinations:", nrow(complete_results), "\n")
  
  # Check CAMK2D across all datasets
  camk2d_all <- complete_results[complete_results$Gene_Symbol == "CAMK2D", ]
  cat("\nCAMK2D results across all datasets:\n")
  if (nrow(camk2d_all) > 0) {
    print(camk2d_all[, c("Dataset", "logFC", "P.Value", "adj.P.Val", "Significant")])
  }
  
} else {
  cat("\nWARNING: No new datasets could be analyzed\n")
}

cat("\n=== Analysis Complete ===\n")