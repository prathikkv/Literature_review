#!/usr/bin/env Rscript
#' Fix GSE57338 Gene Mapping - Reprocess to Get Full Gene Set
#' 
#' GSE57338 currently only has 108 genes due to processing issues
#' This script will reprocess it to get the full ~20K+ genes with proper mapping

library(GEOquery)
library(limma)

# Load functions
source("functions/data_processing.R")
source("functions/camk_definitions.R")

cat("=== Fixing GSE57338 Gene Mapping ===\n\n")

# Check current status
current_file <- "cache/comprehensive/GSE57338_processed.rds"
current_data <- readRDS(current_file)

cat("Current GSE57338 status:\n")
cat("- Success:", current_data$success, "\n")
cat("- Genes:", current_data$n_genes, "\n")
cat("- Samples:", current_data$n_samples, "\n")
cat("- Platform:", current_data$dataset_info$platform, "\n\n")

# The issue is that GSE57338 was processed in a limited way
# Let's try to reprocess it properly with full gene mapping

cat("Attempting to reprocess GSE57338 with full gene mapping...\n")

# Method 1: Try to download fresh data with proper processing
tryCatch({
  cat("Downloading GSE57338 from GEO...\n")
  
  # Download the dataset
  gset <- getGEO("GSE57338", GSEMatrix = TRUE, AnnotGPL = TRUE, getGPL = TRUE)
  
  if (length(gset) > 0) {
    gset <- gset[[1]]
    
    # Get expression data
    expr_data <- exprs(gset)
    pheno_data <- pData(gset)
    feature_data <- fData(gset)
    
    cat("Downloaded data:\n")
    cat("- Expression matrix:", nrow(expr_data), "x", ncol(expr_data), "\n")
    cat("- Feature data:", nrow(feature_data), "columns:", ncol(feature_data), "\n")
    cat("- Phenotype data:", nrow(pheno_data), "columns:", ncol(pheno_data), "\n")
    
    # Look for gene symbol columns in feature data
    cat("Feature data columns:", paste(names(feature_data), collapse = ", "), "\n")
    
    # Common gene symbol column names  
    gene_symbol_cols <- c("Gene symbol", "Gene.Symbol", "GENE_SYMBOL", "gene_symbol", "Symbol", "symbol", "Gene_Symbol")
    gene_col <- NULL
    
    for (col in gene_symbol_cols) {
      if (col %in% names(feature_data) && !all(is.na(feature_data[[col]]))) {
        gene_col <- col
        break
      }
    }
    
    if (!is.null(gene_col)) {
      cat("Found gene symbol column:", gene_col, "\n")
      gene_symbols <- feature_data[[gene_col]]
      
      # Filter out empty/NA gene symbols
      valid_genes <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
      
      if (sum(valid_genes) > 1000) {  # Should have many genes
        # Filter expression matrix to genes with valid symbols
        expr_filtered <- expr_data[valid_genes, ]
        gene_symbols_filtered <- gene_symbols[valid_genes]
        
        # Set rownames to gene symbols
        rownames(expr_filtered) <- gene_symbols_filtered
        
        # Handle duplicate gene symbols by averaging
        if (any(duplicated(gene_symbols_filtered))) {
          cat("Handling", sum(duplicated(gene_symbols_filtered)), "duplicate gene symbols...\n")
          
          unique_genes <- unique(gene_symbols_filtered)
          aggregated_matrix <- matrix(NA, nrow = length(unique_genes), ncol = ncol(expr_filtered))
          rownames(aggregated_matrix) <- unique_genes
          colnames(aggregated_matrix) <- colnames(expr_filtered)
          
          for (gene in unique_genes) {
            gene_indices <- which(gene_symbols_filtered == gene)
            if (length(gene_indices) == 1) {
              aggregated_matrix[gene, ] <- expr_filtered[gene_indices, ]
            } else {
              aggregated_matrix[gene, ] <- colMeans(expr_filtered[gene_indices, , drop = FALSE])
            }
          }
          
          expr_filtered <- aggregated_matrix
        }
        
        cat("Final expression matrix:", nrow(expr_filtered), "genes x", ncol(expr_filtered), "samples\n")
        
        # Create proper dataset structure
        fixed_dataset <- list(
          success = TRUE,
          dataset_id = "GSE57338",
          expression_matrix = expr_filtered,
          phenotype_data = pheno_data,
          feature_data = feature_data,
          dataset_info = list(
            platform = "GPL570",
            tissue = "heart",
            condition = "heart_failure",
            species = "human"
          ),
          n_genes = nrow(expr_filtered),
          n_samples = ncol(expr_filtered),
          download_time = Sys.time()
        )
        
        # Save the fixed dataset
        output_file <- "cache/comprehensive/GSE57338_processed_FIXED.rds"
        saveRDS(fixed_dataset, output_file)
        
        cat("SUCCESS: Fixed GSE57338 saved to:", output_file, "\n")
        cat("- Genes increased from", current_data$n_genes, "to", nrow(expr_filtered), "\n")
        cat("- Samples:", ncol(expr_filtered), "\n")
        
        # Verify CAMK genes are present
        camk_genes <- get_camk_gene_categories()$core
        camk_present <- intersect(rownames(expr_filtered), camk_genes)
        cat("- CAMK genes found:", length(camk_present), "/", length(camk_genes), "\n")
        cat("- CAMK genes:", paste(camk_present, collapse = ", "), "\n")
        
      } else {
        cat("ERROR: Too few valid gene symbols found (", sum(valid_genes), ")\n")
      }
    } else {
      cat("ERROR: No gene symbol column found in feature data\n")
    }
    
  } else {
    cat("ERROR: Could not download GSE57338\n")
  }
  
}, error = function(e) {
  cat("ERROR downloading GSE57338:", e$message, "\n")
})

cat("\n=== GSE57338 Gene Mapping Fix Attempt Complete ===\n")