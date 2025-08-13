#!/usr/bin/env Rscript
#' CAMK-Focused Analysis Pipeline
#' 
#' Performs differential expression analysis focusing exclusively on CAMK family genes

# Load required functions and libraries
source("functions/data_processing.R")
source("functions/analysis.R")
source("enhanced_group_detection.R")
library(limma)
library(ggplot2)
library(pheatmap)

cat("🎯 CAMK-FOCUSED ANALYSIS PIPELINE\n")
cat("==================================\n\n")

# Define CAMK family genes
camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMKK1", "CAMKK2", "CAMK1", "CAMK1D", "CAMK1G", "CAMK4", "CAMKV")

cat("🧬 Target CAMK genes:", length(camk_genes), "\n")
cat("   ", paste(camk_genes, collapse = ", "), "\n\n")

# Load and analyze datasets with proper gene symbols
cache_dir <- "cache/comprehensive"
processed_files <- list.files(cache_dir, pattern = "_processed.rds$", full.names = FALSE)

camk_analysis_results <- list()

for (file in processed_files) {
  dataset_id <- gsub("_processed.rds", "", file)
  
  cat("📊 Analyzing", dataset_id, "\n")
  cat("   Loading dataset...\n")
  
  # Load dataset
  dataset <- readRDS(file.path(cache_dir, file))
  
  if (!dataset$success || is.null(dataset$expression_matrix)) {
    cat("   ❌ Dataset failed or no expression matrix - skipping\n\n")
    next
  }
  
  expr_matrix <- dataset$expression_matrix
  n_genes <- nrow(expr_matrix)
  n_samples <- ncol(expr_matrix)
  
  cat("   Dataset info:", n_genes, "genes x", n_samples, "samples\n")
  
  # Check for CAMK genes in this dataset
  gene_names <- rownames(expr_matrix)
  camk_present <- camk_genes[camk_genes %in% gene_names]
  
  cat("   CAMK genes found:", length(camk_present), "/", length(camk_genes), "\n")
  
  if (length(camk_present) == 0) {
    cat("   ⚠️ No CAMK genes detected - skipping analysis\n\n")
    next
  }
  
  # Print found CAMK genes
  cat("   Present:", paste(camk_present, collapse = ", "), "\n")
  
  # Filter expression matrix to CAMK genes only
  camk_expr_matrix <- expr_matrix[camk_present, , drop = FALSE]
  
  cat("   📈 Filtered to CAMK genes:", nrow(camk_expr_matrix), "genes x", ncol(camk_expr_matrix), "samples\n")
  
  # Try to detect healthy vs disease groups
  detected_groups <- enhanced_auto_detect_groups(dataset)
  
  if (is.null(detected_groups)) {
    cat("   ⚠️ No suitable groups detected - skipping DGE analysis\n\n")
    
    # Still save basic expression data
    camk_analysis_results[[dataset_id]] <- list(
      dataset_id = dataset_id,
      camk_genes_present = camk_present,
      camk_expression_matrix = camk_expr_matrix,
      n_samples = n_samples,
      groups_detected = FALSE
    )
    next
  }
  
  cat("   ✅ Groups detected:", detected_groups$pattern_type, "\n")
  cat("      Healthy:", sum(detected_groups$groups == "Healthy"), "samples\n")
  cat("      Disease:", sum(detected_groups$groups == "Disease"), "samples\n")
  
  # Prepare data for DGE analysis
  if (!is.null(detected_groups$sample_indices)) {
    # Filter expression matrix to samples with group assignments
    camk_expr_filtered <- camk_expr_matrix[, detected_groups$sample_indices, drop = FALSE]
    groups_vector <- detected_groups$groups
  } else {
    camk_expr_filtered <- camk_expr_matrix
    groups_vector <- detected_groups$groups
  }
  
  # Perform CAMK-focused DGE analysis
  tryCatch({
    cat("   🔬 Running CAMK-focused DGE analysis...\n")
    
    # Create design matrix (Disease vs Healthy)
    design <- model.matrix(~ groups_vector)
    colnames(design) <- c("Intercept", "Disease_vs_Healthy")
    
    # Fit linear model
    fit <- lmFit(camk_expr_filtered, design)
    fit <- eBayes(fit)
    
    # Get results
    camk_dge_results <- topTable(fit, coef = "Disease_vs_Healthy", number = Inf, adjust.method = "BH")
    
    # Add gene information
    camk_dge_results$Gene_Symbol <- rownames(camk_dge_results)
    camk_dge_results$Dataset <- dataset_id
    camk_dge_results$Significant <- camk_dge_results$adj.P.Val < 0.05
    camk_dge_results$Regulation <- ifelse(camk_dge_results$logFC > 0, "UP in Disease", "DOWN in Disease")
    
    cat("   ✅ DGE completed!\n")
    cat("      Total CAMK genes analyzed:", nrow(camk_dge_results), "\n")
    cat("      Significant (FDR < 0.05):", sum(camk_dge_results$Significant), "\n")
    
    # Display individual CAMK gene results
    if (nrow(camk_dge_results) > 0) {
      cat("      CAMK gene results:\n")
      
      # Sort by significance
      camk_sorted <- camk_dge_results[order(camk_dge_results$adj.P.Val), ]
      
      for (i in 1:nrow(camk_sorted)) {
        gene <- camk_sorted$Gene_Symbol[i]
        logfc <- round(camk_sorted$logFC[i], 3)
        pval <- camk_sorted$P.Value[i]
        adj_pval <- camk_sorted$adj.P.Val[i]
        sig <- camk_sorted$Significant[i]
        regulation <- camk_sorted$Regulation[i]
        
        sig_marker <- if (adj_pval < 0.001) "***" else if (adj_pval < 0.01) "**" else if (adj_pval < 0.05) "*" else ""
        direction <- if (logfc > 0) "↑" else "↓"
        
        cat(sprintf("        %-8s: %s logFC=%6.3f, FDR=%8.2e %s [%s]\n", 
                    gene, direction, logfc, adj_pval, sig_marker, regulation))
      }
    }
    
    # Store results
    camk_analysis_results[[dataset_id]] <- list(
      dataset_id = dataset_id,
      camk_genes_present = camk_present,
      camk_expression_matrix = camk_expr_filtered,
      groups = detected_groups,
      dge_results = camk_dge_results,
      n_samples_analyzed = ncol(camk_expr_filtered),
      groups_detected = TRUE,
      significant_camk_genes = sum(camk_dge_results$Significant)
    )
    
  }, error = function(e) {
    cat("   ❌ DGE analysis error:", e$message, "\n")
    
    # Store basic info even if DGE failed
    camk_analysis_results[[dataset_id]] <- list(
      dataset_id = dataset_id,
      camk_genes_present = camk_present,
      camk_expression_matrix = camk_expr_filtered,
      groups = detected_groups,
      n_samples_analyzed = ncol(camk_expr_filtered),
      groups_detected = TRUE,
      dge_error = e$message
    )
  })
  
  cat("\n")
}

# Generate summary report
cat("📋 CAMK ANALYSIS SUMMARY\n")
cat("========================\n")

datasets_with_camk <- sum(sapply(camk_analysis_results, function(x) length(x$camk_genes_present) > 0))
datasets_with_dge <- sum(sapply(camk_analysis_results, function(x) !is.null(x$dge_results)))

cat("Datasets processed:", length(camk_analysis_results), "\n")
cat("Datasets with CAMK genes:", datasets_with_camk, "\n")
cat("Datasets with successful DGE:", datasets_with_dge, "\n\n")

# Create combined CAMK results
if (datasets_with_dge > 0) {
  cat("📊 Combined CAMK DGE Results:\n")
  cat("============================\n")
  
  all_camk_results <- do.call(rbind, lapply(camk_analysis_results, function(x) {
    if (!is.null(x$dge_results)) x$dge_results else NULL
  }))
  
  if (!is.null(all_camk_results) && nrow(all_camk_results) > 0) {
    # Summary by gene
    gene_summary <- aggregate(cbind(logFC = all_camk_results$logFC, 
                                   adj.P.Val = all_camk_results$adj.P.Val,
                                   Significant = all_camk_results$Significant), 
                             by = list(Gene = all_camk_results$Gene_Symbol), 
                             FUN = function(x) c(mean = mean(x, na.rm = TRUE), 
                                               count = length(x),
                                               sig_count = sum(x < 0.05, na.rm = TRUE)))
    
    cat("CAMK genes across all datasets:\n")
    for (gene in unique(all_camk_results$Gene_Symbol)) {
      gene_data <- all_camk_results[all_camk_results$Gene_Symbol == gene, ]
      n_datasets <- nrow(gene_data)
      n_significant <- sum(gene_data$Significant)
      avg_logfc <- round(mean(gene_data$logFC), 3)
      
      cat(sprintf("  %-8s: %d datasets, %d significant, avg logFC=%6.3f\n", 
                  gene, n_datasets, n_significant, avg_logfc))
    }
  }
}

# Save results
if (!dir.exists("output")) dir.create("output")

saveRDS(camk_analysis_results, "output/CAMK_focused_analysis_results.rds")

if (exists("all_camk_results") && !is.null(all_camk_results)) {
  write.csv(all_camk_results, "output/CAMK_focused_DGE_all_datasets.csv", row.names = FALSE)
  cat("\n💾 Results saved:\n")
  cat("   • output/CAMK_focused_analysis_results.rds\n")
  cat("   • output/CAMK_focused_DGE_all_datasets.csv\n")
}

cat("\n🎉 CAMK-focused analysis completed!\n")
cat("✨ Focus on CAMK family genes achieved - ready for clinical interpretation\n")