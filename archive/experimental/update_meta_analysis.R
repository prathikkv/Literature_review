#!/usr/bin/env Rscript
#' Update Meta-Analysis to Include All Available Datasets
#' 
#' This script recalculates the meta-analysis with all available datasets

library(metafor)
library(tidyverse)

cat("=== Updating Meta-Analysis with All Datasets ===\n\n")

# Load updated DGE results
dge_file <- "output/CAMK_focused_DGE_all_datasets_UPDATED.csv"
if (!file.exists(dge_file)) {
  cat("ERROR: Updated DGE results not found\n")
  quit(status = 1)
}

dge_results <- read.csv(dge_file, stringsAsFactors = FALSE)
datasets_included <- unique(dge_results$Dataset)

cat("DGE Results Summary:\n")
cat("Total datasets:", length(datasets_included), "\n")
cat("Datasets:", paste(datasets_included, collapse = ", "), "\n")
cat("Total gene-dataset combinations:", nrow(dge_results), "\n\n")

# Get unique genes
unique_genes <- unique(dge_results$Gene_Symbol)
cat("Unique CAMK genes analyzed:", length(unique_genes), "\n")
cat("Genes:", paste(unique_genes, collapse = ", "), "\n\n")

# Function to perform meta-analysis for a single gene
meta_analyze_gene <- function(gene_name) {
  
  # Get data for this gene
  gene_data <- dge_results[dge_results$Gene_Symbol == gene_name, ]
  
  if (nrow(gene_data) < 2) {
    # Need at least 2 studies for meta-analysis
    return(NULL)
  }
  
  # Prepare data for meta-analysis
  logfc_values <- gene_data$logFC
  se_values <- abs(logfc_values / gene_data$t)  # Standard error from t-statistic
  
  # Handle cases where SE calculation fails
  se_values[is.na(se_values) | se_values == 0] <- 0.1
  
  tryCatch({
    # Perform fixed-effects meta-analysis
    meta_result <- rma(yi = logfc_values, sei = se_values, method = "FE")
    
    # Calculate heterogeneity statistics
    i2 <- max(0, (meta_result$QE - meta_result$k + 1) / meta_result$QE * 100)
    het_p <- meta_result$QEp
    
    # Determine regulation direction
    regulation <- ifelse(meta_result$beta > 0, "UP in Disease", "DOWN in Disease")
    
    # Determine significance
    significant <- meta_result$pval < 0.05
    
    # Create result
    result <- data.frame(
      Gene = gene_name,
      N_Studies = nrow(gene_data),
      Combined_logFC = as.numeric(meta_result$beta),
      Combined_SE = as.numeric(meta_result$se),
      Combined_P_Value = as.numeric(meta_result$pval),
      CI_Lower = as.numeric(meta_result$ci.lb),
      CI_Upper = as.numeric(meta_result$ci.ub),
      Heterogeneity_I2 = round(i2, 1),
      Heterogeneity_P = het_p,
      Regulation = regulation,
      Significant = significant,
      Datasets = paste(gene_data$Dataset, collapse = ", "),
      Individual_logFC = paste(round(gene_data$logFC, 3), collapse = ", "),
      Individual_P_Values = paste(round(gene_data$P.Value, 6), collapse = ", "),
      stringsAsFactors = FALSE
    )
    
    return(result)
    
  }, error = function(e) {
    cat("Warning: Meta-analysis failed for", gene_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Perform meta-analysis for all genes
cat("Performing meta-analysis for all genes...\n")

meta_results <- list()
for (gene in unique_genes) {
  gene_result <- meta_analyze_gene(gene)
  if (!is.null(gene_result)) {
    meta_results[[gene]] <- gene_result
    cat("  ", gene, ": Combined logFC =", round(gene_result$Combined_logFC, 3), 
        ", p =", format(gene_result$Combined_P_Value, scientific = TRUE, digits = 3), "\n")
  }
}

# Combine results
if (length(meta_results) > 0) {
  combined_meta <- do.call(rbind, meta_results)
  rownames(combined_meta) <- NULL
  
  # Sort by p-value
  combined_meta <- combined_meta[order(combined_meta$Combined_P_Value), ]
  
  cat("\n=== Meta-Analysis Results Summary ===\n")
  cat("Genes with meta-analysis:", nrow(combined_meta), "\n")
  cat("Significant genes (p < 0.05):", sum(combined_meta$Significant), "\n")
  
  # Save results
  output_file <- "output/CAMK_meta_analysis_summary_UPDATED.csv"
  write.csv(combined_meta, output_file, row.names = FALSE)
  cat("Updated meta-analysis saved to:", output_file, "\n")
  
  # Display significant results
  if (any(combined_meta$Significant)) {
    cat("\nSignificant genes:\n")
    sig_genes <- combined_meta[combined_meta$Significant, ]
    for (i in 1:nrow(sig_genes)) {
      gene <- sig_genes$Gene[i]
      logfc <- round(sig_genes$Combined_logFC[i], 4)
      pval <- sig_genes$Combined_P_Value[i]
      n_studies <- sig_genes$N_Studies[i]
      i2 <- sig_genes$Heterogeneity_I2[i]
      
      cat(sprintf("  %-8s: %4s logFC=%7.4f, p=%8.2e (%d studies, I²=%3.1f%%)\n",
                 gene, ifelse(logfc > 0, "UP", "DOWN"), logfc, pval, n_studies, i2))
    }
  } else {
    cat("\nNo genes reached statistical significance in meta-analysis\n")
  }
  
  # Special focus on CAMK2D
  camk2d_result <- combined_meta[combined_meta$Gene == "CAMK2D", ]
  if (nrow(camk2d_result) > 0) {
    cat("\n=== CAMK2D Meta-Analysis Results ===\n")
    cat("Studies included:", camk2d_result$N_Studies, "\n")
    cat("Datasets:", camk2d_result$Datasets, "\n")
    cat("Combined logFC:", round(camk2d_result$Combined_logFC, 4), "\n")
    cat("95% CI: [", round(camk2d_result$CI_Lower, 4), ",", round(camk2d_result$CI_Upper, 4), "]\n")
    cat("P-value:", format(camk2d_result$Combined_P_Value, scientific = TRUE, digits = 4), "\n")
    cat("Significant:", camk2d_result$Significant, "\n")
    cat("Heterogeneity I²:", camk2d_result$Heterogeneity_I2, "%\n")
    cat("Direction:", camk2d_result$Regulation, "\n")
  }
  
} else {
  cat("ERROR: No meta-analysis results could be generated\n")
}

cat("\n=== Meta-Analysis Update Complete ===\n")