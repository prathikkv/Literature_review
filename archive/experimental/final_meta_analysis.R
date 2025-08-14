#!/usr/bin/env Rscript
#' Final Meta-Analysis with All Available High-Quality Datasets
#' 
#' Performs comprehensive meta-analysis on final DGE results

library(metafor)
library(tidyverse)

cat("=== FINAL META-ANALYSIS ===\n\n")

# Load final DGE results
dge_file <- "output/CAMK_focused_DGE_all_datasets_FINAL.csv"
if (!file.exists(dge_file)) {
  cat("ERROR: Final DGE results not found. Run final_comprehensive_analysis.R first.\n")
  quit(status = 1)
}

dge_results <- read.csv(dge_file, stringsAsFactors = FALSE)
datasets_included <- unique(dge_results$Dataset)

cat("Loading final DGE results...\n")
cat("Datasets included:", length(datasets_included), "-", paste(datasets_included, collapse = ", "), "\n")
cat("Total gene-dataset combinations:", nrow(dge_results), "\n")

# Calculate total samples
sample_counts <- c("GSE57338" = 313, "GSE115574" = 59, "GSE41177" = 38, "GSE79768" = 26)
total_samples <- sum(sample_counts[names(sample_counts) %in% datasets_included])
cat("Total samples in meta-analysis:", total_samples, "\n\n")

# Get unique CAMK genes
unique_genes <- unique(dge_results$Gene_Symbol)
cat("CAMK genes for meta-analysis:", length(unique_genes), "\n")
cat("Genes:", paste(unique_genes, collapse = ", "), "\n\n")

# Function to perform meta-analysis for a single gene
perform_gene_meta_analysis <- function(gene_name) {
  
  # Get data for this gene across all datasets
  gene_data <- dge_results[dge_results$Gene_Symbol == gene_name, ]
  
  if (nrow(gene_data) < 2) {
    return(NULL)  # Need at least 2 studies
  }
  
  # Calculate standard errors from t-statistics
  se_values <- abs(gene_data$logFC / gene_data$t)
  
  # Handle invalid SEs
  se_values[is.na(se_values) | se_values == 0 | !is.finite(se_values)] <- 0.1
  
  tryCatch({
    # Perform fixed-effects meta-analysis
    meta_result <- rma(yi = gene_data$logFC, sei = se_values, method = "FE")
    
    # Calculate heterogeneity
    i2 <- max(0, (meta_result$QE - meta_result$k + 1) / meta_result$QE * 100)
    
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
      Heterogeneity_P = meta_result$QEp,
      Regulation = ifelse(meta_result$beta > 0, "UP in Disease", "DOWN in Disease"),
      Significant = meta_result$pval < 0.05,
      Datasets = paste(gene_data$Dataset, collapse = ", "),
      Individual_logFC = paste(round(gene_data$logFC, 4), collapse = ", "),
      Individual_P_Values = paste(format(gene_data$P.Value, scientific = TRUE, digits = 3), collapse = ", "),
      stringsAsFactors = FALSE
    )
    
    return(result)
    
  }, error = function(e) {
    cat("Warning: Meta-analysis failed for", gene_name, ":", e$message, "\n")
    return(NULL)
  })
}

# Perform meta-analysis for all genes
cat("Performing meta-analysis...\n")

meta_results_list <- list()
for (gene in unique_genes) {
  result <- perform_gene_meta_analysis(gene)
  if (!is.null(result)) {
    meta_results_list[[gene]] <- result
    direction <- if (result$Combined_logFC > 0) "UP" else "DOWN"
    sig_marker <- if (result$Significant) "*" else ""
    cat(sprintf("  %-8s: %4s logFC=%7.4f, p=%8.2e, I²=%5.1f%% %s\n", 
               gene, direction, result$Combined_logFC, result$Combined_P_Value, 
               result$Heterogeneity_I2, sig_marker))
  }
}

# Combine results
if (length(meta_results_list) > 0) {
  final_meta_results <- do.call(rbind, meta_results_list)
  rownames(final_meta_results) <- NULL
  
  # Sort by significance
  final_meta_results <- final_meta_results[order(final_meta_results$Combined_P_Value), ]
  
  cat("\n=== FINAL META-ANALYSIS SUMMARY ===\n")
  cat("Total CAMK genes analyzed:", nrow(final_meta_results), "\n")
  cat("Significant genes (p < 0.05):", sum(final_meta_results$Significant), "\n")
  cat("Datasets contributing:", paste(datasets_included, collapse = ", "), "\n")
  cat("Total samples:", total_samples, "\n")
  
  # Save results
  output_file <- "output/CAMK_meta_analysis_summary_FINAL.csv"
  write.csv(final_meta_results, output_file, row.names = FALSE)
  cat("Final meta-analysis saved to:", output_file, "\n")
  
  # Display significant results
  sig_results <- final_meta_results[final_meta_results$Significant, ]
  if (nrow(sig_results) > 0) {
    cat("\nSIGNIFICANT CAMK GENES:\n")
    for (i in 1:nrow(sig_results)) {
      gene <- sig_results$Gene[i]
      logfc <- round(sig_results$Combined_logFC[i], 4)
      pval <- sig_results$Combined_P_Value[i]
      n_studies <- sig_results$N_Studies[i]
      i2 <- sig_results$Heterogeneity_I2[i]
      ci_lower <- round(sig_results$CI_Lower[i], 4)
      ci_upper <- round(sig_results$CI_Upper[i], 4)
      
      direction <- if (logfc > 0) "UP" else "DOWN"
      
      cat(sprintf("%-2d. %-8s: %4s logFC=%7.4f [%7.4f, %7.4f], p=%8.2e (%d studies, I²=%3.1f%%)\n",
                 i, gene, direction, logfc, ci_lower, ci_upper, pval, n_studies, i2))
    }
  } else {
    cat("\nNo genes reached statistical significance (p < 0.05)\n")
  }
  
  # CAMK2D detailed results
  camk2d_result <- final_meta_results[final_meta_results$Gene == "CAMK2D", ]
  if (nrow(camk2d_result) > 0) {
    cat("\n=== CAMK2D META-ANALYSIS SPOTLIGHT ===\n")
    cat("Studies included:", camk2d_result$N_Studies, "\n")
    cat("Datasets:", camk2d_result$Datasets, "\n")
    cat("Combined effect (logFC):", round(camk2d_result$Combined_logFC, 4), "\n")
    cat("95% Confidence Interval: [", round(camk2d_result$CI_Lower, 4), ",", 
        round(camk2d_result$CI_Upper, 4), "]\n")
    cat("Meta-analysis p-value:", format(camk2d_result$Combined_P_Value, scientific = TRUE, digits = 4), "\n")
    cat("Statistical significance:", if (camk2d_result$Significant) "YES" else "NO", "\n")
    cat("Effect direction:", camk2d_result$Regulation, "\n")
    cat("Between-study heterogeneity (I²):", camk2d_result$Heterogeneity_I2, "%\n")
    cat("Individual study effects:", camk2d_result$Individual_logFC, "\n")
    cat("Individual study p-values:", camk2d_result$Individual_P_Values, "\n")
  }
  
} else {
  cat("ERROR: No meta-analysis results could be generated\n")
  quit(status = 1)
}

cat("\n=== FINAL META-ANALYSIS COMPLETE ===\n")
cat("Analysis based on", total_samples, "samples across", length(datasets_included), "high-quality datasets\n")