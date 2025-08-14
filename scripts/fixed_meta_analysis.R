#!/usr/bin/env Rscript
#' Fixed Meta-Analysis Script
#' 
#' Simplified meta-analysis focusing on core functionality

library(metafor)
library(tidyverse)

cat("═══════════════════════════════════════════════════════════════\n")
cat("🔬 FIXED META-ANALYSIS PIPELINE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Load comprehensive results
dge_results <- read.csv("output/CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv", stringsAsFactors = FALSE)

# Apply quality filters
quality_filtered <- dge_results[abs(dge_results$logFC) <= 0.8, ]
gene_counts <- table(quality_filtered$Gene_Symbol)
genes_multi_dataset <- names(gene_counts)[gene_counts >= 2]
meta_data <- quality_filtered[quality_filtered$Gene_Symbol %in% genes_multi_dataset, ]

cat("📊 Data prepared for meta-analysis:\n")
cat("Genes:", length(genes_multi_dataset), "\n")
cat("Gene-dataset combinations:", nrow(meta_data), "\n\n")

# Perform meta-analysis
meta_results <- list()

for (gene in genes_multi_dataset) {
  gene_data <- meta_data[meta_data$Gene_Symbol == gene, ]
  
  if (nrow(gene_data) >= 2) {
    # Calculate standard errors
    gene_data$SE <- abs(gene_data$logFC / gene_data$t)
    
    tryCatch({
      # Simple fixed-effects meta-analysis
      meta_result <- rma(yi = logFC, sei = SE, data = gene_data, method = "FE")
      
      # Extract results
      meta_results[[gene]] <- data.frame(
        Gene = gene,
        N_Studies = nrow(gene_data),
        Combined_logFC = as.numeric(meta_result$beta),
        Combined_SE = as.numeric(meta_result$se),
        Combined_P_Value = as.numeric(meta_result$pval),
        CI_Lower = as.numeric(meta_result$ci.lb),
        CI_Upper = as.numeric(meta_result$ci.ub),
        Heterogeneity_I2 = if (!is.na(meta_result$I2)) meta_result$I2 else 0,
        Heterogeneity_P = if (!is.na(meta_result$QEp)) meta_result$QEp else 1,
        Datasets = paste(gene_data$Dataset, collapse = ", "),
        Individual_logFC = paste(round(gene_data$logFC, 4), collapse = ", "),
        Individual_P_Values = paste(sprintf("%.3e", gene_data$P.Value), collapse = ", "),
        Regulation = ifelse(meta_result$beta > 0, "UP in Disease", "DOWN in Disease"),
        Significant = meta_result$pval < 0.05,
        stringsAsFactors = FALSE
      )
      
      cat("✅", gene, "meta-analysis complete\n")
      
    }, error = function(e) {
      cat("❌", gene, "failed:", e$message, "\n")
    })
  }
}

# Combine results
if (length(meta_results) > 0) {
  final_results <- do.call(rbind, meta_results)
  rownames(final_results) <- NULL
  
  # Sort by significance
  final_results <- final_results[order(final_results$Combined_P_Value), ]
  
  # Save results
  write.csv(final_results, "output/CAMK_meta_analysis_FINAL.csv", row.names = FALSE)
  
  cat("\n📁 Results saved to: output/CAMK_meta_analysis_FINAL.csv\n\n")
  
  # Summary
  cat("📊 META-ANALYSIS SUMMARY:\n")
  cat("Genes analyzed:", nrow(final_results), "\n")
  cat("Significant genes:", sum(final_results$Significant), "\n\n")
  
  # CAMK2D results
  camk2d <- final_results[final_results$Gene == "CAMK2D", ]
  if (nrow(camk2d) > 0) {
    cat("⭐ CAMK2D RESULTS:\n")
    cat(sprintf("Combined logFC: %.4f (95%% CI: %.4f to %.4f)\n",
                camk2d$Combined_logFC, camk2d$CI_Lower, camk2d$CI_Upper))
    cat(sprintf("P-value: %.2e\n", camk2d$Combined_P_Value))
    cat(sprintf("Significant: %s\n", if (camk2d$Significant) "YES ✅" else "NO"))
    cat(sprintf("Direction: %s\n", camk2d$Regulation))
    cat(sprintf("Datasets: %s\n", camk2d$Datasets))
    cat(sprintf("Individual logFC: %s\n", camk2d$Individual_logFC))
  }
  
  # All significant results
  sig_results <- final_results[final_results$Significant, ]
  if (nrow(sig_results) > 0) {
    cat("\n🏆 SIGNIFICANT GENES:\n")
    for (i in 1:nrow(sig_results)) {
      gene <- sig_results$Gene[i]
      logfc <- sig_results$Combined_logFC[i]
      pval <- sig_results$Combined_P_Value[i]
      n_studies <- sig_results$N_Studies[i]
      direction <- if (logfc > 0) "UP" else "DOWN"
      
      cat(sprintf("  🧬 %-10s: %s (logFC=%.4f, p=%.2e, n=%d)\n",
                  gene, direction, logfc, pval, n_studies))
    }
  }
  
  cat("\n✅ Meta-analysis complete!\n")
  
} else {
  cat("❌ No successful meta-analysis results\n")
}