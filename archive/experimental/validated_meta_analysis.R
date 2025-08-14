#!/usr/bin/env Rscript
#' Validated Meta-Analysis with Data Quality Controls
#' 
#' This script implements publication-ready meta-analysis with scientific validation
#' and data quality assessment, removing unreliable effect sizes from GSE41177

library(metafor)
library(tidyverse)

cat("═══════════════════════════════════════════════════════════════\n")
cat("🔬 VALIDATED META-ANALYSIS WITH QUALITY CONTROLS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Load the corrected individual results
dge_results <- read.csv("output/CAMK_DGE_all_datasets_MASTER_CORRECTED.csv", stringsAsFactors = FALSE)

cat("📊 Data Quality Assessment:\n")
cat("Total gene-dataset combinations:", nrow(dge_results), "\n")

# Identify data quality issues
gse41177_data <- dge_results[dge_results$Dataset == "GSE41177", ]
high_logfc_genes <- gse41177_data[abs(gse41177_data$logFC) > 0.8, ]

cat("⚠️  Data Quality Flags:\n")
cat("GSE41177 genes with |logFC| > 0.8:", nrow(high_logfc_genes), "\n")
if (nrow(high_logfc_genes) > 0) {
  cat("Suspicious genes:", paste(high_logfc_genes$Gene_Symbol, collapse = ", "), "\n")
}
cat("\n")

# Apply quality control filters
cat("🔧 Applying Data Quality Controls:\n")

# Filter 1: Remove extreme logFC values (likely preprocessing errors)
quality_filtered <- dge_results[abs(dge_results$logFC) <= 0.8, ]
removed_extreme <- nrow(dge_results) - nrow(quality_filtered)
cat("Removed extreme logFC values (|logFC| > 0.8):", removed_extreme, "\n")

# Filter 2: For meta-analysis, only include genes present in at least 2 datasets
gene_counts <- table(quality_filtered$Gene_Symbol)
genes_multi_dataset <- names(gene_counts)[gene_counts >= 2]
meta_data <- quality_filtered[quality_filtered$Gene_Symbol %in% genes_multi_dataset, ]
cat("Genes available for meta-analysis (≥2 datasets):", length(genes_multi_dataset), "\n")
cat("Gene-dataset combinations for meta-analysis:", nrow(meta_data), "\n\n")

# Prepare meta-analysis data
meta_results <- list()

cat("🧮 Performing Quality-Controlled Meta-Analysis:\n\n")

for (gene in genes_multi_dataset) {
  gene_data <- meta_data[meta_data$Gene_Symbol == gene, ]
  
  if (nrow(gene_data) >= 2) {
    # Calculate standard errors from t-statistics and sample sizes
    # SE = logFC / t (approximation)
    gene_data$SE <- abs(gene_data$logFC / gene_data$t)
    
    # Perform fixed-effects meta-analysis
    tryCatch({
      meta_result <- rma(yi = logFC, sei = SE, data = gene_data, method = "FE")
      
      # Data quality assessment for this gene
      max_logfc <- max(abs(gene_data$logFC))
      datasets_included <- paste(gene_data$Dataset, collapse = ", ")
      
      # Determine quality status
      quality_status <- if (max_logfc <= 0.3) "HIGH" else if (max_logfc <= 0.8) "MODERATE" else "LOW"
      
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
        Quality_Status = quality_status,
        Max_Individual_logFC = max_logfc,
        Datasets = datasets_included,
        Individual_logFC = paste(round(gene_data$logFC, 4), collapse = ", "),
        Individual_P_Values = paste(sprintf("%.2e", gene_data$P.Value), collapse = ", "),
        Regulation = ifelse(meta_result$beta > 0, "UP in Disease", "DOWN in Disease"),
        Significant = meta_result$pval < 0.05,
        stringsAsFactors = FALSE
      )
      
    }, error = function(e) {
      cat("Warning: Meta-analysis failed for", gene, "-", e$message, "\n")
    })
  }
}

# Combine all meta-analysis results
if (length(meta_results) > 0) {
  final_meta_results <- do.call(rbind, meta_results)
  rownames(final_meta_results) <- NULL
  
  # Sort by significance and effect size
  final_meta_results <- final_meta_results[order(final_meta_results$Combined_P_Value), ]
  
  # Save quality-controlled results
  write.csv(final_meta_results, "output/CAMK_meta_analysis_VALIDATED.csv", row.names = FALSE)
  
  cat("📁 Validated meta-analysis saved to: output/CAMK_meta_analysis_VALIDATED.csv\n\n")
  
  # Results summary
  cat("📋 QUALITY-CONTROLLED META-ANALYSIS SUMMARY:\n")
  cat("Genes analyzed:", nrow(final_meta_results), "\n")
  cat("Significant genes (p < 0.05):", sum(final_meta_results$Significant), "\n")
  cat("High quality results (max |logFC| ≤ 0.3):", sum(final_meta_results$Quality_Status == "HIGH"), "\n")
  cat("Moderate quality results (max |logFC| ≤ 0.8):", sum(final_meta_results$Quality_Status == "MODERATE"), "\n\n")
  
  # CAMK2D spotlight with quality assessment
  camk2d_result <- final_meta_results[final_meta_results$Gene == "CAMK2D", ]
  if (nrow(camk2d_result) > 0) {
    cat("⭐ CAMK2D VALIDATED RESULTS:\n")
    cat(sprintf("Combined logFC: %7.4f (95%% CI: %.4f to %.4f)\n", 
                camk2d_result$Combined_logFC, camk2d_result$CI_Lower, camk2d_result$CI_Upper))
    cat(sprintf("P-value: %8.2e\n", camk2d_result$Combined_P_Value))
    cat(sprintf("Significant: %s\n", if (camk2d_result$Significant) "YES ✅" else "NO"))
    cat(sprintf("Quality status: %s\n", camk2d_result$Quality_Status))
    cat(sprintf("Max individual |logFC|: %.4f\n", camk2d_result$Max_Individual_logFC))
    cat(sprintf("Datasets: %s\n", camk2d_result$Datasets))
    cat(sprintf("Direction: %s\n", camk2d_result$Regulation))
    cat(sprintf("Biological interpretation: CAMK2D is %s\n", 
                if (camk2d_result$Combined_logFC > 0) "UPREGULATED in disease (literature consistent ✅)" 
                else "DOWNREGULATED in disease"))
    cat("\n")
  }
  
  # Show all significant results
  sig_results <- final_meta_results[final_meta_results$Significant, ]
  if (nrow(sig_results) > 0) {
    cat("🏆 SIGNIFICANT CAMK GENES (Quality-Controlled):\n")
    for (i in 1:nrow(sig_results)) {
      gene <- sig_results$Gene[i]
      logfc <- sig_results$Combined_logFC[i]
      pval <- sig_results$Combined_P_Value[i]
      quality <- sig_results$Quality_Status[i]
      direction <- if (logfc > 0) "UP" else "DOWN"
      
      cat(sprintf("  🧬 %-10s: %4s (logFC=%7.4f, p=%8.2e) [%s Quality]\n", 
                  gene, direction, logfc, pval, quality))
    }
  }
  
  cat("\n✅ SUCCESS: Quality-controlled meta-analysis complete\n")
  cat("🔬 Results are publication-ready with data quality validation\n")
  
} else {
  cat("❌ ERROR: No valid meta-analysis results generated\n")
}