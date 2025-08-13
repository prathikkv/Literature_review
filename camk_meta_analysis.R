#!/usr/bin/env Rscript
#' CAMK Meta-Analysis Across All Datasets
#' 
#' Combines CAMK DGE results across all datasets using meta-analysis methods

library(metafor)
library(meta)
library(ggplot2)
library(openxlsx)

cat("🔬 CAMK META-ANALYSIS ACROSS ALL DATASETS\n")
cat("=========================================\n\n")

# Load CAMK-focused DGE results
if (!file.exists("output/CAMK_focused_DGE_all_datasets.csv")) {
  stop("❌ CAMK DGE results not found. Please run camk_focused_analysis.R first.")
}

camk_dge_data <- read.csv("output/CAMK_focused_DGE_all_datasets.csv", stringsAsFactors = FALSE)

cat("📊 Loaded CAMK DGE data:", nrow(camk_dge_data), "results across datasets\n")

# Summary of available data
cat("   Datasets:", length(unique(camk_dge_data$Dataset)), "\n")
cat("   CAMK genes:", length(unique(camk_dge_data$Gene_Symbol)), "\n")
cat("   Total comparisons:", nrow(camk_dge_data), "\n\n")

# Display dataset summary
dataset_summary <- table(camk_dge_data$Dataset)
cat("📋 Results per dataset:\n")
for (dataset in names(dataset_summary)) {
  cat(sprintf("   %-12s: %d CAMK genes\n", dataset, dataset_summary[dataset]))
}
cat("\n")

# Function to perform meta-analysis for each CAMK gene
perform_camk_gene_meta_analysis <- function(gene_data, gene_name) {
  
  if (nrow(gene_data) < 2) {
    cat("   ⚠️", gene_name, "- insufficient data for meta-analysis (n =", nrow(gene_data), ")\n")
    return(NULL)
  }
  
  # Calculate standard errors from t-statistics and sample sizes (approximate)
  gene_data$se <- abs(gene_data$logFC / gene_data$t)
  
  # Remove infinite or missing values
  valid_data <- gene_data[is.finite(gene_data$logFC) & is.finite(gene_data$se) & 
                         gene_data$se > 0 & gene_data$se < 10, ]
  
  if (nrow(valid_data) < 2) {
    cat("   ⚠️", gene_name, "- insufficient valid data after filtering\n")
    return(NULL)
  }
  
  tryCatch({
    # Perform random-effects meta-analysis
    meta_result <- rma(yi = valid_data$logFC, 
                       sei = valid_data$se, 
                       data = valid_data,
                       method = "REML")
    
    # Create summary
    summary_data <- list(
      gene = gene_name,
      n_studies = nrow(valid_data),
      combined_logFC = meta_result$beta[1,1],
      combined_se = meta_result$se,
      combined_pval = meta_result$pval,
      combined_ci_lower = meta_result$ci.lb,
      combined_ci_upper = meta_result$ci.ub,
      heterogeneity_I2 = meta_result$I2,
      heterogeneity_pval = meta_result$QEp,
      tau2 = meta_result$tau2,
      datasets = paste(valid_data$Dataset, collapse = ", "),
      individual_logFC = paste(round(valid_data$logFC, 3), collapse = ", "),
      individual_pvals = paste(signif(valid_data$adj.P.Val, 3), collapse = ", ")
    )
    
    cat(sprintf("   ✅ %-8s: Combined logFC = %6.3f ± %5.3f, P = %8.2e (%d studies)\n",
                gene_name, summary_data$combined_logFC, summary_data$combined_se, 
                summary_data$combined_pval, summary_data$n_studies))
    
    return(summary_data)
    
  }, error = function(e) {
    cat("   ❌", gene_name, "- meta-analysis error:", e$message, "\n")
    return(NULL)
  })
}

# Perform meta-analysis for each CAMK gene
cat("🧬 Performing meta-analysis for each CAMK gene:\n")
cat("===============================================\n")

camk_genes <- unique(camk_dge_data$Gene_Symbol)
meta_results <- list()

for (gene in camk_genes) {
  gene_data <- camk_dge_data[camk_dge_data$Gene_Symbol == gene, ]
  meta_result <- perform_camk_gene_meta_analysis(gene_data, gene)
  
  if (!is.null(meta_result)) {
    meta_results[[gene]] <- meta_result
  }
}

cat(sprintf("\n📊 Meta-analysis completed for %d/%d CAMK genes\n\n", length(meta_results), length(camk_genes)))

# Create comprehensive meta-analysis summary table
if (length(meta_results) > 0) {
  meta_summary_df <- do.call(rbind, lapply(meta_results, function(x) {
    data.frame(
      Gene = x$gene,
      N_Studies = x$n_studies,
      Combined_logFC = round(x$combined_logFC, 4),
      Combined_SE = round(x$combined_se, 4),
      Combined_P_Value = signif(x$combined_pval, 3),
      CI_Lower = round(x$combined_ci_lower, 4),
      CI_Upper = round(x$combined_ci_upper, 4),
      Heterogeneity_I2 = round(x$heterogeneity_I2, 1),
      Heterogeneity_P = signif(x$heterogeneity_pval, 3),
      Regulation = ifelse(x$combined_logFC > 0, "UP in Disease", "DOWN in Disease"),
      Significant = x$combined_pval < 0.05,
      Datasets = x$datasets,
      Individual_logFC = x$individual_logFC,
      Individual_P_Values = x$individual_pvals,
      stringsAsFactors = FALSE
    )
  }))
  
  # Sort by combined p-value
  meta_summary_df <- meta_summary_df[order(meta_summary_df$Combined_P_Value), ]
  
  cat("🎯 CAMK FAMILY META-ANALYSIS RESULTS\n")
  cat("===================================\n\n")
  
  # Display top results
  for (i in 1:nrow(meta_summary_df)) {
    gene <- meta_summary_df$Gene[i]
    logfc <- meta_summary_df$Combined_logFC[i]
    pval <- meta_summary_df$Combined_P_Value[i]
    n_studies <- meta_summary_df$N_Studies[i]
    regulation <- meta_summary_df$Regulation[i]
    significant <- meta_summary_df$Significant[i]
    i2 <- meta_summary_df$Heterogeneity_I2[i]
    
    sig_marker <- if (significant) {
      if (pval < 0.001) "***" else if (pval < 0.01) "**" else "*"
    } else ""
    
    direction <- if (logfc > 0) "↑" else "↓"
    
    cat(sprintf("%-8s: %s Meta-logFC = %6.3f, P = %8.2e %s [%s] (n=%d, I²=%4.1f%%)\n",
                gene, direction, logfc, pval, sig_marker, regulation, n_studies, i2))
  }
  
  # Summary statistics
  significant_genes <- sum(meta_summary_df$Significant)
  upregulated_genes <- sum(meta_summary_df$Combined_logFC > 0 & meta_summary_df$Significant)
  downregulated_genes <- sum(meta_summary_df$Combined_logFC < 0 & meta_summary_df$Significant)
  
  cat(sprintf("\n📈 META-ANALYSIS SUMMARY:\n"))
  cat("========================\n")
  cat(sprintf("CAMK genes analyzed: %d\n", nrow(meta_summary_df)))
  cat(sprintf("Significant genes (P < 0.05): %d (%d%%)\n", significant_genes, 
              round(significant_genes/nrow(meta_summary_df)*100)))
  cat(sprintf("Up-regulated in disease: %d\n", upregulated_genes))
  cat(sprintf("Down-regulated in disease: %d\n", downregulated_genes))
  cat(sprintf("Average studies per gene: %.1f\n", mean(meta_summary_df$N_Studies)))
  cat(sprintf("Average heterogeneity (I²): %.1f%%\n", mean(meta_summary_df$Heterogeneity_I2, na.rm = TRUE)))
  
  # Clinical interpretation
  cat("\n🏥 CLINICAL INTERPRETATION:\n")
  cat("==========================\n")
  
  if (significant_genes > 0) {
    most_significant <- meta_summary_df[1, ]
    cat(sprintf("Most significant CAMK gene: %s (P = %8.2e, %s)\n", 
                most_significant$Gene, most_significant$Combined_P_Value, 
                most_significant$Regulation))
    
    if (most_significant$Combined_logFC > 0) {
      cat("→ This suggests increased CAMK activity in cardiovascular disease\n")
    } else {
      cat("→ This suggests decreased CAMK activity in cardiovascular disease\n")
    }
    
    # Effect size interpretation
    abs_effect <- abs(most_significant$Combined_logFC)
    if (abs_effect > 1) {
      effect_desc <- "Large effect"
    } else if (abs_effect > 0.5) {
      effect_desc <- "Moderate effect"  
    } else {
      effect_desc <- "Small effect"
    }
    
    cat(sprintf("→ Effect size: %s (|logFC| = %.3f)\n", effect_desc, abs_effect))
  } else {
    cat("No CAMK genes showed statistically significant meta-analysis results.\n")
    cat("This may indicate:\n")
    cat("• Heterogeneity between studies\n") 
    cat("• Need for larger sample sizes\n")
    cat("• Disease-specific CAMK regulation patterns\n")
  }
  
} else {
  cat("❌ No successful meta-analyses - insufficient data\n")
}

# Save results
cat("\n💾 Saving results...\n")

# Save meta-analysis summary
if (exists("meta_summary_df")) {
  write.csv(meta_summary_df, "output/CAMK_meta_analysis_summary.csv", row.names = FALSE)
  
  # Create Excel workbook with multiple sheets
  wb <- createWorkbook()
  
  # Sheet 1: Meta-analysis summary
  addWorksheet(wb, "Meta_Analysis_Summary")
  writeData(wb, "Meta_Analysis_Summary", meta_summary_df)
  
  # Sheet 2: Individual study results
  addWorksheet(wb, "Individual_Studies")
  writeData(wb, "Individual_Studies", camk_dge_data)
  
  # Sheet 3: Study information
  study_info <- data.frame(
    Dataset = names(dataset_summary),
    N_CAMK_Genes = as.numeric(dataset_summary),
    Platform = c("Microarray", "RNA-seq", "Microarray", "Microarray", "Microarray")[1:length(dataset_summary)],
    Disease_Type = c("AF", "HF", "AF", "AF", "AF")[1:length(dataset_summary)]
  )
  
  addWorksheet(wb, "Study_Information")  
  writeData(wb, "Study_Information", study_info)
  
  saveWorkbook(wb, "output/CAMK_Comprehensive_Meta_Analysis.xlsx", overwrite = TRUE)
  
  cat("   ✅ output/CAMK_meta_analysis_summary.csv\n")
  cat("   ✅ output/CAMK_Comprehensive_Meta_Analysis.xlsx\n")
}

# Save detailed results
saveRDS(list(
  meta_results = meta_results,
  meta_summary = if(exists("meta_summary_df")) meta_summary_df else NULL,
  raw_data = camk_dge_data
), "output/CAMK_meta_analysis_complete_results.rds")

cat("   ✅ output/CAMK_meta_analysis_complete_results.rds\n")

cat("\n🎉 CAMK META-ANALYSIS COMPLETED!\n")
cat("✨ Comprehensive CAMK family analysis across multiple datasets complete\n")
cat("🏥 Results are ready for clinical interpretation and publication\n")