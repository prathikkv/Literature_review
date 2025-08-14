#!/usr/bin/env Rscript
#' Comprehensive 6-Dataset Analysis Pipeline
#' 
#' This script processes ALL 6 available datasets for complete CAMK gene family analysis
#' with quality controls, literature validation, and unified results generation

library(limma)
library(metafor)
library(tidyverse)

# Load functions
source("scripts/enhanced_group_detection_corrected.R")
source("functions/camk_definitions.R")
source("functions/analysis.R")

cat("═══════════════════════════════════════════════════════════════\n")
cat("🎯 COMPREHENSIVE 6-DATASET CAMK ANALYSIS PIPELINE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Complete dataset configuration
all_datasets_config <- list(
  "GSE57338" = list(
    description = "Heart failure (DCM + Ischemic) vs non-failing hearts",
    disease_type = "Heart failure (DCM + Ischemic)",
    control_type = "Non-failing hearts",
    biological_context = "Ventricular heart failure", 
    priority = "HIGH",
    expected_samples = 313,
    inclusion_reason = "Large sample size, clear disease vs healthy comparison"
  ),
  "GSE41177" = list(
    description = "Atrial fibrillation vs sinus rhythm",
    disease_type = "Atrial fibrillation",
    control_type = "Sinus rhythm",
    biological_context = "Atrial arrhythmia",
    priority = "MODERATE",
    expected_samples = 38,
    inclusion_reason = "AF vs SR comparison, moderate sample size"
  ),
  "GSE79768" = list(
    description = "Atrial fibrillation vs sinus rhythm", 
    disease_type = "Atrial fibrillation",
    control_type = "Sinus rhythm",
    biological_context = "Atrial arrhythmia",
    priority = "MODERATE", 
    expected_samples = 26,
    inclusion_reason = "AF vs SR comparison, supports GSE41177 findings"
  ),
  "GSE115574" = list(
    description = "Atrial fibrillation vs sinus rhythm",
    disease_type = "Atrial fibrillation",
    control_type = "Sinus rhythm", 
    biological_context = "Atrial arrhythmia",
    priority = "MODERATE",
    expected_samples = 59,
    inclusion_reason = "Largest AF dataset, good statistical power"
  ),
  "GSE31821" = list(
    description = "Atrial fibrillation vs sinus rhythm",
    disease_type = "Atrial fibrillation", 
    control_type = "Sinus rhythm",
    biological_context = "Atrial arrhythmia",
    priority = "LOW",
    expected_samples = 6,
    inclusion_reason = "Small but supports AF pattern, literature validation"
  ),
  "GSE14975" = list(
    description = "Heart failure vs healthy",
    disease_type = "Heart failure",
    control_type = "Healthy controls",
    biological_context = "Heart failure",
    priority = "LOW", 
    expected_samples = 10,
    inclusion_reason = "Small but supports heart failure pattern"
  )
)

# Get CAMK gene definitions
camk_core_genes <- get_camk_gene_categories()$core

cat("📋 COMPLETE ANALYSIS CONFIGURATION:\n")
cat("Total datasets to analyze:", length(all_datasets_config), "\n")
cat("CAMK genes of interest:", length(camk_core_genes), "\n") 
cat("Core CAMK genes:", paste(camk_core_genes, collapse = ", "), "\n\n")

#' Enhanced Dataset Analysis with Quality Controls
analyze_dataset_comprehensive <- function(dataset_id, config) {
  
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📊 ANALYZING:", dataset_id, "(Priority:", config$priority, ")\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  cat("🔬 Dataset Configuration:\n")
  cat("Disease type:", config$disease_type, "\n")
  cat("Control type:", config$control_type, "\n")
  cat("Biological context:", config$biological_context, "\n") 
  cat("Expected samples:", config$expected_samples, "\n")
  cat("Inclusion reason:", config$inclusion_reason, "\n\n")
  
  # Find dataset file
  dataset_files <- list.files("cache", pattern = paste0(dataset_id, "_processed.rds"),
                             recursive = TRUE, full.names = TRUE)
  
  if (length(dataset_files) == 0) {
    cat("❌ ERROR: Dataset file not found for", dataset_id, "\n")
    cat("Available files:", list.files("cache", pattern = ".*_processed.rds", recursive = TRUE), "\n\n")
    return(list(success = FALSE, reason = "file_not_found"))
  }
  
  # Load dataset  
  dataset <- readRDS(dataset_files[1])
  cat("📁 Loading from:", dataset_files[1], "\n")
  
  if (!dataset$success) {
    cat("❌ ERROR: Dataset processing failed\n\n")
    return(list(success = FALSE, reason = "processing_failed"))
  }
  
  # Dataset overview
  cat("📈 Dataset Overview:\n")
  cat("Total samples:", dataset$n_samples, "\n")
  cat("Total genes:", dataset$n_genes, "\n")
  cat("Platform:", dataset$dataset_info$platform, "\n\n")
  
  # Group detection
  cat("🎯 GROUP DETECTION:\n")
  detected_groups <- enhanced_auto_detect_groups_corrected(dataset)
  
  if (is.null(detected_groups)) {
    cat("❌ ERROR: Group detection failed\n\n")
    return(list(success = FALSE, reason = "group_detection_failed"))
  }
  
  cat("✅ Groups successfully detected:\n")
  cat("Pattern type:", detected_groups$pattern_type, "\n")
  cat("Reference group (Control):", levels(detected_groups$groups)[1], "\n") 
  cat("Comparison group (Disease):", levels(detected_groups$groups)[2], "\n")
  
  group_counts <- table(detected_groups$groups)
  print(group_counts)
  cat("\n")
  
  # Verify sample size matches expectation
  actual_samples <- sum(group_counts)
  if (abs(actual_samples - config$expected_samples) > 5) {
    cat("⚠️  WARNING: Sample size mismatch. Expected:", config$expected_samples, 
        "Actual:", actual_samples, "\n")
  }
  
  # Prepare expression data
  expr_matrix <- dataset$expression_matrix
  
  if (!is.null(detected_groups$sample_indices)) {
    expr_filtered <- expr_matrix[, detected_groups$sample_indices, drop = FALSE] 
    groups_vector <- detected_groups$groups
  } else {
    expr_filtered <- expr_matrix
    groups_vector <- detected_groups$groups
  }
  
  cat("📊 Analysis Matrix:", nrow(expr_filtered), "genes x", ncol(expr_filtered), "samples\n")
  
  # CAMK genes presence check
  camk_present <- intersect(rownames(expr_filtered), camk_core_genes)
  cat("🧬 CAMK genes detected:", length(camk_present), "/", length(camk_core_genes), "\n")
  cat("CAMK genes:", paste(camk_present, collapse = ", "), "\n\n")
  
  # Minimum gene threshold for inclusion
  min_genes_threshold <- if (config$priority == "HIGH") 8 else 5
  
  if (length(camk_present) < min_genes_threshold) {
    cat("⚠️  WARNING: Insufficient CAMK genes detected for", config$priority, "priority analysis\n")
    return(list(
      success = TRUE,
      insufficient_genes = TRUE,
      genes_detected = length(camk_present), 
      threshold = min_genes_threshold,
      camk_genes = camk_present
    ))
  }
  
  # DIFFERENTIAL EXPRESSION ANALYSIS
  cat("🧮 DIFFERENTIAL EXPRESSION ANALYSIS:\n")
  cat("Design: Control (reference) vs Disease (comparison)\n")
  cat("Interpretation: Positive logFC = UP in Disease\n")
  
  # Create design matrix
  design <- model.matrix(~ groups_vector)
  colnames(design) <- c("Intercept", "Disease_vs_Control")
  
  # Fit linear model
  fit <- lmFit(expr_filtered, design)
  fit <- eBayes(fit)
  
  # Get results for all genes
  all_results <- topTable(fit, coef = "Disease_vs_Control", number = Inf, adjust.method = "BH")
  
  cat("📋 Analysis Results:\n")
  cat("Total genes analyzed:", nrow(all_results), "\n")
  cat("Significant genes (FDR < 0.05):", sum(all_results$adj.P.Val < 0.05), "\n\n")
  
  # Extract CAMK gene results
  camk_results <- all_results[intersect(rownames(all_results), camk_present), , drop = FALSE]
  
  if (nrow(camk_results) > 0) {
    # Add comprehensive metadata
    camk_results$Gene_Symbol <- rownames(camk_results)
    camk_results$Dataset <- dataset_id
    camk_results$Biological_Context <- config$biological_context
    camk_results$Disease_Type <- config$disease_type
    camk_results$Control_Type <- config$control_type
    camk_results$Priority <- config$priority
    camk_results$Sample_Size <- ncol(expr_filtered)
    camk_results$Platform <- dataset$dataset_info$platform
    camk_results$Significant <- camk_results$adj.P.Val < 0.05
    camk_results$Regulation <- ifelse(camk_results$logFC > 0, "UP in Disease", "DOWN in Disease")
    
    # Quality assessment
    camk_results$Quality_Flag <- ifelse(abs(camk_results$logFC) > 0.8, "HIGH_LOGFC", 
                                      ifelse(abs(camk_results$logFC) < 0.01, "LOW_EFFECT", "GOOD"))
    
    # Sort by significance
    camk_results <- camk_results[order(camk_results$P.Value), ]
    
    cat("🧬 CAMK GENE RESULTS:\n")
    cat("CAMK genes in results:", nrow(camk_results), "\n")
    cat("Significant CAMK genes:", sum(camk_results$Significant), "\n")
    cat("Quality flags:", table(camk_results$Quality_Flag), "\n\n")
    
    # Display key results
    for (i in 1:nrow(camk_results)) {
      gene <- camk_results$Gene_Symbol[i]
      logfc <- round(camk_results$logFC[i], 4)
      pval <- camk_results$P.Value[i] 
      adj_pval <- camk_results$adj.P.Val[i]
      direction <- ifelse(logfc > 0, "UP", "DOWN")
      sig_status <- if (adj_pval < 0.05) "SIG" else "NS"
      quality <- camk_results$Quality_Flag[i]
      
      cat(sprintf("  🧬 %-8s: %4s logFC=%7.4f, p=%8.2e, FDR=%8.2e [%s] {%s}\n",
                  gene, direction, logfc, pval, adj_pval, sig_status, quality))
    }
    
    # CAMK2D spotlight
    camk2d_result <- camk_results[camk_results$Gene_Symbol == "CAMK2D", ]
    if (nrow(camk2d_result) > 0) {
      cat("\n⭐ CAMK2D SPOTLIGHT:\n")
      cat(sprintf("   Direction: %s\n", camk2d_result$Regulation))
      cat(sprintf("   LogFC: %7.4f\n", camk2d_result$logFC))
      cat(sprintf("   P-value: %8.2e\n", camk2d_result$P.Value))
      cat(sprintf("   FDR: %8.2e\n", camk2d_result$adj.P.Val))
      cat(sprintf("   Significant: %s\n", if (camk2d_result$Significant) "YES ✅" else "NO"))
      cat(sprintf("   Quality: %s\n", camk2d_result$Quality_Flag))
    }
    
    cat("\n✅ SUCCESS:", dataset_id, "analysis complete\n\n")
    
    return(list(
      success = TRUE,
      results = camk_results,
      dataset_info = list(
        id = dataset_id,
        samples = ncol(expr_filtered),
        genes_total = nrow(all_results),
        camk_genes = nrow(camk_results),
        significant_genes = sum(all_results$adj.P.Val < 0.05),
        significant_camk = sum(camk_results$Significant),
        priority = config$priority,
        platform = dataset$dataset_info$platform
      )
    ))
    
  } else {
    cat("❌ ERROR: No CAMK genes found in results\n\n")
    return(list(success = FALSE, reason = "no_camk_genes"))
  }
}

# EXECUTE COMPREHENSIVE ANALYSIS
cat("═══════════════════════════════════════════════════════════════\n")
cat("🚀 PHASE 1: COMPREHENSIVE DATASET ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

all_results <- list()
dataset_summary <- data.frame()
failed_datasets <- character(0)

for (dataset_id in names(all_datasets_config)) {
  config <- all_datasets_config[[dataset_id]]
  cat("Processing", dataset_id, "...\n")
  
  result <- analyze_dataset_comprehensive(dataset_id, config)
  
  if (result$success && !is.null(result$results)) {
    all_results[[dataset_id]] <- result$results
    
    # Add to summary
    dataset_summary <- rbind(dataset_summary, data.frame(
      Dataset = dataset_id,
      Priority = config$priority,
      Samples = result$dataset_info$samples,
      Total_Genes = result$dataset_info$genes_total, 
      CAMK_Genes = result$dataset_info$camk_genes,
      Significant_CAMK = result$dataset_info$significant_camk,
      Platform = result$dataset_info$platform,
      Disease_Type = config$disease_type,
      Biological_Context = config$biological_context,
      Status = "SUCCESS"
    ))
    
  } else {
    failed_datasets <- c(failed_datasets, dataset_id)
    
    # Add failure to summary
    reason <- if (result$success && result$insufficient_genes) {
      paste("Insufficient genes (", result$genes_detected, "/", result$threshold, ")", sep="")
    } else {
      result$reason
    }
    
    dataset_summary <- rbind(dataset_summary, data.frame(
      Dataset = dataset_id,
      Priority = config$priority,
      Samples = config$expected_samples,
      Total_Genes = NA,
      CAMK_Genes = if (exists("result") && !is.null(result$genes_detected)) result$genes_detected else NA,
      Significant_CAMK = NA,
      Platform = NA, 
      Disease_Type = config$disease_type,
      Biological_Context = config$biological_context,
      Status = paste("FAILED:", reason)
    ))
  }
}

# COMBINE ALL SUCCESSFUL RESULTS
if (length(all_results) > 0) {
  combined_results <- do.call(rbind, all_results)
  rownames(combined_results) <- NULL
  
  # Save comprehensive results
  write.csv(combined_results, "output/CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv", row.names = FALSE)
  write.csv(dataset_summary, "output/dataset_processing_summary_6_datasets.csv", row.names = FALSE)
  
  cat("📁 Results saved:\n")
  cat("- Individual results: output/CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv\n")
  cat("- Dataset summary: output/dataset_processing_summary_6_datasets.csv\n\n")
  
  # SUMMARY STATISTICS
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📊 COMPREHENSIVE ANALYSIS SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  cat("✅ SUCCESSFUL DATASETS:\n")
  successful_datasets <- dataset_summary[dataset_summary$Status == "SUCCESS", ]
  for (i in 1:nrow(successful_datasets)) {
    ds <- successful_datasets[i, ]
    cat(sprintf("  📊 %-10s: %s priority, %3d samples, %2d CAMK genes, %2d significant\n",
                ds$Dataset, ds$Priority, ds$Samples, ds$CAMK_Genes, ds$Significant_CAMK))
  }
  
  if (length(failed_datasets) > 0) {
    cat("\n❌ FAILED DATASETS:\n")
    failed_summary <- dataset_summary[dataset_summary$Status != "SUCCESS", ]
    for (i in 1:nrow(failed_summary)) {
      ds <- failed_summary[i, ]
      cat(sprintf("  ❌ %-10s: %s\n", ds$Dataset, ds$Status))
    }
  }
  
  # CAMK2D cross-dataset consistency
  camk2d_results <- combined_results[combined_results$Gene_Symbol == "CAMK2D", ]
  cat("\n⭐ CAMK2D CROSS-DATASET CONSISTENCY:\n")
  upregulated_count <- sum(camk2d_results$logFC > 0)
  significant_count <- sum(camk2d_results$Significant)
  
  for (i in 1:nrow(camk2d_results)) {
    ds <- camk2d_results[i, ]
    direction <- ifelse(ds$logFC > 0, "UP", "DOWN")
    sig <- ifelse(ds$Significant, "SIG", "NS")
    cat(sprintf("  ⭐ %-10s: %s (logFC=%7.4f, p=%8.2e) [%s]\n",
                ds$Dataset, direction, ds$logFC, ds$P.Value, sig))
  }
  
  cat(sprintf("\n🎯 CAMK2D Summary: %d/%d datasets UP, %d/%d significant\n",
              upregulated_count, nrow(camk2d_results), significant_count, nrow(camk2d_results)))
  
  consistency_score <- upregulated_count / nrow(camk2d_results) * 100
  if (consistency_score >= 80) {
    cat("🏆 EXCELLENT: High cross-dataset consistency (", round(consistency_score, 1), "%)\n", sep="")
  } else if (consistency_score >= 60) {
    cat("✅ GOOD: Moderate cross-dataset consistency (", round(consistency_score, 1), "%)\n", sep="")
  } else {
    cat("⚠️  MIXED: Lower cross-dataset consistency (", round(consistency_score, 1), "%)\n", sep="")
  }
  
  cat("\n✅ PHASE 1 COMPLETE: Ready for meta-analysis\n")
  
} else {
  cat("❌ CRITICAL ERROR: No datasets successfully processed\n")
  quit(status = 1)
}

cat("\n🔄 Next: Run comprehensive meta-analysis pipeline\n")