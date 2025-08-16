#!/usr/bin/env Rscript
#' Step 04: Meta-Analysis
#' 
#' Performs fixed-effects meta-analysis on differential expression results
#' Extracts meta-analysis logic from fixed_meta_analysis.R

source("scripts/utilities/step_interface.R")
library(metafor)
library(tidyverse)

#' Execute Meta-Analysis Step
#'
#' @param step_name Name of this step (should be "step_04_meta_analysis")
#' @param input_data Output from step_03_dge_analysis
#' @param config Full pipeline configuration
#' @param checkpoint_dir Directory for saving checkpoints
#' @return Step result with meta-analysis results
step_04_meta_analysis <- function(step_name, input_data, config, checkpoint_dir = "output/checkpoints") {
  
  cat("🔬 STEP 04: META-ANALYSIS\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Performing fixed-effects meta-analysis...\n\n")
  
  # Validate input from DGE analysis
  validate_step_input(
    input_data = input_data,
    required_fields = c("dge_results", "dge_summary", "analysis_stats"),
    step_name = step_name
  )
  
  dge_results <- input_data$dge_results
  dge_summary <- input_data$dge_summary
  analysis_stats <- input_data$analysis_stats
  
  if (nrow(dge_results) == 0) {
    return(create_step_result(
      success = FALSE,
      error_message = "No DGE results provided for meta-analysis",
      step_name = step_name
    ))
  }
  
  # Get meta-analysis configuration
  meta_config <- config$analysis$meta_analysis
  quality_filters <- meta_config$quality_filters
  target_genes <- config$gene_analysis$gene_family
  
  cat("📋 META-ANALYSIS CONFIGURATION:\n")
  cat("Method:", meta_config$method %||% "fixed_effects", "\n")
  cat("Package:", meta_config$package %||% "metafor", "\n")
  cat("Min datasets per gene:", meta_config$min_datasets %||% 2, "\n")
  cat("Max absolute logFC filter:", quality_filters$max_absolute_logfc %||% 0.8, "\n")
  cat("Significance threshold:", meta_config$significance_threshold %||% 0.05, "\n\n")
  
  cat("📊 INPUT DATA SUMMARY:\n")
  cat("Total DGE results:", nrow(dge_results), "\n")
  cat("Datasets analyzed:", analysis_stats$total_datasets_analyzed, "\n")
  cat("Genes to consider:", length(target_genes), "\n\n")
  
  # Apply quality filters
  cat("🔍 APPLYING QUALITY FILTERS:\n")
  
  # Filter out extreme logFC values if configured
  max_logfc <- quality_filters$max_absolute_logfc %||% 0.8
  quality_filtered <- dge_results[abs(dge_results$logFC) <= max_logfc, ]
  
  filtered_out <- nrow(dge_results) - nrow(quality_filtered)
  if (filtered_out > 0) {
    cat("🚫 Filtered out", filtered_out, "results with |logFC| >", max_logfc, "\n")
  }
  
  # Count genes across datasets
  gene_counts <- table(quality_filtered$Gene_Symbol)
  min_datasets <- meta_config$min_datasets %||% 2
  genes_multi_dataset <- names(gene_counts)[gene_counts >= min_datasets]
  
  cat("📊 Genes with ≥", min_datasets, "datasets:", length(genes_multi_dataset), "\n")
  
  if (length(genes_multi_dataset) == 0) {
    return(create_step_result(
      success = FALSE,
      error_message = "No genes have sufficient datasets for meta-analysis",
      step_name = step_name
    ))
  }
  
  # Prepare data for meta-analysis
  meta_data <- quality_filtered[quality_filtered$Gene_Symbol %in% genes_multi_dataset, ]
  
  cat("📊 Data prepared for meta-analysis:\n")
  cat("Genes:", length(genes_multi_dataset), "\n")
  cat("Gene-dataset combinations:", nrow(meta_data), "\n\n")
  
  # Perform meta-analysis for each gene
  cat("🔬 PERFORMING META-ANALYSIS:\n")
  
  meta_results <- list()
  meta_warnings <- character(0)
  
  for (gene in genes_multi_dataset) {
    gene_data <- meta_data[meta_data$Gene_Symbol == gene, ]
    
    if (nrow(gene_data) >= min_datasets) {
      
      tryCatch({
        # Calculate standard errors from t-statistics
        gene_data$SE <- abs(gene_data$logFC / gene_data$t)
        
        # Perform fixed-effects meta-analysis
        meta_result <- rma(yi = logFC, sei = SE, data = gene_data, method = "FE")
        
        # Extract comprehensive results
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
          Significant = meta_result$pval < (meta_config$significance_threshold %||% 0.05),
          stringsAsFactors = FALSE
        )
        
        cat("✅", gene, "meta-analysis complete\n")
        
        # Check for high heterogeneity
        if (!is.na(meta_result$I2) && meta_result$I2 > 75) {
          warning_msg <- paste(gene, "shows high heterogeneity (I² =", round(meta_result$I2, 1), "%)")
          meta_warnings <- c(meta_warnings, warning_msg)
        }
        
      }, error = function(e) {
        cat("❌", gene, "failed:", e$message, "\n")
        warning_msg <- paste(gene, "meta-analysis failed:", e$message)
        meta_warnings <- c(meta_warnings, warning_msg)
      })
    }
  }
  
  # Combine and process results
  if (length(meta_results) == 0) {
    return(create_step_result(
      success = FALSE,
      error_message = "No successful meta-analysis results generated",
      step_name = step_name,
      warnings = meta_warnings
    ))
  }
  
  final_results <- do.call(rbind, meta_results)
  rownames(final_results) <- NULL
  
  # Sort by significance
  final_results <- final_results[order(final_results$Combined_P_Value), ]
  
  cat("\n📊 META-ANALYSIS RESULTS SUMMARY:\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  genes_analyzed <- nrow(final_results)
  significant_genes <- sum(final_results$Significant)
  significance_threshold <- meta_config$significance_threshold %||% 0.05
  
  cat("Genes successfully analyzed:", genes_analyzed, "\n")
  cat("Significant genes (p <", significance_threshold, "):", significant_genes, "\n")
  
  if (length(meta_warnings) > 0) {
    cat("Meta-analysis warnings:", length(meta_warnings), "\n")
  }
  
  # Primary gene spotlight
  primary_gene <- config$research_target$primary_gene
  primary_gene_results <- final_results[final_results$Gene == primary_gene, ]
  if (nrow(primary_gene_results) > 0) {
    cat("\n⭐", toupper(primary_gene), "META-ANALYSIS RESULTS:\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    
    primary_data <- primary_gene_results[1, ]
    cat(sprintf("Combined logFC: %.4f (95%% CI: %.4f to %.4f)\n",
                primary_data$Combined_logFC, primary_data$CI_Lower, primary_data$CI_Upper))
    cat(sprintf("P-value: %.2e\n", primary_data$Combined_P_Value))
    cat(sprintf("Significant: %s\n", if (primary_data$Significant) "YES ✅" else "NO"))
    cat(sprintf("Direction: %s\n", primary_data$Regulation))
    cat(sprintf("Datasets included: %s\n", primary_data$Datasets))
    cat(sprintf("Individual logFC values: %s\n", primary_data$Individual_logFC))
    
    if (primary_data$Heterogeneity_I2 > 0) {
      cat(sprintf("Heterogeneity (I²): %.1f%% (p = %.3f)\n", primary_data$Heterogeneity_I2, primary_data$Heterogeneity_P))
    }
  }
  
  # Display all significant genes
  sig_results <- final_results[final_results$Significant, ]
  if (nrow(sig_results) > 0) {
    cat("\n🏆 ALL SIGNIFICANT GENES:\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    
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
  
  # Validation against baseline (if configured)
  validation_config <- config$validation$baseline_comparison
  validation_passed <- TRUE
  validation_messages <- character(0)
  
  if (!is.null(validation_config) && validation_config$enabled && nrow(primary_gene_results) > 0) {
    cat("\n🔍 BASELINE VALIDATION:\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    
    baseline <- validation_config$baseline_results
    tolerances <- config$validation$tolerances
    
    # Check primary gene logFC
    primary_gene_lower <- tolower(primary_gene)
    if (!is.null(baseline[[paste0(primary_gene_lower, "_combined_logfc")]])) {
      logfc_diff <- abs(primary_data$Combined_logFC - baseline[[paste0(primary_gene_lower, "_combined_logfc")]])
      logfc_tolerance <- tolerances$logfc_tolerance %||% 1e-4
      
      if (logfc_diff <= logfc_tolerance) {
        cat("✅", primary_gene, "logFC matches baseline (diff =", sprintf("%.6f", logfc_diff), ")\n")
      } else {
        validation_passed <- FALSE
        msg <- paste(primary_gene, "logFC differs from baseline:", sprintf("%.6f", logfc_diff), "> tolerance:", logfc_tolerance)
        validation_messages <- c(validation_messages, msg)
        cat("❌", msg, "\n")
      }
    }
    
    # Check primary gene p-value
    if (!is.null(baseline[[paste0(primary_gene_lower, "_p_value")]])) {
      pval_diff <- abs(primary_data$Combined_P_Value - baseline[[paste0(primary_gene_lower, "_p_value")]])
      pval_tolerance <- tolerances$p_value_tolerance %||% 1e-6
      
      if (pval_diff <= pval_tolerance) {
        cat("✅", primary_gene, "p-value matches baseline (diff =", sprintf("%.2e", pval_diff), ")\n")
      } else {
        validation_passed <- FALSE
        msg <- paste(primary_gene, "p-value differs from baseline:", sprintf("%.2e", pval_diff), "> tolerance:", pval_tolerance)
        validation_messages <- c(validation_messages, msg)
        cat("❌", msg, "\n")
      }
    }
    
    # Check total significant genes
    if (!is.null(baseline$significant_genes)) {
      if (significant_genes == baseline$significant_genes) {
        cat("✅ Significant gene count matches baseline (", significant_genes, ")\n")
      } else {
        validation_passed <- FALSE
        msg <- paste("Significant genes differ from baseline:", significant_genes, "vs", baseline$significant_genes)
        validation_messages <- c(validation_messages, msg)
        cat("❌", msg, "\n")
      }
    }
    
    if (validation_passed) {
      cat("🎉 BASELINE VALIDATION PASSED\n")
    } else {
      cat("⚠️  BASELINE VALIDATION WARNINGS\n")
      meta_warnings <- c(meta_warnings, validation_messages)
    }
  }
  
  # Save results
  output_files <- config$paths$output_files
  
  if (!is.null(output_files$meta_results)) {
    write.csv(final_results, output_files$meta_results, row.names = FALSE)
    cat("\n💾 Meta-analysis results saved to:", output_files$meta_results, "\n")
    
    # Also save to legacy location for compatibility
    if (!is.null(output_files$legacy_meta)) {
      write.csv(final_results, output_files$legacy_meta, row.names = FALSE)
    }
  }
  
  # Create output data structure
  output_data <- list(
    meta_results = final_results,
    meta_summary = list(
      genes_analyzed = genes_analyzed,
      significant_genes = significant_genes,
      significance_threshold = significance_threshold,
      primary_gene_significant = nrow(primary_gene_results) > 0 && primary_gene_results$Significant[1],
      baseline_validation_passed = validation_passed
    ),
    dge_results = dge_results,  # Pass through individual gene results for reporting
    dge_summary = dge_summary,  # Pass through dataset summary for reporting
    analysis_stats = analysis_stats,
    quality_filtering = list(
      original_results = nrow(dge_results),
      after_quality_filter = nrow(quality_filtered),
      genes_with_sufficient_data = length(genes_multi_dataset)
    )
  )
  
  # SAVE OUTPUT FILES
  cat("\n💾 SAVING OUTPUT FILES:\n")
  
  # Ensure output directory exists
  output_dir <- "output/current"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("📁 Created directory:", output_dir, "\n")
  }
  
  # Save meta-analysis results with dynamic filename
  primary_gene_upper <- toupper(primary_gene)
  meta_csv_path <- paste0("output/current/", primary_gene_upper, "_meta_analysis_FINAL.csv")
  write.csv(meta_results, meta_csv_path, row.names = FALSE)
  cat("✅ Meta-analysis results saved:", meta_csv_path, "\n")
  
  # Save DGE summary with dynamic filename
  dge_csv_path <- paste0("output/current/", primary_gene_upper, "_DGE_all_datasets_COMPREHENSIVE.csv")
  write.csv(dge_results, dge_csv_path, row.names = FALSE)
  cat("✅ DGE results saved:", dge_csv_path, "\n")
  
  # Save dataset summary
  dataset_csv_path <- "output/current/dataset_processing_summary_6_datasets.csv"
  write.csv(dge_summary, dataset_csv_path, row.names = FALSE)
  cat("✅ Dataset summary saved:", dataset_csv_path, "\n")
  
  # Validate output
  validate_step_output(
    output_data = output_data,
    required_fields = c("meta_results", "meta_summary"),
    step_name = step_name
  )
  
  cat("\n🎉 META-ANALYSIS COMPLETED SUCCESSFULLY\n")
  cat("Ready for report generation...\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(create_step_result(
    success = TRUE,
    output_data = output_data,
    step_name = step_name,
    warnings = meta_warnings,
    metadata = list(
      genes_analyzed = genes_analyzed,
      significant_genes = significant_genes,
      primary_gene_significant = nrow(primary_gene_results) > 0 && primary_gene_results$Significant[1],
      validation_passed = validation_passed
    )
  ))
}

cat("✅ STEP 04: Meta-Analysis loaded successfully\n\n")