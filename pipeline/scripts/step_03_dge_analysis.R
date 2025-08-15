#!/usr/bin/env Rscript
#' Step 03: Differential Gene Expression Analysis
#' 
#' Performs differential expression analysis on preprocessed datasets
#' Extracts DGE logic from comprehensive_6_dataset_pipeline.R

source("scripts/utilities/step_interface.R")
library(limma)

#' Execute DGE Analysis Step
#'
#' @param step_name Name of this step (should be "step_03_dge_analysis")
#' @param input_data Output from step_02_preprocessing
#' @param config Full pipeline configuration
#' @param checkpoint_dir Directory for saving checkpoints
#' @return Step result with DGE results
step_03_dge_analysis <- function(step_name, input_data, config, checkpoint_dir = "output/checkpoints") {
  
  cat("🧮 STEP 03: DIFFERENTIAL GENE EXPRESSION ANALYSIS\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Performing differential expression analysis...\n\n")
  
  # Validate input from preprocessing
  validate_step_input(
    input_data = input_data,
    required_fields = c("processed_datasets", "preprocessing_summary"),
    step_name = step_name
  )
  
  processed_datasets <- input_data$processed_datasets
  preprocessing_summary <- input_data$preprocessing_summary
  
  if (length(processed_datasets) == 0) {
    return(create_step_result(
      success = FALSE,
      error_message = "No preprocessed datasets provided for DGE analysis",
      step_name = step_name
    ))
  }
  
  # Get analysis configuration
  de_config <- config$analysis$differential_expression
  qc_config <- config$analysis$quality_control
  camk_genes <- config$genes$camk_core_genes
  
  cat("📋 DGE ANALYSIS CONFIGURATION:\n")
  cat("Method:", de_config$method %||% "limma", "\n")
  cat("FDR threshold:", de_config$fdr_threshold %||% 0.05, "\n")
  cat("Adjustment method:", de_config$adjust_method %||% "BH", "\n")
  cat("CAMK genes of interest:", length(camk_genes), "\n")
  cat("Datasets to analyze:", length(processed_datasets), "\n\n")
  
  # Perform DGE analysis on each dataset
  dge_results <- list()
  dge_summary <- data.frame()
  failed_dge <- character(0)
  warnings_list <- character(0)
  
  for (dataset_id in names(processed_datasets)) {
    processed_data <- processed_datasets[[dataset_id]]
    dataset_config <- config$datasets$active_datasets[[dataset_id]]
    
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("🧮 ANALYZING:", dataset_id, "(Priority:", dataset_config$priority, ")\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    
    # Perform DGE analysis for this dataset
    dge_result <- analyze_dataset_dge(dataset_id, processed_data, dataset_config, config)
    
    if (dge_result$success) {
      dge_results[[dataset_id]] <- dge_result$results
      
      # Add to summary
      dge_summary <- rbind(dge_summary, data.frame(
        Dataset = dataset_id,
        Priority = dataset_config$priority,
        Samples = dge_result$sample_count,
        Total_Genes = dge_result$total_genes_analyzed,
        Significant_Genes = dge_result$significant_genes_total,
        CAMK_Genes_Analyzed = dge_result$camk_genes_analyzed,
        CAMK_Significant = dge_result$camk_genes_significant,
        Platform = processed_data$dataset_info$platform %||% "Unknown",
        Disease_Type = dataset_config$disease_type,
        Biological_Context = dataset_config$biological_context,
        Status = "SUCCESS",
        stringsAsFactors = FALSE
      ))
      
      cat("✅ SUCCESS: DGE analysis completed\n")
      cat("📊 Total genes analyzed:", dge_result$total_genes_analyzed, "\n")
      cat("📈 Significant genes (FDR < 0.05):", dge_result$significant_genes_total, "\n")
      cat("🧬 CAMK genes analyzed:", dge_result$camk_genes_analyzed, "\n")
      cat("⭐ CAMK genes significant:", dge_result$camk_genes_significant, "\n")
      
      # Display CAMK2D result if available
      camk2d_result <- dge_result$results[dge_result$results$Gene_Symbol == "CAMK2D", ]
      if (nrow(camk2d_result) > 0) {
        cat("🎯 CAMK2D Result:\n")
        cat("   LogFC:", round(camk2d_result$logFC, 4), "\n")
        cat("   P-value:", sprintf("%.2e", camk2d_result$P.Value), "\n")
        cat("   FDR:", sprintf("%.2e", camk2d_result$adj.P.Val), "\n")
        cat("   Significant:", if (camk2d_result$Significant) "YES ✅" else "NO", "\n")
        cat("   Direction:", camk2d_result$Regulation, "\n")
      }
      
      # Collect warnings
      if (length(dge_result$warnings) > 0) {
        warnings_list <- c(warnings_list, paste(dataset_id, ":", dge_result$warnings))
      }
      
    } else {
      failed_dge <- c(failed_dge, dataset_id)
      
      # Add failure to summary
      dge_summary <- rbind(dge_summary, data.frame(
        Dataset = dataset_id,
        Priority = dataset_config$priority,
        Samples = processed_data$processing_metadata$final_samples %||% NA,
        Total_Genes = NA,
        Significant_Genes = NA,
        CAMK_Genes_Analyzed = NA,
        CAMK_Significant = NA,
        Platform = processed_data$dataset_info$platform %||% "Unknown",
        Disease_Type = dataset_config$disease_type,
        Biological_Context = dataset_config$biological_context,
        Status = paste("FAILED:", dge_result$reason),
        stringsAsFactors = FALSE
      ))
      
      cat("❌ FAILED: DGE analysis failed -", dge_result$reason, "\n")
    }
    cat("\n")
  }
  
  # Combine all successful results
  if (length(dge_results) > 0) {
    combined_results <- do.call(rbind, dge_results)
    rownames(combined_results) <- NULL
  } else {
    combined_results <- data.frame()
  }
  
  # Generate DGE analysis summary
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📊 DIFFERENTIAL EXPRESSION ANALYSIS SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  successful_count <- length(dge_results)
  failed_count <- length(failed_dge)
  
  cat("✅ Successfully analyzed datasets:", successful_count, "\n")
  cat("❌ Failed DGE analysis:", failed_count, "\n")
  
  if (length(warnings_list) > 0) {
    cat("⚠️  Analysis warnings:", length(warnings_list), "\n")
  }
  
  if (successful_count > 0) {
    cat("\n📋 ANALYSIS RESULTS:\n")
    successful_datasets <- dge_summary[dge_summary$Status == "SUCCESS", ]
    total_samples <- sum(successful_datasets$Samples)
    total_significant_camk <- sum(successful_datasets$CAMK_Significant)
    
    for (i in 1:nrow(successful_datasets)) {
      ds <- successful_datasets[i, ]
      cat(sprintf("  🧮 %-10s: %3d samples, %5d genes, %2d CAMK significant\n",
                  ds$Dataset, ds$Samples, ds$Total_Genes, ds$CAMK_Significant))
    }
    
    cat(sprintf("\n📊 COMBINED RESULTS: %d total genes analyzed, %d CAMK significances\n",
                nrow(combined_results), total_significant_camk))
    
    # CAMK2D cross-dataset consistency check
    if (nrow(combined_results) > 0) {
      camk2d_results <- combined_results[combined_results$Gene_Symbol == "CAMK2D", ]
      if (nrow(camk2d_results) > 0) {
        upregulated_count <- sum(camk2d_results$logFC > 0)
        significant_count <- sum(camk2d_results$Significant)
        
        cat("\n⭐ CAMK2D CROSS-DATASET CONSISTENCY:\n")
        for (i in 1:nrow(camk2d_results)) {
          ds <- camk2d_results[i, ]
          direction <- ifelse(ds$logFC > 0, "UP", "DOWN")
          sig <- ifelse(ds$Significant, "SIG", "NS")
          cat(sprintf("  ⭐ %-10s: %s (logFC=%7.4f, p=%8.2e) [%s]\n",
                      ds$Dataset, direction, ds$logFC, ds$P.Value, sig))
        }
        
        consistency_score <- upregulated_count / nrow(camk2d_results) * 100
        cat(sprintf("\n🎯 CAMK2D Summary: %d/%d datasets UP, %d/%d significant\n",
                    upregulated_count, nrow(camk2d_results), significant_count, nrow(camk2d_results)))
        
        if (consistency_score >= 80) {
          cat("🏆 EXCELLENT: High cross-dataset consistency (", round(consistency_score, 1), "%)\n")
        } else if (consistency_score >= 60) {
          cat("✅ GOOD: Moderate cross-dataset consistency (", round(consistency_score, 1), "%)\n")
        } else {
          cat("⚠️  MIXED: Lower cross-dataset consistency (", round(consistency_score, 1), "%)\n")
        }
      }
    }
  }
  
  if (failed_count > 0) {
    cat("\n❌ DGE ANALYSIS FAILURES:\n")
    failed_summary <- dge_summary[dge_summary$Status != "SUCCESS", ]
    for (i in 1:nrow(failed_summary)) {
      ds <- failed_summary[i, ]
      cat(sprintf("  ❌ %-10s: %s\n", ds$Dataset, ds$Status))
    }
  }
  
  # Validate minimum requirements
  min_datasets <- config$validation$success_criteria$min_datasets_processed %||% 3
  
  if (successful_count < min_datasets) {
    return(create_step_result(
      success = FALSE,
      error_message = paste("Insufficient datasets analyzed:", successful_count, "< minimum:", min_datasets),
      step_name = step_name,
      warnings = warnings_list,
      metadata = list(
        datasets_analyzed = successful_count,
        datasets_failed = failed_count
      )
    ))
  }
  
  # Save results to files (for compatibility with existing pipeline)
  output_files <- config$paths$output_files
  
  if (!is.null(output_files$dge_results) && nrow(combined_results) > 0) {
    write.csv(combined_results, output_files$dge_results, row.names = FALSE)
    cat("💾 DGE results saved to:", output_files$dge_results, "\n")
    
    # Also save to legacy location for compatibility
    if (!is.null(output_files$legacy_dge)) {
      write.csv(combined_results, output_files$legacy_dge, row.names = FALSE)
    }
  }
  
  if (!is.null(output_files$dataset_summary)) {
    # Merge with preprocessing summary for complete picture
    complete_summary <- merge(preprocessing_summary, dge_summary, by = "Dataset", all = TRUE, suffixes = c("_preprocessing", "_dge"))
    write.csv(complete_summary, output_files$dataset_summary, row.names = FALSE)
    cat("💾 Dataset summary saved to:", output_files$dataset_summary, "\n")
    
    # Also save to legacy location
    if (!is.null(output_files$legacy_summary)) {
      write.csv(complete_summary, output_files$legacy_summary, row.names = FALSE)
    }
  }
  
  # Create output data structure
  output_data <- list(
    dge_results = combined_results,
    dge_summary = dge_summary,
    preprocessing_summary = preprocessing_summary,
    analysis_stats = list(
      total_datasets_analyzed = successful_count,
      total_genes_in_results = nrow(combined_results),
      total_samples_analyzed = sum(dge_summary$Samples[!is.na(dge_summary$Samples)]),
      camk_genes_analyzed = length(camk_genes),
      camk_significant_total = sum(dge_summary$CAMK_Significant[!is.na(dge_summary$CAMK_Significant)])
    )
  )
  
  # Validate output
  validate_step_output(
    output_data = output_data,
    required_fields = c("dge_results", "dge_summary", "analysis_stats"),
    step_name = step_name
  )
  
  cat("\n🎉 DIFFERENTIAL EXPRESSION ANALYSIS COMPLETED SUCCESSFULLY\n")
  cat("Ready for meta-analysis...\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(create_step_result(
    success = TRUE,
    output_data = output_data,
    step_name = step_name,
    warnings = warnings_list,
    metadata = list(
      datasets_analyzed = successful_count,
      total_genes_analyzed = nrow(combined_results),
      camk_genes_significant = sum(dge_summary$CAMK_Significant[!is.na(dge_summary$CAMK_Significant)])
    )
  ))
}

#' Analyze Dataset DGE
#'
#' Performs differential expression analysis for a single dataset
#' @param dataset_id Dataset identifier
#' @param processed_data Preprocessed dataset
#' @param dataset_config Dataset configuration
#' @param config Full pipeline configuration
#' @return DGE analysis result
analyze_dataset_dge <- function(dataset_id, processed_data, dataset_config, config) {
  
  warnings_generated <- character(0)
  
  # Extract analysis parameters
  de_config <- config$analysis$differential_expression
  qc_config <- config$analysis$quality_control
  camk_genes <- config$genes$camk_core_genes
  
  # Get data
  expr_matrix <- processed_data$expression_matrix
  groups_vector <- processed_data$groups
  
  # Validate data
  if (is.null(expr_matrix) || is.null(groups_vector)) {
    return(list(
      success = FALSE,
      reason = "missing_expression_or_groups"
    ))
  }
  
  if (ncol(expr_matrix) != length(groups_vector)) {
    return(list(
      success = FALSE,
      reason = "dimension_mismatch"
    ))
  }
  
  cat("📊 Analysis Matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
  
  # Check CAMK genes presence
  camk_present <- intersect(rownames(expr_matrix), camk_genes)
  cat("🧬 CAMK genes detected:", length(camk_present), "/", length(camk_genes), "\n")
  
  if (length(camk_present) == 0) {
    return(list(
      success = FALSE,
      reason = "no_camk_genes_detected"
    ))
  }
  
  cat("CAMK genes:", paste(camk_present, collapse = ", "), "\n")
  
  # Perform differential expression analysis
  cat("🧮 DIFFERENTIAL EXPRESSION ANALYSIS:\n")
  cat("Design: Control (reference) vs Disease (comparison)\n")
  cat("Interpretation: Positive logFC = UP in Disease\n")
  
  tryCatch({
    # Create design matrix
    design <- model.matrix(~ groups_vector)
    colnames(design) <- c("Intercept", "Disease_vs_Control")
    
    # Fit linear model
    fit <- lmFit(expr_matrix, design)
    fit <- eBayes(fit)
    
    # Get results for all genes
    all_results <- topTable(fit, coef = "Disease_vs_Control", number = Inf, 
                           adjust.method = de_config$adjust_method %||% "BH")
    
    total_genes <- nrow(all_results)
    fdr_threshold <- de_config$fdr_threshold %||% 0.05
    significant_genes <- sum(all_results$adj.P.Val < fdr_threshold)
    
    cat("📋 Analysis Results:\n")
    cat("Total genes analyzed:", total_genes, "\n")
    cat("Significant genes (FDR <", fdr_threshold, "):", significant_genes, "\n")
    
    # Extract and enhance CAMK gene results
    camk_results <- all_results[intersect(rownames(all_results), camk_present), , drop = FALSE]
    
    if (nrow(camk_results) > 0) {
      # Add comprehensive metadata
      camk_results$Gene_Symbol <- rownames(camk_results)
      camk_results$Dataset <- dataset_id
      camk_results$Biological_Context <- dataset_config$biological_context
      camk_results$Disease_Type <- dataset_config$disease_type
      camk_results$Control_Type <- dataset_config$control_type
      camk_results$Priority <- dataset_config$priority
      camk_results$Sample_Size <- ncol(expr_matrix)
      camk_results$Platform <- processed_data$dataset_info$platform %||% "Unknown"
      camk_results$Significant <- camk_results$adj.P.Val < fdr_threshold
      camk_results$Regulation <- ifelse(camk_results$logFC > 0, "UP in Disease", "DOWN in Disease")
      
      # Quality assessment using configuration thresholds
      high_logfc_threshold <- qc_config$high_logfc_threshold %||% 0.8
      low_effect_threshold <- qc_config$low_effect_threshold %||% 0.01
      
      camk_results$Quality_Flag <- ifelse(abs(camk_results$logFC) > high_logfc_threshold, "HIGH_LOGFC",
                                         ifelse(abs(camk_results$logFC) < low_effect_threshold, "LOW_EFFECT", "GOOD"))
      
      # Sort by significance
      camk_results <- camk_results[order(camk_results$P.Value), ]
      
      camk_significant <- sum(camk_results$Significant)
      
      cat("\n🧬 CAMK GENE RESULTS:\n")
      cat("CAMK genes in results:", nrow(camk_results), "\n")
      cat("Significant CAMK genes:", camk_significant, "\n")
      cat("Quality flags:", table(camk_results$Quality_Flag), "\n")
      
      # Check for quality issues
      high_logfc_count <- sum(camk_results$Quality_Flag == "HIGH_LOGFC")
      if (high_logfc_count > 0) {
        warnings_generated <- c(warnings_generated, 
                              paste("High logFC values detected:", high_logfc_count, "CAMK genes"))
      }
      
      # Display key CAMK results
      for (i in 1:nrow(camk_results)) {
        gene <- camk_results$Gene_Symbol[i]
        logfc <- round(camk_results$logFC[i], 4)
        pval <- camk_results$P.Value[i]
        adj_pval <- camk_results$adj.P.Val[i]
        direction <- ifelse(logfc > 0, "UP", "DOWN")
        sig_status <- if (adj_pval < fdr_threshold) "SIG" else "NS"
        quality <- camk_results$Quality_Flag[i]
        
        cat(sprintf("  🧬 %-8s: %4s logFC=%7.4f, p=%8.2e, FDR=%8.2e [%s] {%s}\n",
                    gene, direction, logfc, pval, adj_pval, sig_status, quality))
      }
      
    } else {
      return(list(
        success = FALSE,
        reason = "no_camk_results_generated"
      ))
    }
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      reason = paste("analysis_error:", e$message)
    ))
  })
  
  return(list(
    success = TRUE,
    results = camk_results,
    sample_count = ncol(expr_matrix),
    total_genes_analyzed = total_genes,
    significant_genes_total = significant_genes,
    camk_genes_analyzed = length(camk_present),
    camk_genes_significant = sum(camk_results$Significant),
    warnings = warnings_generated
  ))
}

cat("✅ STEP 03: Differential Expression Analysis loaded successfully\n\n")