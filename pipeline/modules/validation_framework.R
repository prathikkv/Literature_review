#!/usr/bin/env Rscript
#' Validation Framework Module
#' 
#' Comprehensive validation framework for pipeline consistency
#' Ensures new dynamic pipeline maintains scientific accuracy
#' 
#' @author Claude Code Enhancement Module
#' @version 1.0.0

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(yaml)
  library(testthat)
})

#' Validate DGE Results
#'
#' Compares new DGE results with original/baseline results
#' @param new_results New pipeline DGE results
#' @param original_results Original/baseline DGE results
#' @param tolerance Numeric tolerance for comparisons (default 0.001)
#' @param output_file Optional file to save validation report
#' @return Validation results with pass/fail status
validate_dge_results <- function(new_results, 
                               original_results,
                               tolerance = 0.001,
                               output_file = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║           DGE RESULTS VALIDATION                             ║\n")
  cat("║           Ensuring Scientific Consistency                    ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  validation_results <- list(
    passed = TRUE,
    tests_run = 0,
    tests_passed = 0,
    tests_failed = 0,
    discrepancies = list(),
    summary = character()
  )
  
  # Check if data frames exist
  if (!is.data.frame(new_results) || !is.data.frame(original_results)) {
    validation_results$passed <- FALSE
    validation_results$summary <- "Invalid input: Results must be data frames"
    return(validation_results)
  }
  
  cat("📊 Validating DGE Results\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  # Test 1: Sample size consistency
  cat("1. Sample size validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  if ("Sample_Size" %in% names(new_results) && "Sample_Size" %in% names(original_results)) {
    sample_diff <- abs(sum(new_results$Sample_Size) - sum(original_results$Sample_Size))
    if (sample_diff == 0) {
      cat(" ✅ PASSED\n")
      validation_results$tests_passed <- validation_results$tests_passed + 1
    } else {
      cat(" ❌ FAILED\n")
      validation_results$tests_failed <- validation_results$tests_failed + 1
      validation_results$discrepancies$sample_size <- paste("Difference:", sample_diff)
    }
  }
  
  # Test 2: Gene overlap
  cat("2. Gene overlap validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  new_genes <- unique(new_results$Gene_Symbol)
  orig_genes <- unique(original_results$Gene_Symbol)
  common_genes <- intersect(new_genes, orig_genes)
  
  overlap_pct <- length(common_genes) / length(orig_genes) * 100
  
  if (overlap_pct >= 95) {
    cat(" ✅ PASSED (", round(overlap_pct, 1), "% overlap)\n", sep = "")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED (", round(overlap_pct, 1), "% overlap)\n", sep = "")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$gene_overlap <- paste(overlap_pct, "% overlap")
  }
  
  # Test 3: LogFC consistency
  cat("3. LogFC consistency validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  logfc_comparison <- compare_numeric_columns(
    new_results, original_results, 
    "logFC", common_genes, tolerance
  )
  
  if (logfc_comparison$passed) {
    cat(" ✅ PASSED\n")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED\n")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$logfc <- logfc_comparison$details
  }
  
  # Test 4: P-value consistency
  cat("4. P-value consistency validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  pval_comparison <- compare_numeric_columns(
    new_results, original_results,
    "P.Value", common_genes, tolerance
  )
  
  if (pval_comparison$passed) {
    cat(" ✅ PASSED\n")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED\n")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$pvalue <- pval_comparison$details
  }
  
  # Test 5: Significance consistency
  cat("5. Significance consistency validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  if ("Significant" %in% names(new_results) && "Significant" %in% names(original_results)) {
    sig_comparison <- compare_significance(new_results, original_results, common_genes)
    
    if (sig_comparison$passed) {
      cat(" ✅ PASSED\n")
      validation_results$tests_passed <- validation_results$tests_passed + 1
    } else {
      cat(" ❌ FAILED\n")
      validation_results$tests_failed <- validation_results$tests_failed + 1
      validation_results$discrepancies$significance <- sig_comparison$details
    }
  }
  
  # Overall pass/fail
  validation_results$passed <- validation_results$tests_failed == 0
  
  # Generate summary
  validation_results$summary <- generate_validation_summary(validation_results)
  
  # Save report if requested
  if (!is.null(output_file)) {
    save_validation_report(validation_results, output_file, "DGE Validation")
  }
  
  # Print summary
  print_validation_summary(validation_results, "DGE")
  
  return(validation_results)
}

#' Validate Meta-Analysis Results
#'
#' Validates meta-analysis consistency
#' @param new_meta New meta-analysis results
#' @param original_meta Original meta-analysis results
#' @param tolerance Numeric tolerance
#' @param output_file Optional output file
#' @return Validation results
validate_meta_analysis <- function(new_meta,
                                 original_meta,
                                 tolerance = 0.001,
                                 output_file = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║           META-ANALYSIS VALIDATION                           ║\n")
  cat("║           Verifying Statistical Consistency                  ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  validation_results <- list(
    passed = TRUE,
    tests_run = 0,
    tests_passed = 0,
    tests_failed = 0,
    discrepancies = list(),
    summary = character()
  )
  
  cat("📊 Validating Meta-Analysis Results\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  # Test 1: Number of genes analyzed
  cat("1. Gene count validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  new_genes <- nrow(new_meta)
  orig_genes <- nrow(original_meta)
  
  if (abs(new_genes - orig_genes) <= 2) {  # Allow small difference
    cat(" ✅ PASSED\n")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED\n")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$gene_count <- paste("New:", new_genes, "Original:", orig_genes)
  }
  
  # Test 2: Combined effect sizes
  cat("2. Effect size validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  effect_comparison <- compare_meta_effects(new_meta, original_meta, tolerance)
  
  if (effect_comparison$passed) {
    cat(" ✅ PASSED\n")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED\n")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$effect_sizes <- effect_comparison$details
  }
  
  # Test 3: P-values
  cat("3. Combined p-value validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  pval_comparison <- compare_meta_pvalues(new_meta, original_meta, tolerance)
  
  if (pval_comparison$passed) {
    cat(" ✅ PASSED\n")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED\n")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$pvalues <- pval_comparison$details
  }
  
  # Test 4: Heterogeneity statistics
  cat("4. Heterogeneity validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  if ("Heterogeneity_I2" %in% names(new_meta) && "Heterogeneity_I2" %in% names(original_meta)) {
    het_comparison <- compare_heterogeneity(new_meta, original_meta, tolerance = 5)  # 5% tolerance for I2
    
    if (het_comparison$passed) {
      cat(" ✅ PASSED\n")
      validation_results$tests_passed <- validation_results$tests_passed + 1
    } else {
      cat(" ❌ FAILED\n")
      validation_results$tests_failed <- validation_results$tests_failed + 1
      validation_results$discrepancies$heterogeneity <- het_comparison$details
    }
  }
  
  # Test 5: CAMK2D specific validation
  cat("5. CAMK2D results validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  camk2d_comparison <- validate_camk2d_results(new_meta, original_meta, tolerance)
  
  if (camk2d_comparison$passed) {
    cat(" ✅ PASSED\n")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED\n")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$camk2d <- camk2d_comparison$details
  }
  
  # Overall pass/fail
  validation_results$passed <- validation_results$tests_failed == 0
  
  # Generate summary
  validation_results$summary <- generate_validation_summary(validation_results)
  
  # Save report if requested
  if (!is.null(output_file)) {
    save_validation_report(validation_results, output_file, "Meta-Analysis Validation")
  }
  
  # Print summary
  print_validation_summary(validation_results, "Meta-Analysis")
  
  return(validation_results)
}

#' Validate Report Content
#'
#' Validates HTML report content consistency
#' @param new_html Path to new HTML report
#' @param original_html Path to original HTML report
#' @param output_file Optional output file
#' @return Validation results
validate_report_content <- function(new_html,
                                  original_html,
                                  output_file = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║           REPORT CONTENT VALIDATION                          ║\n")
  cat("║           Verifying Report Consistency                       ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  validation_results <- list(
    passed = TRUE,
    tests_run = 0,
    tests_passed = 0,
    tests_failed = 0,
    discrepancies = list(),
    summary = character()
  )
  
  # Check files exist
  if (!file.exists(new_html) || !file.exists(original_html)) {
    validation_results$passed <- FALSE
    validation_results$summary <- "Report files not found"
    return(validation_results)
  }
  
  cat("📊 Validating Report Content\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  # Test 1: File size comparison
  cat("1. Report size validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  new_size <- file.info(new_html)$size / 1024^2  # MB
  orig_size <- file.info(original_html)$size / 1024^2  # MB
  size_diff_pct <- abs(new_size - orig_size) / orig_size * 100
  
  if (size_diff_pct < 20) {  # Within 20% size difference
    cat(" ✅ PASSED\n")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED\n")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$file_size <- paste(round(size_diff_pct, 1), "% difference")
  }
  
  # Read HTML content
  new_content <- readLines(new_html, warn = FALSE)
  orig_content <- readLines(original_html, warn = FALSE)
  
  # Test 2: Key sections present
  cat("2. Key sections validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  key_sections <- c("Executive Summary", "Meta-Analysis", "CAMK2D", 
                   "Statistical", "Pathway", "Clinical")
  
  sections_found <- 0
  for (section in key_sections) {
    if (any(grepl(section, new_content, ignore.case = TRUE))) {
      sections_found <- sections_found + 1
    }
  }
  
  if (sections_found >= length(key_sections) - 1) {  # Allow one missing
    cat(" ✅ PASSED\n")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED\n")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$sections <- paste(sections_found, "of", length(key_sections), "found")
  }
  
  # Test 3: Figures present
  cat("3. Figure validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  new_figures <- length(grep("<img", new_content, ignore.case = TRUE))
  orig_figures <- length(grep("<img", orig_content, ignore.case = TRUE))
  
  if (abs(new_figures - orig_figures) <= 2) {
    cat(" ✅ PASSED\n")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED\n")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$figures <- paste("New:", new_figures, "Original:", orig_figures)
  }
  
  # Test 4: Tables present
  cat("4. Table validation...")
  validation_results$tests_run <- validation_results$tests_run + 1
  
  new_tables <- length(grep("<table", new_content, ignore.case = TRUE))
  orig_tables <- length(grep("<table", orig_content, ignore.case = TRUE))
  
  if (abs(new_tables - orig_tables) <= 2) {
    cat(" ✅ PASSED\n")
    validation_results$tests_passed <- validation_results$tests_passed + 1
  } else {
    cat(" ❌ FAILED\n")
    validation_results$tests_failed <- validation_results$tests_failed + 1
    validation_results$discrepancies$tables <- paste("New:", new_tables, "Original:", orig_tables)
  }
  
  # Overall pass/fail
  validation_results$passed <- validation_results$tests_failed == 0
  
  # Generate summary
  validation_results$summary <- generate_validation_summary(validation_results)
  
  # Save report if requested
  if (!is.null(output_file)) {
    save_validation_report(validation_results, output_file, "Report Content Validation")
  }
  
  # Print summary
  print_validation_summary(validation_results, "Report Content")
  
  return(validation_results)
}

#' Validate Statistical Significance
#'
#' Validates statistical significance consistency
#' @param new_stats New statistical results
#' @param original_stats Original statistical results
#' @param tolerance P-value tolerance
#' @param output_file Optional output file
#' @return Validation results
validate_statistical_significance <- function(new_stats,
                                            original_stats,
                                            tolerance = 0.001,
                                            output_file = NULL) {
  
  validation_results <- list(
    passed = TRUE,
    tests_run = 0,
    tests_passed = 0,
    tests_failed = 0,
    discrepancies = list()
  )
  
  # Validate p-value distributions
  if ("p_value" %in% names(new_stats) && "p_value" %in% names(original_stats)) {
    validation_results$tests_run <- validation_results$tests_run + 1
    
    p_diff <- abs(new_stats$p_value - original_stats$p_value)
    
    if (all(p_diff < tolerance, na.rm = TRUE)) {
      validation_results$tests_passed <- validation_results$tests_passed + 1
    } else {
      validation_results$tests_failed <- validation_results$tests_failed + 1
      validation_results$discrepancies$p_values <- paste("Max difference:", max(p_diff, na.rm = TRUE))
    }
  }
  
  validation_results$passed <- validation_results$tests_failed == 0
  
  return(validation_results)
}

# Helper functions

#' Compare Numeric Columns
#'
#' Helper to compare numeric columns between datasets
compare_numeric_columns <- function(new_data, orig_data, column, common_genes, tolerance) {
  result <- list(passed = TRUE, details = character())
  
  discrepancies <- character()
  
  for (gene in common_genes) {
    new_val <- new_data[new_data$Gene_Symbol == gene, column]
    orig_val <- orig_data[orig_data$Gene_Symbol == gene, column]
    
    if (length(new_val) > 0 && length(orig_val) > 0) {
      diff <- abs(new_val[1] - orig_val[1])
      
      if (diff > tolerance) {
        discrepancies <- c(discrepancies, 
                         paste(gene, ": diff =", round(diff, 6)))
      }
    }
  }
  
  if (length(discrepancies) > 0) {
    result$passed <- FALSE
    result$details <- paste(length(discrepancies), "genes exceed tolerance")
  }
  
  return(result)
}

#' Compare Significance
#'
#' Helper to compare significance calls
compare_significance <- function(new_data, orig_data, common_genes) {
  result <- list(passed = TRUE, details = character())
  
  mismatches <- 0
  
  for (gene in common_genes) {
    new_sig <- new_data[new_data$Gene_Symbol == gene, "Significant"]
    orig_sig <- orig_data[orig_data$Gene_Symbol == gene, "Significant"]
    
    if (length(new_sig) > 0 && length(orig_sig) > 0) {
      if (new_sig[1] != orig_sig[1]) {
        mismatches <- mismatches + 1
      }
    }
  }
  
  if (mismatches > length(common_genes) * 0.05) {  # Allow 5% mismatch
    result$passed <- FALSE
    result$details <- paste(mismatches, "significance mismatches")
  }
  
  return(result)
}

#' Compare Meta-Analysis Effects
compare_meta_effects <- function(new_meta, orig_meta, tolerance) {
  result <- list(passed = TRUE, details = character())
  
  # Find common genes
  common_genes <- intersect(new_meta$Gene, orig_meta$Gene)
  
  if (length(common_genes) == 0) {
    result$passed <- FALSE
    result$details <- "No common genes found"
    return(result)
  }
  
  # Compare effect sizes
  max_diff <- 0
  for (gene in common_genes) {
    new_effect <- new_meta[new_meta$Gene == gene, "Combined_logFC"]
    orig_effect <- orig_meta[orig_meta$Gene == gene, "Combined_logFC"]
    
    if (length(new_effect) > 0 && length(orig_effect) > 0) {
      diff <- abs(new_effect[1] - orig_effect[1])
      max_diff <- max(max_diff, diff)
    }
  }
  
  if (max_diff > tolerance) {
    result$passed <- FALSE
    result$details <- paste("Max effect difference:", round(max_diff, 6))
  }
  
  return(result)
}

#' Compare Meta-Analysis P-values
compare_meta_pvalues <- function(new_meta, orig_meta, tolerance) {
  result <- list(passed = TRUE, details = character())
  
  common_genes <- intersect(new_meta$Gene, orig_meta$Gene)
  
  max_diff <- 0
  for (gene in common_genes) {
    new_p <- new_meta[new_meta$Gene == gene, "Combined_P_Value"]
    orig_p <- orig_meta[orig_meta$Gene == gene, "Combined_P_Value"]
    
    if (length(new_p) > 0 && length(orig_p) > 0) {
      diff <- abs(new_p[1] - orig_p[1])
      max_diff <- max(max_diff, diff)
    }
  }
  
  if (max_diff > tolerance) {
    result$passed <- FALSE
    result$details <- paste("Max p-value difference:", round(max_diff, 6))
  }
  
  return(result)
}

#' Compare Heterogeneity
compare_heterogeneity <- function(new_meta, orig_meta, tolerance) {
  result <- list(passed = TRUE, details = character())
  
  common_genes <- intersect(new_meta$Gene, orig_meta$Gene)
  
  large_diffs <- 0
  for (gene in common_genes) {
    new_i2 <- new_meta[new_meta$Gene == gene, "Heterogeneity_I2"]
    orig_i2 <- orig_meta[orig_meta$Gene == gene, "Heterogeneity_I2"]
    
    if (length(new_i2) > 0 && length(orig_i2) > 0) {
      diff <- abs(new_i2[1] - orig_i2[1])
      if (diff > tolerance) {
        large_diffs <- large_diffs + 1
      }
    }
  }
  
  if (large_diffs > 0) {
    result$passed <- FALSE
    result$details <- paste(large_diffs, "genes with I² differences >", tolerance, "%")
  }
  
  return(result)
}

#' Validate CAMK2D Results
validate_camk2d_results <- function(new_meta, orig_meta, tolerance) {
  result <- list(passed = TRUE, details = character())
  
  # Check CAMK2D specifically
  new_camk2d <- new_meta[new_meta$Gene == "CAMK2D", ]
  orig_camk2d <- orig_meta[orig_meta$Gene == "CAMK2D", ]
  
  if (nrow(new_camk2d) == 0 || nrow(orig_camk2d) == 0) {
    result$passed <- FALSE
    result$details <- "CAMK2D not found in results"
    return(result)
  }
  
  # Compare key metrics
  logfc_diff <- abs(new_camk2d$Combined_logFC[1] - orig_camk2d$Combined_logFC[1])
  pval_diff <- abs(new_camk2d$Combined_P_Value[1] - orig_camk2d$Combined_P_Value[1])
  
  if (logfc_diff > tolerance || pval_diff > tolerance) {
    result$passed <- FALSE
    result$details <- paste("LogFC diff:", round(logfc_diff, 6), 
                          "P-value diff:", round(pval_diff, 6))
  }
  
  return(result)
}

#' Generate Validation Summary
generate_validation_summary <- function(results) {
  summary <- paste0(
    "Tests run: ", results$tests_run, "\n",
    "Tests passed: ", results$tests_passed, "\n",
    "Tests failed: ", results$tests_failed, "\n",
    "Overall: ", ifelse(results$passed, "PASSED", "FAILED")
  )
  
  if (length(results$discrepancies) > 0) {
    summary <- paste0(summary, "\n\nDiscrepancies found:\n")
    for (name in names(results$discrepancies)) {
      summary <- paste0(summary, "- ", name, ": ", results$discrepancies[[name]], "\n")
    }
  }
  
  return(summary)
}

#' Save Validation Report
save_validation_report <- function(results, output_file, validation_type) {
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  report <- list(
    validation_type = validation_type,
    timestamp = Sys.time(),
    passed = results$passed,
    tests_run = results$tests_run,
    tests_passed = results$tests_passed,
    tests_failed = results$tests_failed,
    discrepancies = results$discrepancies,
    summary = results$summary
  )
  
  saveRDS(report, output_file)
  cat("\n📄 Validation report saved to:", output_file, "\n")
}

#' Print Validation Summary
print_validation_summary <- function(results, type) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📊", type, "VALIDATION SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  cat("Tests run:", results$tests_run, "\n")
  cat("Tests passed:", results$tests_passed, "✅\n")
  cat("Tests failed:", results$tests_failed, if(results$tests_failed > 0) "❌" else "", "\n")
  
  cat("\n🎯 OVERALL RESULT:", if(results$passed) "PASSED ✅" else "FAILED ❌", "\n")
  
  if (length(results$discrepancies) > 0) {
    cat("\n⚠️  DISCREPANCIES:\n")
    for (name in names(results$discrepancies)) {
      cat("  -", name, ":", results$discrepancies[[name]], "\n")
    }
  }
  
  cat("═══════════════════════════════════════════════════════════════\n\n")
}

# Module loading confirmation
cat("✅ Validation Framework Module loaded successfully\n")
cat("   Functions: validate_dge_results(), validate_meta_analysis(),\n")
cat("             validate_report_content(), validate_statistical_significance()\n")
cat("   Version: 1.0.0\n\n")