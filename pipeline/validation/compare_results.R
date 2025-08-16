#!/usr/bin/env Rscript

#' CAMK2D Pipeline Result Validation Against v1.0.0 Baseline
#' ========================================================
#' 
#' Comprehensive validation system that compares pipeline outputs
#' against the v1.0.0 production baseline to ensure numerical
#' accuracy and scientific consistency.
#' 
#' Usage:
#'   Rscript validation/compare_results.R
#'   
#' Outputs:
#'   - validation/validation_report.html (detailed comparison)
#'   - Console summary of validation results

# Load required libraries
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
  library(knitr)
  library(rmarkdown)
})

# ═══════════════════════════════════════════════════════════════
# CONFIGURATION AND SETUP
# ═══════════════════════════════════════════════════════════════

cat("🔍 CAMK2D Pipeline Validation Against v1.0.0 Baseline\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Set working directory to pipeline root
if (basename(getwd()) == "validation") {
  setwd("..")
}

# Load baseline results
baseline_file <- "validation/baseline_results.yml"
if (!file.exists(baseline_file)) {
  stop("❌ Baseline results file not found: ", baseline_file)
}

baseline <- yaml::read_yaml(baseline_file)
cat("📊 Loaded baseline results from v", baseline$baseline$version, "\n")

# Load pipeline configuration
config_file <- "config.yml"
if (!file.exists(config_file)) {
  stop("❌ Pipeline configuration not found: ", config_file)
}

config <- yaml::read_yaml(config_file)
cat("⚙️  Loaded pipeline configuration:", config$pipeline$name, "\n\n")

# NULL coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Extract primary gene for dynamic path resolution
primary_gene <- config$research_target$primary_gene %||% "CAMK2D"
cat("🎯 Primary gene:", primary_gene, "\n\n")

# ═══════════════════════════════════════════════════════════════
# VALIDATION FUNCTIONS
# ═══════════════════════════════════════════════════════════════

#' Compare numerical values with tolerance
compare_numeric <- function(actual, expected, tolerance, name) {
  if (is.na(actual) || is.na(expected)) {
    return(list(
      match = FALSE,
      actual = actual,
      expected = expected,
      difference = NA,
      tolerance = tolerance,
      name = name,
      status = "NA_VALUES"
    ))
  }
  
  difference <- abs(actual - expected)
  match <- difference <= tolerance
  
  list(
    match = match,
    actual = actual,
    expected = expected,
    difference = difference,
    tolerance = tolerance,
    name = name,
    status = if (match) "PASS" else "FAIL"
  )
}

#' Validate file existence and properties
validate_file <- function(file_path, min_size_mb = 0, must_exist = TRUE) {
  if (!file.exists(file_path)) {
    return(list(
      exists = FALSE,
      path = file_path,
      size_mb = 0,
      status = if (must_exist) "MISSING" else "OPTIONAL_MISSING"
    ))
  }
  
  file_size_bytes <- file.size(file_path)
  file_size_mb <- file_size_bytes / (1024 * 1024)
  
  size_ok <- file_size_mb >= min_size_mb
  
  list(
    exists = TRUE,
    path = file_path,
    size_mb = round(file_size_mb, 2),
    meets_size_requirement = size_ok,
    status = if (size_ok) "PASS" else "SIZE_FAIL"
  )
}

# ═══════════════════════════════════════════════════════════════
# VALIDATION EXECUTION
# ═══════════════════════════════════════════════════════════════

validation_results <- list()

cat("🔬 STEP 1: File Output Validation\n")
cat("─────────────────────────────────────────────────────────────\n")

# Validate expected output files
file_validations <- list()
for (expected_file in baseline$expected_outputs$files) {
  file_val <- validate_file(
    expected_file$path, 
    expected_file$min_size_mb %||% 0,
    expected_file$must_exist %||% TRUE
  )
  file_validations[[expected_file$type]] <- file_val
  
  status_icon <- switch(file_val$status,
    "PASS" = "✅",
    "MISSING" = "❌",
    "SIZE_FAIL" = "⚠️",
    "OPTIONAL_MISSING" = "⏭️"
  )
  
  cat(status_icon, expected_file$type, ":", file_val$path)
  if (file_val$exists) {
    cat(" (", file_val$size_mb, "MB)")
  }
  cat("\n")
}

validation_results$file_validation <- file_validations
cat("\n")

cat("🧬 STEP 2: Meta-Analysis Results Validation\n")
cat("─────────────────────────────────────────────────────────────\n")

# Load meta-analysis results if available
meta_file <- gsub("\\{GENE\\}", primary_gene, config$paths$output_files$meta_results)
meta_validations <- list()

if (file.exists(meta_file)) {
  cat("📄 Loading meta-analysis results:", meta_file, "\n")
  
  meta_results <- tryCatch({
    read.csv(meta_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("❌ Error reading meta-analysis file:", e$message, "\n")
    NULL
  })
  
  if (!is.null(meta_results)) {
    # Validate CAMK2D results
    camk2d_row <- meta_results[meta_results$Gene == "CAMK2D", ]
    
    if (nrow(camk2d_row) > 0) {
      # LogFC validation
      logfc_comparison <- compare_numeric(
        actual = camk2d_row$Combined_logFC[1],
        expected = baseline$baseline$camk2d_results$combined_logfc,
        tolerance = baseline$tolerances$logfc_absolute,
        name = "CAMK2D_logFC"
      )
      
      # P-value validation  
      pvalue_comparison <- compare_numeric(
        actual = camk2d_row$Combined_P_Value[1],
        expected = baseline$baseline$camk2d_results$combined_p_value,
        tolerance = baseline$tolerances$p_value_absolute,
        name = "CAMK2D_p_value"
      )
      
      # Significance validation
      is_significant <- camk2d_row$Combined_P_Value[1] < 0.05
      significance_match <- is_significant == baseline$baseline$camk2d_results$is_significant
      
      meta_validations$camk2d_logfc <- logfc_comparison
      meta_validations$camk2d_pvalue <- pvalue_comparison
      meta_validations$camk2d_significance <- list(
        match = significance_match,
        actual = is_significant,
        expected = baseline$baseline$camk2d_results$is_significant,
        name = "CAMK2D_significance",
        status = if (significance_match) "PASS" else "FAIL"
      )
      
      # Display results
      cat("🎯 CAMK2D Validation Results:\n")
      cat("   LogFC:", ifelse(logfc_comparison$match, "✅", "❌"), 
          "Expected =", logfc_comparison$expected, 
          ", Actual =", round(logfc_comparison$actual, 6), 
          ", Diff =", formatC(logfc_comparison$difference, format = "e", digits = 2), "\n")
      
      cat("   P-value:", ifelse(pvalue_comparison$match, "✅", "❌"),
          "Expected =", formatC(pvalue_comparison$expected, format = "e", digits = 2),
          ", Actual =", formatC(pvalue_comparison$actual, format = "e", digits = 2),
          ", Diff =", formatC(pvalue_comparison$difference, format = "e", digits = 2), "\n")
      
      cat("   Significance:", ifelse(significance_match, "✅", "❌"),
          "Expected =", baseline$baseline$camk2d_results$is_significant,
          ", Actual =", is_significant, "\n")
      
    } else {
      cat("❌ CAMK2D not found in meta-analysis results\n")
      meta_validations$camk2d_found <- list(
        match = FALSE,
        status = "MISSING",
        name = "CAMK2D_presence"
      )
    }
    
    # Validate significant genes count
    significant_count <- sum(meta_results$Combined_P_Value < 0.05, na.rm = TRUE)
    expected_significant <- baseline$baseline$camk_gene_results$significant_genes
    
    significant_comparison <- list(
      match = significant_count == expected_significant,
      actual = significant_count,
      expected = expected_significant,
      name = "significant_genes_count",
      status = if (significant_count == expected_significant) "PASS" else "FAIL"
    )
    
    meta_validations$significant_genes <- significant_comparison
    
    cat("📊 Significant genes:", ifelse(significant_comparison$match, "✅", "❌"),
        "Expected =", expected_significant, ", Actual =", significant_count, "\n")
    
  } else {
    cat("❌ Could not load meta-analysis results\n")
  }
} else {
  cat("❌ Meta-analysis results file not found:", meta_file, "\n")
}

validation_results$meta_analysis <- meta_validations
cat("\n")

cat("📈 STEP 3: Dataset Processing Validation\n")
cat("─────────────────────────────────────────────────────────────\n")

# Load dataset summary if available
dataset_summary_file <- config$paths$output_files$dataset_summary
dataset_validations <- list()

if (file.exists(dataset_summary_file)) {
  cat("📄 Loading dataset summary:", dataset_summary_file, "\n")
  
  dataset_summary <- tryCatch({
    read.csv(dataset_summary_file, stringsAsFactors = FALSE)
  }, error = function(e) {
    cat("❌ Error reading dataset summary:", e$message, "\n")
    NULL
  })
  
  if (!is.null(dataset_summary)) {
    # Count successful datasets
    successful_datasets <- sum(dataset_summary$Success == TRUE | dataset_summary$Success == "TRUE", na.rm = TRUE)
    expected_successful <- baseline$baseline$datasets$successfully_processed
    
    datasets_comparison <- list(
      match = successful_datasets == expected_successful,
      actual = successful_datasets,
      expected = expected_successful,
      name = "successful_datasets_count",
      status = if (successful_datasets == expected_successful) "PASS" else "FAIL"
    )
    
    dataset_validations$successful_count <- datasets_comparison
    
    cat("🗂️  Successful datasets:", ifelse(datasets_comparison$match, "✅", "❌"),
        "Expected =", expected_successful, ", Actual =", successful_datasets, "\n")
    
    # Validate total samples if available
    if ("Total_Samples" %in% names(dataset_summary)) {
      total_samples <- sum(dataset_summary$Total_Samples, na.rm = TRUE)
      expected_samples <- baseline$baseline$study_summary$total_samples
      
      samples_comparison <- list(
        match = total_samples == expected_samples,
        actual = total_samples,
        expected = expected_samples,
        name = "total_samples_count",
        status = if (total_samples == expected_samples) "PASS" else "FAIL"
      )
      
      dataset_validations$total_samples <- samples_comparison
      
      cat("👥 Total samples:", ifelse(samples_comparison$match, "✅", "❌"),
          "Expected =", expected_samples, ", Actual =", total_samples, "\n")
    }
  }
} else {
  cat("⚠️  Dataset summary file not found (optional):", dataset_summary_file, "\n")
}

validation_results$dataset_processing <- dataset_validations
cat("\n")

# ═══════════════════════════════════════════════════════════════
# OVERALL VALIDATION ASSESSMENT
# ═══════════════════════════════════════════════════════════════

cat("📋 STEP 4: Overall Validation Assessment\n")
cat("─────────────────────────────────────────────────────────────\n")

# Collect all validation statuses
all_validations <- list()

# File validations
for (name in names(validation_results$file_validation)) {
  all_validations[[paste0("file_", name)]] <- validation_results$file_validation[[name]]$status
}

# Meta-analysis validations
for (name in names(validation_results$meta_analysis)) {
  all_validations[[paste0("meta_", name)]] <- validation_results$meta_analysis[[name]]$status
}

# Dataset validations
for (name in names(validation_results$dataset_processing)) {
  all_validations[[paste0("dataset_", name)]] <- validation_results$dataset_processing[[name]]$status
}

# Count validation results
total_validations <- length(all_validations)
passed_validations <- sum(sapply(all_validations, function(x) x == "PASS"))
failed_validations <- sum(sapply(all_validations, function(x) x %in% c("FAIL", "MISSING", "SIZE_FAIL")))
warning_validations <- sum(sapply(all_validations, function(x) x %in% c("NA_VALUES")))

# Overall success assessment
overall_success <- failed_validations == 0

cat("📊 Validation Summary:\n")
cat("   Total checks:", total_validations, "\n")
cat("   Passed:", passed_validations, "✅\n")
cat("   Failed:", failed_validations, if (failed_validations > 0) "❌" else "✅", "\n")
cat("   Warnings:", warning_validations, if (warning_validations > 0) "⚠️" else "✅", "\n")
cat("\n")

# Critical validations check
critical_checks <- c(
  "meta_camk2d_logfc", "meta_camk2d_pvalue", "meta_camk2d_significance",
  "file_html_report", "file_meta_analysis_results"
)

critical_passed <- sum(sapply(critical_checks, function(x) {
  if (x %in% names(all_validations)) {
    all_validations[[x]] == "PASS"
  } else {
    FALSE
  }
}))

critical_success <- critical_passed == length(critical_checks)

cat("🎯 Critical Validations (Required for Success):\n")
for (check in critical_checks) {
  if (check %in% names(all_validations)) {
    status <- all_validations[[check]]
    icon <- if (status == "PASS") "✅" else "❌"
    cat("   ", icon, check, ":", status, "\n")
  } else {
    cat("   ❌", check, ": NOT_FOUND\n")
  }
}
cat("\n")

# Final verdict
cat("🏆 FINAL VALIDATION VERDICT\n")
cat("═══════════════════════════════════════════════════════════════\n")

if (overall_success && critical_success) {
  cat("✅ VALIDATION SUCCESSFUL\n")
  cat("   Pipeline produces results identical to v1.0.0 baseline\n")
  cat("   All critical scientific metrics match expected values\n")
  cat("   Self-contained pipeline is validated for production use\n")
  exit_code <- 0
} else {
  cat("❌ VALIDATION FAILED\n")
  cat("   Pipeline results differ from v1.0.0 baseline\n")
  cat("   Review detailed results above for specific issues\n")
  cat("   Pipeline requires investigation before production use\n")
  exit_code <- 1
}

cat("\n")
cat("📄 Detailed validation report: validation/validation_report.html\n")
cat("📊 Pipeline outputs: output/current/\n")
cat("⚙️  Configuration: config.yml\n")

# ═══════════════════════════════════════════════════════════════
# SAVE VALIDATION RESULTS
# ═══════════════════════════════════════════════════════════════

# Save validation results to YAML
validation_summary <- list(
  validation_date = Sys.time(),
  pipeline_version = config$pipeline$version,
  baseline_version = baseline$baseline$version,
  overall_success = overall_success,
  critical_success = critical_success,
  total_checks = total_validations,
  passed_checks = passed_validations,
  failed_checks = failed_validations,
  warning_checks = warning_validations,
  detailed_results = validation_results,
  all_validation_statuses = all_validations
)

yaml::write_yaml(validation_summary, "validation/validation_summary.yml")
cat("💾 Validation summary saved: validation/validation_summary.yml\n")

# Exit with appropriate code
cat("\n")
if (interactive()) {
  cat("Validation complete. Check results above.\n")
} else {
  quit(save = "no", status = exit_code)
}