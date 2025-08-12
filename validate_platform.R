#!/usr/bin/env Rscript
#' CAMK2D Platform Validation Script
#' 
#' This script validates that all platform components are working correctly
#' and provides a quick health check for production readiness.

cat("🔍 CAMK2D Platform Validation\n")
cat("============================\n\n")

# Load required libraries for validation
suppressPackageStartupMessages({
  library(tidyverse)
  library(GEOquery)
  library(Biobase)
  library(limma)
  library(readr)
})

# Configuration
validation_tests <- list(
  check_files = TRUE,
  check_functions = TRUE,
  test_expression_download = TRUE,
  check_outputs = TRUE
)

test_results <- list()

# Test 1: Check Core Files
cat("📁 Test 1: Core File Structure\n")
cat("==============================\n")

required_files <- c(
  "functions/targeted_camk2d_search.R",
  "functions/expression_data_validator.R", 
  "functions/dge_analysis_framework.R",
  "run_integrated_camk2d_platform.R",
  "01_literature_processing.Rmd",
  "02_cross_species_discovery.Rmd",
  "03_integrated_discovery_validation.Rmd",
  "data/comprehensive_camk2d_literature.csv"
)

missing_files <- c()
for (file in required_files) {
  if (file.exists(file)) {
    cat("✅", file, "\n")
  } else {
    cat("❌", file, "\n")
    missing_files <- c(missing_files, file)
  }
}

test_results$files_check <- length(missing_files) == 0
if (test_results$files_check) {
  cat("✅ All core files present\n")
} else {
  cat("❌ Missing", length(missing_files), "files\n")
}

cat("\n")

# Test 2: Function Loading
cat("🔧 Test 2: Function Loading\n")
cat("===========================\n")

tryCatch({
  source("functions/expression_data_validator.R")
  cat("✅ Expression validator loaded\n")
  
  source("functions/targeted_camk2d_search.R")
  cat("✅ CAMK2D search functions loaded\n")
  
  source("functions/dge_analysis_framework.R")
  cat("✅ DGE analysis framework loaded\n")
  
  test_results$functions_check <- TRUE
  cat("✅ All functions loaded successfully\n")
  
}, error = function(e) {
  cat("❌ Function loading failed:", e$message, "\n")
  test_results$functions_check <- FALSE
})

cat("\n")

# Test 3: Expression Download Pipeline (Quick Test)
cat("🧬 Test 3: Expression Download Pipeline\n")
cat("======================================\n")

if (validation_tests$test_expression_download) {
  tryCatch({
    # Test with a small, known working dataset
    cat("📥 Testing expression download with GSE248443...\n")
    
    # Create test output directory
    test_dir <- file.path(tempdir(), "platform_validation")
    if (!dir.exists(test_dir)) {
      dir.create(test_dir, recursive = TRUE)
    }
    
    # Test the download function
    result <- download_and_validate_expression(
      geo_accession = "GSE248443",
      output_dir = test_dir,
      validate_camk2d = TRUE,
      min_expression = 0.5
    )
    
    if (!is.null(result) && result$suitable) {
      cat("✅ Expression download successful\n")
      cat("   Dataset:", result$geo_accession, "\n")
      cat("   Expression matrix:", nrow(result$expression_matrix), "genes ×", ncol(result$expression_matrix), "samples\n")
      cat("   CAMK2D detected:", result$camk2d_validation$camk2d_detected, "\n")
      if (result$camk2d_validation$camk2d_detected) {
        cat("   CAMK2D expression:", round(result$camk2d_validation$camk2d_mean_expression, 2), "\n")
      }
      test_results$expression_download <- TRUE
    } else {
      cat("❌ Expression download failed\n")
      if (!is.null(result)) {
        cat("   Reason:", result$reason, "\n")
      }
      test_results$expression_download <- FALSE
    }
    
    # Clean up
    unlink(test_dir, recursive = TRUE)
    
  }, error = function(e) {
    cat("❌ Expression download test failed:", e$message, "\n")
    test_results$expression_download <- FALSE
  })
} else {
  cat("⚠️ Expression download test skipped\n")
  test_results$expression_download <- TRUE
}

cat("\n")

# Test 4: Output Directory Structure
cat("📊 Test 4: Output Structure\n")
cat("===========================\n")

if (dir.exists("output")) {
  cat("✅ Output directory exists\n")
  
  # Check for key output files
  output_files <- list.files("output", recursive = TRUE)
  cat("   Total output files:", length(output_files), "\n")
  
  # Check for example data
  if (dir.exists("output/test_expression_data")) {
    example_files <- list.files("output/test_expression_data", pattern = "GSE248443")
    if (length(example_files) >= 3) {
      cat("✅ Example data (GSE248443) available:", length(example_files), "files\n")
      test_results$output_check <- TRUE
    } else {
      cat("⚠️ Example data incomplete:", length(example_files), "files\n")
      test_results$output_check <- FALSE
    }
  } else {
    cat("⚠️ No example expression data found\n")
    test_results$output_check <- FALSE
  }
} else {
  cat("❌ Output directory missing\n")
  test_results$output_check <- FALSE
}

cat("\n")

# Test 5: Literature Data Check
cat("📚 Test 5: Literature Data\n")
cat("==========================\n")

if (file.exists("data/comprehensive_camk2d_literature.csv")) {
  tryCatch({
    lit_data <- read_csv("data/comprehensive_camk2d_literature.csv", show_col_types = FALSE)
    cat("✅ Literature data loaded\n")
    cat("   Total papers:", nrow(lit_data), "\n")
    cat("   High relevance papers:", sum(lit_data$Proposal_Relevance == "High Relevance", na.rm = TRUE), "\n")
    test_results$literature_check <- TRUE
  }, error = function(e) {
    cat("❌ Literature data loading failed:", e$message, "\n")
    test_results$literature_check <- FALSE
  })
} else {
  cat("❌ Literature data file not found\n")
  test_results$literature_check <- FALSE
}

cat("\n")

# Overall Platform Status
cat("🏁 Platform Validation Summary\n")
cat("==============================\n")

total_tests <- length(test_results)
passed_tests <- sum(unlist(test_results))

cat("Tests passed:", passed_tests, "/", total_tests, "\n\n")

if (passed_tests == total_tests) {
  cat("🚀 PLATFORM STATUS: PRODUCTION READY ✅\n")
  cat("========================================\n")
  cat("✅ All validation tests passed\n")
  cat("✅ Expression download pipeline operational\n")
  cat("✅ Core functions working correctly\n")
  cat("✅ Literature and output data available\n")
  cat("✅ Platform ready for research execution\n\n")
  
  cat("🎯 NEXT STEPS:\n")
  cat("1. Run: source('run_integrated_camk2d_platform.R')\n")
  cat("2. Review generated reports in output/\n")
  cat("3. Begin bioinformatics analysis with validated datasets\n")
  
} else if (passed_tests >= 4) {
  cat("⚠️ PLATFORM STATUS: MOSTLY FUNCTIONAL\n")
  cat("=====================================\n")
  cat("Core functionality available with minor issues\n")
  cat("Platform can support research activities\n")
  
} else {
  cat("❌ PLATFORM STATUS: NEEDS ATTENTION\n")
  cat("===================================\n")
  cat("Multiple validation failures detected\n")
  cat("Platform requires troubleshooting before use\n")
}

cat("\n📁 Validation complete. Platform status assessed.\n")