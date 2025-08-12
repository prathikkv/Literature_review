#!/usr/bin/env Rscript
#' CAMK2D Research Platform - Simple Execution Script
#' 
#' This script provides a simplified interface to run the complete 
#' CAMK2D research platform with production settings.
#' 
#' Usage:
#'   Rscript run_platform.R
#'   source("run_platform.R")

cat("🚀 CAMK2D Research Platform\n")
cat("===========================\n")
cat("Production-ready bioinformatics pipeline for CAMK2D research\n\n")

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(rmarkdown)
  library(yaml)
})

# Load configuration
config_file <- "config.yml"
if (file.exists(config_file)) {
  cat("📋 Loading configuration from", config_file, "\n")
  config <- yaml::read_yaml(config_file)
} else {
  cat("⚠️ Configuration file not found, using defaults\n")
  # Default configuration
  config <- list(
    research = list(focus_area = "both", species = "human"),
    datasets = list(max_datasets = 12, min_samples = 10, min_camk2d_expression = 2),
    expression = list(validation_enabled = TRUE, validate_camk2d = TRUE),
    output = list(directory = "output", generate_reports = TRUE),
    performance = list(verbose_logging = TRUE)
  )
}

# Display configuration
cat("\n📊 Platform Configuration:\n")
cat("   Research focus:", config$research$focus_area, "\n")
cat("   Max datasets:", config$datasets$max_datasets, "\n")
cat("   Expression validation:", config$expression$validation_enabled, "\n")
cat("   Output directory:", config$output$directory, "\n\n")

# Ensure output directory exists
if (!dir.exists(config$output$directory)) {
  dir.create(config$output$directory, recursive = TRUE)
  cat("📁 Created output directory:", config$output$directory, "\n")
}

# Start pipeline execution
start_time <- Sys.time()
cat("⏱️ Pipeline started at:", format(start_time), "\n\n")

# Step 1: Literature Analysis
cat("📚 Step 1: Literature Analysis\n")
cat("==============================\n")
tryCatch({
  render(
    "01_literature_processing.Rmd",
    output_dir = config$output$directory,
    quiet = !config$performance$verbose_logging
  )
  cat("✅ Literature analysis completed\n\n")
}, error = function(e) {
  cat("❌ Literature analysis failed:", e$message, "\n\n")
})

# Step 2: Dataset Discovery
cat("🔍 Step 2: Dataset Discovery\n")
cat("============================\n")
tryCatch({
  render(
    "02_cross_species_discovery.Rmd",
    output_dir = config$output$directory,
    quiet = !config$performance$verbose_logging
  )
  cat("✅ Dataset discovery completed\n\n")
}, error = function(e) {
  cat("❌ Dataset discovery failed:", e$message, "\n\n")
})

# Step 3: Integrated Discovery & Validation
cat("🔬 Step 3: Integrated Discovery & Validation\n")
cat("===========================================\n")
tryCatch({
  render(
    "03_integrated_discovery_validation.Rmd",
    params = list(
      focus_area = config$research$focus_area,
      max_datasets = config$datasets$max_datasets,
      min_samples = config$datasets$min_samples,
      min_camk2d_expression = config$datasets$min_camk2d_expression,
      generate_analysis_ready = TRUE,
      output_dir = config$output$directory
    ),
    output_dir = config$output$directory,
    quiet = !config$performance$verbose_logging
  )
  cat("✅ Integrated validation completed\n\n")
}, error = function(e) {
  cat("❌ Integrated validation failed:", e$message, "\n\n")
})

# Step 4: Platform Summary
cat("📊 Step 4: Platform Summary\n")
cat("===========================\n")

# Check outputs
output_files <- list.files(config$output$directory, recursive = TRUE)
expression_files <- list.files(
  file.path(config$output$directory, "expression_data"), 
  pattern = "_expression\\.csv$"
)

cat("📁 Output files generated:", length(output_files), "\n")
cat("📊 Expression datasets:", length(expression_files), "\n")

# Check for key outputs
key_outputs <- c(
  "Literature" = any(grepl("Literature.*\\.csv$", output_files)),
  "Discovery" = any(grepl("Discovery.*\\.csv$", output_files)),
  "Expression" = length(expression_files) > 0,
  "Reports" = any(grepl("\\.html$", output_files))
)

cat("\n🎯 Pipeline Components:\n")
for (component in names(key_outputs)) {
  status <- if (key_outputs[[component]]) "✅" else "❌"
  cat("   ", status, component, "\n")
}

# Execution time
end_time <- Sys.time()
execution_time <- difftime(end_time, start_time, units = "mins")

cat("\n⏱️ Execution completed at:", format(end_time), "\n")
cat("🕒 Total execution time:", round(execution_time, 1), "minutes\n")

# Final status
successful_components <- sum(key_outputs)
total_components <- length(key_outputs)

if (successful_components == total_components) {
  cat("\n🏆 PLATFORM EXECUTION SUCCESSFUL\n")
  cat("===================================\n")
  cat("All pipeline components completed successfully!\n")
  cat("Results are available in:", config$output$directory, "\n")
} else {
  cat("\n⚠️ PLATFORM EXECUTION PARTIAL\n")
  cat("================================\n")
  cat("Some components may have failed. Check the logs above.\n")
  cat(successful_components, "of", total_components, "components successful.\n")
}

# Cleanup (optional)
if (exists("config") && 
    !is.null(config$performance) && 
    !is.null(config$performance$cleanup_temp_files) && 
    config$performance$cleanup_temp_files) {
  cat("\n🧹 Cleaning up temporary files...\n")
  # Clean R temporary files
  temp_files <- list.files(tempdir(), full.names = TRUE)
  unlink(temp_files[grepl("^file", basename(temp_files))], recursive = TRUE)
  cat("✅ Cleanup completed\n")
}

cat("\n🎯 CAMK2D Platform Execution Complete\n")
cat("=====================================\n")