#!/usr/bin/env Rscript
#' Complete Pipeline Runner with Real Analysis Execution
#' 
#' Executes the full CAMK2D analysis pipeline with configurable features
#' 
#' Usage:
#'   Rscript run_pipeline_complete.R --quick  # 2-3 minutes (core analysis)
#'   Rscript run_pipeline_complete.R          # 5 minutes (standard)
#'   Rscript run_pipeline_complete.R --full   # 10+ minutes (all features)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
execution_mode <- ifelse("--quick" %in% args, "quick",
                  ifelse("--full" %in% args, "full", "standard"))

# Start timer
pipeline_start <- Sys.time()

# Load required libraries
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
})

# INITIALIZATION - Load config first for dynamic banner
cat("🔍 Initializing pipeline configuration...\n")
config <- yaml::read_yaml("config.yml")

# Get gene name for dynamic banner
primary_gene <- config$research_target$primary_gene %||% "UNKNOWN"

# Display banner with gene name
cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║      ", toupper(primary_gene), "PIPELINE - COMPLETE EXECUTION", paste(rep(" ", max(0, 18-nchar(primary_gene))), collapse=""), "║\n")

# Show mode and time estimate
if (execution_mode == "quick") {
  cat("║        Mode: QUICK | Estimated Time: 2-3 minutes             ║\n")
} else if (execution_mode == "full") {
  cat("║        Mode: FULL | Estimated Time: 10+ minutes              ║\n")
} else {
  cat("║        Mode: STANDARD | Estimated Time: 5 minutes            ║\n")
}

cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Adjust config based on mode
if (execution_mode == "quick") {
  # Disable enhanced features for speed
  config$dynamic_features$enabled <- FALSE
  config$dynamic_features$auto_download <- FALSE
  config$dynamic_features$dataset_discovery <- FALSE
  config$dynamic_features$pathway_analysis <- FALSE
  config$dynamic_features$gene_family_discovery <- FALSE
  config$dynamic_features$literature_mining <- FALSE
  config$dynamic_features$drug_target_prediction <- FALSE
  cat("⚡ Quick mode: Enhanced features disabled for speed\n")
  
} else if (execution_mode == "full") {
  # Enable everything
  config$dynamic_features$enabled <- TRUE
  cat("🔥 Full mode: All features enabled\n")
  
} else {
  # Standard mode - selective features
  config$dynamic_features$enabled <- TRUE
  config$dynamic_features$auto_download <- FALSE
  config$dynamic_features$dataset_discovery <- FALSE
  config$dynamic_features$gene_family_discovery <- TRUE
  cat("📊 Standard mode: Core features + gene family discovery\n")
}

# Write updated config
yaml::write_yaml(config, "config.yml")
cat("✅ Configuration updated for", toupper(execution_mode), "mode\n\n")

# Load the real pipeline orchestrator
cat("🚀 Loading pipeline orchestrator...\n")
source("scripts/pipeline_orchestrator.R")

# Execute the real pipeline
cat("🔄 Starting CAMK2D analysis pipeline...\n")
cat("═══════════════════════════════════════════════════════════════\n")

pipeline_result <- execute_dynamic_pipeline(
  config_file = "config.yml",
  force_rerun = FALSE  # Use checkpoints when available
)

cat("═══════════════════════════════════════════════════════════════\n")

# Generate interactive documentation
cat("🌐 GENERATING INTERACTIVE DOCUMENTATION\n")
cat("═══════════════════════════════════════════════════════════════\n")

tryCatch({
  # Source from current working directory 
  if (file.exists("generate_interactive_documentation.R")) {
    source("generate_interactive_documentation.R")
  } else {
    stop("Cannot find generate_interactive_documentation.R in current directory: ", getwd())
  }
  
  cat("📄 Converting 70+ flowcharts to interactive HTML...\n")
  
  # Ensure output directory exists
  output_dir <- "output/current"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("📁 Created output directory:", output_dir, "\n")
  }
  
  doc_result <- generate_interactive_documentation(
    input_file = "Technical_Documentation_CAMK2D_Pipeline.md",
    output_file = "output/current/Interactive_Technical_Documentation.html",
    title = "CAMK2D Pipeline - Interactive Technical Documentation",
    include_results_summary = TRUE
  )
  
  if (doc_result) {
    cat("✅ Interactive documentation generated successfully\n")
  } else {
    cat("⚠️  Documentation generation had issues but completed\n")
  }
}, error = function(e) {
  cat("❌ Documentation generation failed:", e$message, "\n")
})

# Calculate total execution time
pipeline_end <- Sys.time()
total_duration <- difftime(pipeline_end, pipeline_start, units = "secs")
total_duration_min <- round(total_duration / 60, 2)

# Check generated files
analysis_report <- paste0("output/current/", toupper(primary_gene), "_Analysis_Report.html")
doc_report <- "output/current/Interactive_Technical_Documentation.html"

cat("═══════════════════════════════════════════════════════════════\n")
cat("📊 COMPLETE PIPELINE EXECUTION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Report status
if (file.exists(analysis_report)) {
  size_kb <- round(file.info(analysis_report)$size / 1024, 1)
  cat("✅ Analysis Report: Generated (", size_kb, "KB)\n")
  cat("   📄", analysis_report, "\n")
} else {
  cat("❌ Analysis Report: Not found\n")
}

if (file.exists(doc_report)) {
  size_kb <- round(file.info(doc_report)$size / 1024, 1)
  cat("✅ Interactive Documentation: Generated (", size_kb, "KB)\n")
  cat("   🌐", doc_report, "\n")
} else {
  cat("❌ Interactive Documentation: Not generated\n")
}

cat("\n📊 Final Execution Summary:\n")
cat("   Mode:", toupper(execution_mode), "\n")
cat("   Pipeline Success:", if (pipeline_result$success) "YES ✅" else "NO ❌", "\n")
cat("   Steps Completed:", pipeline_result$steps_completed, "\n")
cat("   Steps Failed:", pipeline_result$steps_failed, "\n")
cat("   Total Time:", total_duration_min, "minutes\n")
cat("   Pipeline Time:", pipeline_result$execution_time_minutes, "minutes\n")

if (length(pipeline_result$warnings) > 0) {
  cat("   Warnings:", length(pipeline_result$warnings), "\n")
}

if (pipeline_result$success) {
  cat("\n🎉 PIPELINE EXECUTION COMPLETED SUCCESSFULLY!\n")
  
  if (execution_mode == "quick") {
    cat("\n⚡ Quick mode complete! Enhanced features were skipped for speed.\n")
    cat("   For full analysis with all features, run: Rscript run_pipeline_complete.R --full\n")
  }
  
  cat("\n🎯 Next Steps:\n")
  cat("   1. Open the Analysis Report in your browser\n")
  cat("   2. Explore the Interactive Documentation with 70+ flowcharts\n")
  cat("   3. Review detailed results in output/current/\n")
  
} else {
  cat("\n❌ PIPELINE EXECUTION FAILED\n")
  cat("Check output/logs/ for detailed error information\n")
}

cat("\n✨ Complete pipeline execution finished!\n\n")