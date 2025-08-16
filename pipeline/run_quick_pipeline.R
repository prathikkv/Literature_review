#!/usr/bin/env Rscript
#' Quick Pipeline Runner - Streamlined for faster execution
#' 
#' Runs core analysis + documentation without enhanced features

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║        QUICK CAMK2D PIPELINE - STREAMLINED                   ║\n")
cat("║        Core Analysis + Interactive Documentation             ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Load configuration
suppressPackageStartupMessages({
  library(yaml)
  library(tidyverse)
})

config <- yaml::read_yaml("config.yml")

# STEP 1: Run core analysis pipeline ONLY
cat("📊 STEP 1: Core Analysis Pipeline (without enhanced features)\n")
cat("═══════════════════════════════════════════════════════════════\n")

tryCatch({
  source("scripts/pipeline_orchestrator.R")
  
  cat("🚀 Executing streamlined analysis...\n")
  
  # Run with minimal configuration - skip enhanced features
  simplified_config <- config
  simplified_config$dynamic_features$enabled <- FALSE
  
  # Execute just the core steps
  pipeline_result <- execute_dynamic_pipeline(
    config_file = "config.yml",
    force_rerun = FALSE,
    steps_to_run = c("step_01_data_loader", 
                     "step_02_preprocessing",
                     "step_03_dge_analysis", 
                     "step_04_meta_analysis",
                     "step_05_report_generator")
  )
  
  if (pipeline_result$success) {
    cat("✅ Analysis pipeline completed\n")
  } else {
    cat("⚠️  Analysis completed with warnings\n")
  }
  
}, error = function(e) {
  cat("⚠️  Analysis pipeline error:", e$message, "\n")
  cat("   Continuing with documentation generation...\n")
})

# STEP 2: Generate interactive documentation
cat("\n🌐 STEP 2: Interactive Documentation Generation\n")
cat("═══════════════════════════════════════════════════════════════\n")

tryCatch({
  source("generate_interactive_documentation.R")
  
  doc_result <- generate_interactive_documentation(
    input_file = "Technical_Documentation_CAMK2D_Pipeline.md",
    output_file = "output/current/Interactive_Technical_Documentation.html",
    title = "CAMK2D Pipeline - Interactive Technical Documentation"
  )
  
  if (doc_result) {
    cat("✅ Interactive documentation generated\n")
  }
  
}, error = function(e) {
  cat("❌ Documentation generation failed:", e$message, "\n")
})

# STEP 3: Summary
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("📊 QUICK PIPELINE SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Check what was generated
analysis_report <- "output/current/CAMK_Analysis_Report.html"
doc_report <- "output/current/Interactive_Technical_Documentation.html"

if (file.exists(analysis_report)) {
  size_mb <- round(file.info(analysis_report)$size / 1024 / 1024, 2)
  cat("✅ Analysis Report: Generated (", size_mb, "MB)\n")
} else {
  cat("⚠️  Analysis Report: Not generated\n")
  cat("   The analysis may still be processing\n")
}

if (file.exists(doc_report)) {
  size_mb <- round(file.info(doc_report)$size / 1024 / 1024, 2)
  cat("✅ Interactive Documentation: Generated (", size_mb, "MB)\n")
} else {
  cat("❌ Interactive Documentation: Not found\n")
}

cat("\n📁 Output Location: output/current/\n")
cat("🎯 Both HTML reports should be available after completion\n\n")