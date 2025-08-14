#!/usr/bin/env Rscript
#' Repository Structure Validation
#' 
#' This script validates that the cleaned repository structure is working correctly

cat("═══════════════════════════════════════════════════════════════\n")
cat("🔍 REPOSITORY STRUCTURE VALIDATION\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Check key files exist
key_files <- c(
  "reports/CAMK_Analysis_Professional_Report.Rmd",
  "reports/CAMK_Professional_Analysis_Report.html", 
  "scripts/core/comprehensive_6_dataset_pipeline.R",
  "scripts/core/fixed_meta_analysis.R",
  "output/current/CAMK_meta_analysis_FINAL.csv",
  "output/current/CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv",
  "cache/microarray/GSE57338_processed.rds",
  "README.md"
)

cat("📋 Checking key files:\n")
all_good <- TRUE
for (file in key_files) {
  if (file.exists(file)) {
    cat("✅", file, "\n")
  } else {
    cat("❌", file, "- MISSING\n")
    all_good <- FALSE
  }
}

# Check results are loadable
cat("\n📊 Validating key results:\n")
tryCatch({
  meta_results <- read.csv("output/current/CAMK_meta_analysis_FINAL.csv")
  camk2d <- meta_results[meta_results$Gene == "CAMK2D", ]
  
  if (nrow(camk2d) == 1) {
    cat("✅ CAMK2D meta-analysis result found\n")
    cat(sprintf("   logFC: %.4f, p-value: %.2e\n", camk2d$Combined_logFC, camk2d$Combined_P_Value))
  } else {
    cat("❌ CAMK2D result not found\n")
    all_good <- FALSE
  }
}, error = function(e) {
  cat("❌ Error loading meta-analysis results:", e$message, "\n")
  all_good <- FALSE
})

# Check dataset results
tryCatch({
  individual_results <- read.csv("output/current/CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv")
  datasets <- unique(individual_results$Dataset)
  cat("✅ Individual dataset results loaded\n")
  cat("   Datasets:", paste(datasets, collapse = ", "), "\n")
  cat("   Total gene-dataset combinations:", nrow(individual_results), "\n")
}, error = function(e) {
  cat("❌ Error loading individual results:", e$message, "\n")
  all_good <- FALSE
})

# Check archive structure
cat("\n📁 Checking archive organization:\n")
archive_dirs <- c("archive/experimental", "archive/old_reports", "archive/old_analysis")
for (dir in archive_dirs) {
  if (dir.exists(dir)) {
    file_count <- length(list.files(dir, recursive = TRUE))
    cat("✅", dir, "- ", file_count, "files archived\n")
  } else {
    cat("❌", dir, "- Missing\n")
    all_good <- FALSE
  }
}

# Summary
cat("\n═══════════════════════════════════════════════════════════════\n")
if (all_good) {
  cat("🏆 VALIDATION SUCCESSFUL - Repository structure is clean and functional\n")
} else {
  cat("❌ VALIDATION FAILED - Some issues detected\n")
}
cat("═══════════════════════════════════════════════════════════════\n")

# Quick stats
cat("\n📊 Repository cleanup summary:\n")
core_scripts <- length(list.files("scripts/core", pattern = "\\.R$"))
util_scripts <- length(list.files("scripts/utilities", pattern = "\\.R$")) 
archived_scripts <- length(list.files("archive/experimental", pattern = "\\.R$"))
archived_reports <- length(list.files("archive/old_reports"))

cat("Core scripts:", core_scripts, "\n")
cat("Utility scripts:", util_scripts, "\n") 
cat("Archived scripts:", archived_scripts, "\n")
cat("Archived reports:", archived_reports, "\n")
cat("Total reduction: ~", archived_scripts + archived_reports, "files moved to archive\n")

# Validation complete