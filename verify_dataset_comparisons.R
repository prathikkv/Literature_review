#!/usr/bin/env Rscript
#' Verify Dataset Comparisons - Identify True Healthy vs Disease Studies
#' 
#' Examines each dataset to determine what comparisons are actually being made

cat("🔍 VERIFYING DATASET COMPARISONS\n")
cat("===============================\n\n")

# Load the CAMK analysis results to see current comparisons
camk_results <- readRDS("output/CAMK_focused_analysis_results.rds")

cat("📊 Current datasets in CAMK analysis:\n")

for (dataset_id in names(camk_results)) {
  
  cat("\n📋 Dataset:", dataset_id, "\n")
  cat("   ", paste(rep("=", nchar(dataset_id) + 9), collapse = ""), "\n")
  
  dataset_info <- camk_results[[dataset_id]]
  
  # Check if groups were detected
  if (is.null(dataset_info$groups)) {
    cat("   ❌ No groups detected\n")
    next
  }
  
  # Display group information
  if (is.list(dataset_info$groups)) {
    groups <- dataset_info$groups$groups
    pattern_type <- dataset_info$groups$pattern_type
    column_used <- dataset_info$groups$column
    
    cat("   🔬 Group detection successful\n")
    cat("   📊 Pattern type:", pattern_type, "\n")
    cat("   📝 Column used:", column_used, "\n")
    
    if (length(groups) > 0) {
      group_table <- table(groups)
      cat("   👥 Groups found:\n")
      for (group_name in names(group_table)) {
        cat(sprintf("      %-12s: %3d samples\n", group_name, group_table[group_name]))
      }
      
      # Determine comparison type
      comparison_type <- "Unknown"
      clinical_relevance <- "❓ Unknown"
      
      if (pattern_type == "healthy_vs_disease" || pattern_type == "healthy_vs_disease_general") {
        comparison_type <- "✅ HEALTHY vs DISEASE"
        clinical_relevance <- "🎯 HIGH - Direct disease relevance"
      } else if (pattern_type == "AF_SR") {
        comparison_type <- "⚠️ AF vs SR"
        clinical_relevance <- "🟡 MODERATE - Disease subtype comparison"
      } else if (pattern_type == "traditional") {
        # Need to examine the group names to determine
        group_names <- names(group_table)
        if (any(grepl("healthy|control|normal|non.failing", group_names, ignore.case = TRUE))) {
          comparison_type <- "✅ HEALTHY vs DISEASE"
          clinical_relevance <- "🎯 HIGH - Direct disease relevance"
        } else if (any(grepl("AF|SR", group_names, ignore.case = TRUE))) {
          comparison_type <- "⚠️ AF vs SR"
          clinical_relevance <- "🟡 MODERATE - Disease subtype comparison"
        } else {
          comparison_type <- "❌ DISEASE vs DISEASE"
          clinical_relevance <- "🔴 LOW - Not healthy vs disease"
        }
      }
      
      cat("   🔍 Comparison type:", comparison_type, "\n")
      cat("   📈 Clinical relevance:", clinical_relevance, "\n")
      
      # Check sample size adequacy
      min_group_size <- min(group_table)
      total_samples <- sum(group_table)
      
      if (total_samples >= 50 && min_group_size >= 10) {
        power_assessment <- "🟢 EXCELLENT - Good statistical power"
      } else if (total_samples >= 20 && min_group_size >= 5) {
        power_assessment <- "🟡 MODERATE - Adequate for analysis"
      } else {
        power_assessment <- "🔴 LOW - Limited statistical power"
      }
      
      cat("   📊 Statistical power:", power_assessment, "\n")
      
      # CAMK gene availability
      if (!is.null(dataset_info$camk_genes_present)) {
        camk_count <- length(dataset_info$camk_genes_present)
        cat(sprintf("   🧬 CAMK genes: %d/11 available\n", camk_count))
      }
      
    }
  }
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("🎯 RECOMMENDATIONS FOR PROPER CAMK DYSREGULATION ANALYSIS:\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Categorize datasets
healthy_vs_disease <- c()
disease_subtype <- c()  
insufficient_data <- c()

for (dataset_id in names(camk_results)) {
  dataset_info <- camk_results[[dataset_id]]
  
  if (!is.null(dataset_info$groups) && is.list(dataset_info$groups)) {
    pattern_type <- dataset_info$groups$pattern_type
    groups <- dataset_info$groups$groups
    group_table <- table(groups)
    
    # Check if this is healthy vs disease
    if (pattern_type == "healthy_vs_disease" || pattern_type == "healthy_vs_disease_general" ||
        (pattern_type == "traditional" && any(grepl("healthy|control|normal|non.failing", names(group_table), ignore.case = TRUE)))) {
      
      if (sum(group_table) >= 20 && min(group_table) >= 5) {
        healthy_vs_disease <- c(healthy_vs_disease, dataset_id)
      } else {
        insufficient_data <- c(insufficient_data, dataset_id)
      }
      
    } else if (pattern_type == "AF_SR" || any(grepl("AF|SR", names(group_table), ignore.case = TRUE))) {
      disease_subtype <- c(disease_subtype, dataset_id)
    } else {
      insufficient_data <- c(insufficient_data, dataset_id)
    }
  } else {
    insufficient_data <- c(insufficient_data, dataset_id)
  }
}

cat("🎯 PRIMARY ANALYSIS (Healthy vs Disease):\n")
if (length(healthy_vs_disease) > 0) {
  cat("   ✅ Recommended datasets:\n")
  for (ds in healthy_vs_disease) {
    dataset_info <- camk_results[[ds]]
    if (!is.null(dataset_info$groups) && is.list(dataset_info$groups)) {
      groups <- dataset_info$groups$groups
      group_table <- table(groups)
      total <- sum(group_table)
      cat(sprintf("      • %-12s: %d samples (%s)\n", ds, total, paste(names(group_table), collapse = " vs ")))
    }
  }
} else {
  cat("   ❌ No healthy vs disease datasets identified!\n")
  cat("   🔧 Need to check GSE57338 - it should have healthy controls\n")
}

cat("\n🟡 SECONDARY ANALYSIS (Disease Subtypes):\n")
if (length(disease_subtype) > 0) {
  cat("   ⚠️ Disease subtype datasets (separate analysis):\n")
  for (ds in disease_subtype) {
    dataset_info <- camk_results[[ds]]
    if (!is.null(dataset_info$groups) && is.list(dataset_info$groups)) {
      groups <- dataset_info$groups$groups
      group_table <- table(groups)
      total <- sum(group_table)
      cat(sprintf("      • %-12s: %d samples (%s)\n", ds, total, paste(names(group_table), collapse = " vs ")))
    }
  }
} else {
  cat("   No disease subtype datasets\n")
}

cat("\n❌ EXCLUDED FROM MAIN ANALYSIS:\n")
if (length(insufficient_data) > 0) {
  cat("   Datasets with insufficient data or unclear groups:\n")
  for (ds in insufficient_data) {
    cat("      •", ds, "\n")
  }
} else {
  cat("   None\n")
}

cat("\n💡 NEXT STEPS:\n")
cat("1. Verify GSE57338 (313 samples) has healthy controls - it should be our PRIMARY dataset\n")
cat("2. Re-run analysis focusing ONLY on healthy vs disease comparisons\n")
cat("3. Create separate meta-analysis excluding disease vs disease studies\n")
cat("4. Interpret results in context of cardiovascular disease pathophysiology\n")

cat("\n✅ Dataset verification completed!\n")