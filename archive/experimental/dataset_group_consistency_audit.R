#!/usr/bin/env Rscript
#' Dataset Group Consistency Audit
#' 
#' Critical investigation: Are all datasets using consistent Disease vs Control comparisons?
#' This addresses the concerning CAMK2D contradiction across datasets

library(limma)
source("scripts/enhanced_group_detection.R")

cat("=== DATASET GROUP CONSISTENCY AUDIT ===\n\n")
cat("INVESTIGATING: Why only 1/3 datasets show CAMK2D upregulation\n")
cat("HYPOTHESIS: Datasets may have flipped group assignments or inconsistent methodology\n\n")

# List of datasets to audit
datasets <- c("GSE57338", "GSE41177", "GSE79768")
audit_results <- list()

for (dataset_id in datasets) {
  cat("=====================================\n")
  cat("AUDITING:", dataset_id, "\n") 
  cat("=====================================\n")
  
  # Find and load dataset
  dataset_files <- list.files("cache", pattern = paste0(dataset_id, "_processed.rds"), 
                             recursive = TRUE, full.names = TRUE)
  
  if (length(dataset_files) == 0) {
    cat("ERROR: Dataset file not found\n\n")
    next
  }
  
  dataset <- readRDS(dataset_files[1])
  
  if (!dataset$success) {
    cat("ERROR: Dataset processing failed\n\n")
    next
  }
  
  # Basic dataset info
  cat("DATASET INFO:\n")
  cat("- Samples:", dataset$n_samples, "\n")
  cat("- Genes:", dataset$n_genes, "\n")
  cat("- Platform:", dataset$dataset_info$platform, "\n")
  cat("- File:", dataset_files[1], "\n\n")
  
  # Examine phenotype data
  pheno_data <- dataset$phenotype_data
  cat("PHENOTYPE DATA ANALYSIS:\n")
  cat("- Phenotype columns:", ncol(pheno_data), "\n")
  
  # Key columns to examine
  key_columns <- c("source_name_ch1", "title", "description", "characteristics_ch1", 
                   "disease status:ch1", "heart failure:ch1")
  
  available_key_cols <- key_columns[key_columns %in% colnames(pheno_data)]
  cat("- Available key columns:", paste(available_key_cols, collapse = ", "), "\n\n")
  
  # Examine source_name_ch1 (most informative)
  if ("source_name_ch1" %in% colnames(pheno_data)) {
    source_names <- pheno_data$source_name_ch1
    unique_sources <- unique(source_names)
    
    cat("SOURCE NAME ANALYSIS:\n")
    cat("- Unique source types:", length(unique_sources), "\n")
    for (i in 1:length(unique_sources)) {
      count <- sum(source_names == unique_sources[i])
      cat(sprintf("  %d. %s (%d samples)\n", i, unique_sources[i], count))
    }
    cat("\n")
    
    # Disease/Control identification
    disease_keywords <- c("cardiomyopathy", "CMP", "heart failure", "HF", "ischemic", "MI", "infarct", "disease", "pathological")
    control_keywords <- c("control", "normal", "non-failing", "healthy", "donor")
    
    disease_samples <- rep(FALSE, length(source_names))
    control_samples <- rep(FALSE, length(source_names))
    
    for (keyword in disease_keywords) {
      disease_samples <- disease_samples | grepl(keyword, source_names, ignore.case = TRUE)
    }
    
    for (keyword in control_keywords) {
      control_samples <- control_samples | grepl(keyword, source_names, ignore.case = TRUE)
    }
    
    cat("DISEASE/CONTROL CLASSIFICATION:\n")
    cat("- Disease samples:", sum(disease_samples), "\n")
    cat("- Control samples:", sum(control_samples), "\n")
    cat("- Unclassified samples:", sum(!disease_samples & !control_samples), "\n")
    cat("- Overlap (both disease & control keywords):", sum(disease_samples & control_samples), "\n\n")
    
    if (sum(disease_samples & control_samples) > 0) {
      cat("WARNING: Some samples match both disease and control keywords!\n")
      overlapping_sources <- unique(source_names[disease_samples & control_samples])
      for (source in overlapping_sources) {
        cat("  - Overlapping:", source, "\n")
      }
      cat("\n")
    }
    
  } else {
    cat("WARNING: source_name_ch1 not available\n\n")
  }
  
  # Run enhanced group detection
  cat("ENHANCED GROUP DETECTION ANALYSIS:\n")
  detected_groups <- enhanced_auto_detect_groups(dataset)
  
  if (!is.null(detected_groups)) {
    cat("- Pattern detected:", detected_groups$pattern_type, "\n")
    cat("- Groups found:", length(unique(detected_groups$groups)), "\n")
    
    group_table <- table(detected_groups$groups)
    print(group_table)
    
    # Check factor levels (important for contrast direction)
    if (is.factor(detected_groups$groups)) {
      cat("- Factor levels:", paste(levels(detected_groups$groups), collapse = ", "), "\n")
      cat("- Reference level (first):", levels(detected_groups$groups)[1], "\n")
    }
    
    # Examine which samples are in each group
    group_names <- names(group_table)
    for (group_name in group_names) {
      group_indices <- which(detected_groups$groups == group_name)
      cat(sprintf("- %s group samples:\n", group_name))
      
      if ("source_name_ch1" %in% colnames(pheno_data)) {
        group_sources <- unique(pheno_data$source_name_ch1[group_indices])
        for (source in group_sources) {
          source_count <- sum(pheno_data$source_name_ch1[group_indices] == source)
          cat(sprintf("    %s (%d samples)\n", source, source_count))
        }
      }
      cat("\n")
    }
    
  } else {
    cat("ERROR: Group detection failed\n\n")
  }
  
  # Test limma design matrix construction
  cat("LIMMA DESIGN MATRIX TEST:\n")
  if (!is.null(detected_groups)) {
    # Simulate design matrix creation
    groups_vector <- detected_groups$groups
    
    if (is.character(groups_vector)) {
      groups_vector <- factor(groups_vector)
    }
    
    design <- model.matrix(~ groups_vector)
    cat("- Design matrix dimensions:", nrow(design), "x", ncol(design), "\n")
    cat("- Column names:", paste(colnames(design), collapse = ", "), "\n")
    
    # Check contrast interpretation
    if (ncol(design) >= 2) {
      contrast_col <- colnames(design)[2]
      cat("- Contrast column:", contrast_col, "\n")
      
      # Determine what positive logFC means
      ref_level <- levels(groups_vector)[1]
      alt_level <- levels(groups_vector)[2]
      
      cat(sprintf("- Positive logFC means: %s > %s\n", alt_level, ref_level))
      cat(sprintf("- Negative logFC means: %s < %s\n", alt_level, ref_level))
      
      # Is this the expected direction?
      disease_in_alt <- grepl("disease|cardiomyopathy|CMP|failure|HF|ischemic", alt_level, ignore.case = TRUE)
      control_in_ref <- grepl("control|normal|non-failing|healthy|donor", ref_level, ignore.case = TRUE)
      
      correct_direction <- disease_in_alt && control_in_ref
      cat("- Expected direction (Disease > Control):", correct_direction, "\n")
      
      if (!correct_direction) {
        cat("⚠️  WARNING: Contrast direction may be FLIPPED!\n")
        cat("   Expected: Disease vs Control (positive = UP in disease)\n")
        cat(sprintf("   Actual: %s vs %s\n", alt_level, ref_level))
      }
    }
    
  } else {
    cat("Cannot create design matrix - no groups detected\n")
  }
  
  cat("\n")
  
  # Store audit results
  audit_results[[dataset_id]] <- list(
    dataset_id = dataset_id,
    success = !is.null(detected_groups),
    groups = if (!is.null(detected_groups)) detected_groups$groups else NULL,
    pattern_type = if (!is.null(detected_groups)) detected_groups$pattern_type else NA,
    phenotype_data = pheno_data
  )
}

# SUMMARY ANALYSIS
cat("=========================================\n")
cat("CROSS-DATASET CONSISTENCY ANALYSIS\n")
cat("=========================================\n\n")

# Compare group detection patterns
cat("GROUP DETECTION PATTERNS:\n")
for (dataset_id in names(audit_results)) {
  result <- audit_results[[dataset_id]]
  if (result$success) {
    cat(sprintf("- %s: %s\n", dataset_id, result$pattern_type))
  } else {
    cat(sprintf("- %s: FAILED\n", dataset_id))
  }
}
cat("\n")

# Check for consistent biological hypothesis
cat("BIOLOGICAL HYPOTHESIS CONSISTENCY CHECK:\n")

# Load the DGE results to see current directions
dge_results <- read.csv("output/CAMK_focused_DGE_all_datasets_CORRECTED.csv")
camk2d_results <- dge_results[dge_results$Gene_Symbol == "CAMK2D", ]

cat("Current CAMK2D results:\n")
for (i in 1:nrow(camk2d_results)) {
  dataset <- camk2d_results$Dataset[i]
  logfc <- round(camk2d_results$logFC[i], 4)
  direction <- ifelse(logfc > 0, "UP", "DOWN")
  cat(sprintf("- %s: %s in Disease (logFC = %7.4f)\n", dataset, direction, logfc))
}
cat("\n")

# Literature expectation vs results
cat("LITERATURE EXPECTATION vs ACTUAL RESULTS:\n")
cat("Expected: CAMK2D UPREGULATED in cardiovascular disease\n")
upregulated_count <- sum(camk2d_results$logFC > 0)
downregulated_count <- sum(camk2d_results$logFC < 0)

cat(sprintf("Actual: %d datasets UP, %d datasets DOWN\n", upregulated_count, downregulated_count))

if (upregulated_count < downregulated_count) {
  cat("🚨 CRITICAL ISSUE: Majority of datasets contradict literature expectation!\n")
  cat("   This suggests potential methodological inconsistency\n")
} else if (upregulated_count == downregulated_count) {
  cat("⚠️  MIXED RESULTS: Equal split between UP and DOWN\n")
  cat("   This suggests methodological issues or biological complexity\n")
} else {
  cat("✅ MOSTLY CONSISTENT: Majority align with literature expectation\n")
}

cat("\n=== AUDIT COMPLETE ===\n")
cat("RECOMMENDATION: Review individual dataset group assignments for methodological corrections\n")