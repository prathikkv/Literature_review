#!/usr/bin/env Rscript
#' Step 02: Preprocessing
#' 
#' Handles group detection, quality control, and data preparation for analysis
#' Extracts group detection logic from enhanced_group_detection_corrected.R

source("scripts/utilities/step_interface.R")

#' Execute Preprocessing Step
#'
#' @param step_name Name of this step (should be "step_02_preprocessing")
#' @param input_data Output from step_01_data_loader
#' @param config Full pipeline configuration
#' @param checkpoint_dir Directory for saving checkpoints
#' @return Step result with processed datasets
step_02_preprocessing <- function(step_name, input_data, config, checkpoint_dir = "output/checkpoints") {
  
  cat("🔍 STEP 02: PREPROCESSING\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Processing datasets for group detection and quality control...\n\n")
  
  # Validate input from data loader
  validate_step_input(
    input_data = input_data,
    required_fields = c("loaded_datasets", "dataset_summary"),
    step_name = step_name
  )
  
  loaded_datasets <- input_data$loaded_datasets
  original_summary <- input_data$dataset_summary
  
  if (length(loaded_datasets) == 0) {
    return(create_step_result(
      success = FALSE,
      error_message = "No datasets provided for preprocessing",
      step_name = step_name
    ))
  }
  
  cat("📋 PREPROCESSING CONFIGURATION:\n")
  cat("Datasets to process:", length(loaded_datasets), "\n")
  cat("Group detection method: Enhanced biological logic\n")
  cat("Quality control: Enabled\n\n")
  
  # Process each dataset
  processed_datasets <- list()
  preprocessing_summary <- data.frame()
  failed_preprocessing <- character(0)
  warnings_list <- character(0)
  
  for (dataset_id in names(loaded_datasets)) {
    dataset <- loaded_datasets[[dataset_id]]
    dataset_config <- config$datasets$active_datasets[[dataset_id]]
    
    cat("═══════════════════════════════════════════════════════════════\n")
    cat("🔍 PREPROCESSING:", dataset_id, "\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    
    # Perform preprocessing for this dataset
    preprocessing_result <- preprocess_single_dataset(dataset_id, dataset, dataset_config, config)
    
    if (preprocessing_result$success) {
      processed_datasets[[dataset_id]] <- preprocessing_result$processed_data
      
      # Add to summary
      preprocessing_summary <- rbind(preprocessing_summary, data.frame(
        Dataset = dataset_id,
        Priority = dataset_config$priority,
        Original_Samples = dataset$n_samples %||% ncol(dataset$expression_matrix),
        Processed_Samples = preprocessing_result$final_sample_count,
        Control_Samples = preprocessing_result$group_counts[["Control"]] %||% 0,
        Disease_Samples = preprocessing_result$group_counts[["Disease"]] %||% 0,
        Total_Genes = nrow(dataset$expression_matrix),
        CAMK_Genes = preprocessing_result$camk_genes_detected,
        Group_Detection = preprocessing_result$group_detection_method,
        QC_Flags = preprocessing_result$qc_flags %||% 0,
        Status = "SUCCESS",
        stringsAsFactors = FALSE
      ))
      
      cat("✅ SUCCESS: Preprocessing completed\n")
      cat("📊 Groups detected:", preprocessing_result$group_detection_method, "\n")
      cat("👥 Sample distribution:", 
          "Control =", preprocessing_result$group_counts[["Control"]] %||% 0,
          "| Disease =", preprocessing_result$group_counts[["Disease"]] %||% 0, "\n")
      cat("🧬 CAMK genes detected:", preprocessing_result$camk_genes_detected, "\n")
      
      # Collect warnings
      if (length(preprocessing_result$warnings) > 0) {
        warnings_list <- c(warnings_list, paste(dataset_id, ":", preprocessing_result$warnings))
      }
      
    } else {
      failed_preprocessing <- c(failed_preprocessing, dataset_id)
      
      # Add failure to summary
      preprocessing_summary <- rbind(preprocessing_summary, data.frame(
        Dataset = dataset_id,
        Priority = dataset_config$priority,
        Original_Samples = dataset$n_samples %||% ncol(dataset$expression_matrix),
        Processed_Samples = NA,
        Control_Samples = NA,
        Disease_Samples = NA,
        Total_Genes = nrow(dataset$expression_matrix),
        CAMK_Genes = NA,
        Group_Detection = "FAILED",
        QC_Flags = NA,
        Status = paste("FAILED:", preprocessing_result$reason),
        stringsAsFactors = FALSE
      ))
      
      cat("❌ FAILED: Preprocessing failed -", preprocessing_result$reason, "\n")
    }
    cat("\n")
  }
  
  # Generate preprocessing summary
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📊 PREPROCESSING SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  successful_count <- length(processed_datasets)
  failed_count <- length(failed_preprocessing)
  
  cat("✅ Successfully preprocessed:", successful_count, "\n")
  cat("❌ Failed preprocessing:", failed_count, "\n")
  
  if (length(warnings_list) > 0) {
    cat("⚠️  Warnings generated:", length(warnings_list), "\n")
  }
  
  if (successful_count > 0) {
    cat("\n📋 PREPROCESSED DATASETS:\n")
    successful_datasets <- preprocessing_summary[preprocessing_summary$Status == "SUCCESS", ]
    total_samples <- sum(successful_datasets$Processed_Samples)
    total_camk_detections <- sum(successful_datasets$CAMK_Genes)
    
    for (i in 1:nrow(successful_datasets)) {
      ds <- successful_datasets[i, ]
      cat(sprintf("  🔍 %-10s: %3d samples (%2d control, %2d disease), %2d CAMK genes\n",
                  ds$Dataset, ds$Processed_Samples, ds$Control_Samples, ds$Disease_Samples, ds$CAMK_Genes))
    }
    
    cat(sprintf("\n📊 TOTALS: %d samples across %d datasets, %d total CAMK gene detections\n",
                total_samples, successful_count, total_camk_detections))
  }
  
  if (failed_count > 0) {
    cat("\n❌ PREPROCESSING FAILURES:\n")
    failed_summary <- preprocessing_summary[preprocessing_summary$Status != "SUCCESS", ]
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
      error_message = paste("Insufficient datasets preprocessed:", successful_count, "< minimum:", min_datasets),
      step_name = step_name,
      warnings = warnings_list,
      metadata = list(
        datasets_processed = successful_count,
        datasets_failed = failed_count,
        preprocessing_summary = preprocessing_summary
      )
    ))
  }
  
  # Create output data structure
  output_data <- list(
    processed_datasets = processed_datasets,
    preprocessing_summary = preprocessing_summary,
    original_summary = original_summary,
    processing_stats = list(
      total_input_datasets = length(loaded_datasets),
      successfully_processed = successful_count,
      failed_processing = failed_count,
      total_samples = sum(preprocessing_summary$Processed_Samples[!is.na(preprocessing_summary$Processed_Samples)]),
      total_camk_detections = sum(preprocessing_summary$CAMK_Genes[!is.na(preprocessing_summary$CAMK_Genes)])
    )
  )
  
  # Validate output
  validate_step_output(
    output_data = output_data,
    required_fields = c("processed_datasets", "preprocessing_summary", "processing_stats"),
    step_name = step_name
  )
  
  cat("\n🎉 PREPROCESSING COMPLETED SUCCESSFULLY\n")
  cat("Ready for differential expression analysis...\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(create_step_result(
    success = TRUE,
    output_data = output_data,
    step_name = step_name,
    warnings = warnings_list,
    metadata = list(
      datasets_processed = successful_count,
      total_samples = sum(preprocessing_summary$Processed_Samples[!is.na(preprocessing_summary$Processed_Samples)]),
      quality_warnings = length(warnings_list)
    )
  ))
}

#' Preprocess Single Dataset
#'
#' Performs group detection and quality control for a single dataset
#' @param dataset_id Dataset identifier
#' @param dataset Loaded dataset object
#' @param dataset_config Configuration for this dataset
#' @param config Full pipeline configuration
#' @return List with preprocessing results
preprocess_single_dataset <- function(dataset_id, dataset, dataset_config, config) {
  
  warnings_generated <- character(0)
  
  # Perform group detection
  cat("🎯 GROUP DETECTION:\n")
  
  group_detection_result <- enhanced_auto_detect_groups_config_driven(dataset, config)
  
  if (is.null(group_detection_result) || !group_detection_result$success) {
    return(list(
      success = FALSE,
      reason = "group_detection_failed",
      warnings = warnings_generated
    ))
  }
  
  cat("✅ Groups successfully detected:\n")
  cat("Pattern type:", group_detection_result$pattern_type, "\n")
  cat("Reference group (Control):", levels(group_detection_result$groups)[1], "\n")
  cat("Comparison group (Disease):", levels(group_detection_result$groups)[2], "\n")
  
  group_counts <- table(group_detection_result$groups)
  print(group_counts)
  
  # Validate group sizes
  min_control <- config$analysis$group_detection$validation$min_control_samples %||% 3
  min_disease <- config$analysis$group_detection$validation$min_disease_samples %||% 3
  
  if (group_counts[["Control"]] < min_control) {
    warnings_generated <- c(warnings_generated, paste("Low control sample count:", group_counts[["Control"]]))
  }
  
  if (group_counts[["Disease"]] < min_disease) {
    warnings_generated <- c(warnings_generated, paste("Low disease sample count:", group_counts[["Disease"]]))
  }
  
  # Check for sample imbalance
  max_imbalance_ratio <- config$analysis$quality_control$max_sample_imbalance_ratio %||% 10
  imbalance_ratio <- max(group_counts) / min(group_counts)
  
  if (imbalance_ratio > max_imbalance_ratio) {
    warnings_generated <- c(warnings_generated, paste("Sample imbalance detected. Ratio:", round(imbalance_ratio, 2)))
  }
  
  # Prepare expression data
  expr_matrix <- dataset$expression_matrix
  
  if (!is.null(group_detection_result$sample_indices)) {
    expr_filtered <- expr_matrix[, group_detection_result$sample_indices, drop = FALSE]
    groups_vector <- group_detection_result$groups
  } else {
    expr_filtered <- expr_matrix
    groups_vector <- group_detection_result$groups
  }
  
  cat("📊 Analysis Matrix:", nrow(expr_filtered), "genes x", ncol(expr_filtered), "samples\n")
  
  # CAMK genes detection
  camk_genes <- config$genes$camk_core_genes
  camk_present <- intersect(rownames(expr_filtered), camk_genes)
  
  cat("🧬 CAMK genes detected:", length(camk_present), "/", length(camk_genes), "\n")
  cat("CAMK genes:", paste(camk_present, collapse = ", "), "\n")
  
  # Validate minimum CAMK genes
  priority <- dataset_config$priority
  min_genes_threshold <- config$genes$minimum_genes[[priority]] %||% 5
  
  if (length(camk_present) < min_genes_threshold) {
    warnings_generated <- c(warnings_generated, 
                          paste("Insufficient CAMK genes detected:", length(camk_present), "<", min_genes_threshold))
  }
  
  # Quality control checks
  qc_flags <- 0
  
  # Check for extreme expression values (potential artifacts)
  if (any(is.infinite(expr_filtered)) || any(is.na(expr_filtered))) {
    qc_flags <- qc_flags + 1
    warnings_generated <- c(warnings_generated, "Infinite or NA values detected in expression data")
  }
  
  # Create processed dataset object
  processed_data <- list(
    dataset_id = dataset_id,
    expression_matrix = expr_filtered,
    groups = groups_vector,
    group_counts = as.list(group_counts),
    camk_genes_present = camk_present,
    phenotype_data = dataset$phenotype_data,
    dataset_info = dataset$dataset_info,
    group_detection_info = group_detection_result,
    qc_flags = qc_flags,
    processing_metadata = list(
      preprocessed_at = Sys.time(),
      sample_filtering = !is.null(group_detection_result$sample_indices),
      original_samples = ncol(expr_matrix),
      final_samples = ncol(expr_filtered)
    )
  )
  
  return(list(
    success = TRUE,
    processed_data = processed_data,
    final_sample_count = ncol(expr_filtered),
    group_counts = group_counts,
    camk_genes_detected = length(camk_present),
    group_detection_method = group_detection_result$pattern_type,
    qc_flags = qc_flags,
    warnings = warnings_generated
  ))
}

#' Enhanced Auto-Detect Groups (Configuration-Driven)
#'
#' Configuration-driven version of enhanced group detection
#' @param dataset Processed dataset object
#' @param config Pipeline configuration
#' @return Group detection result
enhanced_auto_detect_groups_config_driven <- function(dataset, config) {
  
  if (!dataset$success || is.null(dataset$phenotype_data)) {
    return(list(success = FALSE, reason = "invalid_dataset"))
  }
  
  pheno_data <- dataset$phenotype_data
  n_samples <- nrow(pheno_data)
  
  cat("SEARCH: Analyzing phenotype data for healthy vs disease patterns...\n")
  cat("   Samples:", n_samples, "\n")
  cat("   Columns:", ncol(pheno_data), "\n")
  
  # Get configuration for group detection
  group_config <- config$analysis$group_detection
  priority_columns <- group_config$priority_columns
  control_patterns <- group_config$control_patterns
  disease_patterns <- group_config$disease_patterns
  validation_config <- group_config$validation
  
  # Try each column in priority order
  for (col_name in priority_columns) {
    if (col_name %in% colnames(pheno_data)) {
      cat("   Checking column:", col_name, "\n")
      
      source_info <- pheno_data[[col_name]]
      unique_values <- unique(source_info)
      cat("     Unique values:", length(unique_values), "\n")
      
      # Classify each unique value
      control_samples <- rep(FALSE, n_samples)
      disease_samples <- rep(FALSE, n_samples)
      
      for (pattern in control_patterns) {
        control_samples <- control_samples | grepl(pattern, source_info, ignore.case = TRUE)
      }
      
      for (pattern in disease_patterns) {
        disease_samples <- disease_samples | grepl(pattern, source_info, ignore.case = TRUE)
      }
      
      # Check for valid group separation
      n_control <- sum(control_samples & !disease_samples)  # Pure control
      n_disease <- sum(disease_samples & !control_samples)  # Pure disease
      n_overlap <- sum(control_samples & disease_samples)   # Ambiguous
      n_unclassified <- sum(!control_samples & !disease_samples)  # Neither
      
      cat(sprintf("     Control samples: %d, Disease samples: %d, Overlap: %d, Unclassified: %d\n",
                  n_control, n_disease, n_overlap, n_unclassified))
      
      # Validate against configuration
      min_control <- validation_config$min_control_samples
      min_disease <- validation_config$min_disease_samples
      max_overlap <- validation_config$max_overlap_samples
      max_unclassified_percent <- validation_config$max_unclassified_percent
      
      if (n_control >= min_control && n_disease >= min_disease && 
          n_overlap <= max_overlap && n_unclassified <= n_samples * max_unclassified_percent) {
        
        # Create groups with biological reference logic
        groups <- rep(NA, n_samples)
        groups[control_samples] <- "Control"
        groups[disease_samples] <- "Disease"
        
        # Handle unclassified samples
        if (n_unclassified > 0) {
          unclassified_indices <- which(!control_samples & !disease_samples)
          
          # Try secondary classification
          for (i in unclassified_indices) {
            sample_info <- source_info[i]
            
            if (grepl("normal|control|wild|WT|baseline", sample_info, ignore.case = TRUE)) {
              groups[i] <- "Control"
            } else if (grepl("treatment|disease|pathology|failure|dysfunction", sample_info, ignore.case = TRUE)) {
              groups[i] <- "Disease"
            }
          }
        }
        
        # Remove any remaining NA values
        valid_samples <- !is.na(groups)
        
        if (sum(valid_samples) >= min_control + min_disease) {
          # Create factor with proper reference level
          groups_factor <- factor(groups[valid_samples], levels = c("Control", "Disease"))
          
          sample_indices <- if (sum(valid_samples) < n_samples) which(valid_samples) else NULL
          
          cat("     SUCCESS: HEALTHY VS DISEASE pattern detected!\n")
          cat("       Control samples:", sum(groups_factor == "Control"), "\n")
          cat("       Disease samples:", sum(groups_factor == "Disease"), "\n")
          
          return(list(
            success = TRUE,
            groups = groups_factor,
            sample_indices = sample_indices,
            pattern_type = "healthy_vs_disease",
            column_used = col_name,
            n_control = sum(groups_factor == "Control"),
            n_disease = sum(groups_factor == "Disease")
          ))
        }
      }
    }
  }
  
  cat("❌ No valid group pattern found\n")
  return(list(success = FALSE, reason = "no_valid_groups"))
}

cat("✅ STEP 02: Preprocessing loaded successfully\n\n")