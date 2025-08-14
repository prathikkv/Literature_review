#!/usr/bin/env Rscript
#' Step 01: Data Loader
#' 
#' Loads and validates datasets from cache based on pipeline configuration
#' Replaces hardcoded dataset loading with dynamic, configuration-driven approach

source("scripts/utilities/step_interface.R")

#' Execute Data Loading Step
#'
#' @param step_name Name of this step (should be "step_01_data_loader")
#' @param input_data Not used (first step in pipeline)
#' @param config Full pipeline configuration
#' @param checkpoint_dir Directory for saving checkpoints
#' @return Step result with loaded datasets
step_01_data_loader <- function(step_name, input_data = NULL, config, checkpoint_dir = "output/checkpoints") {
  
  cat("📊 STEP 01: DATA LOADER\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Loading datasets from configuration...\n\n")
  
  # Validate that this is the first step (no input expected)
  if (!is.null(input_data)) {
    warning("Data loader is first step but received input data - ignoring")
  }
  
  # Get dataset configuration
  if (is.null(config$datasets$active_datasets)) {
    return(create_step_result(
      success = FALSE,
      error_message = "No active datasets configured",
      step_name = step_name
    ))
  }
  
  datasets_config <- config$datasets$active_datasets
  loaded_datasets <- list()
  failed_datasets <- character(0)
  dataset_summary <- data.frame()
  
  cat("📋 DATASET LOADING CONFIGURATION:\n")
  cat("Total datasets configured:", length(datasets_config), "\n")
  
  # Process each configured dataset
  for (dataset_id in names(datasets_config)) {
    dataset_config <- datasets_config[[dataset_id]]
    
    cat("\n═══════════════════════════════════════════════════════════════\n")
    cat("📊 LOADING:", dataset_id, "(Priority:", dataset_config$priority, ")\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    
    # Display dataset configuration
    cat("🔬 Dataset Configuration:\n")
    cat("Disease type:", dataset_config$disease_type, "\n")
    cat("Control type:", dataset_config$control_type, "\n") 
    cat("Biological context:", dataset_config$biological_context, "\n")
    cat("Expected samples:", dataset_config$expected_samples, "\n")
    cat("Platform:", dataset_config$platform %||% "Unknown", "\n")
    cat("Inclusion reason:", dataset_config$inclusion_reason, "\n\n")
    
    # Attempt to load dataset
    dataset_result <- load_single_dataset(dataset_id, dataset_config, config)
    
    if (dataset_result$success) {
      loaded_datasets[[dataset_id]] <- dataset_result$dataset
      
      # Add to summary
      dataset_summary <- rbind(dataset_summary, data.frame(
        Dataset = dataset_id,
        Priority = dataset_config$priority,
        Samples = dataset_result$n_samples,
        Total_Genes = dataset_result$n_genes,
        Platform = dataset_result$platform,
        Disease_Type = dataset_config$disease_type,
        Biological_Context = dataset_config$biological_context,
        Status = "SUCCESS",
        Cache_File = dataset_result$cache_file,
        stringsAsFactors = FALSE
      ))
      
      cat("✅ SUCCESS:", dataset_id, "loaded successfully\n")
      cat("📈 Samples:", dataset_result$n_samples, "| Genes:", dataset_result$n_genes, "\n")
      
    } else {
      failed_datasets <- c(failed_datasets, dataset_id)
      
      # Check if dataset is optional
      is_optional <- !is.null(dataset_config$status) && dataset_config$status == "optional"
      
      # Add failure to summary
      dataset_summary <- rbind(dataset_summary, data.frame(
        Dataset = dataset_id,
        Priority = dataset_config$priority,
        Samples = NA,
        Total_Genes = NA,
        Platform = dataset_config$platform %||% NA,
        Disease_Type = dataset_config$disease_type,
        Biological_Context = dataset_config$biological_context,
        Status = if (is_optional) paste("OPTIONAL FAILED:", dataset_result$reason) else paste("FAILED:", dataset_result$reason),
        Cache_File = dataset_result$attempted_file %||% NA,
        stringsAsFactors = FALSE
      ))
      
      if (is_optional) {
        cat("ℹ️  OPTIONAL dataset failed:", dataset_id, "-", dataset_result$reason, "\n")
      } else {
        cat("❌ FAILED:", dataset_id, "-", dataset_result$reason, "\n")
      }
    }
  }
  
  # Generate loading summary
  cat("\n═══════════════════════════════════════════════════════════════\n")
  cat("📊 DATA LOADING SUMMARY\n") 
  cat("═══════════════════════════════════════════════════════════════\n")
  
  successful_count <- length(loaded_datasets)
  failed_count <- length(failed_datasets)
  optional_failures <- sum(grepl("OPTIONAL FAILED", dataset_summary$Status))
  
  cat("✅ Successfully loaded datasets:", successful_count, "\n")
  cat("❌ Failed datasets:", failed_count, "\n")
  cat("ℹ️  Optional failures:", optional_failures, "\n")
  
  if (successful_count > 0) {
    cat("\n📋 LOADED DATASETS:\n")
    successful_datasets <- dataset_summary[dataset_summary$Status == "SUCCESS", ]
    for (i in 1:nrow(successful_datasets)) {
      ds <- successful_datasets[i, ]
      cat(sprintf("  📊 %-10s: %s priority, %3d samples, %5d genes\n",
                  ds$Dataset, ds$Priority, ds$Samples, ds$Total_Genes))
    }
  }
  
  if (failed_count > optional_failures) {
    cat("\n❌ CRITICAL FAILURES:\n")
    critical_failures <- dataset_summary[grepl("^FAILED:", dataset_summary$Status), ]
    for (i in 1:nrow(critical_failures)) {
      ds <- critical_failures[i, ]
      cat(sprintf("  ❌ %-10s: %s\n", ds$Dataset, ds$Status))
    }
  }
  
  # Validate minimum dataset requirements
  min_datasets <- config$validation$success_criteria$min_datasets_processed %||% 3
  
  if (successful_count < min_datasets) {
    return(create_step_result(
      success = FALSE,
      error_message = paste("Insufficient datasets loaded:", successful_count, "< minimum required:", min_datasets),
      step_name = step_name,
      metadata = list(
        datasets_loaded = successful_count,
        datasets_failed = failed_count,
        dataset_summary = dataset_summary
      )
    ))
  }
  
  # Create output data structure
  output_data <- list(
    loaded_datasets = loaded_datasets,
    dataset_summary = dataset_summary,
    loading_stats = list(
      total_configured = length(datasets_config),
      successfully_loaded = successful_count,
      failed_datasets = failed_count,
      optional_failures = optional_failures
    )
  )
  
  # Validate output
  validate_step_output(
    output_data = output_data,
    required_fields = c("loaded_datasets", "dataset_summary", "loading_stats"),
    step_name = step_name
  )
  
  cat("\n🎉 DATA LOADING COMPLETED SUCCESSFULLY\n")
  cat("Ready for preprocessing step...\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(create_step_result(
    success = TRUE,
    output_data = output_data,
    step_name = step_name,
    metadata = list(
      datasets_loaded = successful_count,
      datasets_configured = length(datasets_config),
      total_samples = sum(dataset_summary$Samples[!is.na(dataset_summary$Samples)]),
      platforms = unique(dataset_summary$Platform[!is.na(dataset_summary$Platform)])
    )
  ))
}

#' Load Single Dataset
#'
#' Loads a single dataset from cache with error handling
#' @param dataset_id Dataset identifier
#' @param dataset_config Configuration for this dataset
#' @param config Full pipeline configuration
#' @return List with success status and dataset info
load_single_dataset <- function(dataset_id, dataset_config, config) {
  
  # Determine cache file location
  cache_file <- dataset_config$cache_location
  
  # Try fallback location if primary not found
  if (!file.exists(cache_file) && !is.null(dataset_config$fallback_location)) {
    cache_file <- dataset_config$fallback_location
  }
  
  # Check if file exists
  if (!file.exists(cache_file)) {
    # Try auto-discovery if enabled
    if (config$datasets$discovery$auto_detect %||% FALSE) {
      discovered_file <- discover_dataset_file(dataset_id, config)
      if (!is.null(discovered_file)) {
        cache_file <- discovered_file
        cat("🔍 AUTO-DISCOVERED cache file:", cache_file, "\n")
      }
    }
  }
  
  if (!file.exists(cache_file)) {
    return(list(
      success = FALSE,
      reason = "cache_file_not_found",
      attempted_file = cache_file
    ))
  }
  
  cat("📁 Loading from:", cache_file, "\n")
  
  # Load dataset
  tryCatch({
    dataset <- readRDS(cache_file)
    
    # Validate dataset structure
    if (!validate_dataset_structure(dataset, dataset_id)) {
      return(list(
        success = FALSE,
        reason = "invalid_dataset_structure",
        attempted_file = cache_file
      ))
    }
    
    # Check if dataset processing was successful
    if (!dataset$success) {
      return(list(
        success = FALSE,
        reason = "dataset_processing_failed",
        attempted_file = cache_file
      ))
    }
    
    # Extract dataset information
    n_samples <- dataset$n_samples %||% 
                 (if (!is.null(dataset$expression_matrix)) ncol(dataset$expression_matrix) else 0)
    n_genes <- dataset$n_genes %||%
               (if (!is.null(dataset$expression_matrix)) nrow(dataset$expression_matrix) else 0)
    platform <- dataset$dataset_info$platform %||% "Unknown"
    
    # Validate sample count against expected
    expected_samples <- dataset_config$expected_samples
    if (!is.null(expected_samples) && abs(n_samples - expected_samples) > 5) {
      cat("⚠️  WARNING: Sample count mismatch. Expected:", expected_samples, "Actual:", n_samples, "\n")
    }
    
    return(list(
      success = TRUE,
      dataset = dataset,
      n_samples = n_samples,
      n_genes = n_genes,
      platform = platform,
      cache_file = cache_file
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      reason = paste("load_error:", e$message),
      attempted_file = cache_file
    ))
  })
}

#' Discover Dataset File
#'
#' Attempts to find dataset file using auto-discovery
#' @param dataset_id Dataset identifier
#' @param config Pipeline configuration
#' @return Path to discovered file or NULL
discover_dataset_file <- function(dataset_id, config) {
  
  # Get discovery configuration
  cache_dirs <- config$datasets$discovery$cache_directories %||% c("cache")
  file_pattern <- config$datasets$discovery$file_pattern %||% "*_processed.rds"
  
  # Search for files matching dataset ID
  for (cache_dir in cache_dirs) {
    if (dir.exists(cache_dir)) {
      # Look for files matching the dataset ID
      pattern <- paste0(dataset_id, "_processed.rds")
      files <- list.files(cache_dir, pattern = pattern, recursive = TRUE, full.names = TRUE)
      
      if (length(files) > 0) {
        # Return first match (could be enhanced with selection logic)
        return(files[1])
      }
    }
  }
  
  return(NULL)
}

#' Validate Dataset Structure
#'
#' Validates that loaded dataset has required structure
#' @param dataset Loaded dataset object
#' @param dataset_id Dataset identifier for error reporting
#' @return TRUE if valid, FALSE otherwise
validate_dataset_structure <- function(dataset, dataset_id) {
  
  # Check basic structure
  required_fields <- c("success", "expression_matrix", "phenotype_data")
  
  for (field in required_fields) {
    if (is.null(dataset[[field]])) {
      cat("❌ Dataset", dataset_id, "missing required field:", field, "\n")
      return(FALSE)
    }
  }
  
  # Check expression matrix
  if (is.null(dataset$expression_matrix) || !is.matrix(dataset$expression_matrix)) {
    cat("❌ Dataset", dataset_id, "has invalid expression matrix\n")
    return(FALSE)
  }
  
  # Check phenotype data
  if (is.null(dataset$phenotype_data) || !is.data.frame(dataset$phenotype_data)) {
    cat("❌ Dataset", dataset_id, "has invalid phenotype data\n")
    return(FALSE)
  }
  
  # Check dimensions match
  if (ncol(dataset$expression_matrix) != nrow(dataset$phenotype_data)) {
    cat("❌ Dataset", dataset_id, "dimension mismatch: expression vs phenotype\n")
    return(FALSE)
  }
  
  return(TRUE)
}

# Test function for standalone execution
test_data_loader <- function() {
  cat("🧪 TESTING DATA LOADER STEP\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  # Load configuration
  source("scripts/utilities/config_validator.R")
  config_result <- load_and_validate_config("config_dynamic_pipeline.yml", strict_mode = FALSE)
  
  if (!config_result$validation$success) {
    cat("❌ Configuration validation failed\n")
    return(FALSE)
  }
  
  # Execute step
  result <- step_01_data_loader(
    step_name = "step_01_data_loader",
    input_data = NULL,
    config = config_result$config
  )
  
  # Print result
  print(result)
  
  return(result$success)
}

# Auto-test if run as script
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) == 0) {
  test_data_loader()
}

cat("✅ STEP 01: Data Loader loaded successfully\n\n")