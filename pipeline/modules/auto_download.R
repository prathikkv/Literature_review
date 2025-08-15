#!/usr/bin/env Rscript
#' Auto-Download GEO Datasets Module
#' 
#' Automatically downloads missing GEO datasets from config.yml
#' Non-disruptive enhancement to existing pipeline
#' 
#' @author Claude Code Enhancement Module
#' @version 1.0.0

# Load required libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(yaml)
  library(tidyverse)
})

#' Auto-Download GEO Datasets
#'
#' Reads dataset list from config.yml and downloads missing datasets
#' @param config_file Path to configuration file
#' @param cache_dir Cache directory for storing datasets
#' @param force_download Force re-download even if cached
#' @return List with download status and results
auto_download_geo_datasets <- function(config_file = "config.yml",
                                      cache_dir = "cache",
                                      force_download = FALSE) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║           AUTO-DOWNLOAD GEO DATASETS MODULE                  ║\n")
  cat("║           Non-Disruptive Pipeline Enhancement                ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Initialize results tracking
  download_results <- list(
    success = TRUE,
    datasets_checked = 0,
    datasets_downloaded = 0,
    datasets_cached = 0,
    datasets_failed = 0,
    details = list(),
    log = character()
  )
  
  # Load configuration
  tryCatch({
    config <- yaml::read_yaml(config_file)
    cat("✅ Configuration loaded from:", config_file, "\n")
  }, error = function(e) {
    cat("❌ Failed to load config:", e$message, "\n")
    download_results$success <- FALSE
    return(download_results)
  })
  
  # Extract dataset list
  datasets <- names(config$datasets$active_datasets)
  if (length(datasets) == 0) {
    cat("⚠️  No datasets found in configuration\n")
    return(download_results)
  }
  
  cat("📊 Found", length(datasets), "datasets in configuration\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  # Create cache directories if needed
  cache_paths <- c(
    file.path(cache_dir, "microarray"),
    file.path(cache_dir, "comprehensive"),
    file.path(cache_dir, "downloads")
  )
  
  for (path in cache_paths) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE)
      cat("📁 Created cache directory:", path, "\n")
    }
  }
  
  # Process each dataset
  for (gse_id in datasets) {
    download_results$datasets_checked <- download_results$datasets_checked + 1
    
    cat("\n🔍 Processing:", gse_id, "\n")
    
    # Get dataset configuration
    dataset_config <- config$datasets$active_datasets[[gse_id]]
    
    # Check cache locations
    cache_locations <- c(
      dataset_config$cache_location,
      dataset_config$fallback_location,
      file.path(cache_dir, "microarray", paste0(gse_id, "_processed.rds")),
      file.path(cache_dir, "comprehensive", paste0(gse_id, "_processed.rds"))
    )
    
    # Remove NA values
    cache_locations <- cache_locations[!is.na(cache_locations)]
    
    # Check if dataset exists in cache
    dataset_exists <- FALSE
    existing_file <- NULL
    
    for (cache_file in cache_locations) {
      if (file.exists(cache_file) && !force_download) {
        dataset_exists <- TRUE
        existing_file <- cache_file
        break
      }
    }
    
    if (dataset_exists) {
      cat("✅ Dataset already cached:", existing_file, "\n")
      download_results$datasets_cached <- download_results$datasets_cached + 1
      download_results$details[[gse_id]] <- list(
        status = "cached",
        location = existing_file
      )
      next
    }
    
    # Download dataset
    cat("📥 Downloading dataset from GEO...\n")
    
    download_status <- download_geo_dataset(
      gse_id = gse_id,
      platform = dataset_config$platform,
      cache_dir = cache_dir,
      expected_samples = dataset_config$expected_samples
    )
    
    if (download_status$success) {
      cat("✅ Successfully downloaded:", gse_id, "\n")
      download_results$datasets_downloaded <- download_results$datasets_downloaded + 1
      download_results$details[[gse_id]] <- download_status
    } else {
      cat("❌ Failed to download:", gse_id, "-", download_status$error, "\n")
      download_results$datasets_failed <- download_results$datasets_failed + 1
      download_results$details[[gse_id]] <- download_status
    }
  }
  
  # Generate summary
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📊 DOWNLOAD SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Datasets checked:", download_results$datasets_checked, "\n")
  cat("Already cached:", download_results$datasets_cached, "\n")
  cat("Downloaded:", download_results$datasets_downloaded, "\n")
  cat("Failed:", download_results$datasets_failed, "\n")
  
  # Save download log
  log_file <- file.path(cache_dir, "download_log.txt")
  write_download_log(download_results, log_file)
  cat("\n💾 Download log saved to:", log_file, "\n")
  
  return(download_results)
}

#' Download Individual GEO Dataset
#'
#' Downloads and processes a single GEO dataset
#' @param gse_id GEO Series ID (e.g., "GSE57338")
#' @param platform GPL platform ID
#' @param cache_dir Cache directory
#' @param expected_samples Expected number of samples for validation
#' @return List with download status and file location
download_geo_dataset <- function(gse_id, platform, cache_dir, expected_samples = NULL) {
  
  result <- list(
    success = FALSE,
    gse_id = gse_id,
    location = NULL,
    samples = 0,
    platform_detected = NULL,
    error = NULL
  )
  
  # Set download directory
  download_dir <- file.path(cache_dir, "downloads")
  
  tryCatch({
    # Download from GEO
    cat("  ⏳ Connecting to GEO database...\n")
    
    # Try to get GSE object
    gse <- getGEO(gse_id, 
                  destdir = download_dir,
                  GSEMatrix = TRUE,
                  AnnotGPL = FALSE)
    
    # Handle list return (multiple platforms)
    if (is.list(gse) && length(gse) > 0) {
      # Try to find matching platform
      platform_match <- NULL
      for (i in seq_along(gse)) {
        gse_platform <- annotation(gse[[i]])
        if (gse_platform == platform || is.null(platform)) {
          platform_match <- i
          break
        }
      }
      
      if (!is.null(platform_match)) {
        gse <- gse[[platform_match]]
      } else {
        gse <- gse[[1]]  # Use first platform if no match
      }
    }
    
    # Extract expression matrix and phenotype data
    expr_matrix <- exprs(gse)
    pheno_data <- pData(gse)
    
    # Validate sample count
    actual_samples <- ncol(expr_matrix)
    cat("  📊 Samples found:", actual_samples, "\n")
    
    if (!is.null(expected_samples) && abs(actual_samples - expected_samples) > 5) {
      cat("  ⚠️  Warning: Sample count mismatch (expected:", expected_samples, ")\n")
    }
    
    # Detect platform
    platform_detected <- annotation(gse)
    cat("  🔬 Platform detected:", platform_detected, "\n")
    
    # Create processed dataset object
    processed_dataset <- list(
      success = TRUE,
      gse_id = gse_id,
      expression_matrix = expr_matrix,
      phenotype_data = pheno_data,
      dataset_info = list(
        platform = platform_detected,
        samples = actual_samples,
        genes = nrow(expr_matrix),
        download_date = Sys.Date()
      )
    )
    
    # Save to cache
    output_file <- file.path(cache_dir, "microarray", 
                            paste0(gse_id, "_processed.rds"))
    
    saveRDS(processed_dataset, output_file)
    cat("  💾 Saved to:", output_file, "\n")
    
    # Update result
    result$success <- TRUE
    result$location <- output_file
    result$samples <- actual_samples
    result$platform_detected <- platform_detected
    
  }, error = function(e) {
    result$error <- e$message
    cat("  ❌ Download error:", e$message, "\n")
  })
  
  return(result)
}

#' Write Download Log
#'
#' Saves download results to a log file
#' @param results Download results object
#' @param log_file Path to log file
write_download_log <- function(results, log_file) {
  
  log_content <- character()
  
  # Add header
  log_content <- c(log_content,
    "GEO DATASET AUTO-DOWNLOAD LOG",
    paste("Date:", Sys.Date()),
    paste("Time:", format(Sys.time(), "%H:%M:%S")),
    "",
    "SUMMARY:",
    paste("  Datasets checked:", results$datasets_checked),
    paste("  Already cached:", results$datasets_cached),
    paste("  Downloaded:", results$datasets_downloaded),
    paste("  Failed:", results$datasets_failed),
    "",
    "DETAILS:"
  )
  
  # Add dataset details
  for (gse_id in names(results$details)) {
    detail <- results$details[[gse_id]]
    log_content <- c(log_content,
      paste0("  ", gse_id, ":"),
      paste0("    Status: ", detail$status %||% "unknown"),
      paste0("    Location: ", detail$location %||% "N/A")
    )
    
    if (!is.null(detail$error)) {
      log_content <- c(log_content,
        paste0("    Error: ", detail$error)
      )
    }
  }
  
  # Write to file
  writeLines(log_content, log_file)
}

#' Validate Download Integrity
#'
#' Checks if downloaded dataset is valid and complete
#' @param file_path Path to downloaded RDS file
#' @return TRUE if valid, FALSE otherwise
validate_download_integrity <- function(file_path) {
  
  if (!file.exists(file_path)) {
    return(FALSE)
  }
  
  tryCatch({
    data <- readRDS(file_path)
    
    # Check required fields
    required_fields <- c("expression_matrix", "phenotype_data")
    for (field in required_fields) {
      if (is.null(data[[field]])) {
        return(FALSE)
      }
    }
    
    # Check data dimensions
    if (nrow(data$expression_matrix) == 0 || ncol(data$expression_matrix) == 0) {
      return(FALSE)
    }
    
    return(TRUE)
    
  }, error = function(e) {
    return(FALSE)
  })
}

# NULL coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Test function for module
test_auto_download <- function() {
  cat("🧪 Testing Auto-Download Module\n")
  
  # Test with small dataset
  test_result <- auto_download_geo_datasets(
    config_file = "config.yml",
    cache_dir = "cache",
    force_download = FALSE
  )
  
  if (test_result$success) {
    cat("✅ Auto-download module test passed\n")
  } else {
    cat("❌ Auto-download module test failed\n")
  }
  
  return(test_result)
}

cat("✅ Auto-Download Module loaded successfully\n")
cat("   Functions: auto_download_geo_datasets(), download_geo_dataset()\n")
cat("   Usage: source('modules/auto_download.R')\n")
cat("         auto_download_geo_datasets()\n")