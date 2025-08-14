#!/usr/bin/env Rscript
#' Multi-Modal Dataset Download Module
#' 
#' Specialized functions for downloading multi-modal datasets

#' Download Multi-Modal Datasets
#' 
#' @param target_datasets Character vector of dataset IDs
#' @param cache_dirs Named list of cache directories by modality
#' @param config Configuration from config.yml
#' @param max_retries Maximum retry attempts
#' @param timeout_seconds Timeout for downloads
#' @param force_fresh Force fresh downloads
#' @param verify_integrity Verify data integrity
download_multimodal_datasets <- function(target_datasets, 
                                        cache_dirs,
                                        config,
                                        max_retries = 3,
                                        timeout_seconds = 1200,
                                        force_fresh = FALSE,
                                        verify_integrity = TRUE) {
  
  cat("PROCESS: Multi-Modal Dataset Download System\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  # Create all cache directories
  for (dir_path in cache_dirs) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE)
      cat("CREATED: Cache directory", dir_path, "\n")
    }
  }
  
  download_results <- list()
  
  # Get dataset configurations from config
  all_datasets <- list()
  
  # Add microarray datasets
  if (!is.null(config$datasets$microarray$datasets)) {
    for (dataset_id in names(config$datasets$microarray$datasets)) {
      all_datasets[[dataset_id]] <- config$datasets$microarray$datasets[[dataset_id]]
      all_datasets[[dataset_id]]$modality <- "microarray"
      all_datasets[[dataset_id]]$platform <- config$datasets$microarray$platform
    }
  }
  
  # Add RNA-seq datasets
  if (!is.null(config$datasets$rna_seq$datasets)) {
    for (dataset_id in names(config$datasets$rna_seq$datasets)) {
      all_datasets[[dataset_id]] <- config$datasets$rna_seq$datasets[[dataset_id]]
      all_datasets[[dataset_id]]$modality <- "rna_seq"
    }
  }
  
  # Add single-cell datasets
  if (!is.null(config$datasets$single_cell$datasets)) {
    for (dataset_id in names(config$datasets$single_cell$datasets)) {
      all_datasets[[dataset_id]] <- config$datasets$single_cell$datasets[[dataset_id]]
      all_datasets[[dataset_id]]$modality <- "single_cell"
    }
  }
  
  total_datasets <- length(target_datasets)
  cat("DATA: Target Datasets:", total_datasets, "\n")
  cat("MODALITIES: Available datasets:", length(all_datasets), "\n\n")
  
  # Download each dataset
  for (i in seq_along(target_datasets)) {
    dataset_id <- target_datasets[i]
    cat("\\nDATA: Processing", i, "of", total_datasets, ":", dataset_id, "\n")
    cat(paste(rep("=", 40), collapse = ""), "\n")
    
    # Check if dataset is configured
    if (!dataset_id %in% names(all_datasets)) {
      cat("WARNING: Dataset", dataset_id, "not found in configuration\n")
      download_results[[dataset_id]] <- list(success = FALSE, error = "Dataset not configured")
      next
    }
    
    dataset_info <- all_datasets[[dataset_id]]
    modality <- dataset_info$modality
    
    # Determine cache directory based on modality
    if (modality == "microarray") {
      cache_dir <- cache_dirs$microarray
    } else if (modality == "rna_seq") {
      cache_dir <- cache_dirs$rna_seq
    } else if (modality == "single_cell") {
      cache_dir <- cache_dirs$single_cell
    } else {
      cache_dir <- cache_dirs$comprehensive
    }
    
    cache_file <- file.path(cache_dir, paste0(dataset_id, "_processed.rds"))
    
    # Check if cached and not forcing fresh
    if (file.exists(cache_file) && !force_fresh) {
      cat("CACHED: Loading existing data for", dataset_id, "\n")
      tryCatch({
        cached_data <- readRDS(cache_file)
        download_results[[dataset_id]] <- cached_data
        cat("SUCCESS: Loaded cached data (", cached_data$n_samples, "samples )\n")
        next
      }, error = function(e) {
        cat("WARNING: Failed to load cache, will re-download\n")
      })
    }
    
    # Download the dataset
    cat("DOWNLOAD: Attempting download for", dataset_id, "(", modality, ")\n")
    
    download_result <- download_single_dataset(
      dataset_id = dataset_id,
      dataset_info = dataset_info,
      cache_file = cache_file,
      max_retries = max_retries,
      timeout_seconds = timeout_seconds
    )
    
    download_results[[dataset_id]] <- download_result
    
    if (download_result$success) {
      cat("SUCCESS: Downloaded", dataset_id, "(", download_result$n_samples, "samples )\n")
    } else {
      cat("FAILED: Download failed for", dataset_id, "\n")
      cat("ERROR:", download_result$error, "\n")
    }
  }
  
  # Summary
  successful <- sum(sapply(download_results, function(x) x$success == TRUE))
  cat("\\nSUMMARY: Download Results\n")
  cat("=========================\n")
  cat("Total datasets attempted:", total_datasets, "\n")
  cat("Successful downloads:", successful, "\n")
  cat("Failed downloads:", total_datasets - successful, "\n")
  
  return(download_results)
}

#' Download Single Dataset
#' 
#' Download a single dataset with retry logic
download_single_dataset <- function(dataset_id, dataset_info, cache_file, max_retries = 3, timeout_seconds = 1200) {
  
  # Load required libraries
  suppressPackageStartupMessages({
    library(GEOquery)
    library(Biobase)
  })
  
  for (attempt in 1:max_retries) {
    cat("  Attempt", attempt, "of", max_retries, "\n")
    
    result <- tryCatch({
      # Set timeout
      old_timeout <- getOption("timeout")
      options(timeout = timeout_seconds)
      
      # Download using GEOquery
      cat("  Downloading from GEO...\n")
      gse <- getGEO(dataset_id, GSEMatrix = TRUE, getGPL = FALSE)
      
      if (length(gse) == 0) {
        stop("No data returned from GEO")
      }
      
      # Get the first (and usually only) series matrix
      eset <- gse[[1]]
      
      # Extract expression data
      expression_matrix <- exprs(eset)
      phenotype_data <- pData(eset)
      feature_data <- fData(eset)
      
      # Basic validation
      if (nrow(expression_matrix) == 0) {
        stop("Empty expression matrix")
      }
      
      if (ncol(expression_matrix) == 0) {
        stop("No samples in expression matrix")
      }
      
      # Create result object
      processed_data <- list(
        success = TRUE,
        dataset_id = dataset_id,
        expression_matrix = expression_matrix,
        phenotype_data = as.data.frame(phenotype_data),
        feature_data = as.data.frame(feature_data),
        dataset_info = dataset_info,
        n_genes = nrow(expression_matrix),
        n_samples = ncol(expression_matrix),
        download_time = Sys.time()
      )
      
      # Save to cache
      saveRDS(processed_data, cache_file)
      cat("  Saved to cache:", basename(cache_file), "\n")
      
      # Restore timeout
      options(timeout = old_timeout)
      
      return(processed_data)
      
    }, error = function(e) {
      options(timeout = old_timeout)
      cat("  ERROR on attempt", attempt, ":", e$message, "\n")
      if (attempt < max_retries) {
        cat("  Retrying in 5 seconds...\n")
        Sys.sleep(5)
      }
      return(list(success = FALSE, error = e$message))
    })
    
    if (result$success) {
      return(result)
    }
  }
  
  # All attempts failed
  return(list(success = FALSE, error = paste("Failed after", max_retries, "attempts")))
}

cat("SUCCESS: Multi-Modal Download Module loaded successfully\n")
cat("FUNCTIONS: download_multimodal_datasets(), download_single_dataset()\n")