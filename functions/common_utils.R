#!/usr/bin/env Rscript
#' Common Utility Functions
#' 
#' Shared utilities to eliminate code redundancy across CAMK analysis scripts

#' Load Required Libraries
#' 
#' Centralized library loading to ensure consistency
#' @param analysis_type Character: "basic", "visualization", "meta_analysis", "network"
load_required_libraries <- function(analysis_type = "basic") {
  
  # Core libraries for all analyses
  basic_libs <- c("tidyverse", "limma", "ggplot2")
  
  # Additional libraries by analysis type
  additional_libs <- switch(analysis_type,
    "basic" = c(),
    "visualization" = c("pheatmap", "corrplot", "RColorBrewer"),
    "meta_analysis" = c("metafor", "meta"),
    "network" = c("WGCNA", "igraph", "cluster"),
    c()  # Default empty
  )
  
  all_libs <- unique(c(basic_libs, additional_libs))
  
  # Load libraries with error handling
  for (lib in all_libs) {
    tryCatch({
      suppressPackageStartupMessages(library(lib, character.only = TRUE))
    }, error = function(e) {
      warning(paste("Could not load library:", lib))
    })
  }
  
  invisible(all_libs)
}

#' Initialize Output Directory
#' 
#' Standardized output directory creation
#' @param output_dir Character, path to output directory
#' @param create_subdirs Character vector of subdirectories to create
ensure_output_directory <- function(output_dir = "output", create_subdirs = NULL) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("SAVED: Created output directory:", output_dir, "\n")
  }
  
  # Create subdirectories if specified
  if (!is.null(create_subdirs)) {
    for (subdir in create_subdirs) {
      full_path <- file.path(output_dir, subdir)
      if (!dir.exists(full_path)) {
        dir.create(full_path, recursive = TRUE)
        cat("SAVED: Created subdirectory:", full_path, "\n")
      }
    }
  }
  
  invisible(output_dir)
}

#' Save Analysis Results
#' 
#' Standardized result saving with consistent naming
#' @param results Object to save
#' @param filename Character, base filename (without extension)
#' @param output_dir Character, output directory
#' @param save_csv Logical, whether to save CSV version if applicable
#' @param add_timestamp Logical, whether to add timestamp to filename
save_analysis_results <- function(results, filename, output_dir = "output", 
                                 save_csv = FALSE, add_timestamp = FALSE) {
  
  ensure_output_directory(output_dir)
  
  # Add timestamp if requested
  if (add_timestamp) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    filename <- paste0(filename, "_", timestamp)
  }
  
  # Save RDS file
  rds_path <- file.path(output_dir, paste0(filename, ".rds"))
  saveRDS(results, rds_path)
  cat("SAVED: Saved RDS:", rds_path, "\n")
  
  # Save CSV if requested and data is data.frame-like
  if (save_csv && is.data.frame(results)) {
    csv_path <- file.path(output_dir, paste0(filename, ".csv"))
    write.csv(results, csv_path, row.names = FALSE)
    cat("SAVED: Saved CSV:", csv_path, "\n")
  } else if (save_csv && is.list(results)) {
    # Try to save first data.frame element as CSV
    df_elements <- results[sapply(results, is.data.frame)]
    if (length(df_elements) > 0) {
      csv_path <- file.path(output_dir, paste0(filename, "_main.csv"))
      write.csv(df_elements[[1]], csv_path, row.names = FALSE)
      cat("SAVED: Saved main results as CSV:", csv_path, "\n")
    }
  }
  
  invisible(rds_path)
}

#' Standardized Analysis Header
#' 
#' Consistent header formatting across analysis scripts
#' @param title Character, analysis title
#' @param subtitle Character, optional subtitle
#' @param methodology Character, methodology description
print_analysis_header <- function(title, subtitle = NULL, methodology = NULL) {
  
  # Main title
  cat("TARGET:", toupper(title), "\n")
  cat(paste(rep("=", nchar(title) + 3), collapse = ""), "\n")
  
  # Subtitle if provided
  if (!is.null(subtitle)) {
    cat("DATA:", subtitle, "\n")
  }
  
  # Methodology note if provided
  if (!is.null(methodology)) {
    cat("METHOD: Methodology:", methodology, "\n")
  }
  
  cat("\n")
}

#' Format P-values
#' 
#' Consistent p-value formatting
#' @param p_values Numeric vector of p-values
#' @param threshold Numeric, threshold for scientific notation
format_pvalues <- function(p_values, threshold = 0.001) {
  ifelse(p_values < threshold, 
         format(p_values, scientific = TRUE, digits = 2),
         format(round(p_values, 3), nsmall = 3))
}

#' Calculate Effect Size Categories
#' 
#' Standardized effect size interpretation
#' @param effect_sizes Numeric vector of effect sizes (absolute log fold changes)
#' @return Character vector of effect size categories
categorize_effect_sizes <- function(effect_sizes) {
  abs_effects <- abs(effect_sizes)
  
  ifelse(abs_effects >= 1.0, "Large", 
         ifelse(abs_effects >= 0.5, "Moderate",
                ifelse(abs_effects >= 0.2, "Small", "Minimal")))
}

#' Null-safe Operator
#' 
#' Utility operator for handling NULL values
#' @param x First value
#' @param y Default value if x is NULL
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Progress Reporter
#' 
#' Standardized progress reporting
#' @param current Current step number
#' @param total Total steps
#' @param message Character, progress message
report_progress <- function(current, total, message = NULL) {
  percentage <- round(current / total * 100)
  progress_bar <- paste0("[", 
                        paste(rep("=", floor(percentage / 5)), collapse = ""),
                        paste(rep(" ", 20 - floor(percentage / 5)), collapse = ""),
                        "]")
  
  if (!is.null(message)) {
    cat("RESULTS: Progress", progress_bar, paste0(percentage, "%"), "-", message, "\n")
  } else {
    cat("RESULTS: Progress", progress_bar, paste0(percentage, "%\n"))
  }
}

#' Validate Dataset Structure
#' 
#' Check if dataset has required components for analysis
#' @param dataset List, dataset object to validate
#' @return List with validation results
validate_dataset_structure <- function(dataset) {
  
  required_components <- c("expression_matrix", "phenotype_data")
  optional_components <- c("feature_data", "platform_info")
  
  results <- list(
    is_valid = TRUE,
    missing_required = character(),
    available_components = names(dataset),
    sample_count = 0,
    gene_count = 0
  )
  
  # Check required components
  for (component in required_components) {
    if (!component %in% names(dataset) || is.null(dataset[[component]])) {
      results$missing_required <- c(results$missing_required, component)
      results$is_valid <- FALSE
    }
  }
  
  # Get dataset dimensions if valid
  if (results$is_valid && !is.null(dataset$expression_matrix)) {
    results$gene_count <- nrow(dataset$expression_matrix)
    results$sample_count <- ncol(dataset$expression_matrix)
  }
  
  results
}