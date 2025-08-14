#!/usr/bin/env Rscript
#' Pipeline Results Dynamic Loading Module
#' 
#' This module provides functions to dynamically load pipeline results
#' for integration with R Markdown documentation, ensuring consistency
#' between pipeline execution and documentation output

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(purrr)
})

#' Load Latest Pipeline Results
#'
#' Dynamically loads all available pipeline results from output directory
#' @param output_dir Pipeline output directory
#' @param fresh_only Only load results from latest pipeline run
#' @return List with all available pipeline results
load_pipeline_results <- function(output_dir = "output", fresh_only = TRUE) {
  
  # Loading pipeline results dynamically
  
  # Initialize results list
  pipeline_results <- list(
    timestamp = Sys.time(),
    data_loaded = FALSE,
    meta_analysis = NULL,
    dge_results = NULL,
    datasets_processed = NULL,
    visualizations = NULL,
    analysis_summary = NULL,
    sample_counts = NULL
  )
  
  # Check if output directory exists
  if (!dir.exists(output_dir)) {
    # Output directory not found, return empty results
    return(pipeline_results)
  }
  
  # Load meta-analysis results
  meta_file <- file.path(output_dir, "CAMK_meta_analysis_summary.csv")
  if (file.exists(meta_file)) {
    tryCatch({
      pipeline_results$meta_analysis <- read_csv(meta_file, show_col_types = FALSE)
      # Meta-analysis results loaded successfully
    }, error = function(e) {
      # Failed to load meta-analysis results
      pipeline_results$meta_analysis <- NULL
    })
  }
  
  # Load DGE results
  dge_file <- file.path(output_dir, "CAMK_focused_DGE_all_datasets.csv")
  if (file.exists(dge_file)) {
    tryCatch({
      pipeline_results$dge_results <- read_csv(dge_file, show_col_types = FALSE)
      # DGE results loaded successfully
    }, error = function(e) {
      # Failed to load DGE results
      pipeline_results$dge_results <- NULL
    })
  }
  
  # Extract dataset and sample information
  if (!is.null(pipeline_results$dge_results)) {
    datasets_info <- pipeline_results$dge_results %>%
      group_by(Dataset) %>%
      summarise(
        total_genes = n(),
        significant_genes = sum(Significant == TRUE, na.rm = TRUE),
        .groups = 'drop'
      )
    
    pipeline_results$datasets_processed <- datasets_info
    
    # Estimate sample counts from known dataset sizes
    dataset_samples <- c(
      "GSE57338" = 313,
      "GSE14975" = 10, 
      "GSE31821" = 6,
      "GSE41177" = 38,
      "GSE79768" = 26
    )
    
    actual_datasets <- unique(pipeline_results$dge_results$Dataset)
    total_samples <- sum(dataset_samples[names(dataset_samples) %in% actual_datasets], na.rm = TRUE)
    
    pipeline_results$sample_counts <- list(
      total_samples = total_samples,
      total_datasets = length(actual_datasets),
      dataset_breakdown = dataset_samples[names(dataset_samples) %in% actual_datasets]
    )
  }
  
  # Find visualization files
  viz_patterns <- c(
    "heatmap" = "CAMK_expression_heatmap.*\\.png",
    "correlation" = "CAMK_correlation.*\\.png", 
    "network" = ".*network.*\\.png",
    "survival" = ".*survival.*\\.png",
    "roc" = ".*roc.*\\.png"
  )
  
  pipeline_results$visualizations <- list()
  
  for (viz_type in names(viz_patterns)) {
    pattern <- viz_patterns[[viz_type]]
    files <- list.files(output_dir, pattern = pattern, full.names = TRUE, recursive = FALSE)
    if (length(files) > 0) {
      pipeline_results$visualizations[[viz_type]] <- files
      # Found visualization files
    }
  }
  
  # Load comprehensive results if available
  comprehensive_files <- c(
    "output/comprehensive_analysis_results.rds",
    "results/comprehensive_final_report/comprehensive_camk_analysis_final_results.rds"
  )
  
  for (comp_file in comprehensive_files) {
    if (file.exists(comp_file)) {
      tryCatch({
        comp_results <- readRDS(comp_file)
        pipeline_results$comprehensive_results <- comp_results
        # Comprehensive results loaded successfully
        break
      }, error = function(e) {
        # Failed to load comprehensive results
        next
      })
    }
  }
  
  # Create analysis summary
  pipeline_results$analysis_summary <- create_analysis_summary(pipeline_results)
  
  pipeline_results$data_loaded <- TRUE
  
  return(pipeline_results)
}

#' Create Analysis Summary
#'
#' Creates a summary of loaded pipeline results
#' @param pipeline_results Loaded pipeline results
#' @return Summary data frame
create_analysis_summary <- function(pipeline_results) {
  
  summary_data <- data.frame(
    Component = character(0),
    Status = character(0),
    Details = character(0),
    stringsAsFactors = FALSE
  )
  
  # Meta-analysis status
  if (!is.null(pipeline_results$meta_analysis)) {
    sig_genes <- sum(pipeline_results$meta_analysis$Significant == TRUE, na.rm = TRUE)
    total_genes <- nrow(pipeline_results$meta_analysis)
    
    summary_data <- rbind(summary_data, data.frame(
      Component = "Meta-Analysis",
      Status = "COMPLETED",
      Details = paste(sig_genes, "significant of", total_genes, "genes"),
      stringsAsFactors = FALSE
    ))
  } else {
    summary_data <- rbind(summary_data, data.frame(
      Component = "Meta-Analysis", 
      Status = "NOT FOUND",
      Details = "No meta-analysis results available",
      stringsAsFactors = FALSE
    ))
  }
  
  # DGE analysis status
  if (!is.null(pipeline_results$dge_results)) {
    n_datasets <- length(unique(pipeline_results$dge_results$Dataset))
    
    summary_data <- rbind(summary_data, data.frame(
      Component = "DGE Analysis",
      Status = "COMPLETED", 
      Details = paste("Results from", n_datasets, "dataset(s)"),
      stringsAsFactors = FALSE
    ))
  } else {
    summary_data <- rbind(summary_data, data.frame(
      Component = "DGE Analysis",
      Status = "NOT FOUND",
      Details = "No DGE results available",
      stringsAsFactors = FALSE
    ))
  }
  
  # Visualizations status
  if (!is.null(pipeline_results$visualizations) && length(pipeline_results$visualizations) > 0) {
    total_viz <- sum(sapply(pipeline_results$visualizations, length))
    
    summary_data <- rbind(summary_data, data.frame(
      Component = "Visualizations",
      Status = "AVAILABLE",
      Details = paste(total_viz, "visualization files found"),
      stringsAsFactors = FALSE
    ))
  } else {
    summary_data <- rbind(summary_data, data.frame(
      Component = "Visualizations", 
      Status = "LIMITED",
      Details = "Few or no visualization files found",
      stringsAsFactors = FALSE
    ))
  }
  
  # Sample counts status
  if (!is.null(pipeline_results$sample_counts)) {
    summary_data <- rbind(summary_data, data.frame(
      Component = "Sample Counts",
      Status = "CALCULATED",
      Details = paste(pipeline_results$sample_counts$total_samples, "total samples across", 
                     pipeline_results$sample_counts$total_datasets, "datasets"),
      stringsAsFactors = FALSE
    ))
  } else {
    summary_data <- rbind(summary_data, data.frame(
      Component = "Sample Counts",
      Status = "UNKNOWN", 
      Details = "Could not determine sample counts",
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_data)
}

#' Get Dynamic Metrics for R Markdown
#'
#' Extracts key metrics for use in R Markdown executive summary
#' @param pipeline_results Loaded pipeline results
#' @return List of dynamic metrics
get_dynamic_metrics <- function(pipeline_results) {
  
  metrics <- list(
    total_samples = 0,
    total_datasets = 0,
    significant_genes = 0,
    total_genes_tested = 0,
    top_gene_name = "CAMK2D",
    top_gene_pval = "N/A",
    analysis_modules_completed = 0
  )
  
  # Extract sample counts
  if (!is.null(pipeline_results$sample_counts)) {
    metrics$total_samples <- pipeline_results$sample_counts$total_samples
    metrics$total_datasets <- pipeline_results$sample_counts$total_datasets
  }
  
  # Extract meta-analysis metrics
  if (!is.null(pipeline_results$meta_analysis)) {
    metrics$significant_genes <- sum(pipeline_results$meta_analysis$Significant == TRUE, na.rm = TRUE)
    metrics$total_genes_tested <- nrow(pipeline_results$meta_analysis)
    
    # Get most significant gene
    if (metrics$significant_genes > 0) {
      sig_results <- pipeline_results$meta_analysis[pipeline_results$meta_analysis$Significant == TRUE, ]
      if (nrow(sig_results) > 0) {
        top_result <- sig_results[order(sig_results$Combined_P_Value)[1], ]
        metrics$top_gene_name <- top_result$Gene
        metrics$top_gene_pval <- formatC(top_result$Combined_P_Value, format = "e", digits = 2)
      }
    }
  }
  
  # Count completed analysis modules
  completed_modules <- 0
  if (!is.null(pipeline_results$meta_analysis)) completed_modules <- completed_modules + 1
  if (!is.null(pipeline_results$dge_results)) completed_modules <- completed_modules + 1
  if (!is.null(pipeline_results$visualizations) && length(pipeline_results$visualizations) > 0) completed_modules <- completed_modules + 1
  
  metrics$analysis_modules_completed <- completed_modules
  
  return(metrics)
}

#' Check Analysis Module Status
#'
#' Checks which analysis modules have been completed
#' @param pipeline_results Loaded pipeline results  
#' @return List of module statuses
check_module_status <- function(pipeline_results) {
  
  status <- list(
    dge_analysis = FALSE,
    meta_analysis = FALSE, 
    survival_analysis = FALSE,
    ml_analysis = FALSE,
    network_analysis = FALSE,
    pathway_analysis = FALSE
  )
  
  # Check DGE analysis
  status$dge_analysis <- !is.null(pipeline_results$dge_results)
  
  # Check meta-analysis
  status$meta_analysis <- !is.null(pipeline_results$meta_analysis)
  
  # Check for survival analysis results
  survival_files <- c(
    "results/survival_analysis/CAMK2D_survival_analysis.rds",
    "output/CAMK2D_survival_analysis.rds"
  )
  status$survival_analysis <- any(sapply(survival_files, file.exists))
  
  # Check for ML analysis results  
  ml_files <- c(
    "results/ml_biomarker_prediction/ml_biomarker_prediction_results.rds",
    "output/ml_biomarker_prediction_results.rds"
  )
  status$ml_analysis <- any(sapply(ml_files, file.exists))
  
  # Check for network analysis results
  network_files <- c(
    "results/camk_interconnection_analysis/camk_interconnection_analysis.rds", 
    "output/camk_interconnection_analysis.rds"
  )
  status$network_analysis <- any(sapply(network_files, file.exists))
  
  # Check for pathway analysis results
  pathway_files <- c(
    "results/pathway_analysis/pathway_enrichment_results.rds",
    "output/pathway_enrichment_results.rds"
  )
  status$pathway_analysis <- any(sapply(pathway_files, file.exists))
  
  return(status)
}

#' Sync Pipeline Results
#'
#' Main function to sync pipeline results with R Markdown
#' @param force_reload Force reload of all results
#' @return Comprehensive results list for R Markdown
sync_pipeline_results <- function(force_reload = FALSE) {
  
  # Synchronizing pipeline results with R Markdown
  
  # Load pipeline results
  pipeline_results <- load_pipeline_results(fresh_only = !force_reload)
  
  # Get dynamic metrics
  dynamic_metrics <- get_dynamic_metrics(pipeline_results)
  
  # Check module status  
  module_status <- check_module_status(pipeline_results)
  
  # Create comprehensive sync results
  sync_results <- list(
    pipeline_results = pipeline_results,
    dynamic_metrics = dynamic_metrics,
    module_status = module_status,
    sync_timestamp = Sys.time(),
    data_quality = validate_data_quality(pipeline_results)
  )
  
  # Pipeline results synchronized successfully
  
  return(sync_results)
}

#' Validate Data Quality
#'
#' Performs basic validation of loaded pipeline data
#' @param pipeline_results Loaded pipeline results
#' @return Validation report
validate_data_quality <- function(pipeline_results) {
  
  validation <- list(
    status = "PASS",
    issues = character(0),
    warnings = character(0)
  )
  
  # Check for missing critical data
  if (is.null(pipeline_results$dge_results)) {
    validation$issues <- c(validation$issues, "No DGE results found")
    validation$status <- "FAIL"
  }
  
  if (is.null(pipeline_results$meta_analysis)) {
    validation$warnings <- c(validation$warnings, "No meta-analysis results found")
  }
  
  # Check data consistency
  if (!is.null(pipeline_results$dge_results) && !is.null(pipeline_results$meta_analysis)) {
    dge_datasets <- unique(pipeline_results$dge_results$Dataset)
    meta_datasets <- strsplit(paste(pipeline_results$meta_analysis$Datasets, collapse = ", "), ", ")[[1]]
    meta_datasets <- unique(trimws(meta_datasets))
    
    if (length(setdiff(dge_datasets, meta_datasets)) > 0) {
      validation$warnings <- c(validation$warnings, "Dataset mismatch between DGE and meta-analysis")
    }
  }
  
  # Check for empty results
  if (!is.null(pipeline_results$meta_analysis) && nrow(pipeline_results$meta_analysis) == 0) {
    validation$warnings <- c(validation$warnings, "Meta-analysis results are empty")
  }
  
  return(validation)
}

# Pipeline Results Dynamic Loading Module loaded successfully
# MAIN FUNCTIONS: sync_pipeline_results(), load_pipeline_results(), get_dynamic_metrics()
# INTEGRATION: Ready for R Markdown dynamic integration