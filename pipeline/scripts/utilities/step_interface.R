#!/usr/bin/env Rscript
#' Standard Step Interface Framework
#' 
#' Provides standardized contracts and utilities for all pipeline steps
#' Enables modular, testable, and maintainable pipeline components

library(tidyverse)
library(yaml)

#' Standard Step Interface
#' 
#' All pipeline steps must implement this interface for consistency
#' @param step_name Unique identifier for this step
#' @param input_data Input data from previous step (or NULL for first step)
#' @param config Full pipeline configuration
#' @param checkpoint_dir Directory for saving checkpoints
#' @param step_params Step-specific parameters (optional)
#' @return Standardized step result object
execute_step <- function(step_name, input_data = NULL, config, checkpoint_dir = "output/checkpoints", step_params = list()) {
  
  # This is the interface definition - actual steps will override this
  stop("execute_step must be implemented by each pipeline step")
}

#' Create Standard Step Result
#' 
#' Factory function for creating consistent step results
#' @param success TRUE if step completed successfully
#' @param output_data Step output data
#' @param metadata Step execution metadata
#' @param next_steps Optional list of conditional next steps
#' @param checkpoint_path Path to saved checkpoint (if applicable)
#' @param step_name Name of the step
#' @param execution_time Time taken to execute step
#' @param error_message Error message if step failed
#' @return Standardized step result object
create_step_result <- function(success = TRUE, 
                              output_data = NULL, 
                              metadata = list(),
                              next_steps = NULL,
                              checkpoint_path = NULL,
                              step_name = NULL,
                              execution_time = NULL,
                              error_message = NULL,
                              warnings = character(0)) {
  
  # Add standard metadata
  result_metadata <- list(
    step_name = step_name,
    timestamp = Sys.time(),
    execution_time = execution_time,
    success = success,
    data_summary = if (!is.null(output_data)) summarize_data(output_data) else NULL
  )
  
  # Merge with provided metadata
  result_metadata <- c(result_metadata, metadata)
  
  # Create result object
  result <- list(
    success = success,
    output_data = output_data,
    metadata = result_metadata,
    next_steps = next_steps,
    checkpoint_path = checkpoint_path,
    error_message = error_message,
    warnings = warnings,
    
    # Interface version for compatibility checking
    interface_version = "2.0.0"
  )
  
  # Add result class for method dispatch
  class(result) <- c("pipeline_step_result", "list")
  
  return(result)
}

#' Validate Step Input
#' 
#' Validates that input data meets step requirements
#' @param input_data Input data to validate
#' @param required_fields List of required fields in input data
#' @param step_name Name of step performing validation
#' @return TRUE if valid, stops execution if invalid
validate_step_input <- function(input_data, required_fields = NULL, step_name = "Unknown") {
  
  cat("🔍 VALIDATING INPUT for", step_name, "\n")
  
  # Check if input is needed
  if (is.null(input_data) && !is.null(required_fields)) {
    stop(paste("Step", step_name, "requires input data but none provided"))
  }
  
  # Check required fields
  if (!is.null(input_data) && !is.null(required_fields)) {
    if (is.list(input_data)) {
      missing_fields <- setdiff(required_fields, names(input_data))
      if (length(missing_fields) > 0) {
        stop(paste("Step", step_name, "missing required input fields:", paste(missing_fields, collapse = ", ")))
      }
    }
  }
  
  cat("✅ Input validation passed\n")
  return(TRUE)
}

#' Validate Step Output
#' 
#' Validates that output data meets step requirements
#' @param output_data Output data to validate
#' @param required_fields List of required fields in output data
#' @param step_name Name of step performing validation
#' @return TRUE if valid, stops execution if invalid
validate_step_output <- function(output_data, required_fields = NULL, step_name = "Unknown") {
  
  cat("🔍 VALIDATING OUTPUT for", step_name, "\n")
  
  # Check if output is provided
  if (is.null(output_data)) {
    stop(paste("Step", step_name, "produced no output data"))
  }
  
  # Check required fields
  if (!is.null(required_fields)) {
    if (is.list(output_data)) {
      missing_fields <- setdiff(required_fields, names(output_data))
      if (length(missing_fields) > 0) {
        stop(paste("Step", step_name, "missing required output fields:", paste(missing_fields, collapse = ", ")))
      }
    }
  }
  
  cat("✅ Output validation passed\n")
  return(TRUE)
}

#' Save Step Checkpoint
#' 
#' Saves step results for resumable execution
#' @param step_result Step result object
#' @param checkpoint_dir Directory for checkpoints
#' @param compress Use compression for checkpoint files
#' @return Path to saved checkpoint
save_step_checkpoint <- function(step_result, checkpoint_dir = "output/checkpoints", compress = TRUE) {
  
  if (!dir.exists(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE)
  }
  
  # Create checkpoint filename with timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  step_name <- step_result$metadata$step_name %||% "unknown_step"
  
  checkpoint_filename <- paste0(step_name, "_", timestamp, ".rds")
  checkpoint_path <- file.path(checkpoint_dir, checkpoint_filename)
  
  # Save checkpoint
  if (compress) {
    saveRDS(step_result, checkpoint_path, compress = TRUE)
  } else {
    saveRDS(step_result, checkpoint_path, compress = FALSE)
  }
  
  cat("💾 CHECKPOINT saved:", checkpoint_path, "\n")
  return(checkpoint_path)
}

#' Load Step Checkpoint
#' 
#' Loads previously saved step checkpoint
#' @param checkpoint_path Path to checkpoint file
#' @return Loaded step result object
load_step_checkpoint <- function(checkpoint_path) {
  
  if (!file.exists(checkpoint_path)) {
    stop(paste("Checkpoint file not found:", checkpoint_path))
  }
  
  step_result <- readRDS(checkpoint_path)
  
  # Validate checkpoint structure
  if (!inherits(step_result, "pipeline_step_result")) {
    stop("Invalid checkpoint file format")
  }
  
  cat("💾 CHECKPOINT loaded:", checkpoint_path, "\n")
  return(step_result)
}

#' Find Latest Checkpoint
#' 
#' Finds the most recent checkpoint for a given step
#' @param step_name Name of step to find checkpoint for
#' @param checkpoint_dir Directory containing checkpoints
#' @return Path to latest checkpoint or NULL if none found
find_latest_checkpoint <- function(step_name, checkpoint_dir = "output/checkpoints") {
  
  if (!dir.exists(checkpoint_dir)) {
    return(NULL)
  }
  
  # Find all checkpoints for this step
  pattern <- paste0("^", step_name, "_.*\\.rds$")
  checkpoint_files <- list.files(checkpoint_dir, pattern = pattern, full.names = TRUE)
  
  if (length(checkpoint_files) == 0) {
    return(NULL)
  }
  
  # Get file modification times
  file_times <- file.mtime(checkpoint_files)
  
  # Return most recent
  latest_file <- checkpoint_files[which.max(file_times)]
  
  cat("🔍 FOUND latest checkpoint for", step_name, ":", basename(latest_file), "\n")
  return(latest_file)
}

#' Summarize Data Object
#' 
#' Creates a summary of data object for metadata
#' @param data Data object to summarize
#' @return Summary list
summarize_data <- function(data) {
  
  if (is.null(data)) {
    return(list(type = "NULL", size = 0))
  }
  
  summary <- list(
    type = class(data)[1],
    size = object.size(data)
  )
  
  if (is.data.frame(data)) {
    summary$rows <- nrow(data)
    summary$cols <- ncol(data)
    summary$column_names <- names(data)
  } else if (is.list(data)) {
    summary$length <- length(data)
    summary$names <- names(data)
  } else if (is.vector(data)) {
    summary$length <- length(data)
  }
  
  return(summary)
}

#' Execute Step with Error Handling
#' 
#' Wrapper that provides comprehensive error handling for step execution
#' @param step_function Function implementing the step
#' @param step_name Name of the step
#' @param input_data Input data for step
#' @param config Pipeline configuration
#' @param checkpoint_dir Checkpoint directory
#' @param save_checkpoint Whether to save checkpoint on success
#' @param max_retries Maximum retry attempts
#' @return Step result object
execute_step_safely <- function(step_function, 
                                step_name, 
                                input_data = NULL, 
                                config, 
                                checkpoint_dir = "output/checkpoints",
                                save_checkpoint = TRUE,
                                max_retries = 3) {
  
  start_time <- Sys.time()
  attempt <- 1
  last_error <- NULL
  
  while (attempt <= max_retries) {
    cat("🚀 EXECUTING", step_name, "- Attempt", attempt, "\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    
    tryCatch({
      # Execute the step
      result <- step_function(
        step_name = step_name,
        input_data = input_data,
        config = config,
        checkpoint_dir = checkpoint_dir
      )
      
      # Validate result structure
      if (!inherits(result, "pipeline_step_result")) {
        stop("Step must return a pipeline_step_result object")
      }
      
      # Calculate execution time
      end_time <- Sys.time()
      execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
      result$metadata$execution_time <- execution_time
      
      # Save checkpoint if successful and requested
      if (result$success && save_checkpoint) {
        checkpoint_path <- save_step_checkpoint(result, checkpoint_dir)
        result$checkpoint_path <- checkpoint_path
      }
      
      cat("✅ STEP COMPLETED:", step_name, "\n")
      cat("⏱️  Execution time:", round(execution_time, 2), "seconds\n")
      cat("═══════════════════════════════════════════════════════════════\n\n")
      
      return(result)
      
    }, error = function(e) {
      last_error <- e
      cat("❌ STEP FAILED:", step_name, "- Attempt", attempt, "\n")
      cat("Error:", e$message, "\n")
      
      if (attempt < max_retries) {
        cat("🔄 Retrying in 5 seconds...\n")
        Sys.sleep(5)
        attempt <<- attempt + 1
      } else {
        cat("💥 All retry attempts exhausted\n")
        cat("═══════════════════════════════════════════════════════════════\n\n")
      }
    })
  }
  
  # If we get here, all attempts failed
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  result <- create_step_result(
    success = FALSE,
    output_data = NULL,
    step_name = step_name,
    execution_time = execution_time,
    error_message = last_error$message,
    metadata = list(
      attempts = max_retries,
      final_error = last_error$message
    )
  )
  
  return(result)
}

#' Print Step Result
#' 
#' Pretty print method for step results
#' @param result Step result object
print.pipeline_step_result <- function(result) {
  cat("🔧 PIPELINE STEP RESULT\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Step:", result$metadata$step_name %||% "Unknown", "\n")
  cat("Success:", if (result$success) "✅ YES" else "❌ NO", "\n")
  cat("Timestamp:", as.character(result$metadata$timestamp), "\n")
  
  if (!is.null(result$metadata$execution_time)) {
    cat("Execution time:", round(result$metadata$execution_time, 2), "seconds\n")
  }
  
  if (!is.null(result$error_message)) {
    cat("Error:", result$error_message, "\n")
  }
  
  if (length(result$warnings) > 0) {
    cat("Warnings:", length(result$warnings), "\n")
  }
  
  if (!is.null(result$metadata$data_summary)) {
    cat("Output data:", result$metadata$data_summary$type, "\n")
  }
  
  if (!is.null(result$checkpoint_path)) {
    cat("Checkpoint:", basename(result$checkpoint_path), "\n")
  }
  
  cat("═══════════════════════════════════════════════════════════════\n")
}

# Utility operator for NULL coalescing
`%||%` <- function(x, y) if (is.null(x)) y else x

cat("✅ STEP INTERFACE FRAMEWORK loaded successfully\n")
cat("Interface version: 2.0.0\n")
cat("Available functions: execute_step_safely, create_step_result, validate_step_input/output\n")
cat("Checkpoint functions: save/load/find_latest_checkpoint\n\n")