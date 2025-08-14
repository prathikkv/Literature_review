#!/usr/bin/env Rscript
#' Pipeline Orchestrator
#' 
#' Main controller for the dynamic CAMK2D analysis pipeline
#' Manages step execution, dependencies, checkpointing, and error handling

# Load required utilities
source("scripts/utilities/step_interface.R")
source("scripts/utilities/config_validator.R")

# Load all pipeline steps
source("scripts/pipeline/step_01_data_loader.R")
source("scripts/pipeline/step_02_preprocessing.R")
source("scripts/pipeline/step_03_dge_analysis.R")
source("scripts/pipeline/step_04_meta_analysis.R")
source("scripts/pipeline/step_05_report_generator.R")

#' Execute Dynamic Pipeline
#'
#' Main orchestrator function that manages the entire pipeline execution
#' @param config_file Path to configuration file
#' @param resume_from_checkpoint Resume from a specific checkpoint (step name)
#' @param force_rerun Force rerun of all steps (ignore checkpoints)
#' @param steps_to_run Specific steps to run (NULL for all steps)
#' @param checkpoint_dir Directory for checkpoints
#' @return Pipeline execution result
execute_dynamic_pipeline <- function(config_file = "config_dynamic_pipeline.yml",
                                    resume_from_checkpoint = NULL,
                                    force_rerun = FALSE,
                                    steps_to_run = NULL,
                                    checkpoint_dir = "output/checkpoints") {
  
  pipeline_start_time <- Sys.time()
  
  cat("🚀 CAMK2D DYNAMIC PIPELINE ORCHESTRATOR\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Pipeline execution started at:", as.character(pipeline_start_time), "\n")
  cat("Configuration file:", config_file, "\n")
  cat("Checkpoint directory:", checkpoint_dir, "\n")
  if (!is.null(resume_from_checkpoint)) {
    cat("Resume from checkpoint:", resume_from_checkpoint, "\n")
  }
  if (force_rerun) {
    cat("Force rerun: TRUE (ignoring checkpoints)\n")
  }
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Step 1: Load and validate configuration
  cat("🔧 PHASE 1: CONFIGURATION VALIDATION\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  config_result <- load_and_validate_config(config_file, strict_mode = FALSE)
  
  if (!config_result$validation$success) {
    stop("❌ Configuration validation failed. Please fix errors and try again.")
  }
  
  config <- config_result$config
  cat("✅ Configuration loaded and validated successfully\n\n")
  
  # Step 2: Initialize pipeline execution plan
  cat("📋 PHASE 2: PIPELINE EXECUTION PLANNING\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  execution_plan <- create_execution_plan(config, steps_to_run, resume_from_checkpoint, force_rerun, checkpoint_dir)
  
  cat("Pipeline steps configured:", length(config$pipeline$steps), "\n")
  cat("Steps to execute:", length(execution_plan$steps_to_execute), "\n")
  cat("Steps to execute:", paste(execution_plan$steps_to_execute, collapse = ", "), "\n")
  
  if (execution_plan$resume_mode) {
    cat("Resume mode: Resuming from", execution_plan$resume_from_step, "\n")
  }
  
  cat("✅ Execution plan created successfully\n\n")
  
  # Step 3: Execute pipeline steps
  cat("🔄 PHASE 3: PIPELINE EXECUTION\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  execution_result <- execute_pipeline_steps(execution_plan, config, checkpoint_dir)
  
  # Step 4: Generate pipeline summary
  pipeline_end_time <- Sys.time()
  total_execution_time <- as.numeric(difftime(pipeline_end_time, pipeline_start_time, units = "secs"))
  
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📊 PIPELINE EXECUTION SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  cat("🕐 Total execution time:", round(total_execution_time / 60, 2), "minutes\n")
  cat("📊 Steps planned:", length(execution_plan$steps_to_execute), "\n")
  cat("✅ Steps completed successfully:", execution_result$successful_steps, "\n")
  cat("❌ Steps failed:", execution_result$failed_steps, "\n")
  cat("⚠️  Total warnings:", length(execution_result$all_warnings), "\n")
  
  if (execution_result$pipeline_success) {
    cat("\n🎉 PIPELINE EXECUTION COMPLETED SUCCESSFULLY!\n")
    
    # Display key results if available
    if (!is.null(execution_result$final_result) && !is.null(execution_result$final_result$output_data$meta_summary)) {
      meta_summary <- execution_result$final_result$output_data$meta_summary
      cat("\n📊 KEY SCIENTIFIC RESULTS:\n")
      cat("  Genes analyzed:", meta_summary$genes_analyzed %||% "N/A", "\n")
      cat("  Significant genes:", meta_summary$significant_genes %||% "N/A", "\n")
      cat("  CAMK2D significant:", if (meta_summary$camk2d_significant %||% FALSE) "YES ✅" else "NO", "\n")
      
      if (!is.null(execution_result$final_result$output_data$report_path)) {
        cat("  Report generated:", execution_result$final_result$output_data$report_path, "\n")
      }
    }
    
    cat("\n🚀 PIPELINE IS READY FOR SCIENTIFIC USE\n")
    
  } else {
    cat("\n❌ PIPELINE EXECUTION FAILED\n")
    cat("Last successful step:", execution_result$last_successful_step %||% "None", "\n")
    cat("Failed step:", execution_result$failed_step %||% "Unknown", "\n")
    cat("Error message:", execution_result$error_message %||% "Unknown error", "\n")
  }
  
  if (length(execution_result$all_warnings) > 0) {
    cat("\n⚠️  PIPELINE WARNINGS:\n")
    for (warning in execution_result$all_warnings) {
      cat("  ⚠️ ", warning, "\n")
    }
  }
  
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(list(
    success = execution_result$pipeline_success,
    execution_time_minutes = round(total_execution_time / 60, 2),
    steps_completed = execution_result$successful_steps,
    steps_failed = execution_result$failed_steps,
    warnings = execution_result$all_warnings,
    final_result = execution_result$final_result,
    execution_log = execution_result$step_results
  ))
}

#' Create Execution Plan
#'
#' Analyzes configuration and creates step execution plan
#' @param config Pipeline configuration
#' @param steps_to_run Specific steps to run
#' @param resume_from_checkpoint Resume point
#' @param force_rerun Force rerun flag
#' @param checkpoint_dir Checkpoint directory
#' @return Execution plan
create_execution_plan <- function(config, steps_to_run, resume_from_checkpoint, force_rerun, checkpoint_dir) {
  
  # Get all configured steps
  all_steps <- config$pipeline$steps
  dependencies <- config$pipeline$dependencies
  
  # Determine which steps to execute
  if (!is.null(steps_to_run)) {
    # Specific steps requested
    steps_to_execute <- intersect(steps_to_run, all_steps)
    
    # Add dependencies if needed
    for (step in steps_to_execute) {
      step_deps <- dependencies[[step]]
      if (!is.null(step_deps) && length(step_deps) > 0) {
        steps_to_execute <- unique(c(step_deps, steps_to_execute))
      }
    }
    
  } else if (!is.null(resume_from_checkpoint)) {
    # Resume from a specific step
    resume_index <- which(all_steps == resume_from_checkpoint)
    if (length(resume_index) == 0) {
      stop(paste("Resume checkpoint not found:", resume_from_checkpoint))
    }
    
    steps_to_execute <- all_steps[resume_index:length(all_steps)]
    
  } else {
    # Execute all steps
    steps_to_execute <- all_steps
  }
  
  # Check for available checkpoints (unless force_rerun)
  checkpoint_info <- list()
  if (!force_rerun) {
    for (step in steps_to_execute) {
      latest_checkpoint <- find_latest_checkpoint(step, checkpoint_dir)
      if (!is.null(latest_checkpoint)) {
        checkpoint_info[[step]] <- latest_checkpoint
      }
    }
  }
  
  return(list(
    steps_to_execute = steps_to_execute,
    all_steps = all_steps,
    dependencies = dependencies,
    checkpoint_info = checkpoint_info,
    resume_mode = !is.null(resume_from_checkpoint),
    resume_from_step = resume_from_checkpoint,
    force_rerun = force_rerun
  ))
}

#' Execute Pipeline Steps
#'
#' Executes the planned pipeline steps with error handling
#' @param execution_plan Execution plan
#' @param config Pipeline configuration
#' @param checkpoint_dir Checkpoint directory
#' @return Execution result
execute_pipeline_steps <- function(execution_plan, config, checkpoint_dir) {
  
  steps_to_execute <- execution_plan$steps_to_execute
  dependencies <- execution_plan$dependencies
  checkpoint_info <- execution_plan$checkpoint_info
  
  # Initialize execution tracking
  step_results <- list()
  current_data <- NULL
  successful_steps <- 0
  failed_steps <- 0
  all_warnings <- character(0)
  
  # Step execution mapping
  step_functions <- list(
    "step_01_data_loader" = step_01_data_loader,
    "step_02_preprocessing" = step_02_preprocessing,
    "step_03_dge_analysis" = step_03_dge_analysis,
    "step_04_meta_analysis" = step_04_meta_analysis,
    "step_05_report_generator" = step_05_report_generator
  )
  
  # Execute steps in order
  for (step_name in steps_to_execute) {
    
    cat("\n🔄 EXECUTING STEP:", step_name, "\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    
    step_start_time <- Sys.time()
    
    # Check if we can resume from checkpoint
    if (!execution_plan$force_rerun && step_name %in% names(checkpoint_info)) {
      cat("💾 ATTEMPTING CHECKPOINT RESUME\n")
      
      tryCatch({
        checkpoint_result <- load_step_checkpoint(checkpoint_info[[step_name]])
        
        if (checkpoint_result$success) {
          cat("✅ RESUMED FROM CHECKPOINT:", basename(checkpoint_info[[step_name]]), "\n")
          
          # Validate checkpoint data for next step
          if (validate_checkpoint_for_next_step(checkpoint_result, step_name, dependencies)) {
            step_results[[step_name]] <- checkpoint_result
            current_data <- checkpoint_result$output_data
            successful_steps <- successful_steps + 1
            
            # Collect warnings from checkpoint
            if (length(checkpoint_result$warnings) > 0) {
              all_warnings <- c(all_warnings, checkpoint_result$warnings)
            }
            
            next  # Skip to next step
          } else {
            cat("⚠️  CHECKPOINT DATA INVALID - RERUNNING STEP\n")
          }
        }
      }, error = function(e) {
        cat("❌ CHECKPOINT LOAD FAILED:", e$message, "- RERUNNING STEP\n")
      })
    }
    
    # Get step function
    step_function <- step_functions[[step_name]]
    if (is.null(step_function)) {
      cat("❌ STEP FUNCTION NOT FOUND:", step_name, "\n")
      
      return(list(
        pipeline_success = FALSE,
        successful_steps = successful_steps,
        failed_steps = failed_steps + 1,
        last_successful_step = if (successful_steps > 0) names(step_results)[successful_steps] else NULL,
        failed_step = step_name,
        error_message = paste("Step function not found:", step_name),
        all_warnings = all_warnings,
        step_results = step_results
      ))
    }
    
    # Execute step with comprehensive error handling
    step_result <- execute_step_safely(
      step_function = step_function,
      step_name = step_name,
      input_data = current_data,
      config = config,
      checkpoint_dir = checkpoint_dir,
      save_checkpoint = config$execution$checkpoints$enabled %||% TRUE,
      max_retries = config$execution$error_handling$max_retries %||% 3
    )
    
    step_end_time <- Sys.time()
    step_execution_time <- as.numeric(difftime(step_end_time, step_start_time, units = "secs"))
    
    # Process step result
    if (step_result$success) {
      cat("✅ STEP COMPLETED SUCCESSFULLY:", step_name, "\n")
      cat("⏱️  Step execution time:", round(step_execution_time, 2), "seconds\n")
      
      step_results[[step_name]] <- step_result
      current_data <- step_result$output_data
      successful_steps <- successful_steps + 1
      
      # Collect warnings
      if (length(step_result$warnings) > 0) {
        all_warnings <- c(all_warnings, step_result$warnings)
        cat("⚠️  Step warnings:", length(step_result$warnings), "\n")
      }
      
    } else {
      cat("❌ STEP FAILED:", step_name, "\n")
      cat("❌ Error:", step_result$error_message, "\n")
      
      failed_steps <- failed_steps + 1
      
      # Collect any warnings from failed step
      if (length(step_result$warnings) > 0) {
        all_warnings <- c(all_warnings, step_result$warnings)
      }
      
      # Check if we should continue or fail fast
      fail_fast <- config$execution$error_handling$fail_fast %||% FALSE
      
      if (fail_fast) {
        cat("💥 FAIL FAST MODE - STOPPING PIPELINE\n")
        
        return(list(
          pipeline_success = FALSE,
          successful_steps = successful_steps,
          failed_steps = failed_steps,
          last_successful_step = if (successful_steps > 0) names(step_results)[successful_steps] else NULL,
          failed_step = step_name,
          error_message = step_result$error_message,
          all_warnings = all_warnings,
          step_results = step_results,
          final_result = NULL
        ))
      } else {
        cat("🔄 CONTINUE ON ERROR MODE - ATTEMPTING NEXT STEP\n")
        # Continue with NULL data (next step will need to handle gracefully)
        current_data <- NULL
      }
    }
    
    cat("═══════════════════════════════════════════════════════════════\n")
  }
  
  # Determine overall success
  pipeline_success <- (failed_steps == 0) && (successful_steps > 0)
  final_result <- if (length(step_results) > 0) step_results[[length(step_results)]] else NULL
  
  return(list(
    pipeline_success = pipeline_success,
    successful_steps = successful_steps,
    failed_steps = failed_steps,
    last_successful_step = if (successful_steps > 0) names(step_results)[successful_steps] else NULL,
    failed_step = if (failed_steps > 0) "Multiple or Unknown" else NULL,
    error_message = if (!pipeline_success) "Pipeline had step failures" else NULL,
    all_warnings = all_warnings,
    step_results = step_results,
    final_result = final_result
  ))
}

#' Validate Checkpoint for Next Step
#'
#' Validates that checkpoint data is suitable for the next pipeline step
#' @param checkpoint_result Loaded checkpoint
#' @param current_step Current step name
#' @param dependencies Step dependencies
#' @return TRUE if valid
validate_checkpoint_for_next_step <- function(checkpoint_result, current_step, dependencies) {
  
  # Basic validation
  if (is.null(checkpoint_result$output_data)) {
    return(FALSE)
  }
  
  # Check interface version compatibility
  if (is.null(checkpoint_result$interface_version) || checkpoint_result$interface_version != "2.0.0") {
    cat("⚠️  CHECKPOINT INTERFACE VERSION MISMATCH\n")
    return(FALSE)
  }
  
  # Step-specific validation
  if (current_step == "step_01_data_loader") {
    required_fields <- c("loaded_datasets", "dataset_summary", "loading_stats")
  } else if (current_step == "step_02_preprocessing") {
    required_fields <- c("processed_datasets", "preprocessing_summary", "processing_stats")
  } else if (current_step == "step_03_dge_analysis") {
    required_fields <- c("dge_results", "dge_summary", "analysis_stats")
  } else if (current_step == "step_04_meta_analysis") {
    required_fields <- c("meta_results", "meta_summary")
  } else if (current_step == "step_05_report_generator") {
    required_fields <- c("report_path", "report_info")
  } else {
    return(TRUE)  # Unknown step, assume valid
  }
  
  # Check required fields
  for (field in required_fields) {
    if (is.null(checkpoint_result$output_data[[field]])) {
      cat("⚠️  CHECKPOINT MISSING REQUIRED FIELD:", field, "\n")
      return(FALSE)
    }
  }
  
  return(TRUE)
}

#' Quick Pipeline Test
#'
#' Quick test function for development
#' @param config_file Configuration file path
test_dynamic_pipeline <- function(config_file = "config_dynamic_pipeline.yml") {
  cat("🧪 TESTING DYNAMIC PIPELINE\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  result <- execute_dynamic_pipeline(
    config_file = config_file,
    force_rerun = FALSE  # Use checkpoints if available
  )
  
  cat("\n🧪 TEST RESULT:", if (result$success) "SUCCESS ✅" else "FAILED ❌", "\n")
  
  return(result)
}

# Auto-execute if run as script
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("🚀 AUTO-EXECUTION MODE\n\n")
    result <- execute_dynamic_pipeline()
  } else {
    config_file <- args[1]
    cat("🚀 EXECUTION MODE with config:", config_file, "\n\n")
    result <- execute_dynamic_pipeline(config_file = config_file)
  }
  
  # Exit with appropriate code
  if (result$success) {
    cat("🎉 Pipeline execution completed successfully!\n")
    quit(status = 0)
  } else {
    cat("❌ Pipeline execution failed!\n")
    quit(status = 1)
  }
}

cat("✅ PIPELINE ORCHESTRATOR loaded successfully\n")
cat("🚀 Ready to execute: execute_dynamic_pipeline()\n")
cat("🧪 Quick test: test_dynamic_pipeline()\n\n")