#!/usr/bin/env Rscript
#' Configuration Validation System
#' 
#' Comprehensive validation framework for dynamic pipeline configuration
#' Ensures all required parameters are present and valid before pipeline execution

library(yaml)
library(tidyverse)

#' Load and Validate Pipeline Configuration
#'
#' @param config_file Path to YAML configuration file
#' @param strict_mode Enable strict validation (fail on warnings)
#' @return List with validated configuration and validation report
load_and_validate_config <- function(config_file = "config.yml", strict_mode = TRUE) {
  
  cat("🔍 CONFIGURATION VALIDATION\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Configuration file:", config_file, "\n")
  cat("Strict mode:", strict_mode, "\n\n")
  
  # Initialize validation report
  validation_report <- list(
    errors = character(0),
    warnings = character(0),
    info = character(0),
    success = FALSE
  )
  
  # Load configuration file
  tryCatch({
    config <- yaml::read_yaml(config_file)
    validation_report$info <- c(validation_report$info, "✅ Configuration file loaded successfully")
  }, error = function(e) {
    validation_report$errors <- c(validation_report$errors, paste("❌ Failed to load config file:", e$message))
    return(list(config = NULL, validation = validation_report))
  })
  
  # SECTION 1: Pipeline Architecture Validation
  cat("📋 Validating Pipeline Architecture...\n")
  
  if (is.null(config$pipeline)) {
    validation_report$errors <- c(validation_report$errors, "❌ Missing 'pipeline' section")
  } else {
    # Validate required pipeline fields
    required_pipeline_fields <- c("name", "version", "execution_mode", "steps", "dependencies")
    for (field in required_pipeline_fields) {
      if (is.null(config$pipeline[[field]])) {
        validation_report$errors <- c(validation_report$errors, paste("❌ Missing pipeline field:", field))
      }
    }
    
    # Validate execution mode
    if (!is.null(config$pipeline$execution_mode)) {
      if (!config$pipeline$execution_mode %in% c("dynamic", "legacy")) {
        validation_report$errors <- c(validation_report$errors, "❌ Invalid execution_mode. Must be 'dynamic' or 'legacy'")
      }
    }
    
    # Validate step dependencies
    if (!is.null(config$pipeline$steps) && !is.null(config$pipeline$dependencies)) {
      for (step in config$pipeline$steps) {
        if (is.null(config$pipeline$dependencies[[step]])) {
          validation_report$warnings <- c(validation_report$warnings, paste("⚠️  No dependencies specified for step:", step))
        }
      }
    }
    
    validation_report$info <- c(validation_report$info, "✅ Pipeline architecture validated")
  }
  
  # SECTION 2: Dataset Configuration Validation
  cat("📊 Validating Dataset Configuration...\n")
  
  if (is.null(config$datasets)) {
    validation_report$errors <- c(validation_report$errors, "❌ Missing 'datasets' section")
  } else {
    if (is.null(config$datasets$active_datasets)) {
      validation_report$errors <- c(validation_report$errors, "❌ No active datasets configured")
    } else {
      # Validate each dataset
      required_dataset_fields <- c("description", "disease_type", "control_type", "biological_context", "priority", "expected_samples")
      
      for (dataset_id in names(config$datasets$active_datasets)) {
        dataset_config <- config$datasets$active_datasets[[dataset_id]]
        
        for (field in required_dataset_fields) {
          if (is.null(dataset_config[[field]])) {
            validation_report$errors <- c(validation_report$errors, 
                                        paste("❌ Dataset", dataset_id, "missing field:", field))
          }
        }
        
        # Validate priority levels
        if (!is.null(dataset_config$priority)) {
          if (!dataset_config$priority %in% c("HIGH", "MODERATE", "LOW")) {
            validation_report$warnings <- c(validation_report$warnings,
                                          paste("⚠️  Dataset", dataset_id, "has invalid priority:", dataset_config$priority))
          }
        }
        
        # Check if cache files exist
        if (!is.null(dataset_config$cache_location)) {
          if (file.exists(dataset_config$cache_location)) {
            validation_report$info <- c(validation_report$info, paste("✅ Dataset", dataset_id, "cache file found"))
          } else {
            if (!is.null(dataset_config$status) && dataset_config$status == "optional") {
              validation_report$info <- c(validation_report$info, paste("ℹ️  Dataset", dataset_id, "marked optional, cache file not found"))
            } else {
              validation_report$warnings <- c(validation_report$warnings,
                                            paste("⚠️  Dataset", dataset_id, "cache file not found:", dataset_config$cache_location))
            }
          }
        }
      }
      
      # Count available datasets
      available_datasets <- length(config$datasets$active_datasets)
      validation_report$info <- c(validation_report$info, paste("✅", available_datasets, "datasets configured"))
    }
  }
  
  # SECTION 3: Gene Configuration Validation
  cat("🧬 Validating Gene Configuration...\n")
  
  if (is.null(config$genes)) {
    validation_report$errors <- c(validation_report$errors, "❌ Missing 'genes' section")
  } else {
    # Validate CAMK genes
    if (is.null(config$genes$camk_core_genes)) {
      validation_report$errors <- c(validation_report$errors, "❌ Missing CAMK core genes list")
    } else {
      expected_camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMKK1", "CAMKK2", 
                              "CAMK1", "CAMK1D", "CAMK1G", "CAMK4", "CAMKV")
      
      missing_genes <- setdiff(expected_camk_genes, config$genes$camk_core_genes)
      if (length(missing_genes) > 0) {
        validation_report$warnings <- c(validation_report$warnings, 
                                       paste("⚠️  Missing expected CAMK genes:", paste(missing_genes, collapse = ", ")))
      }
      
      if (!"CAMK2D" %in% config$genes$camk_core_genes) {
        validation_report$errors <- c(validation_report$errors, "❌ Primary gene CAMK2D not in core genes list")
      }
      
      validation_report$info <- c(validation_report$info, paste("✅", length(config$genes$camk_core_genes), "CAMK genes configured"))
    }
    
    # Validate minimum gene thresholds
    if (!is.null(config$genes$minimum_genes)) {
      required_priorities <- c("HIGH", "MODERATE", "LOW")
      for (priority in required_priorities) {
        if (is.null(config$genes$minimum_genes[[priority]])) {
          validation_report$warnings <- c(validation_report$warnings,
                                        paste("⚠️  Missing minimum genes threshold for priority:", priority))
        }
      }
    }
  }
  
  # SECTION 4: Analysis Parameters Validation
  cat("🧮 Validating Analysis Parameters...\n")
  
  if (is.null(config$analysis)) {
    validation_report$errors <- c(validation_report$errors, "❌ Missing 'analysis' section")
  } else {
    # Validate differential expression settings
    if (!is.null(config$analysis$differential_expression)) {
      de_config <- config$analysis$differential_expression
      
      # Check FDR threshold
      if (!is.null(de_config$fdr_threshold)) {
        if (de_config$fdr_threshold <= 0 || de_config$fdr_threshold >= 1) {
          validation_report$errors <- c(validation_report$errors, 
                                       "❌ FDR threshold must be between 0 and 1")
        }
      }
      
      # Check method
      if (!is.null(de_config$method)) {
        if (!de_config$method %in% c("limma", "DESeq2", "edgeR")) {
          validation_report$warnings <- c(validation_report$warnings,
                                        paste("⚠️  Unusual DE method:", de_config$method))
        }
      }
    }
    
    # Validate quality control thresholds
    if (!is.null(config$analysis$quality_control)) {
      qc_config <- config$analysis$quality_control
      
      if (!is.null(qc_config$high_logfc_threshold)) {
        if (qc_config$high_logfc_threshold <= 0) {
          validation_report$errors <- c(validation_report$errors,
                                       "❌ High logFC threshold must be positive")
        }
      }
    }
    
    validation_report$info <- c(validation_report$info, "✅ Analysis parameters validated")
  }
  
  # SECTION 5: File Paths Validation
  cat("📁 Validating File Paths...\n")
  
  if (is.null(config$paths)) {
    validation_report$errors <- c(validation_report$errors, "❌ Missing 'paths' section")
  } else {
    # Check if key directories exist
    key_dirs <- c("cache_root", "functions_dir", "scripts_dir")
    for (dir_key in key_dirs) {
      if (!is.null(config$paths[[dir_key]])) {
        if (dir.exists(config$paths[[dir_key]])) {
          validation_report$info <- c(validation_report$info, paste("✅ Directory found:", dir_key))
        } else {
          validation_report$warnings <- c(validation_report$warnings,
                                        paste("⚠️  Directory not found:", config$paths[[dir_key]]))
        }
      }
    }
    
    # Ensure output directories can be created
    if (!is.null(config$paths$output)) {
      for (output_dir in config$paths$output) {
        if (!dir.exists(output_dir)) {
          tryCatch({
            dir.create(output_dir, recursive = TRUE)
            validation_report$info <- c(validation_report$info, paste("✅ Created output directory:", output_dir))
          }, error = function(e) {
            validation_report$errors <- c(validation_report$errors,
                                        paste("❌ Cannot create output directory:", output_dir))
          })
        }
      }
    }
  }
  
  # SECTION 6: Final Validation Summary
  cat("\n📋 VALIDATION SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  # Print results
  if (length(validation_report$errors) > 0) {
    cat("❌ ERRORS FOUND:\n")
    for (error in validation_report$errors) {
      cat("  ", error, "\n")
    }
  }
  
  if (length(validation_report$warnings) > 0) {
    cat("\n⚠️  WARNINGS:\n")
    for (warning in validation_report$warnings) {
      cat("  ", warning, "\n")
    }
  }
  
  if (length(validation_report$info) > 0) {
    cat("\n✅ VALIDATION SUCCESS:\n")
    for (info in validation_report$info) {
      cat("  ", info, "\n")
    }
  }
  
  # Determine final validation status
  has_errors <- length(validation_report$errors) > 0
  has_warnings <- length(validation_report$warnings) > 0
  
  if (has_errors) {
    validation_report$success <- FALSE
    cat("\n❌ CONFIGURATION VALIDATION FAILED\n")
    cat("Please fix the errors above before running the pipeline.\n")
  } else if (has_warnings && strict_mode) {
    validation_report$success <- FALSE
    cat("\n⚠️  CONFIGURATION VALIDATION FAILED (Strict Mode)\n")
    cat("Please resolve warnings or disable strict mode.\n")
  } else {
    validation_report$success <- TRUE
    if (has_warnings) {
      cat("\n✅ CONFIGURATION VALIDATION PASSED (with warnings)\n")
    } else {
      cat("\n🎉 CONFIGURATION VALIDATION PASSED\n")
    }
    cat("Configuration is ready for pipeline execution.\n")
  }
  
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(list(
    config = config,
    validation = validation_report
  ))
}

#' Quick Configuration Check
#' 
#' Lightweight validation for development use
#' @param config_file Configuration file path
#' @return TRUE if config is valid, FALSE otherwise
validate_config_quick <- function(config_file = "config.yml") {
  result <- load_and_validate_config(config_file, strict_mode = FALSE)
  return(result$validation$success)
}

#' Configuration Summary
#'
#' Print a summary of the loaded configuration
#' @param config Loaded configuration object
print_config_summary <- function(config) {
  cat("🔧 CONFIGURATION SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  if (!is.null(config$pipeline)) {
    cat("Pipeline:", config$pipeline$name, "v", config$pipeline$version, "\n")
    cat("Execution mode:", config$pipeline$execution_mode, "\n")
    cat("Steps configured:", length(config$pipeline$steps), "\n")
  }
  
  if (!is.null(config$datasets$active_datasets)) {
    cat("Active datasets:", length(config$datasets$active_datasets), "\n")
    high_priority <- sum(sapply(config$datasets$active_datasets, function(x) x$priority == "HIGH"))
    cat("High priority datasets:", high_priority, "\n")
  }
  
  if (!is.null(config$genes$camk_core_genes)) {
    cat("CAMK genes:", length(config$genes$camk_core_genes), "\n")
  }
  
  if (!is.null(config$analysis$differential_expression$fdr_threshold)) {
    cat("FDR threshold:", config$analysis$differential_expression$fdr_threshold, "\n")
  }
  
  cat("═══════════════════════════════════════════════════════════════\n\n")
}

# Auto-validate if run as script (skip if called from pipeline)
if (!interactive() && length(commandArgs(trailingOnly = TRUE)) == 0 && 
    !exists("PIPELINE_EXECUTION_MODE", envir = globalenv())) {
  cat("🚀 AUTO-VALIDATION MODE\n\n")
  result <- load_and_validate_config()
  
  if (result$validation$success) {
    print_config_summary(result$config)
    cat("✅ Configuration validation completed successfully!\n")
    quit(status = 0)
  } else {
    cat("❌ Configuration validation failed!\n")
    quit(status = 1)
  }
}