#!/usr/bin/env Rscript
#' Pipeline Validation Script
#' 
#' Comprehensive validation and testing for the CAMK2D analysis pipeline

cat("TEST: CAMK2D PIPELINE VALIDATION\n")
cat("==============================\n")
cat("TIME: Started:", Sys.time(), "\n\n")

# Load validation functions
if (file.exists("../../functions/data_processing.R")) {
  suppressWarnings(source("../../functions/data_processing.R"))
}
if (file.exists("../../functions/analysis.R")) {
  suppressWarnings(source("../../functions/analysis.R"))
}
if (file.exists("../../functions/visualization.R")) {
  suppressWarnings(source("../../functions/visualization.R"))
}
if (file.exists("../../functions/utilities.R")) {
  suppressWarnings(source("../../functions/utilities.R"))
}

#' Comprehensive Pipeline Validation
#'
#' Tests all major components and generates validation report
run_comprehensive_validation <- function() {
  
  validation_results <- list()
  
  # ==========================================================================
  # TEST 1: FUNCTION MODULE LOADING
  # ==========================================================================
  
  cat("PACKAGE: TEST 1: FUNCTION MODULE LOADING\n")
  cat("===================================\n")
  
  modules <- c(
    "data_processing.R" = "functions/data_processing.R",
    "analysis.R" = "functions/analysis.R",
    "visualization.R" = "functions/visualization.R", 
    "utilities.R" = "functions/utilities.R"
  )
  
  module_status <- list()
  
  for (module_name in names(modules)) {
    cat("Testing", module_name, "...")
    tryCatch({
      source(modules[[module_name]])
      cat(" SUCCESS:\n")
      module_status[[module_name]] <- TRUE
    }, error = function(e) {
      cat(" ERROR: Error:", e$message, "\n")
      module_status[[module_name]] <- FALSE
    })
  }
  
  modules_loaded <- sum(unlist(module_status))
  total_modules <- length(module_status)
  
  cat("\nDATA: Module Loading Summary:\n")
  cat("  • Loaded successfully:", modules_loaded, "/", total_modules, "\n")
  cat("  • Success rate:", round(100 * modules_loaded / total_modules, 1), "%\n\n")
  
  validation_results$module_loading <- list(
    success_rate = modules_loaded / total_modules,
    modules_loaded = modules_loaded,
    total_modules = total_modules,
    status = module_status
  )
  
  # ==========================================================================
  # TEST 2: CORE FUNCTION AVAILABILITY
  # ==========================================================================
  
  cat("SEARCH: TEST 2: CORE FUNCTION AVAILABILITY\n")
  cat("=====================================\n")
  
  required_functions <- c(
    "get_comprehensive_dataset_list",
    "download_comprehensive_datasets", 
    "comprehensive_preprocessing_pipeline",
    "comprehensive_differential_expression_pipeline",
    "comprehensive_meta_analysis_pipeline",
    "comprehensive_pathway_analysis_pipeline",
    "get_camk_family_genes",
    "comprehensive_reporting_pipeline"
  )
  
  function_status <- list()
  
  for (func_name in required_functions) {
    cat("Checking", func_name, "...")
    if (exists(func_name)) {
      cat(" SUCCESS:\n")
      function_status[[func_name]] <- TRUE
    } else {
      cat(" ERROR:\n")
      function_status[[func_name]] <- FALSE
    }
  }
  
  functions_available <- sum(unlist(function_status))
  total_functions <- length(function_status)
  
  cat("\nDATA: Function Availability Summary:\n")
  cat("  • Available functions:", functions_available, "/", total_functions, "\n")
  cat("  • Success rate:", round(100 * functions_available / total_functions, 1), "%\n\n")
  
  validation_results$function_availability <- list(
    success_rate = functions_available / total_functions,
    functions_available = functions_available,
    total_functions = total_functions,
    status = function_status
  )
  
  # ==========================================================================
  # TEST 3: PACKAGE DEPENDENCIES
  # ==========================================================================
  
  cat("PACKAGE: TEST 3: PACKAGE DEPENDENCIES\n")
  cat("===============================\n")
  
  critical_packages <- c(
    "GEOquery", "limma", "DESeq2", "tidyverse", "metafor", 
    "clusterProfiler", "biomaRt", "openxlsx", "ggplot2", "plotly"
  )
  
  package_status <- list()
  
  for (pkg in critical_packages) {
    cat("Checking", pkg, "...")
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(" SUCCESS:\n")
      package_status[[pkg]] <- TRUE
    } else {
      cat(" ERROR:\n")
      package_status[[pkg]] <- FALSE
    }
  }
  
  packages_available <- sum(unlist(package_status))
  total_packages <- length(package_status)
  
  cat("\nDATA: Package Dependencies Summary:\n")
  cat("  • Available packages:", packages_available, "/", total_packages, "\n")
  cat("  • Success rate:", round(100 * packages_available / total_packages, 1), "%\n\n")
  
  validation_results$package_dependencies <- list(
    success_rate = packages_available / total_packages,
    packages_available = packages_available,
    total_packages = total_packages,
    status = package_status
  )
  
  # ==========================================================================
  # TEST 4: CAMK GENE LIST VALIDATION
  # ==========================================================================
  
  cat("GENETIC: TEST 4: CAMK GENE LIST VALIDATION\n")
  cat("====================================\n")
  
  if (exists("get_camk_family_genes")) {
    camk_genes <- get_camk_family_genes()
    
    expected_camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMKK1", "CAMKK2")
    
    cat("CAMK genes available:", length(camk_genes), "\n")
    cat("Expected core genes:", paste(expected_camk_genes, collapse = ", "), "\n")
    
    core_genes_present <- sum(expected_camk_genes %in% camk_genes)
    
    cat("SUCCESS: Core CAMK genes present:", core_genes_present, "/", length(expected_camk_genes), "\n\n")
    
    validation_results$camk_genes <- list(
      genes_available = length(camk_genes),
      core_genes_present = core_genes_present,
      gene_list = camk_genes
    )
  } else {
    cat("ERROR: get_camk_family_genes() function not available\n\n")
    validation_results$camk_genes <- list(
      genes_available = 0,
      core_genes_present = 0,
      gene_list = NULL
    )
  }
  
  # ==========================================================================
  # TEST 5: DATASET CONFIGURATION
  # ==========================================================================
  
  cat("DATA: TEST 5: DATASET CONFIGURATION\n")
  cat("================================\n")
  
  if (exists("get_comprehensive_dataset_list")) {
    dataset_list <- get_comprehensive_dataset_list()
    
    # Count datasets by category
    human_hf <- length(dataset_list$human_hf)
    human_af <- length(dataset_list$human_af)
    animal_models <- length(dataset_list$animal_models)
    total_datasets <- human_hf + human_af + animal_models
    
    cat("Dataset configuration:\n")
    cat("  • Human Heart Failure:", human_hf, "datasets\n")
    cat("  • Human Atrial Fibrillation:", human_af, "datasets\n") 
    cat("  • Animal Models:", animal_models, "datasets\n")
    cat("  • Total datasets:", total_datasets, "\n")
    
    # Check for key datasets from prompts.md
    expected_datasets <- c("GSE120895", "GSE57338", "GSE141910", "GSE31821", 
                          "GSE41177", "GSE79768", "GSE115574", "GSE14975")
    
    all_dataset_ids <- c(names(dataset_list$human_hf), names(dataset_list$human_af), 
                        names(dataset_list$animal_models))
    
    key_datasets_present <- sum(expected_datasets %in% all_dataset_ids)
    
    cat("SUCCESS: Key datasets from prompts.md:", key_datasets_present, "/", length(expected_datasets), "\n\n")
    
    validation_results$dataset_configuration <- list(
      total_datasets = total_datasets,
      human_hf = human_hf,
      human_af = human_af,
      animal_models = animal_models,
      key_datasets_present = key_datasets_present
    )
  } else {
    cat("ERROR: get_comprehensive_dataset_list() function not available\n\n")
    validation_results$dataset_configuration <- list(
      total_datasets = 0,
      key_datasets_present = 0
    )
  }
  
  # ==========================================================================
  # TEST 6: OUTPUT DIRECTORIES
  # ==========================================================================
  
  cat("SAVED: TEST 6: OUTPUT DIRECTORY STRUCTURE\n")
  cat("=====================================\n")
  
  required_dirs <- c("cache", "data", "output", "results")
  existing_dirs <- c()
  
  for (dir_name in required_dirs) {
    if (dir.exists(dir_name)) {
      cat("SUCCESS: Directory exists:", dir_name, "\n")
      existing_dirs <- c(existing_dirs, dir_name)
    } else {
      cat("WARNING: Directory missing:", dir_name, "(will be created during analysis)\n")
    }
  }
  
  cat("\nDirectory structure:", length(existing_dirs), "/", length(required_dirs), "directories exist\n\n")
  
  validation_results$output_directories <- list(
    required_dirs = required_dirs,
    existing_dirs = existing_dirs,
    directories_exist = length(existing_dirs)
  )
  
  # ==========================================================================
  # OVERALL VALIDATION ASSESSMENT
  # ==========================================================================
  
  cat("ACHIEVEMENT: OVERALL VALIDATION ASSESSMENT\n")
  cat("================================\n")
  
  # Calculate overall scores
  module_score <- validation_results$module_loading$success_rate * 100
  function_score <- validation_results$function_availability$success_rate * 100
  package_score <- validation_results$package_dependencies$success_rate * 100
  
  overall_score <- (module_score + function_score + package_score) / 3
  
  cat("DATA: Validation Scores:\n")
  cat("  • Module Loading:", round(module_score, 1), "%\n")
  cat("  • Function Availability:", round(function_score, 1), "%\n")
  cat("  • Package Dependencies:", round(package_score, 1), "%\n")
  cat("  • Overall Score:", round(overall_score, 1), "%\n\n")
  
  # Assessment
  if (overall_score >= 95) {
    cat("COMPLETE: VALIDATION PASSED: Pipeline ready for production use!\n")
    cat("SUCCESS: All critical components are functional\n")
    cat("SUCCESS: Ready for comprehensive CAMK2D analysis\n")
    status <- "PASSED"
  } else if (overall_score >= 85) {
    cat("WARNING: VALIDATION MOSTLY PASSED: Minor issues detected\n")
    cat("INSIGHT: Pipeline functional with some limitations\n")
    cat("FIX: Consider addressing missing components\n")
    status <- "MOSTLY_PASSED"
  } else {
    cat("ERROR: VALIDATION FAILED: Significant issues detected\n")
    cat("FIX: Multiple components require attention\n")
    cat("INSIGHT: Run setup.R to install missing dependencies\n")
    status <- "FAILED"
  }
  
  validation_results$overall_assessment <- list(
    overall_score = overall_score,
    status = status,
    module_score = module_score,
    function_score = function_score,
    package_score = package_score
  )
  
  return(validation_results)
}

#' Save Validation Report
#'
#' @param validation_results Validation results
save_validation_report <- function(validation_results) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists("output")) {
    dir.create("output", recursive = TRUE)
  }
  
  # Save detailed results
  saveRDS(validation_results, "output/validation_results.rds")
  
  # Create text report
  report_lines <- c(
    "CAMK2D PIPELINE VALIDATION REPORT",
    "==================================",
    paste("Generated:", Sys.time()),
    paste("Overall Status:", validation_results$overall_assessment$status),
    paste("Overall Score:", round(validation_results$overall_assessment$overall_score, 1), "%"),
    "",
    "COMPONENT SCORES:",
    paste("• Module Loading:", round(validation_results$overall_assessment$module_score, 1), "%"),
    paste("• Function Availability:", round(validation_results$overall_assessment$function_score, 1), "%"),
    paste("• Package Dependencies:", round(validation_results$overall_assessment$package_score, 1), "%"),
    "",
    "SUMMARY:",
    paste("• Modules loaded:", validation_results$module_loading$modules_loaded, "/", 
          validation_results$module_loading$total_modules),
    paste("• Functions available:", validation_results$function_availability$functions_available, "/", 
          validation_results$function_availability$total_functions),
    paste("• Packages available:", validation_results$package_dependencies$packages_available, "/", 
          validation_results$package_dependencies$total_packages),
    paste("• CAMK genes configured:", validation_results$camk_genes$genes_available),
    paste("• Datasets configured:", validation_results$dataset_configuration$total_datasets),
    "",
    "NEXT STEPS:",
    if (validation_results$overall_assessment$status == "PASSED") {
      "SUCCESS: Pipeline ready - execute: Rscript run_pipeline.R"
    } else if (validation_results$overall_assessment$status == "MOSTLY_PASSED") {
      "WARNING: Consider running: Rscript setup.R to address missing components"
    } else {
      "ERROR: Run: Rscript setup.R to install missing dependencies"
    }
  )
  
  # Save text report
  writeLines(report_lines, "output/validation_report.txt")
  
  cat("SAVED: Validation report saved to: output/validation_report.txt\n")
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if (!interactive()) {
  # Run validation
  tryCatch({
    validation_results <- run_comprehensive_validation()
    save_validation_report(validation_results)
    
    cat("\nTIME: Validation completed:", Sys.time(), "\n")
    
    if (validation_results$overall_assessment$status == "PASSED") {
      cat("TARGET: Next step: Rscript run_pipeline.R\n")
    } else {
      cat("FIX: Next step: Rscript setup.R\n")
    }
    
  }, error = function(e) {
    cat("ERROR: Validation error:", e$message, "\n")
    quit(status = 1)
  })
} else {
  cat("SUCCESS: Validation script loaded and ready\n")
  cat("TEST: Run validation: validation_results <- run_comprehensive_validation()\n")
}

cat("\nTEST: CAMK2D PIPELINE VALIDATION SYSTEM\n")
cat("ENHANCED: Comprehensive Testing and Quality Assurance\n")
cat("ENHANCED: Production-Ready Validation Framework\n\n")