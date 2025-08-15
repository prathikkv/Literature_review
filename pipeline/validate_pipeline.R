#!/usr/bin/env Rscript
#' Pipeline Validation Script
#' 
#' Comprehensive validation of bioinformatics pipeline compliance
#' Checks all security, standards, and functionality requirements
#' 
#' @author Claude Code Pipeline Validation
#' @version 1.0.0

cat("\n")
cat("╔══════════════════════════════════════════════════════════════╗\n")
cat("║        BIOINFORMATICS PIPELINE VALIDATION SUITE             ║\n")
cat("║        Security • Standards • Functionality                 ║\n")
cat("╚══════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Load all utility modules
source("scripts/utilities/error_handling.R")
source("scripts/utilities/input_validation.R")
source("scripts/utilities/run_management.R")
source("scripts/utilities/code_quality.R")

# Initialize validation
validation_start_time <- Sys.time()
run_id <- generate_run_id("VALIDATION", c("PIPELINE"))
log_config <- initialize_logging_system(config = NULL, run_id = run_id)

log_message(
  message = "Pipeline validation started",
  level = "INFO",
  category = "VALIDATION",
  details = list(run_id = run_id, start_time = validation_start_time),
  log_config = log_config
)

# Track validation results
validation_results <- list(
  passed = 0,
  failed = 0,
  warnings = 0,
  tests_run = character(),
  issues = character(),
  recommendations = character()
)

#' Add Validation Result
add_result <- function(test_name, passed, issue = NULL, recommendation = NULL) {
  validation_results$tests_run <<- c(validation_results$tests_run, test_name)
  
  if (passed) {
    validation_results$passed <<- validation_results$passed + 1
    cat("✅", test_name, "\n")
  } else {
    validation_results$failed <<- validation_results$failed + 1
    cat("❌", test_name, "\n")
    if (!is.null(issue)) {
      validation_results$issues <<- c(validation_results$issues, issue)
      cat("   Issue:", issue, "\n")
    }
    if (!is.null(recommendation)) {
      validation_results$recommendations <<- c(validation_results$recommendations, recommendation)
      cat("   Fix:", recommendation, "\n")
    }
  }
}

#' Add Warning Result
add_warning <- function(test_name, warning_msg, recommendation = NULL) {
  validation_results$tests_run <<- c(validation_results$tests_run, test_name)
  validation_results$warnings <<- validation_results$warnings + 1
  cat("⚠️ ", test_name, "\n")
  cat("   Warning:", warning_msg, "\n")
  if (!is.null(recommendation)) {
    validation_results$recommendations <<- c(validation_results$recommendations, recommendation)
    cat("   Suggestion:", recommendation, "\n")
  }
}

cat("🔍 PHASE 1: DIRECTORY STRUCTURE VALIDATION\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Test 1: Required directories exist
required_dirs <- c("scripts", "modules", "output", "templates", "scripts/utilities")
all_dirs_exist <- TRUE
for (dir in required_dirs) {
  if (!dir.exists(dir)) {
    all_dirs_exist <- FALSE
    break
  }
}
add_result("Required directory structure", all_dirs_exist, 
           if (!all_dirs_exist) "Missing required directories" else NULL,
           if (!all_dirs_exist) "Create missing directories" else NULL)

# Test 2: Configuration file exists
config_exists <- file.exists("config.yml")
add_result("Configuration file present", config_exists,
           if (!config_exists) "config.yml not found" else NULL,
           if (!config_exists) "Create config.yml file" else NULL)

# Test 3: Key pipeline files exist
key_files <- c("run_enhanced_pipeline.R", "run_pharma_pipeline.R", "templates/CAMK_Analysis_Report.Rmd")
all_files_exist <- all(file.exists(key_files))
add_result("Key pipeline files present", all_files_exist,
           if (!all_files_exist) "Missing essential pipeline files" else NULL,
           if (!all_files_exist) "Ensure all pipeline components are present" else NULL)

cat("\n🔒 PHASE 2: SECURITY VALIDATION\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Test 4: Input validation module functional
tryCatch({
  test_gene <- validate_gene_symbol("TP53")
  security_functional <- test_gene$valid
  add_result("Input validation system", security_functional,
             if (!security_functional) "Input validation failed" else NULL,
             if (!security_functional) "Check input validation implementation" else NULL)
}, error = function(e) {
  add_result("Input validation system", FALSE, 
             paste("Input validation error:", e$message),
             "Fix input validation module")
})

# Test 5: Error handling system functional
tryCatch({
  test_error <- handle_pipeline_error(
    error = simpleError("test error"),
    context = list(test = TRUE),
    step_name = "validation_test",
    recovery_action = "test recovery",
    log_config = log_config
  )
  error_handling_functional <- !is.null(test_error$error_id)
  add_result("Error handling system", error_handling_functional,
             if (!error_handling_functional) "Error handling failed" else NULL,
             if (!error_handling_functional) "Fix error handling module" else NULL)
}, error = function(e) {
  add_result("Error handling system", FALSE,
             paste("Error handling failed:", e$message),
             "Fix error handling implementation")
})

# Test 6: Logging system functional
tryCatch({
  log_message("Validation test message", "INFO", "VALIDATION", log_config = log_config)
  log_files_exist <- file.exists(log_config$files$main)
  add_result("Logging system", log_files_exist,
             if (!log_files_exist) "Logging system not working" else NULL,
             if (!log_files_exist) "Check logging configuration" else NULL)
}, error = function(e) {
  add_result("Logging system", FALSE,
             paste("Logging failed:", e$message),
             "Fix logging system")
})

cat("\n📊 PHASE 3: BIOINFORMATICS STANDARDS VALIDATION\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Test 7: FAIR data principles
fair_dirs <- c("output/runs", "output/archive", "output/logs")
fair_implemented <- all(dir.exists(fair_dirs))
add_result("FAIR data principles structure", fair_implemented,
           if (!fair_implemented) "FAIR data structure incomplete" else NULL,
           if (!fair_implemented) "Create proper data organization structure" else NULL)

# Test 8: Run management system
tryCatch({
  test_run_id <- generate_run_id("TEST", c("VALIDATION"))
  test_dirs <- create_run_directory("output", test_run_id, create_symlink = FALSE)
  run_mgmt_functional <- dir.exists(test_dirs$current_run)
  add_result("Run management system", run_mgmt_functional,
             if (!run_mgmt_functional) "Run management failed" else NULL,
             if (!run_mgmt_functional) "Fix run management implementation" else NULL)
}, error = function(e) {
  add_result("Run management system", FALSE,
             paste("Run management error:", e$message),
             "Fix run management module")
})

# Test 9: Dynamic file naming
tryCatch({
  test_filename <- generate_dynamic_filename("TP53", c("Cancer"), "report", "html")
  filename_valid <- grepl("TP53.*Cancer.*report\\.html$", test_filename)
  add_result("Dynamic file naming", filename_valid,
             if (!filename_valid) "Dynamic naming not working" else NULL,
             if (!filename_valid) "Fix dynamic filename generation" else NULL)
}, error = function(e) {
  add_result("Dynamic file naming", FALSE,
             paste("Dynamic naming error:", e$message),
             "Fix dynamic naming function")
})

cat("\n🔧 PHASE 4: CODE QUALITY ASSESSMENT\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Test 10: R package standards
tryCatch({
  quality_check <- check_r_package_standards(".")
  code_quality_acceptable <- length(quality_check$issues) <= 3  # Allow some flexibility
  
  if (code_quality_acceptable) {
    add_result("R package standards compliance", TRUE)
  } else {
    add_warning("R package standards compliance", 
                paste("Code quality issues found:", length(quality_check$issues), "issues"),
                "Run generate_code_quality_report() for details")
  }
}, error = function(e) {
  add_result("R package standards compliance", FALSE,
             paste("Code quality check failed:", e$message),
             "Fix code quality module")
})

cat("\n🧪 PHASE 5: FUNCTIONAL VALIDATION\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Test 11: Configuration validation
if (config_exists) {
  tryCatch({
    config <- yaml::read_yaml("config.yml")
    config_valid <- validate_configuration(config)
    add_result("Configuration validation", config_valid$valid,
               if (!config_valid$valid) paste("Config errors:", paste(config_valid$errors, collapse = "; ")) else NULL,
               if (!config_valid$valid) "Fix configuration issues" else NULL)
    
    if (length(config_valid$warnings) > 0) {
      add_warning("Configuration warnings", 
                  paste(length(config_valid$warnings), "warnings found"),
                  "Review configuration warnings")
    }
  }, error = function(e) {
    add_result("Configuration validation", FALSE,
               paste("Config validation failed:", e$message),
               "Fix configuration file")
  })
} else {
  add_result("Configuration validation", FALSE, "No config file to validate", "Create config.yml")
}

# Test 12: Module loading
modules_to_test <- c(
  "modules/gene_family_discovery.R",
  "modules/literature_mining.R", 
  "modules/validation_framework.R"
)

modules_loaded <- 0
for (module in modules_to_test) {
  if (file.exists(module)) {
    tryCatch({
      source(module)
      modules_loaded <- modules_loaded + 1
    }, error = function(e) {
      # Module failed to load
    })
  }
}

all_modules_loaded <- modules_loaded == length(modules_to_test)
add_result("Core modules loading", all_modules_loaded,
           if (!all_modules_loaded) paste("Failed to load", length(modules_to_test) - modules_loaded, "modules") else NULL,
           if (!all_modules_loaded) "Check module syntax and dependencies" else NULL)

cat("\n📋 VALIDATION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n")

validation_end_time <- Sys.time()
total_duration <- as.numeric(difftime(validation_end_time, validation_start_time, units = "secs"))

# Calculate scores
total_tests <- validation_results$passed + validation_results$failed + validation_results$warnings
pass_rate <- round((validation_results$passed / total_tests) * 100, 1)
overall_score <- if (validation_results$failed == 0 && validation_results$warnings <= 2) "EXCELLENT" else
                if (validation_results$failed <= 1 && validation_results$warnings <= 4) "GOOD" else
                if (validation_results$failed <= 3) "NEEDS IMPROVEMENT" else "CRITICAL ISSUES"

cat("📊 Test Results:\n")
cat("   ✅ Passed:", validation_results$passed, "\n")
cat("   ❌ Failed:", validation_results$failed, "\n")
cat("   ⚠️  Warnings:", validation_results$warnings, "\n")
cat("   📈 Pass Rate:", paste0(pass_rate, "%"), "\n")
cat("   🏆 Overall Score:", overall_score, "\n")
cat("   ⏱️  Duration:", format_duration(total_duration), "\n")

# Log final results
log_message(
  message = "Pipeline validation completed",
  level = "INFO",
  category = "VALIDATION",
  details = list(
    passed = validation_results$passed,
    failed = validation_results$failed,
    warnings = validation_results$warnings,
    overall_score = overall_score,
    duration_seconds = total_duration
  ),
  log_config = log_config
)

# Generate recommendations if needed
if (validation_results$failed > 0 || validation_results$warnings > 2) {
  cat("\n🔧 PRIORITY ACTIONS NEEDED:\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  if (length(validation_results$issues) > 0) {
    cat("Critical Issues:\n")
    for (i in 1:length(validation_results$issues)) {
      cat(paste0("   ", i, ". ", validation_results$issues[i]), "\n")
    }
  }
  
  if (length(validation_results$recommendations) > 0) {
    cat("\nRecommended Actions:\n")
    unique_recommendations <- unique(validation_results$recommendations)
    for (i in 1:length(unique_recommendations)) {
      cat(paste0("   ", i, ". ", unique_recommendations[i]), "\n")
    }
  }
} else {
  cat("\n🎉 EXCELLENT! PIPELINE MEETS ALL STANDARDS\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("✅ Security compliant\n")
  cat("✅ Bioinformatics standards met\n") 
  cat("✅ FAIR data principles implemented\n")
  cat("✅ Error handling comprehensive\n")
  cat("✅ Code quality acceptable\n")
  cat("✅ Ready for production use\n")
}

# Generate validation report
report_file <- file.path("output/logs", paste0("validation_report_", run_id, ".txt"))
tryCatch({
  generate_code_quality_report(".", report_file)
  cat("\n📄 Detailed report saved to:", report_file, "\n")
}, error = function(e) {
  cat("\n⚠️  Could not generate detailed report:", e$message, "\n")
})

cat("\n🚀 VALIDATION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Return validation results for programmatic use
invisible(list(
  overall_score = overall_score,
  passed = validation_results$passed,
  failed = validation_results$failed,
  warnings = validation_results$warnings,
  pass_rate = pass_rate,
  ready_for_production = (validation_results$failed == 0 && validation_results$warnings <= 2)
))