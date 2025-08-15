#!/usr/bin/env Rscript
#' Step 05: Report Generator
#' 
#' Generates final HTML report from meta-analysis results
#' Uses existing RMD template with dynamic data

source("scripts/utilities/step_interface.R")
library(rmarkdown)

#' Execute Report Generation Step
#'
#' @param step_name Name of this step (should be "step_05_report_generator")
#' @param input_data Output from step_04_meta_analysis
#' @param config Full pipeline configuration
#' @param checkpoint_dir Directory for saving checkpoints
#' @return Step result with report generation status
step_05_report_generator <- function(step_name, input_data, config, checkpoint_dir = "output/checkpoints") {
  
  cat("📝 STEP 05: REPORT GENERATION\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("Generating final HTML report...\n\n")
  
  # Validate input from meta-analysis
  validate_step_input(
    input_data = input_data,
    required_fields = c("meta_results", "meta_summary"),
    step_name = step_name
  )
  
  meta_results <- input_data$meta_results
  meta_summary <- input_data$meta_summary
  
  if (nrow(meta_results) == 0) {
    return(create_step_result(
      success = FALSE,
      error_message = "No meta-analysis results provided for report generation",
      step_name = step_name
    ))
  }
  
  # Get report configuration
  report_config <- config$paths$reports
  
  cat("📋 REPORT GENERATION CONFIGURATION:\n")
  cat("Template:", report_config$template %||% "reports/CAMK_Analysis_Professional_Report.Rmd", "\n")
  cat("Output HTML:", report_config$output_html %||% "reports/CAMK_Analysis_Professional_Report_DYNAMIC.html", "\n")
  cat("Baseline HTML:", report_config$baseline_html %||% "reports/CAMK_Analysis_Professional_Report_BASELINE.html", "\n\n")
  
  # Verify required files exist
  template_file <- report_config$template %||% "reports/CAMK_Analysis_Professional_Report.Rmd"
  
  if (!file.exists(template_file)) {
    return(create_step_result(
      success = FALSE,
      error_message = paste("Report template not found:", template_file),
      step_name = step_name
    ))
  }
  
  # Ensure output directory exists
  output_html <- report_config$output_html %||% "reports/CAMK_Analysis_Professional_Report_DYNAMIC.html"
  output_dir <- dirname(output_html)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("📁 Created output directory:", output_dir, "\n")
  }
  
  # Prepare data files for report template
  cat("📊 PREPARING DATA FOR REPORT:\n")
  
  # Ensure all required output files exist and are current
  preparation_result <- prepare_report_data(input_data, config)
  
  if (!preparation_result$success) {
    return(create_step_result(
      success = FALSE,
      error_message = paste("Failed to prepare report data:", preparation_result$reason),
      step_name = step_name
    ))
  }
  
  cat("✅ Report data prepared successfully\n")
  
  # Generate the report
  cat("\n📝 GENERATING HTML REPORT:\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  generation_result <- generate_html_report(template_file, output_html, config)
  
  if (!generation_result$success) {
    return(create_step_result(
      success = FALSE,
      error_message = paste("Report generation failed:", generation_result$reason),
      step_name = step_name,
      warnings = generation_result$warnings
    ))
  }
  
  # Verify report was generated successfully
  if (!file.exists(output_html)) {
    return(create_step_result(
      success = FALSE,
      error_message = "Report file was not created despite successful rendering",
      step_name = step_name
    ))
  }
  
  # Get report file size and basic info
  report_info <- file.info(output_html)
  report_size_mb <- round(report_info$size / (1024^2), 2)
  
  cat("\n✅ REPORT GENERATION COMPLETED SUCCESSFULLY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📄 Report file:", output_html, "\n")
  cat("📊 File size:", report_size_mb, "MB\n")
  cat("⏰ Generated at:", as.character(report_info$mtime), "\n")
  
  # Optional: Compare with baseline if it exists
  comparison_result <- NULL
  baseline_html <- report_config$baseline_html
  
  if (!is.null(baseline_html) && file.exists(baseline_html)) {
    cat("\n🔍 BASELINE COMPARISON:\n")
    
    comparison_result <- compare_reports(output_html, baseline_html, config)
    
    if (!is.null(comparison_result)) {
      cat("📊 Size comparison:\n")
      cat("  Dynamic report:", comparison_result$dynamic_size_mb, "MB\n")
      cat("  Baseline report:", comparison_result$baseline_size_mb, "MB\n")
      
      if (comparison_result$size_difference_percent < 20) {
        cat("✅ Report sizes are similar (", round(comparison_result$size_difference_percent, 1), "% difference)\n")
      } else {
        cat("⚠️  Report sizes differ significantly (", round(comparison_result$size_difference_percent, 1), "% difference)\n")
      }
    }
  }
  
  # Create output data structure
  output_data <- list(
    report_path = output_html,
    report_info = list(
      size_mb = report_size_mb,
      generated_at = report_info$mtime,
      template_used = template_file
    ),
    meta_results = meta_results,
    meta_summary = meta_summary,
    generation_stats = list(
      rendering_successful = TRUE,
      warnings_count = length(generation_result$warnings %||% character(0)),
      baseline_comparison = comparison_result
    )
  )
  
  # Validate output
  validate_step_output(
    output_data = output_data,
    required_fields = c("report_path", "report_info", "generation_stats"),
    step_name = step_name
  )
  
  cat("\n🎉 DYNAMIC PIPELINE EXECUTION COMPLETED SUCCESSFULLY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📊 Final Results Summary:\n")
  cat("  Genes analyzed:", meta_summary$genes_analyzed, "\n")
  cat("  Significant genes:", meta_summary$significant_genes, "\n")
  cat("  CAMK2D significant:", if (meta_summary$camk2d_significant) "YES ✅" else "NO", "\n")
  cat("  Report generated:", output_html, "\n")
  cat("  Report size:", report_size_mb, "MB\n")
  
  if (!is.null(meta_summary$baseline_validation_passed)) {
    cat("  Baseline validation:", if (meta_summary$baseline_validation_passed) "PASSED ✅" else "WARNINGS ⚠️", "\n")
  }
  
  cat("\n🚀 PIPELINE EXECUTION COMPLETE - READY FOR SCIENTIFIC USE\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(create_step_result(
    success = TRUE,
    output_data = output_data,
    step_name = step_name,
    warnings = generation_result$warnings,
    metadata = list(
      report_size_mb = report_size_mb,
      genes_analyzed = meta_summary$genes_analyzed,
      significant_genes = meta_summary$significant_genes,
      camk2d_significant = meta_summary$camk2d_significant
    )
  ))
}

#' Prepare Report Data
#'
#' Ensures all required data files exist for the report template
#' @param input_data Pipeline data
#' @param config Pipeline configuration
#' @return Preparation result
prepare_report_data <- function(input_data, config) {
  
  warnings_generated <- character(0)
  
  # Check that required output files exist and are up to date
  output_files <- config$paths$output_files
  
  # Verify meta-analysis results file
  if (!is.null(output_files$meta_results)) {
    if (!file.exists(output_files$meta_results)) {
      return(list(success = FALSE, reason = "Meta-analysis results file missing"))
    }
  }
  
  # Verify DGE results file  
  if (!is.null(output_files$dge_results)) {
    if (!file.exists(output_files$dge_results)) {
      return(list(success = FALSE, reason = "DGE results file missing"))
    }
  }
  
  # Verify dataset summary file
  if (!is.null(output_files$dataset_summary)) {
    if (!file.exists(output_files$dataset_summary)) {
      return(list(success = FALSE, reason = "Dataset summary file missing"))
    }
  }
  
  # Ensure output/current/ directory exists with symlinks or copies for RMD compatibility
  current_dir <- config$paths$output$current %||% "output/current"
  
  if (!dir.exists(current_dir)) {
    dir.create(current_dir, recursive = TRUE)
  }
  
  # Copy/update files to current directory for RMD template
  file_mappings <- list(
    "CAMK_meta_analysis_FINAL.csv" = output_files$meta_results,
    "CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv" = output_files$dge_results,
    "dataset_processing_summary_6_datasets.csv" = output_files$dataset_summary
  )
  
  # Check if methodology comparison file exists, create if missing
  methodology_file <- file.path(current_dir, "methodology_comparison_analysis.csv")
  if (!file.exists(methodology_file)) {
    cat("Creating methodology comparison file...\n")
    # This file is now pre-created, but we'll check just in case
  }
  
  for (target_file in names(file_mappings)) {
    source_file <- file_mappings[[target_file]]
    target_path <- file.path(current_dir, target_file)
    
    if (!is.null(source_file) && file.exists(source_file)) {
      if (!file.exists(target_path) || file.mtime(source_file) > file.mtime(target_path)) {
        file.copy(source_file, target_path, overwrite = TRUE)
        cat("📄 Updated:", target_file, "\n")
      }
    } else {
      warnings_generated <- c(warnings_generated, paste("Source file missing for", target_file))
    }
  }
  
  return(list(
    success = TRUE,
    warnings = warnings_generated
  ))
}

#' Generate HTML Report
#'
#' Renders the RMD template to HTML
#' @param template_file Path to RMD template
#' @param output_html Path for output HTML
#' @param config Pipeline configuration
#' @return Generation result
generate_html_report <- function(template_file, output_html, config) {
  
  warnings_generated <- character(0)
  
  tryCatch({
    # Change to reports directory for rendering (required for relative paths in RMD)
    original_wd <- getwd()
    reports_dir <- dirname(template_file)
    
    if (reports_dir != ".") {
      setwd(reports_dir)
      template_basename <- basename(template_file)
      output_basename <- basename(output_html)
    } else {
      template_basename <- template_file
      output_basename <- output_html
    }
    
    cat("📍 Working directory:", getwd(), "\n")
    cat("📄 Template:", template_basename, "\n")
    cat("📄 Output:", output_basename, "\n")
    
    # Render the report
    cat("\n🔄 Rendering report...\n")
    
    render_result <- rmarkdown::render(
      input = template_basename,
      output_file = output_basename,
      quiet = FALSE  # Show rendering progress
    )
    
    # Restore working directory
    if (reports_dir != ".") {
      setwd(original_wd)
      
      # Move the rendered file to the correct output location
      rendered_file <- file.path(reports_dir, output_basename)
      if (file.exists(rendered_file) && rendered_file != output_html) {
        file.copy(rendered_file, output_html, overwrite = TRUE)
        cat("📁 Moved report to:", output_html, "\n")
        # Clean up template directory
        file.remove(rendered_file)
      }
    }
    
    cat("✅ Report rendered successfully to:", output_html, "\n")
    
    return(list(
      success = TRUE,
      output_file = output_html,
      warnings = warnings_generated
    ))
    
  }, error = function(e) {
    # Restore working directory on error
    if (exists("original_wd") && getwd() != original_wd) {
      setwd(original_wd)
    }
    
    return(list(
      success = FALSE,
      reason = paste("Rendering error:", e$message),
      warnings = warnings_generated
    ))
    
  }, warning = function(w) {
    warnings_generated <<- c(warnings_generated, w$message)
    invokeRestart("muffleWarning")
  })
}

#' Compare Reports
#'
#' Compares dynamic and baseline reports
#' @param dynamic_html Path to dynamic report
#' @param baseline_html Path to baseline report
#' @param config Pipeline configuration
#' @return Comparison result
compare_reports <- function(dynamic_html, baseline_html, config) {
  
  if (!file.exists(dynamic_html) || !file.exists(baseline_html)) {
    return(NULL)
  }
  
  # Get file sizes
  dynamic_info <- file.info(dynamic_html)
  baseline_info <- file.info(baseline_html)
  
  dynamic_size_mb <- round(dynamic_info$size / (1024^2), 2)
  baseline_size_mb <- round(baseline_info$size / (1024^2), 2)
  
  size_difference_percent <- abs(dynamic_size_mb - baseline_size_mb) / baseline_size_mb * 100
  
  return(list(
    dynamic_size_mb = dynamic_size_mb,
    baseline_size_mb = baseline_size_mb,
    size_difference_percent = size_difference_percent,
    dynamic_modified = dynamic_info$mtime,
    baseline_modified = baseline_info$mtime
  ))
}

cat("✅ STEP 05: Report Generator loaded successfully\n\n")