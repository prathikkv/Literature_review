#!/usr/bin/env Rscript
#' Step 05: Report Generator
#' 
#' Generates final HTML report from meta-analysis results
#' Uses existing RMD template with dynamic data

source("scripts/utilities/step_interface.R")
source("scripts/utilities/run_management.R")
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
  
  # SIMPLE REPORT GENERATION (V1.0.0 COMPATIBLE)
  primary_gene <- config$research_target$primary_gene %||% "UNKNOWN"
  diseases <- config$research_target$diseases %||% c("Unknown")
  
  # Use existing output/current directory structure (no timestamped directories)
  current_output_dir <- "output/current"
  if (!dir.exists(current_output_dir)) {
    dir.create(current_output_dir, recursive = TRUE)
  }
  
  # Update config paths to use current directory (replace {GENE} placeholders)
  config$paths$output_files$dge_results <- gsub("\\{GENE\\}", primary_gene, config$paths$output_files$dge_results)
  config$paths$output_files$meta_results <- gsub("\\{GENE\\}", primary_gene, config$paths$output_files$meta_results)
  config$paths$reports$analysis_report$output_html <- gsub("\\{GENE\\}", primary_gene, config$paths$reports$analysis_report$output_html)
  
  # Get report configuration (now with dynamic paths)
  report_config <- config$paths$reports$analysis_report
  
  cat("📋 REPORT GENERATION CONFIGURATION:\n")
  cat("🎯 Gene:", primary_gene, "\n")
  cat("🏥 Diseases:", paste(diseases, collapse = ", "), "\n")
  cat("📁 Output Directory:", current_output_dir, "\n")
  cat("📄 Template:", report_config$template %||% "templates/Gene_Analysis_Report.Rmd", "\n")
  cat("📊 Report File:", report_config$output_html, "\n\n")
  
  # Verify required files exist
  template_file <- report_config$template %||% "templates/Gene_Analysis_Report.Rmd"
  
  if (!file.exists(template_file)) {
    return(create_step_result(
      success = FALSE,
      error_message = paste("Report template not found:", template_file),
      step_name = step_name
    ))
  }
  
  # Ensure output directory exists
  # Replace {GENE} placeholder with actual gene name
  output_html_template <- report_config$output_html %||% "reports/{GENE}_Analysis_Report.html"
  output_html <- gsub("\\{GENE\\}", toupper(primary_gene), output_html_template)
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
  cat("  CAMK2D significant:", if (!is.null(meta_summary$camk2d_significant) && meta_summary$camk2d_significant) "YES ✅" else "NO", "\n")
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
  
  # Replace {GENE} placeholder with actual gene name
  primary_gene <- config$research_target$primary_gene %||% config$gene_analysis$primary_gene
  
  # Create required data files from pipeline input data
  if (!is.null(input_data$meta_results) && nrow(input_data$meta_results) > 0) {
    meta_file <- gsub("\\{GENE\\}", toupper(primary_gene), output_files$meta_results)
    cat("📄 Creating meta-analysis results file:", meta_file, "\n")
    write.csv(input_data$meta_results, meta_file, row.names = FALSE)
  } else {
    return(list(success = FALSE, reason = "Meta-analysis results data missing from input"))
  }
  
  # Create DGE results file from individual DGE results data
  dge_data <- NULL
  if (!is.null(input_data$dge_results) && nrow(input_data$dge_results) > 0) {
    dge_data <- input_data$dge_results
  } else if (!is.null(input_data$individual_results) && nrow(input_data$individual_results) > 0) {
    dge_data <- input_data$individual_results
  } else if (!is.null(input_data$dge_summary) && nrow(input_data$dge_summary) > 0) {
    dge_data <- input_data$dge_summary
  }
  
  if (!is.null(dge_data)) {
    dge_file <- gsub("\\{GENE\\}", toupper(primary_gene), output_files$dge_results)
    cat("📄 Creating DGE results file:", dge_file, "\n")
    write.csv(dge_data, dge_file, row.names = FALSE)
  } else {
    return(list(success = FALSE, reason = "DGE results data missing from input (checked: dge_summary, individual_results, dge_results)"))
  }
  
  # Create dataset summary file from available data
  dataset_summary_data <- NULL
  if (!is.null(input_data$dataset_summary) && nrow(input_data$dataset_summary) > 0) {
    dataset_summary_data <- input_data$dataset_summary
  } else if (!is.null(input_data$dge_summary) && nrow(input_data$dge_summary) > 0) {
    dataset_summary_data <- input_data$dge_summary
  } else if (!is.null(input_data$analysis_stats$dataset_summary) && nrow(input_data$analysis_stats$dataset_summary) > 0) {
    dataset_summary_data <- input_data$analysis_stats$dataset_summary
  }
  
  if (!is.null(dataset_summary_data)) {
    # Create both the configured filename and the template-expected filename
    summary_file <- output_files$dataset_summary
    cat("📄 Creating dataset summary file:", summary_file, "\n")
    write.csv(dataset_summary_data, summary_file, row.names = FALSE)
    
    # Also create the file with the name expected by the template
    template_summary_file <- "output/current/dataset_processing_summary_6_datasets.csv"
    write.csv(dataset_summary_data, template_summary_file, row.names = FALSE)
    cat("📄 Created template-compatible dataset summary file:", template_summary_file, "\n")
  } else {
    warnings_generated <- c(warnings_generated, "Dataset summary data missing from input")
  }
  
  # Ensure output/current/ directory exists 
  current_dir <- config$paths$output$current %||% "output/current"
  
  if (!dir.exists(current_dir)) {
    dir.create(current_dir, recursive = TRUE)
  }
  
  cat("✅ All required data files have been created successfully\n")
  
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
    # Convert to absolute paths to avoid working directory issues
    template_file <- normalizePath(template_file, mustWork = TRUE)
    output_dir <- dirname(output_html)
    output_basename <- basename(output_html)
    
    # Ensure output directory exists
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    output_html <- file.path(normalizePath(output_dir, mustWork = TRUE), output_basename)
    
    cat("📍 Working directory:", getwd(), "\n")
    cat("📄 Template (absolute):", template_file, "\n")
    cat("📄 Output (absolute):", output_html, "\n")
    
    # Render the report using absolute paths - no working directory change needed
    cat("\n🔄 Rendering report...\n")
    
    # Suppress warnings during rendering to avoid restart issues
    render_result <- suppressWarnings({
      rmarkdown::render(
        input = template_file,
        output_file = output_html,
        quiet = FALSE  # Show rendering progress
      )
    })
    
    # Verify the output file was created
    if (!file.exists(output_html)) {
      stop("Report file was not created at expected location: ", output_html)
    }
    
    cat("✅ Report rendered successfully to:", output_html, "\n")
    
    return(list(
      success = TRUE,
      output_file = output_html,
      warnings = warnings_generated
    ))
    
  }, error = function(e) {
    return(list(
      success = FALSE,
      reason = paste("Rendering error:", e$message),
      warnings = warnings_generated
    ))
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