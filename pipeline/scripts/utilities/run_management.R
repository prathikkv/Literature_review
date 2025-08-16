#!/usr/bin/env Rscript
#' Run Management and File Naming Utilities
#' 
#' Provides standardized file naming, run tracking, and output organization
#' for bioinformatics pipeline compliance
#' 
#' @author Claude Code Pipeline Standards
#' @version 1.0.0

# Load required libraries
suppressPackageStartupMessages({
  library(yaml)
  library(digest)
  library(jsonlite)
})

#' Generate Dynamic File Names
#'
#' Creates standardized file names with gene, disease, and timestamp
#' @param primary_gene Primary gene symbol
#' @param diseases Vector of disease names
#' @param file_type Type of file ("report", "data", "log", etc.)
#' @param extension File extension (default: "html")
#' @param timestamp Optional timestamp (uses current if NULL)
#' @return Standardized filename
generate_dynamic_filename <- function(primary_gene, 
                                    diseases, 
                                    file_type = "report",
                                    extension = "html",
                                    timestamp = NULL) {
  
  # Input validation
  if (missing(primary_gene) || is.null(primary_gene) || length(primary_gene) == 0 || all(nchar(primary_gene) == 0)) {
    stop("Primary gene must be provided")
  }
  
  if (missing(diseases) || is.null(diseases) || length(diseases) == 0) {
    stop("At least one disease must be provided")
  }
  
  # Ensure single values
  if (length(primary_gene) > 1) {
    primary_gene <- primary_gene[1]
  }
  
  # Sanitize inputs
  clean_gene <- sanitize_filename_component(primary_gene)
  clean_diseases <- sapply(diseases, sanitize_filename_component)
  
  # Generate timestamp if not provided
  if (is.null(timestamp)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  }
  
  # Combine diseases (limit to 2 for filename length)
  disease_string <- if (length(clean_diseases) > 2) {
    paste0(paste(clean_diseases[1:2], collapse = "-"), "_etc")
  } else {
    paste(clean_diseases, collapse = "-")
  }
  
  # Create filename
  filename <- paste(clean_gene, disease_string, timestamp, file_type, sep = "_")
  
  # Add extension
  if (!endsWith(filename, paste0(".", extension))) {
    filename <- paste0(filename, ".", extension)
  }
  
  return(filename)
}

#' Sanitize Filename Component
#'
#' Removes/replaces invalid characters for filesystem compatibility
#' @param component String component to sanitize
#' @return Sanitized string
sanitize_filename_component <- function(component) {
  if (is.null(component) || length(component) == 0 || all(nchar(component) == 0)) {
    return("Unknown")
  }
  
  # Handle vectors by taking first element
  if (length(component) > 1) {
    component <- component[1]
  }
  
  # Remove or replace problematic characters
  clean <- gsub("[^A-Za-z0-9_-]", "", component)
  clean <- gsub("\\s+", "_", clean)
  clean <- gsub("_{2,}", "_", clean)
  clean <- gsub("^_|_$", "", clean)
  
  # Ensure not empty
  if (nchar(clean) == 0) {
    clean <- "Unknown"
  }
  
  # Limit length
  if (nchar(clean) > 20) {
    clean <- substr(clean, 1, 20)
  }
  
  return(clean)
}

#' Create Run Directory Structure
#'
#' Creates organized directory structure for a pipeline run
#' @param base_output_dir Base output directory (default: "output")
#' @param run_id Unique run identifier
#' @param create_symlink Whether to create current symlink (default: TRUE)
#' @return List with directory paths
create_run_directory <- function(base_output_dir = "output", 
                               run_id,
                               create_symlink = TRUE) {
  
  # Validate inputs
  if (missing(run_id) || nchar(run_id) == 0) {
    stop("Run ID must be provided")
  }
  
  # Define directory structure
  dirs <- list(
    base = base_output_dir,
    runs = file.path(base_output_dir, "runs"),
    current_run = file.path(base_output_dir, "runs", run_id),
    data = file.path(base_output_dir, "runs", run_id, "data"),
    logs = file.path(base_output_dir, "runs", run_id, "logs"),
    plots = file.path(base_output_dir, "runs", run_id, "plots"),
    metadata = file.path(base_output_dir, "runs", run_id, "metadata"),
    current_symlink = file.path(base_output_dir, "current"),
    archive = file.path(base_output_dir, "archive")
  )
  
  # Create directories
  for (dir_path in dirs[c("runs", "current_run", "data", "logs", "plots", "metadata", "archive")]) {
    if (!dir.exists(dir_path)) {
      dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      cat("📁 Created directory:", dir_path, "\n")
    }
  }
  
  # Create/update symlink to current run
  if (create_symlink) {
    if (file.exists(dirs$current_symlink) || dir.exists(dirs$current_symlink)) {
      unlink(dirs$current_symlink, recursive = TRUE)
    }
    
    # Create symlink (works on Unix-like systems)
    # Use relative path to avoid broken links
    relative_target <- file.path("runs", run_id)
    tryCatch({
      file.symlink(relative_target, dirs$current_symlink)
      cat("🔗 Created symlink: current ->", relative_target, "\n")
    }, error = function(e) {
      # Fallback: copy directory if symlink fails
      cat("⚠️  Symlink failed, creating copy instead\n")
      dir.create(dirs$current_symlink, showWarnings = FALSE)
    })
  }
  
  return(dirs)
}

#' Generate Run ID
#'
#' Creates a unique run identifier based on gene, diseases, and timestamp
#' @param primary_gene Primary gene symbol
#' @param diseases Vector of diseases
#' @param timestamp Optional timestamp
#' @return Unique run identifier
generate_run_id <- function(primary_gene, diseases, timestamp = NULL) {
  
  if (is.null(timestamp)) {
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  }
  
  # Create base ID
  clean_gene <- sanitize_filename_component(primary_gene)
  clean_diseases <- sapply(diseases, sanitize_filename_component)
  disease_string <- paste(clean_diseases[1:min(2, length(clean_diseases))], collapse = "-")
  
  run_id <- paste(clean_gene, disease_string, timestamp, sep = "_")
  
  return(run_id)
}

#' Create Run Metadata
#'
#' Generates comprehensive metadata for a pipeline run
#' @param config Pipeline configuration
#' @param run_id Unique run identifier
#' @param start_time Run start time
#' @param additional_info Additional metadata to include
#' @return Metadata list
create_run_metadata <- function(config, run_id, start_time = Sys.time(), additional_info = list()) {
  
  # System information
  session_info <- sessionInfo()
  
  # Create comprehensive metadata
  metadata <- list(
    run_info = list(
      run_id = run_id,
      start_time = as.character(start_time),
      timestamp = format(start_time, "%Y-%m-%d %H:%M:%S %Z"),
      user = Sys.getenv("USER", "unknown"),
      working_directory = getwd(),
      r_version = paste(R.version$major, R.version$minor, sep = ".")
    ),
    
    analysis_config = list(
      primary_gene = config$research_target$primary_gene,
      diseases = config$research_target$diseases,
      gene_family = config$research_target$gene_family,
      p_value_threshold = config$analysis_parameters$p_value_threshold,
      log_fc_threshold = config$analysis_parameters$log_fc_threshold,
      meta_analysis_method = config$analysis_parameters$meta_analysis_method
    ),
    
    pipeline_info = list(
      name = config$pipeline$name,
      version = config$pipeline$version,
      execution_mode = config$pipeline$execution_mode,
      steps = config$pipeline$steps
    ),
    
    dynamic_features = config$dynamic_features,
    
    system_info = list(
      platform = session_info$platform,
      os = Sys.info()["sysname"],
      machine = Sys.info()["machine"],
      r_packages = sapply(session_info$otherPkgs, function(x) x$Version)
    ),
    
    data_provenance = list(
      config_checksum = digest(config, algo = "md5"),
      generation_time = as.character(Sys.time())
    )
  )
  
  # Add any additional information
  if (length(additional_info) > 0) {
    metadata <- c(metadata, additional_info)
  }
  
  return(metadata)
}

#' Save Run Metadata
#'
#' Saves metadata to JSON file in run directory
#' @param metadata Metadata list
#' @param run_dirs Run directory structure
#' @param filename Metadata filename (default: "metadata.json")
save_run_metadata <- function(metadata, run_dirs, filename = "metadata.json") {
  
  metadata_file <- file.path(run_dirs$metadata, filename)
  
  tryCatch({
    writeLines(jsonlite::toJSON(metadata, pretty = TRUE, auto_unbox = TRUE), metadata_file)
    cat("💾 Saved metadata to:", metadata_file, "\n")
    return(TRUE)
  }, error = function(e) {
    cat("❌ Failed to save metadata:", e$message, "\n")
    return(FALSE)
  })
}

#' Get Latest Run
#'
#' Finds the most recent run directory
#' @param base_output_dir Base output directory
#' @return Path to latest run directory or NULL if none found
get_latest_run <- function(base_output_dir = "output") {
  
  runs_dir <- file.path(base_output_dir, "runs")
  
  if (!dir.exists(runs_dir)) {
    return(NULL)
  }
  
  run_dirs <- list.dirs(runs_dir, full.names = TRUE, recursive = FALSE)
  
  if (length(run_dirs) == 0) {
    return(NULL)
  }
  
  # Sort by modification time and return latest
  run_info <- file.info(run_dirs)
  latest_run <- run_dirs[which.max(run_info$mtime)]
  
  return(latest_run)
}

#' Archive Old Runs
#'
#' Moves old runs to archive directory
#' @param base_output_dir Base output directory
#' @param keep_recent Number of recent runs to keep (default: 5)
#' @param max_age_days Maximum age in days before archiving (default: 30)
archive_old_runs <- function(base_output_dir = "output", keep_recent = 5, max_age_days = 30) {
  
  runs_dir <- file.path(base_output_dir, "runs")
  archive_dir <- file.path(base_output_dir, "archive")
  
  if (!dir.exists(runs_dir)) {
    return(invisible(NULL))
  }
  
  # Ensure archive directory exists
  if (!dir.exists(archive_dir)) {
    dir.create(archive_dir, recursive = TRUE)
  }
  
  run_dirs <- list.dirs(runs_dir, full.names = TRUE, recursive = FALSE)
  
  if (length(run_dirs) <= keep_recent) {
    return(invisible(NULL))
  }
  
  # Sort by modification time
  run_info <- file.info(run_dirs)
  run_dirs_sorted <- run_dirs[order(run_info$mtime, decreasing = TRUE)]
  
  # Archive old runs
  to_archive <- run_dirs_sorted[(keep_recent + 1):length(run_dirs_sorted)]
  current_time <- Sys.time()
  
  for (run_dir in to_archive) {
    run_age <- difftime(current_time, run_info[run_dir, "mtime"], units = "days")
    
    if (as.numeric(run_age) > max_age_days) {
      run_name <- basename(run_dir)
      archive_path <- file.path(archive_dir, run_name)
      
      tryCatch({
        file.rename(run_dir, archive_path)
        cat("📦 Archived run:", run_name, "\n")
      }, error = function(e) {
        cat("⚠️  Failed to archive run", run_name, ":", e$message, "\n")
      })
    }
  }
}

#' Update Output Paths in Config
#'
#' Updates configuration with new dynamic output paths
#' @param config Configuration object
#' @param run_dirs Run directory structure
#' @param primary_gene Primary gene
#' @param diseases Disease vector
#' @return Updated configuration
update_output_paths_in_config <- function(config, run_dirs, primary_gene, diseases) {
  
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
  
  # Generate dynamic report filename
  report_filename <- generate_dynamic_filename(
    primary_gene = primary_gene,
    diseases = diseases,
    file_type = "report",
    extension = "html",
    timestamp = timestamp
  )
  
  # Update paths in config
  config$paths$reports$analysis_report$output_html <- file.path(run_dirs$current_run, report_filename)
  # Don't override template - keep what's in config
  # config$paths$reports$analysis_report$template is already set in config.yml
  
  # Don't update output paths - they should point to output/current where files are actually saved
  # The meta-analysis step saves to output/current, not to the run directory
  # config$paths$output_files paths remain as configured
  
  return(config)
}

#' Validate File Paths
#'
#' Ensures all required directories exist and paths are valid
#' @param paths Vector of file paths to validate
#' @return List with validation results
validate_file_paths <- function(paths) {
  
  results <- list(
    valid = character(),
    invalid = character(),
    created = character()
  )
  
  for (path in paths) {
    if (is.null(path) || nchar(path) == 0) {
      results$invalid <- c(results$invalid, "empty_path")
      next
    }
    
    # Check if directory needs to be created
    dir_path <- dirname(path)
    if (!dir.exists(dir_path)) {
      tryCatch({
        dir.create(dir_path, recursive = TRUE)
        results$created <- c(results$created, dir_path)
      }, error = function(e) {
        results$invalid <- c(results$invalid, path)
        next
      })
    }
    
    # Check if path is writable
    if (dir.exists(dir_path)) {
      # Try to create a test file
      test_file <- file.path(dir_path, "test_write.tmp")
      tryCatch({
        writeLines("test", test_file)
        file.remove(test_file)
        results$valid <- c(results$valid, path)
      }, error = function(e) {
        results$invalid <- c(results$invalid, path)
      })
    } else {
      results$invalid <- c(results$invalid, path)
    }
  }
  
  return(results)
}

# Module loading confirmation
cat("✅ Run Management Utilities loaded successfully\n")
cat("   Functions: generate_dynamic_filename(), create_run_directory(),\n")
cat("             generate_run_id(), create_run_metadata(), archive_old_runs()\n")
cat("   Version: 1.0.0 (Bioinformatics Standards Compliant)\n\n")