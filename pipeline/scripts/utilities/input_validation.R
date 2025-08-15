#!/usr/bin/env Rscript
#' Input Validation and Security Module
#' 
#' Comprehensive input validation for bioinformatics pipeline security
#' Prevents code injection, validates biological data, ensures data integrity
#' 
#' @author Claude Code Security Module
#' @version 1.0.0

# Load required libraries
suppressPackageStartupMessages({
  library(stringr)
})

#' Validate Gene Symbol
#'
#' Validates gene symbols against HGNC standards and common patterns
#' @param gene_symbol Gene symbol to validate
#' @param allow_aliases Whether to allow common aliases (default: TRUE)
#' @return List with validation result and cleaned gene symbol
validate_gene_symbol <- function(gene_symbol, allow_aliases = TRUE) {
  
  result <- list(
    valid = FALSE,
    cleaned = NA,
    warnings = character(),
    errors = character()
  )
  
  # Basic input checks
  if (is.null(gene_symbol) || is.na(gene_symbol) || nchar(trimws(gene_symbol)) == 0) {
    result$errors <- c(result$errors, "Gene symbol cannot be empty")
    return(result)
  }
  
  # Clean input
  clean_gene <- trimws(toupper(gene_symbol))
  
  # Remove dangerous characters
  if (grepl("[;<>&|`$(){}\\[\\]\"'\\\\]", clean_gene)) {
    result$errors <- c(result$errors, "Gene symbol contains invalid characters")
    return(result)
  }
  
  # Check length constraints
  if (nchar(clean_gene) > 20) {
    result$errors <- c(result$errors, "Gene symbol too long (max 20 characters)")
    return(result)
  }
  
  if (nchar(clean_gene) < 2) {
    result$errors <- c(result$errors, "Gene symbol too short (min 2 characters)")
    return(result)
  }
  
  # HGNC pattern validation
  # Most human genes follow pattern: Letters + optional numbers
  if (!grepl("^[A-Z][A-Z0-9]*$", clean_gene)) {
    result$warnings <- c(result$warnings, "Gene symbol doesn't match standard HGNC pattern")
  }
  
  # Check for common invalid patterns
  invalid_patterns <- c(
    "^[0-9]+$",     # All numbers
    "^TEST$",       # Test strings
    "^NULL$",       # Null values
    "^NA$",         # Missing values
    "^UNKNOWN$"     # Unknown values
  )
  
  for (pattern in invalid_patterns) {
    if (grepl(pattern, clean_gene)) {
      result$errors <- c(result$errors, paste("Invalid gene symbol pattern:", clean_gene))
      return(result)
    }
  }
  
  # Validation passed
  result$valid <- TRUE
  result$cleaned <- clean_gene
  
  # Add informational warnings for unusual but valid symbols
  if (nchar(clean_gene) > 10) {
    result$warnings <- c(result$warnings, "Unusually long gene symbol")
  }
  
  if (grepl("[0-9]", clean_gene) && !grepl("^[A-Z]+[0-9]+[A-Z]*$", clean_gene)) {
    result$warnings <- c(result$warnings, "Unusual number placement in gene symbol")
  }
  
  return(result)
}

#' Validate Disease Names
#'
#' Validates disease names for safety and standardization
#' @param diseases Vector of disease names
#' @param max_diseases Maximum number of diseases allowed (default: 5)
#' @return List with validation results
validate_disease_names <- function(diseases, max_diseases = 5) {
  
  result <- list(
    valid = FALSE,
    cleaned = character(),
    warnings = character(),
    errors = character()
  )
  
  # Basic input checks
  if (is.null(diseases) || length(diseases) == 0) {
    result$errors <- c(result$errors, "At least one disease must be provided")
    return(result)
  }
  
  if (length(diseases) > max_diseases) {
    result$errors <- c(result$errors, paste("Too many diseases (max", max_diseases, ")"))
    return(result)
  }
  
  cleaned_diseases <- character()
  
  for (i in seq_along(diseases)) {
    disease <- diseases[i]
    
    # Basic checks
    if (is.null(disease) || is.na(disease) || nchar(trimws(disease)) == 0) {
      result$errors <- c(result$errors, paste("Disease", i, "is empty"))
      next
    }
    
    # Clean input
    clean_disease <- trimws(disease)
    
    # Security checks - remove dangerous characters
    if (grepl("[;<>&|`$(){}\\[\\]\"'\\\\]", clean_disease)) {
      result$errors <- c(result$errors, paste("Disease name contains invalid characters:", clean_disease))
      next
    }
    
    # Length checks
    if (nchar(clean_disease) > 100) {
      result$errors <- c(result$errors, paste("Disease name too long:", clean_disease))
      next
    }
    
    if (nchar(clean_disease) < 3) {
      result$errors <- c(result$errors, paste("Disease name too short:", clean_disease))
      next
    }
    
    # Pattern validation - allow letters, numbers, spaces, hyphens
    if (!grepl("^[A-Za-z0-9\\s\\-\\.]+$", clean_disease)) {
      result$warnings <- c(result$warnings, paste("Unusual characters in disease name:", clean_disease))
    }
    
    # Standardize capitalization
    clean_disease <- tools::toTitleCase(tolower(clean_disease))
    
    cleaned_diseases <- c(cleaned_diseases, clean_disease)
  }
  
  if (length(result$errors) == 0) {
    result$valid <- TRUE
    result$cleaned <- cleaned_diseases
  }
  
  return(result)
}

#' Validate File Paths
#'
#' Validates file paths for security and accessibility
#' @param paths Vector of file paths to validate
#' @param check_write_access Whether to check write permissions (default: TRUE)
#' @param allow_creation Whether to allow directory creation (default: TRUE)
#' @return List with validation results
validate_file_paths <- function(paths, check_write_access = TRUE, allow_creation = TRUE) {
  
  result <- list(
    valid = TRUE,
    secure_paths = character(),
    warnings = character(),
    errors = character(),
    created_dirs = character()
  )
  
  if (length(paths) == 0) {
    result$errors <- c(result$errors, "No paths provided")
    result$valid <- FALSE
    return(result)
  }
  
  for (path in paths) {
    # Basic checks
    if (is.null(path) || is.na(path) || nchar(trimws(path)) == 0) {
      result$errors <- c(result$errors, "Empty path provided")
      result$valid <- FALSE
      next
    }
    
    clean_path <- trimws(path)
    
    # Security checks - prevent path traversal attacks
    if (grepl("\\.\\.", clean_path)) {
      result$errors <- c(result$errors, paste("Path traversal detected:", clean_path))
      result$valid <- FALSE
      next
    }
    
    # Check for dangerous characters
    if (grepl("[;<>&|`$\\\\]", clean_path)) {
      result$errors <- c(result$errors, paste("Dangerous characters in path:", clean_path))
      result$valid <- FALSE
      next
    }
    
    # Ensure absolute paths or relative to working directory
    if (!grepl("^(/|[A-Za-z]:|\\./|output/|cache/|temp/)", clean_path)) {
      result$warnings <- c(result$warnings, paste("Unusual path format:", clean_path))
    }
    
    # Check directory existence and creation
    dir_path <- dirname(clean_path)
    
    if (!dir.exists(dir_path)) {
      if (allow_creation) {
        tryCatch({
          dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
          result$created_dirs <- c(result$created_dirs, dir_path)
        }, error = function(e) {
          result$errors <- c(result$errors, paste("Cannot create directory:", dir_path))
          result$valid <- FALSE
        })
      } else {
        result$errors <- c(result$errors, paste("Directory does not exist:", dir_path))
        result$valid <- FALSE
        next
      }
    }
    
    # Check write access if requested
    if (check_write_access && dir.exists(dir_path)) {
      test_file <- file.path(dir_path, paste0("test_", Sys.getpid(), ".tmp"))
      tryCatch({
        writeLines("test", test_file)
        file.remove(test_file)
      }, error = function(e) {
        result$errors <- c(result$errors, paste("No write access:", dir_path))
        result$valid <- FALSE
      })
    }
    
    if (result$valid) {
      result$secure_paths <- c(result$secure_paths, clean_path)
    }
  }
  
  return(result)
}

#' Validate Configuration Parameters
#'
#' Validates pipeline configuration for security and biological validity
#' @param config Configuration object
#' @return List with validation results
validate_configuration <- function(config) {
  
  result <- list(
    valid = TRUE,
    warnings = character(),
    errors = character(),
    security_issues = character()
  )
  
  # Validate research target
  if (!is.null(config$research_target)) {
    # Validate primary gene
    if (!is.null(config$research_target$primary_gene)) {
      gene_validation <- validate_gene_symbol(config$research_target$primary_gene)
      if (!gene_validation$valid) {
        result$errors <- c(result$errors, paste("Invalid primary gene:", paste(gene_validation$errors, collapse = "; ")))
        result$valid <- FALSE
      }
      result$warnings <- c(result$warnings, gene_validation$warnings)
    }
    
    # Validate diseases
    if (!is.null(config$research_target$diseases)) {
      disease_validation <- validate_disease_names(config$research_target$diseases)
      if (!disease_validation$valid) {
        result$errors <- c(result$errors, paste("Invalid diseases:", paste(disease_validation$errors, collapse = "; ")))
        result$valid <- FALSE
      }
      result$warnings <- c(result$warnings, disease_validation$warnings)
    }
    
    # Validate gene family
    if (!is.null(config$research_target$gene_family) && length(config$research_target$gene_family) > 0) {
      for (gene in config$research_target$gene_family) {
        gene_validation <- validate_gene_symbol(gene)
        if (!gene_validation$valid) {
          result$warnings <- c(result$warnings, paste("Invalid gene family member:", gene))
        }
      }
    }
  }
  
  # Validate analysis parameters
  if (!is.null(config$analysis_parameters)) {
    params <- config$analysis_parameters
    
    # P-value threshold
    if (!is.null(params$p_value_threshold)) {
      if (!is.numeric(params$p_value_threshold) || params$p_value_threshold <= 0 || params$p_value_threshold > 1) {
        result$errors <- c(result$errors, "P-value threshold must be between 0 and 1")
        result$valid <- FALSE
      }
    }
    
    # Log FC threshold
    if (!is.null(params$log_fc_threshold)) {
      if (!is.numeric(params$log_fc_threshold) || params$log_fc_threshold < 0) {
        result$errors <- c(result$errors, "Log FC threshold must be non-negative")
        result$valid <- FALSE
      }
    }
    
    # Sample size
    if (!is.null(params$min_sample_size)) {
      if (!is.numeric(params$min_sample_size) || params$min_sample_size < 5 || params$min_sample_size > 10000) {
        result$errors <- c(result$errors, "Minimum sample size must be between 5 and 10000")
        result$valid <- FALSE
      }
    }
  }
  
  # Check for potential security issues in paths
  if (!is.null(config$paths)) {
    all_paths <- unlist(config$paths, recursive = TRUE)
    if (length(all_paths) > 0) {
      path_validation <- validate_file_paths(all_paths)
      if (length(path_validation$invalid) > 0) {
        result$security_issues <- c(result$security_issues, 
                                   paste("Invalid paths:", paste(path_validation$invalid, collapse = ", ")))
        result$valid <- FALSE
      }
    }
  }
  
  return(result)
}

#' Sanitize User Input
#'
#' General purpose input sanitization function
#' @param input Input string to sanitize
#' @param max_length Maximum allowed length (default: 200)
#' @param allow_special_chars Whether to allow special characters (default: FALSE)
#' @return Sanitized string
sanitize_user_input <- function(input, max_length = 200, allow_special_chars = FALSE) {
  
  if (is.null(input) || is.na(input)) {
    return("")
  }
  
  # Convert to string and trim
  clean_input <- trimws(as.character(input))
  
  # Length check
  if (nchar(clean_input) > max_length) {
    clean_input <- substr(clean_input, 1, max_length)
  }
  
  # Remove dangerous characters unless explicitly allowed
  if (!allow_special_chars) {
    clean_input <- gsub("[;<>&|`$(){}\\[\\]\"'\\\\]", "", clean_input)
  }
  
  return(clean_input)
}

#' Check Data Integrity
#'
#' Validates data integrity using checksums
#' @param file_path Path to file to check
#' @param expected_checksum Expected MD5 checksum (optional)
#' @return List with integrity check results
check_data_integrity <- function(file_path, expected_checksum = NULL) {
  
  result <- list(
    valid = FALSE,
    checksum = NA,
    file_exists = FALSE,
    file_readable = FALSE,
    warnings = character(),
    errors = character()
  )
  
  # Check file existence
  if (!file.exists(file_path)) {
    result$errors <- c(result$errors, paste("File does not exist:", file_path))
    return(result)
  }
  
  result$file_exists <- TRUE
  
  # Check readability
  tryCatch({
    file_info <- file.info(file_path)
    if (file_info$size == 0) {
      result$warnings <- c(result$warnings, "File is empty")
    }
    result$file_readable <- TRUE
  }, error = function(e) {
    result$errors <- c(result$errors, paste("Cannot read file:", e$message))
    return(result)
  })
  
  # Calculate checksum
  tryCatch({
    if (requireNamespace("digest", quietly = TRUE)) {
      result$checksum <- digest::digest(file_path, file = TRUE, algo = "md5")
    } else {
      result$warnings <- c(result$warnings, "digest package not available for checksum")
    }
  }, error = function(e) {
    result$warnings <- c(result$warnings, paste("Checksum calculation failed:", e$message))
  })
  
  # Compare with expected checksum if provided
  if (!is.null(expected_checksum) && !is.na(result$checksum)) {
    if (result$checksum == expected_checksum) {
      result$valid <- TRUE
    } else {
      result$errors <- c(result$errors, "Checksum mismatch - data integrity compromised")
    }
  } else {
    result$valid <- result$file_readable
  }
  
  return(result)
}

#' Log Security Event
#'
#' Logs security-related events for audit trail
#' @param event_type Type of security event
#' @param details Event details
#' @param severity Severity level ("INFO", "WARNING", "ERROR")
#' @param log_file Log file path (default: "output/logs/security.log")
log_security_event <- function(event_type, details, severity = "INFO", log_file = "output/logs/security.log") {
  
  tryCatch({
    # Ensure log directory exists
    log_dir <- dirname(log_file)
    if (!dir.exists(log_dir)) {
      dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Create log entry
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    log_entry <- paste(timestamp, severity, event_type, details, sep = " | ")
    
    # Append to log file
    cat(log_entry, "\n", file = log_file, append = TRUE)
    
  }, error = function(e) {
    # Fall back to console if logging fails
    cat("SECURITY LOG:", event_type, "-", details, "\n")
  })
}

# Module loading confirmation
cat("✅ Input Validation and Security Module loaded successfully\n")
cat("   Functions: validate_gene_symbol(), validate_disease_names(),\n")
cat("             validate_file_paths(), validate_configuration(),\n")
cat("             sanitize_user_input(), check_data_integrity()\n")
cat("   Version: 1.0.0 (Security Compliant)\n\n")