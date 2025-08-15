#!/usr/bin/env Rscript
#' Pipeline Core Functions
#' 
#' Essential functions extracted from the main functions/ directory
#' for self-contained pipeline execution

# ═══════════════════════════════════════════════════════════════
# CAMK GENE DEFINITIONS
# ═══════════════════════════════════════════════════════════════

#' Get CAMK Gene Categories
#' 
#' Returns comprehensive CAMK gene family definitions
#' @return List with core CAMK genes
get_camk_gene_categories <- function() {
  list(
    core = c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMKK1", "CAMKK2", 
             "CAMK1", "CAMK1D", "CAMK1G", "CAMK4", "CAMKV"),
    primary_focus = "CAMK2D",
    camkii_subfamily = c("CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G"),
    camki_subfamily = c("CAMK1", "CAMK1D", "CAMK1G"),
    camkk_subfamily = c("CAMKK1", "CAMKK2"),
    other = c("CAMK4", "CAMKV")
  )
}

#' Get CAMK Core Genes
#' 
#' Simplified function for getting core CAMK genes
#' @return Vector of CAMK gene symbols
get_camk_core_genes <- function() {
  c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMKK1", "CAMKK2", 
    "CAMK1", "CAMK1D", "CAMK1G", "CAMK4", "CAMKV")
}

# ═══════════════════════════════════════════════════════════════
# DATA VALIDATION FUNCTIONS
# ═══════════════════════════════════════════════════════════════

#' Validate Dataset Structure
#' 
#' Checks that dataset has required fields and proper structure
#' @param dataset Dataset object to validate
#' @param dataset_id Dataset identifier for error reporting
#' @return TRUE if valid, FALSE otherwise
validate_dataset_structure <- function(dataset, dataset_id) {
  
  required_fields <- c("success", "expression_matrix", "phenotype_data")
  
  for (field in required_fields) {
    if (is.null(dataset[[field]])) {
      cat("❌ Dataset", dataset_id, "missing required field:", field, "\n")
      return(FALSE)
    }
  }
  
  if (!is.matrix(dataset$expression_matrix)) {
    cat("❌ Dataset", dataset_id, "expression_matrix is not a matrix\n")
    return(FALSE)
  }
  
  if (!is.data.frame(dataset$phenotype_data)) {
    cat("❌ Dataset", dataset_id, "phenotype_data is not a data.frame\n")
    return(FALSE)
  }
  
  if (ncol(dataset$expression_matrix) != nrow(dataset$phenotype_data)) {
    cat("❌ Dataset", dataset_id, "dimension mismatch\n")
    return(FALSE)
  }
  
  return(TRUE)
}

#' Validate Expression Data Quality
#' 
#' Checks expression data for common quality issues
#' @param expr_matrix Expression matrix
#' @param dataset_id Dataset identifier
#' @return List with validation results
validate_expression_quality <- function(expr_matrix, dataset_id) {
  
  issues <- character(0)
  warnings <- character(0)
  
  # Check for missing values
  if (any(is.na(expr_matrix))) {
    na_count <- sum(is.na(expr_matrix))
    na_percent <- na_count / length(expr_matrix) * 100
    if (na_percent > 5) {
      issues <- c(issues, paste("High NA percentage:", round(na_percent, 2), "%"))
    } else {
      warnings <- c(warnings, paste("NA values present:", na_count))
    }
  }
  
  # Check for infinite values
  if (any(is.infinite(expr_matrix))) {
    inf_count <- sum(is.infinite(expr_matrix))
    issues <- c(issues, paste("Infinite values present:", inf_count))
  }
  
  # Check expression range
  expr_range <- range(expr_matrix, na.rm = TRUE)
  if (expr_range[1] < -50 || expr_range[2] > 50) {
    warnings <- c(warnings, paste("Unusual expression range:", 
                                 round(expr_range[1], 2), "to", round(expr_range[2], 2)))
  }
  
  # Check for all-zero genes
  zero_genes <- rowSums(expr_matrix == 0, na.rm = TRUE) == ncol(expr_matrix)
  if (any(zero_genes)) {
    warnings <- c(warnings, paste("All-zero genes detected:", sum(zero_genes)))
  }
  
  return(list(
    valid = length(issues) == 0,
    issues = issues,
    warnings = warnings,
    na_count = sum(is.na(expr_matrix)),
    infinite_count = sum(is.infinite(expr_matrix)),
    expression_range = expr_range
  ))
}

# ═══════════════════════════════════════════════════════════════
# ANALYSIS HELPER FUNCTIONS
# ═══════════════════════════════════════════════════════════════

#' Create Design Matrix
#' 
#' Creates design matrix for limma analysis
#' @param groups_vector Factor vector of sample groups
#' @return Design matrix
create_design_matrix <- function(groups_vector) {
  
  if (!is.factor(groups_vector)) {
    groups_vector <- factor(groups_vector)
  }
  
  # Ensure Control is reference level
  if ("Control" %in% levels(groups_vector)) {
    groups_vector <- relevel(groups_vector, ref = "Control")
  }
  
  design <- model.matrix(~ groups_vector)
  colnames(design) <- c("Intercept", "Disease_vs_Control")
  
  return(design)
}

#' Apply Quality Flags
#' 
#' Assigns quality flags based on effect sizes and other criteria
#' @param results DGE results data frame
#' @param high_logfc_threshold Threshold for high logFC flag
#' @param low_effect_threshold Threshold for low effect flag
#' @return Updated results with Quality_Flag column
apply_quality_flags <- function(results, high_logfc_threshold = 0.8, low_effect_threshold = 0.01) {
  
  results$Quality_Flag <- ifelse(
    abs(results$logFC) > high_logfc_threshold, "HIGH_LOGFC",
    ifelse(abs(results$logFC) < low_effect_threshold, "LOW_EFFECT", "GOOD")
  )
  
  return(results)
}

#' Calculate Standard Errors for Meta-Analysis
#' 
#' Calculates standard errors from t-statistics for meta-analysis
#' @param gene_data Data frame with logFC and t columns
#' @return Updated data frame with SE column
calculate_standard_errors <- function(gene_data) {
  
  gene_data$SE <- abs(gene_data$logFC / gene_data$t)
  
  # Handle edge cases
  gene_data$SE[is.na(gene_data$SE)] <- 1
  gene_data$SE[is.infinite(gene_data$SE)] <- 1
  gene_data$SE[gene_data$SE == 0] <- 0.001
  
  return(gene_data)
}

# ═══════════════════════════════════════════════════════════════
# UTILITY FUNCTIONS
# ═══════════════════════════════════════════════════════════════

#' Check Package Installation
#' 
#' Checks if required packages are installed and installs if missing
#' @param packages Vector of package names
#' @param install_missing Whether to install missing packages
#' @return Vector of missing packages
check_required_packages <- function(packages, install_missing = TRUE) {
  
  missing_packages <- character(0)
  
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
      
      if (install_missing) {
        cat("📦 Installing missing package:", pkg, "\n")
        
        # Try CRAN first
        tryCatch({
          install.packages(pkg, repos = "https://cran.r-project.org")
        }, error = function(e) {
          # Try Bioconductor for bioinformatics packages
          if (pkg %in% c("limma", "DESeq2", "edgeR", "metafor")) {
            tryCatch({
              if (!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager")
              }
              BiocManager::install(pkg)
            }, error = function(e2) {
              cat("❌ Failed to install", pkg, ":", e2$message, "\n")
            })
          } else {
            cat("❌ Failed to install", pkg, ":", e$message, "\n")
          }
        })
      }
    }
  }
  
  return(missing_packages)
}

#' Load Required Libraries
#' 
#' Loads all required libraries with error handling
#' @param packages Vector of package names
#' @return TRUE if all packages loaded successfully
load_required_libraries <- function(packages) {
  
  success_count <- 0
  
  for (pkg in packages) {
    tryCatch({
      suppressPackageStartupMessages(library(pkg, character.only = TRUE))
      success_count <- success_count + 1
    }, error = function(e) {
      cat("⚠️  Warning: Could not load package", pkg, "\n")
    })
  }
  
  if (success_count == length(packages)) {
    cat("✅ All", length(packages), "packages loaded successfully\n")
    return(TRUE)
  } else {
    cat("⚠️ ", (length(packages) - success_count), "packages failed to load\n")
    return(FALSE)
  }
}

#' Format Scientific Numbers
#' 
#' Formats numbers in scientific notation for reports
#' @param x Numeric vector
#' @param digits Number of significant digits
#' @return Character vector of formatted numbers
format_scientific <- function(x, digits = 3) {
  formatC(x, format = "e", digits = digits)
}

#' Format P-values
#' 
#' Formats p-values with appropriate precision
#' @param p_values Numeric vector of p-values
#' @return Character vector of formatted p-values
format_p_values <- function(p_values) {
  ifelse(p_values < 0.001, 
         formatC(p_values, format = "e", digits = 2),
         formatC(p_values, format = "f", digits = 3))
}

#' Create Summary Statistics
#' 
#' Creates summary statistics for a dataset
#' @param dataset Dataset object
#' @return List with summary statistics
create_dataset_summary <- function(dataset) {
  
  if (!dataset$success || is.null(dataset$expression_matrix)) {
    return(list(
      success = FALSE,
      reason = "Invalid dataset"
    ))
  }
  
  expr_matrix <- dataset$expression_matrix
  
  summary_stats <- list(
    success = TRUE,
    n_samples = ncol(expr_matrix),
    n_genes = nrow(expr_matrix),
    platform = dataset$dataset_info$platform %||% "Unknown",
    expression_range = range(expr_matrix, na.rm = TRUE),
    na_count = sum(is.na(expr_matrix)),
    has_phenotype = !is.null(dataset$phenotype_data)
  )
  
  if (!is.null(dataset$phenotype_data)) {
    summary_stats$phenotype_samples <- nrow(dataset$phenotype_data)
    summary_stats$phenotype_columns <- ncol(dataset$phenotype_data)
  }
  
  return(summary_stats)
}

# ═══════════════════════════════════════════════════════════════
# FILE I/O FUNCTIONS
# ═══════════════════════════════════════════════════════════════

#' Safe File Write
#' 
#' Writes files with error handling and backup
#' @param data Data to write
#' @param file_path Output file path
#' @param backup Whether to create backup if file exists
#' @return TRUE if successful
safe_file_write <- function(data, file_path, backup = TRUE) {
  
  # Create directory if it doesn't exist
  dir_path <- dirname(file_path)
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  
  # Create backup if file exists
  if (backup && file.exists(file_path)) {
    backup_path <- paste0(file_path, ".backup.", format(Sys.time(), "%Y%m%d_%H%M%S"))
    file.copy(file_path, backup_path)
  }
  
  # Write file
  tryCatch({
    if (is.data.frame(data)) {
      write.csv(data, file_path, row.names = FALSE)
    } else {
      saveRDS(data, file_path)
    }
    cat("💾 File saved:", file_path, "\n")
    return(TRUE)
  }, error = function(e) {
    cat("❌ Failed to save file:", file_path, "-", e$message, "\n")
    return(FALSE)
  })
}

#' Safe File Read
#' 
#' Reads files with error handling
#' @param file_path File path to read
#' @param file_type Type of file ("csv", "rds", "auto")
#' @return Data or NULL if failed
safe_file_read <- function(file_path, file_type = "auto") {
  
  if (!file.exists(file_path)) {
    cat("❌ File not found:", file_path, "\n")
    return(NULL)
  }
  
  # Determine file type
  if (file_type == "auto") {
    file_ext <- tools::file_ext(file_path)
    file_type <- switch(tolower(file_ext),
                       "csv" = "csv",
                       "rds" = "rds",
                       "csv")  # Default to csv
  }
  
  tryCatch({
    if (file_type == "csv") {
      data <- read.csv(file_path, stringsAsFactors = FALSE)
    } else if (file_type == "rds") {
      data <- readRDS(file_path)
    } else {
      stop("Unsupported file type:", file_type)
    }
    
    cat("📄 File loaded:", file_path, "\n")
    return(data)
  }, error = function(e) {
    cat("❌ Failed to read file:", file_path, "-", e$message, "\n")
    return(NULL)
  })
}

# NULL coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

cat("✅ Pipeline core functions loaded successfully\n")
cat("Available functions: CAMK gene definitions, data validation, analysis helpers, utilities\n\n")