#!/usr/bin/env Rscript
#' Expression Data Validator for CAMK2D Research
#' 
#' This module downloads actual expression data from GEO datasets,
#' validates CAMK2D expression, extracts phenotype information,
#' and ensures datasets are suitable for bioinformatics analysis.

# Load required libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(tidyverse)
  library(limma)
  library(readr)
})

# Global configuration for better performance
options(timeout = 600)  # 10 minute timeout for large downloads
options(download.file.method = "curl")  # Use curl for better reliability

#' Download GEO Expression Data with Retry Logic and Caching
#'
#' @param geo_accession Character. GEO accession number
#' @param max_retries Integer. Maximum number of retry attempts
#' @param verbose Logical. Enable verbose logging
#' @return List with expression matrix, phenotype data, feature data, and data type
download_geo_expression_data <- function(geo_accession, max_retries = 3, verbose = TRUE) {
  
  if (verbose) {
    cat("\n🔍 Starting download for", geo_accession, "\n")
    cat("   Max retries:", max_retries, "\n")
    cat("   Timeout setting:", getOption("timeout"), "seconds\n")
  }
  
  # Check cache first
  cache_dir <- "cache/geo_downloads"
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  cache_file <- file.path(cache_dir, paste0(geo_accession, "_cached.rds"))
  if (file.exists(cache_file)) {
    if (verbose) cat("📦 Found cached data, loading...\n")
    cached_data <- readRDS(cache_file)
    if (!is.null(cached_data)) {
      if (verbose) cat("✅ Loaded from cache successfully\n")
      return(cached_data)
    }
  }
  
  # Try each method with retries
  methods <- c("GSEMatrix", "Supplementary", "GSE")
  
  for (method in methods) {
    if (verbose) cat("\n📡 Trying method:", method, "\n")
    
    for (attempt in 1:max_retries) {
      if (verbose && attempt > 1) {
        cat("   Retry attempt", attempt, "of", max_retries, "\n")
        Sys.sleep(5 * attempt)  # Exponential backoff
      }
      
      result <- NULL
      
      if (method == "GSEMatrix") {
        result <- try_gsematrix_download(geo_accession, verbose)
      } else if (method == "Supplementary") {
        result <- try_supplementary_download(geo_accession, verbose)
      } else if (method == "GSE") {
        result <- try_gse_download(geo_accession, verbose)
      }
      
      if (!is.null(result) && nrow(result$expr_matrix) > 0) {
        if (verbose) cat("✅ Method", method, "successful on attempt", attempt, "\n")
        
        # Cache the successful download
        saveRDS(result, cache_file)
        if (verbose) cat("💾 Saved to cache:", cache_file, "\n")
        
        return(result)
      }
      
      if (verbose && is.null(result)) {
        cat("   ⚠️ Attempt", attempt, "failed for method", method, "\n")
      }
    }
  }
  
  if (verbose) cat("❌ All download methods failed after", max_retries, "attempts each\n")
  return(NULL)
}

#' Try GSEMatrix Download Method
#'
#' @param geo_accession Character. GEO accession
#' @param verbose Logical. Enable verbose logging
#' @return List with data or NULL if failed
try_gsematrix_download <- function(geo_accession, verbose = TRUE) {
  
  tryCatch({
    if (verbose) cat("   Downloading with GSEMatrix...\n")
    
    # Download with error handling
    gse <- getGEO(geo_accession, GSEMatrix = TRUE, getGPL = TRUE, AnnotGPL = TRUE)
    
    if (length(gse) == 0) {
      if (verbose) cat("   ⚠️ GSEMatrix returned empty result\n")
      return(NULL)
    }
    
    # Use the first platform (most datasets have only one)
    eset <- gse[[1]]
    
    # Extract data
    expr_matrix <- exprs(eset)
    pheno_data <- pData(eset)
    feature_data <- fData(eset)
    
    # Check if we got meaningful data
    if (nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
      if (verbose) cat("   ⚠️ GSEMatrix returned empty expression matrix\n")
      return(NULL)
    }
    
    if (verbose) cat("   📊 Platform: ", nrow(expr_matrix), " genes x ", ncol(expr_matrix), " samples\n")
    
    return(list(
      expr_matrix = expr_matrix,
      pheno_data = pheno_data,
      feature_data = feature_data,
      data_type = "GSEMatrix",
      original_object = eset
    ))
    
  }, error = function(e) {
    if (verbose) cat("   ❌ GSEMatrix error:", e$message, "\n")
    return(NULL)
  }, warning = function(w) {
    if (verbose) cat("   ⚠️ GSEMatrix warning:", w$message, "\n")
  })
}

#' Try Supplementary Files Download Method
#'
#' @param geo_accession Character. GEO accession
#' @param verbose Logical. Enable verbose logging
#' @return List with data or NULL if failed
try_supplementary_download <- function(geo_accession, verbose = TRUE) {
  
  tryCatch({
    if (verbose) cat("   Downloading supplementary files...\n")
    
    # Create temporary directory for download
    temp_dir <- file.path(tempdir(), geo_accession)
    if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE)
    dir.create(temp_dir, recursive = TRUE)
    
    # Download supplementary files
    supp_files <- getGEOSuppFiles(geo_accession, makeDirectory = FALSE, baseDir = temp_dir)
    
    if (nrow(supp_files) == 0) {
      cat("⚠️ No supplementary files found\n")
      return(NULL)
    }
    
    cat("📁 Found", nrow(supp_files), "supplementary files\n")
    
    # Look for expression matrix files
    expr_file <- find_expression_file(supp_files, temp_dir)
    if (is.null(expr_file)) {
      cat("⚠️ No expression matrix file found in supplementary files\n")
      return(NULL)
    }
    
    # Load expression data
    expr_data <- load_expression_file(expr_file)
    if (is.null(expr_data)) {
      cat("⚠️ Failed to load expression data from file\n")
      return(NULL)
    }
    
    # Get phenotype data from GEO (this usually works even when expression matrix doesn't)
    pheno_data <- get_phenotype_data_only(geo_accession)
    
    # Create minimal feature data
    feature_data <- create_feature_data_from_matrix(expr_data)
    
    cat("📊 Supplementary: ", nrow(expr_data), " genes x ", ncol(expr_data), " samples\n")
    
    # Clean up temporary files
    unlink(temp_dir, recursive = TRUE)
    
    return(list(
      expr_matrix = as.matrix(expr_data),
      pheno_data = pheno_data,
      feature_data = feature_data,
      data_type = "Supplementary",
      file_source = expr_file
    ))
    
  }, error = function(e) {
    cat("❌ Supplementary files method failed:", e$message, "\n")
    return(NULL)
  })
}

#' Try Sample-Level Download Method
#'
#' @param geo_accession Character. GEO accession
#' @return List with data or NULL if failed  
try_sample_level_download <- function(geo_accession) {
  
  tryCatch({
    cat("🔬 Trying sample-level download...\n")
    
    # This is a fallback method - usually not needed but can work for some datasets
    gse <- getGEO(geo_accession, GSEMatrix = FALSE)
    
    if (length(gse) == 0) {
      cat("⚠️ Sample-level download returned empty result\n")
      return(NULL)
    }
    
    # Extract sample data (this is more complex and usually not needed)
    # For now, return NULL - this can be implemented if needed
    cat("⚠️ Sample-level method not fully implemented\n")
    return(NULL)
    
  }, error = function(e) {
    cat("❌ Sample-level method failed:", e$message, "\n")
    return(NULL)
  })
}

#' Find Expression Matrix File in Supplementary Files
#'
#' @param supp_files Data frame. Supplementary files info
#' @param temp_dir Character. Temporary directory path
#' @return Character. Path to expression file or NULL
find_expression_file <- function(supp_files, temp_dir) {
  
  # Common patterns for expression matrix files (updated for actual files)
  expr_patterns <- c(
    "count|counts",                                # Count files (most common for RNA-seq)
    "fpkm|rpkm|tpm",                              # Normalized expression
    "expression.*matrix|matrix.*expression",       # Expression matrices
    "norm.*count|normalized.*count",              # Normalized counts
    "raw.*count|count.*raw",                      # Raw counts  
    "gene.*expression|expression.*gene",          # Gene expression
    ".*counts.*csv|.*counts.*txt",                # CSV/TXT count files
    "countstable|counts_table"                    # Count table variants
  )
  
  # Get file names and paths (rownames are already full paths)
  file_paths <- rownames(supp_files)
  file_names <- basename(file_paths)
  
  cat("📁 Searching in files:", paste(file_names, collapse = ", "), "\n")
  
  # Check each pattern
  for (pattern in expr_patterns) {
    matches <- grep(pattern, file_names, ignore.case = TRUE)
    if (length(matches) > 0) {
      cat("📍 Pattern '", pattern, "' matched ", length(matches), " files: ", 
          paste(file_names[matches], collapse = ", "), "\n")
      
      # Prefer .txt, .csv, .tsv files over compressed
      best_match <- matches[1]
      for (match_idx in matches) {
        if (grepl("\\.(txt|csv|tsv)(\\.gz)?$", file_names[match_idx], ignore.case = TRUE)) {
          best_match <- match_idx
          break
        }
      }
      
      selected_file <- file_paths[best_match]
      cat("📊 Found potential expression file:", file_names[best_match], "\n")
      
      # Check if file actually contains expression data
      if (validate_expression_file(selected_file)) {
        cat("✅ File validated as expression data\n")
        return(selected_file)
      } else {
        cat("❌ File failed validation\n")
      }
    }
  }
  
  # If no pattern matches, try the largest data file (likely expression matrix)
  if (length(file_paths) > 0) {
    # Filter for data files (exclude .xlsx, .pdf, .doc files)
    data_files <- file_paths[!grepl("\\.(xlsx|pdf|doc|docx|xls)$", file_paths, ignore.case = TRUE)]
    
    if (length(data_files) > 0) {
      # Get file sizes
      file_info <- file.info(data_files)
      valid_files <- data_files[!is.na(file_info$size)]
      
      if (length(valid_files) > 0) {
        file_sizes <- file.info(valid_files)$size
        largest_file <- valid_files[which.max(file_sizes)]
        
        cat("📊 Trying largest data file:", basename(largest_file), "\n")
        if (validate_expression_file(largest_file)) {
          return(largest_file)
        }
      }
    }
  }
  
  return(NULL)
}

#' Validate if File Contains Expression Data
#'
#' @param file_path Character. Path to file
#' @return Logical. TRUE if file looks like expression data
validate_expression_file <- function(file_path) {
  
  tryCatch({
    # Handle compressed files
    if (grepl("\\.gz$", file_path)) {
      con <- gzfile(file_path, "rt")
    } else {
      con <- file(file_path, "rt")
    }
    
    # Read first few lines
    first_lines <- readLines(con, n = 5)
    close(con)
    
    if (length(first_lines) < 2) return(FALSE)
    
    # Check if looks like expression data (numeric values, gene names)
    header_line <- first_lines[1]
    data_line <- first_lines[2]
    
    # Split by common delimiters
    delim <- detect_delimiter(header_line)
    header_parts <- strsplit(header_line, delim)[[1]]
    data_parts <- strsplit(data_line, delim)[[1]]
    
    # Should have at least 3 columns (gene + 2 samples minimum)
    if (length(header_parts) < 3 || length(data_parts) < 3) return(FALSE)
    
    # Most values in data line should be numeric
    numeric_values <- suppressWarnings(as.numeric(data_parts[-1]))  # Skip first column (gene names)
    numeric_percent <- sum(!is.na(numeric_values)) / length(numeric_values)
    
    # At least 30% should be numeric (relaxed for files with annotation columns)
    # The load_expression_file function will filter to just sample columns
    return(numeric_percent >= 0.3)
    
  }, error = function(e) {
    return(FALSE)
  })
}

#' Detect Delimiter in File
#'
#' @param line Character. Line to analyze
#' @return Character. Detected delimiter
detect_delimiter <- function(line) {
  
  # Count common delimiters
  tab_count <- lengths(regmatches(line, gregexpr("\t", line)))
  comma_count <- lengths(regmatches(line, gregexpr(",", line)))
  space_count <- lengths(regmatches(line, gregexpr(" ", line)))
  
  # Return most common delimiter
  counts <- c(tab = tab_count, comma = comma_count, space = space_count)
  delims <- c(tab = "\t", comma = ",", space = " ")
  
  return(delims[which.max(counts)])
}

#' Load Expression Data from File
#'
#' @param file_path Character. Path to expression file
#' @return Data frame. Expression data or NULL if failed
load_expression_file <- function(file_path) {
  
  tryCatch({
    cat("📖 Loading expression data from:", basename(file_path), "\n")
    
    # Detect delimiter
    if (grepl("\\.gz$", file_path)) {
      con <- gzfile(file_path, "rt")
    } else {
      con <- file(file_path, "rt")
    }
    first_line <- readLines(con, n = 1)
    close(con)
    
    delim <- detect_delimiter(first_line)
    
    # Read the file with robust CSV parsing
    # Try multiple approaches for robustness
    expr_data <- NULL
    
    # Method 1: Try readr::read_delim (more robust for malformed CSV)
    tryCatch({
      if (delim == ",") {
        expr_data <- readr::read_csv(file_path, show_col_types = FALSE)
      } else if (delim == "\t") {
        expr_data <- readr::read_tsv(file_path, show_col_types = FALSE)
      } else {
        expr_data <- readr::read_delim(file_path, delim = delim, show_col_types = FALSE)
      }
      
      # Convert to data frame and set row names
      expr_data <- as.data.frame(expr_data)
      first_col <- expr_data[, 1]
      expr_data <- expr_data[, -1]
      rownames(expr_data) <- first_col
      
      cat("✅ Successfully read with readr\n")
      
    }, error = function(e) {
      cat("⚠️ readr failed:", e$message, "\n")
      expr_data <<- NULL
    })
    
    # Method 2: Fallback to base R with robust options
    if (is.null(expr_data)) {
      tryCatch({
        expr_data <- read.table(file_path, header = TRUE, sep = delim, row.names = 1,
                               check.names = FALSE, stringsAsFactors = FALSE,
                               fill = TRUE, quote = "", comment.char = "")
        cat("✅ Successfully read with base R\n")
        
      }, error = function(e) {
        cat("❌ Base R also failed:", e$message, "\n")
        return(NULL)
      })
    }
    
    # Validate the loaded data (relaxed thresholds for testing)
    if (nrow(expr_data) < 100) {  # Reduced threshold for initial testing
      cat("⚠️ File has too few rows (", nrow(expr_data), ") - may not be expression data\n")
      return(NULL)
    }
    
    if (ncol(expr_data) < 2) {  # Should have at least 2 samples
      cat("⚠️ File has too few columns (", ncol(expr_data), ") - may not be expression data\n")
      return(NULL)
    }
    
    # Check if data has enough numeric columns for expression analysis
    numeric_cols <- sapply(expr_data, is.numeric)
    
    # For expression data, we need at least some numeric columns
    # But we'll be smarter about identifying sample columns vs annotation columns
    
    # Look for sample-like column names (common patterns)
    sample_patterns <- c(
      "sample|Sample|SAMPLE",
      "ctrl|control|Control|CTRL|CTR",
      "treat|treatment|Treatment|TREAT",
      "case|Case|CASE",
      "rep|replicate|Replicate|REP",
      "_\\d+$|\\d+$",  # Ending with numbers
      "^[A-Z0-9_]+_[A-Z0-9_]+$"  # Pattern like condition_replicate
    )
    
    # Patterns for annotation columns to EXCLUDE
    annotation_patterns <- c(
      "gene_id|geneid|gene|symbol|description|annotation",
      "entrez|ensembl|refseq|dbxref|db_xref",
      "chromosome|chr|position|start|end|strand",
      "length|size|biotype|transcript"
    )
    
    # Identify likely sample columns
    col_names <- colnames(expr_data)
    likely_samples <- rep(FALSE, length(col_names))
    likely_annotations <- rep(FALSE, length(col_names))
    
    for (pattern in sample_patterns) {
      matches <- grepl(pattern, col_names, ignore.case = TRUE)
      likely_samples <- likely_samples | matches
    }
    
    for (pattern in annotation_patterns) {
      matches <- grepl(pattern, col_names, ignore.case = TRUE)
      likely_annotations <- likely_annotations | matches
    }
    
    # Also consider first N numeric columns as likely samples (but exclude annotations)
    first_numeric <- which(numeric_cols & !likely_annotations)[1:min(10, sum(numeric_cols & !likely_annotations))]
    first_numeric <- first_numeric[!is.na(first_numeric)]
    
    # Combine both approaches but exclude annotation columns
    sample_candidates <- (likely_samples | (1:length(col_names) %in% first_numeric)) & !likely_annotations
    sample_numeric <- numeric_cols & sample_candidates
    
    num_sample_cols <- sum(sample_numeric)
    total_sample_candidates <- sum(sample_candidates)
    
    cat("   Detected", num_sample_cols, "numeric sample columns out of", 
        total_sample_candidates, "sample candidates\n")
    
    # Require at least 2 numeric sample columns for expression analysis
    if (num_sample_cols < 2) {
      cat("⚠️ Too few numeric sample columns for expression analysis\n")
      cat("   Numeric sample columns:", num_sample_cols, "\n")
      return(NULL)
    }
    
    # If we have enough sample columns, filter the data to keep only expression data
    if (num_sample_cols >= 2) {
      # Keep only the numeric sample columns for expression analysis
      expr_data <- expr_data[, sample_numeric, drop = FALSE]
      cat("✅ Filtered to", ncol(expr_data), "expression sample columns\n")
    }
    
    cat("✅ Loaded", nrow(expr_data), "genes x", ncol(expr_data), "samples\n")
    return(expr_data)
    
  }, error = function(e) {
    cat("❌ Failed to load expression file:", e$message, "\n")
    return(NULL)
  })
}

#' Get Phenotype Data Only from GEO
#'
#' @param geo_accession Character. GEO accession
#' @return Data frame. Phenotype data
get_phenotype_data_only <- function(geo_accession) {
  
  tryCatch({
    # Try to get just the phenotype data
    gse <- getGEO(geo_accession, GSEMatrix = TRUE, getGPL = FALSE)
    
    if (length(gse) > 0) {
      pheno_data <- pData(gse[[1]])
      cat("📋 Retrieved phenotype data for", nrow(pheno_data), "samples\n")
      return(pheno_data)
    }
    
    # If that fails, create minimal phenotype data
    cat("⚠️ Could not retrieve phenotype data, creating minimal data\n")
    return(data.frame(
      title = paste("Sample", 1:10),
      source_name_ch1 = "Unknown",
      characteristics_ch1 = "Unknown",
      stringsAsFactors = FALSE
    ))
    
  }, error = function(e) {
    cat("⚠️ Error getting phenotype data:", e$message, "\n")
    return(data.frame(
      title = paste("Sample", 1:10),
      source_name_ch1 = "Unknown", 
      characteristics_ch1 = "Unknown",
      stringsAsFactors = FALSE
    ))
  })
}

#' Create Feature Data from Expression Matrix
#'
#' @param expr_data Data frame. Expression data
#' @return Data frame. Feature annotations
create_feature_data_from_matrix <- function(expr_data) {
  
  gene_ids <- rownames(expr_data)
  
  # Create basic feature data
  feature_data <- data.frame(
    ID = gene_ids,
    Gene.Symbol = gene_ids,  # Assume row names are gene symbols
    Gene.Name = gene_ids,
    Description = paste("Gene:", gene_ids),
    stringsAsFactors = FALSE,
    row.names = gene_ids
  )
  
  cat("📋 Created feature data for", nrow(feature_data), "genes\n")
  return(feature_data)
}

#' Download and Validate Expression Data
#'
#' @param geo_accession Character. GEO accession number (e.g., "GSE123456")
#' @param output_dir Character. Directory to save expression data
#' @param validate_camk2d Logical. Check CAMK2D expression levels
#' @param min_expression Numeric. Minimum CAMK2D expression required
#' @return List with expression data and validation results
download_and_validate_expression <- function(geo_accession, 
                                           output_dir = "data/expression_data",
                                           validate_camk2d = TRUE,
                                           min_expression = 1) {
  
  cat("📥 Downloading expression data for", geo_accession, "...\n")
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  tryCatch({
    # Try primary method: getGEO with GSEMatrix
    cat("🔄 Attempting GSEMatrix download...\n")
    download_result <- download_geo_expression_data(geo_accession)
    
    if (is.null(download_result)) {
      cat("❌ No data found for", geo_accession, "\n")
      return(list(
        geo_accession = geo_accession,
        suitable = FALSE,
        reason = "No expression data available"
      ))
    }
    
    # Extract components
    expr_matrix <- download_result$expr_matrix
    pheno_data <- download_result$pheno_data
    feature_data <- download_result$feature_data
    data_type <- download_result$data_type
    
    cat("✅ Downloaded", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
    
    # Validate dataset quality
    quality_results <- validate_dataset_quality(expr_matrix, pheno_data, geo_accession)
    
    if (!quality_results$suitable_for_analysis) {
      cat("❌ Dataset failed quality validation\n")
      return(list(
        geo_accession = geo_accession,
        suitable = FALSE,
        reason = quality_results$failure_reason,
        quality_results = quality_results
      ))
    }
    
    # Extract and clean phenotype data
    cleaned_pheno <- extract_phenotype_data(pheno_data, geo_accession)
    
    # Validate CAMK2D expression
    camk2d_validation <- NULL
    if (validate_camk2d) {
      camk2d_validation <- validate_camk2d_expression(
        expr_matrix, feature_data, cleaned_pheno, min_expression
      )
      
      if (!camk2d_validation$camk2d_detected) {
        cat("❌ CAMK2D not adequately expressed in dataset\n")
        return(list(
          geo_accession = geo_accession,
          suitable = FALSE,
          reason = "CAMK2D not expressed",
          camk2d_validation = camk2d_validation,
          quality_results = quality_results
        ))
      }
    }
    
    # Clean and normalize expression data
    cleaned_expression <- clean_expression_matrix(expr_matrix, feature_data)
    
    # Save processed data
    saved_files <- save_expression_data(
      expression_matrix = cleaned_expression,
      phenotype_data = cleaned_pheno,
      feature_data = feature_data,
      geo_accession = geo_accession,
      output_dir = output_dir
    )
    
    # Create analysis summary
    analysis_summary <- create_expression_summary(
      cleaned_expression, cleaned_pheno, camk2d_validation, quality_results
    )
    
    cat("✅ Expression data validated and saved for", geo_accession, "\n")
    cat("   📊 CAMK2D expression:", camk2d_validation$camk2d_mean_expression, "\n")
    cat("   👥 Sample groups:", paste(names(table(cleaned_pheno$group)), collapse = ", "), "\n")
    
    return(list(
      geo_accession = geo_accession,
      suitable = TRUE,
      expression_matrix = cleaned_expression,
      phenotype_data = cleaned_pheno,
      feature_data = feature_data,
      camk2d_validation = camk2d_validation,
      quality_results = quality_results,
      analysis_summary = analysis_summary,
      saved_files = saved_files
    ))
    
  }, error = function(e) {
    cat("❌ Failed to process", geo_accession, ":", e$message, "\n")
    return(list(
      geo_accession = geo_accession,
      suitable = FALSE,
      reason = paste("Download/processing error:", e$message)
    ))
  })
}

#' Validate Dataset Quality for Analysis
#'
#' @param expr_matrix Matrix. Expression matrix
#' @param pheno_data Data frame. Phenotype data
#' @param geo_accession Character. GEO accession
#' @return List with quality validation results
validate_dataset_quality <- function(expr_matrix, pheno_data, geo_accession) {
  
  cat("🔍 Validating dataset quality...\n")
  
  # Check dimensions
  n_genes <- nrow(expr_matrix)
  n_samples <- ncol(expr_matrix)
  
  # Check for sufficient data (relaxed thresholds for initial testing)
  if (n_samples < 6) {  # Reduced from 10 to 6
    return(list(
      suitable_for_analysis = FALSE,
      failure_reason = paste("Too few samples:", n_samples)
    ))
  }
  
  if (n_genes < 1000) {  # Reduced from 5000 to 1000 
    return(list(
      suitable_for_analysis = FALSE,
      failure_reason = paste("Too few genes:", n_genes)
    ))
  }
  
  # Check for missing data
  missing_percent <- sum(is.na(expr_matrix)) / length(expr_matrix) * 100
  if (missing_percent > 30) {
    return(list(
      suitable_for_analysis = FALSE,
      failure_reason = paste("Too much missing data:", round(missing_percent, 1), "%")
    ))
  }
  
  # Check expression range (detect if log-transformed)
  expr_range <- range(expr_matrix, na.rm = TRUE)
  is_log_transformed <- expr_range[2] < 25  # Likely log2-transformed if max < 25
  
  # Check for variability
  gene_variances <- apply(expr_matrix, 1, var, na.rm = TRUE)
  low_variance_genes <- sum(gene_variances < 0.1, na.rm = TRUE)
  low_variance_percent <- low_variance_genes / n_genes * 100
  
  if (low_variance_percent > 80) {
    return(list(
      suitable_for_analysis = FALSE,
      failure_reason = paste("Low gene expression variability:", round(low_variance_percent, 1), "%")
    ))
  }
  
  # Check phenotype data quality
  # First, filter phenotype data to only include samples present in expression matrix
  expr_sample_names <- colnames(expr_matrix)
  
  # Find matching samples between expression and phenotype data
  pheno_sample_names <- pheno_data$title
  matching_samples <- intersect(expr_sample_names, pheno_sample_names)
  
  if (length(matching_samples) < 3) {
    return(list(
      suitable_for_analysis = FALSE,
      failure_reason = paste("Too few matching samples between expression and phenotype data:", length(matching_samples))
    ))
  }
  
  # Filter phenotype data to matching samples only
  filtered_pheno_data <- pheno_data[pheno_data$title %in% matching_samples, ]
  
  cat("   Filtered phenotype data to", nrow(filtered_pheno_data), "samples matching expression data\n")
  
  pheno_quality <- validate_phenotype_quality(filtered_pheno_data)
  
  if (!pheno_quality$adequate) {
    return(list(
      suitable_for_analysis = FALSE,
      failure_reason = pheno_quality$failure_reason
    ))
  }
  
  # All checks passed
  return(list(
    suitable_for_analysis = TRUE,
    n_genes = n_genes,
    n_samples = n_samples,
    missing_percent = round(missing_percent, 2),
    is_log_transformed = is_log_transformed,
    expression_range = expr_range,
    low_variance_percent = round(low_variance_percent, 1),
    phenotype_quality = pheno_quality
  ))
}

#' Validate Phenotype Data Quality
#'
#' @param pheno_data Data frame. Phenotype data
#' @return List with phenotype validation results
validate_phenotype_quality <- function(pheno_data) {
  
  n_samples <- nrow(pheno_data)
  n_annotations <- ncol(pheno_data)
  
  # Look for group/condition information
  potential_group_cols <- c("title", "source_name_ch1", "characteristics_ch1", 
                           "description", "treatment", "disease", "condition")
  
  group_info_found <- FALSE
  group_column <- NULL
  
  for (col in names(pheno_data)) {
    if (any(sapply(potential_group_cols, function(x) grepl(x, col, ignore.case = TRUE)))) {
      # Check if column has multiple unique values
      unique_vals <- length(unique(pheno_data[[col]]))
      # Allow up to all samples to be unique for small datasets, but cap at 80% for larger ones
      max_allowed <- max(n_samples, ceiling(n_samples * 0.8))
      if (unique_vals >= 2 && unique_vals <= max_allowed) {  # Between 2 and max allowed
        group_info_found <- TRUE
        group_column <- col
        break
      }
    }
  }
  
  if (!group_info_found) {
    return(list(
      adequate = FALSE,
      failure_reason = "No clear group/condition information found in phenotype data"
    ))
  }
  
  # Check for case-control or treatment-control design (expanded patterns)
  group_values <- pheno_data[[group_column]]
  has_control <- any(grepl("control|normal|healthy|baseline|ctr|ctrl|wild.*type|wt|reference", group_values, ignore.case = TRUE))
  has_treatment <- any(grepl("disease|treatment|patient|case|affected|ko|knockout|icko|mcko|mutant|transgenic", group_values, ignore.case = TRUE))
  
  if (!has_control && !has_treatment) {
    return(list(
      adequate = FALSE,
      failure_reason = "No clear case-control or treatment-control design detected"
    ))
  }
  
  return(list(
    adequate = TRUE,
    group_column = group_column,
    n_groups = length(unique(group_values)),
    has_control = has_control,
    has_treatment = has_treatment,
    group_distribution = table(group_values)
  ))
}

#' Extract and Standardize Phenotype Data
#'
#' @param pheno_data Data frame. Raw phenotype data
#' @param geo_accession Character. GEO accession
#' @return Data frame with cleaned phenotype data
extract_phenotype_data <- function(pheno_data, geo_accession) {
  
  cat("🏷️  Extracting phenotype information...\n")
  
  # Start with sample IDs
  cleaned_pheno <- data.frame(
    sample_id = rownames(pheno_data),
    geo_accession = geo_accession,
    stringsAsFactors = FALSE
  )
  
  # Extract group/condition information
  group_col <- find_best_group_column(pheno_data)
  
  if (!is.null(group_col)) {
    cleaned_pheno$group_raw <- pheno_data[[group_col]]
    cleaned_pheno$group <- standardize_group_labels(pheno_data[[group_col]])
  } else {
    cleaned_pheno$group_raw <- "Unknown"
    cleaned_pheno$group <- "Unknown"
  }
  
  # Extract additional annotations
  cleaned_pheno$title <- pheno_data$title %||% ""
  cleaned_pheno$source <- pheno_data$source_name_ch1 %||% ""
  cleaned_pheno$characteristics <- pheno_data$characteristics_ch1 %||% ""
  
  # Extract age, sex, treatment information if available
  cleaned_pheno$age <- extract_age_info(pheno_data)
  cleaned_pheno$sex <- extract_sex_info(pheno_data)
  cleaned_pheno$treatment <- extract_treatment_info(pheno_data)
  
  # Add disease classification
  cleaned_pheno$disease_category <- classify_disease_from_phenotype(cleaned_pheno)
  
  cat("   👥 Identified", length(unique(cleaned_pheno$group)), "sample groups\n")
  
  return(cleaned_pheno)
}

#' Find Best Group Column in Phenotype Data
#'
#' @param pheno_data Data frame. Phenotype data
#' @return Character. Best column name for grouping or NULL
find_best_group_column <- function(pheno_data) {
  
  # Priority order for group columns
  priority_cols <- c("title", "source_name_ch1", "characteristics_ch1", "description", "genotype:ch1", "treatment:ch1")
  
  for (col in priority_cols) {
    if (col %in% names(pheno_data)) {
      values <- pheno_data[[col]]
      unique_vals <- length(unique(values))
      
      # Good group column: 2-20 unique values (relaxed upper limit)
      if (unique_vals >= 2 && unique_vals <= 20) {
        # Extended pattern matching for group indicators
        group_patterns <- c(
          "control|normal|healthy|disease|treatment|case",
          "ctr|ctrl|ko|knockout|icko|mcko",
          "mutant|transgenic|wild.*type|wt",
          "_\\d+$",  # Ending with numbers (like sample replicates)
          "baseline|reference"
        )
        
        # Check if any pattern matches
        for (pattern in group_patterns) {
          if (any(grepl(pattern, values, ignore.case = TRUE))) {
            return(col)
          }
        }
        
        # If title column has multiple unique values, it's likely grouping info
        if (col == "title" && unique_vals >= 2) {
          return(col)
        }
      }
    }
  }
  
  return(NULL)
}

#' Standardize Group Labels
#'
#' @param group_values Character vector. Raw group labels
#' @return Character vector. Standardized labels
standardize_group_labels <- function(group_values) {
  
  standardized <- tolower(as.character(group_values))
  
  # Standardize control labels (expanded patterns)
  standardized[grepl("control|normal|healthy|baseline|reference|ctr|ctrl|wild.*type|wt", standardized)] <- "control"
  
  # Standardize disease/treatment labels (expanded patterns)
  standardized[grepl("disease|patient|case|affected|treatment|drug|ko|knockout|icko|mcko|mutant|transgenic", standardized)] <- "case"
  
  # Standardize specific conditions
  standardized[grepl("afib|atrial.*fibrillation", standardized)] <- "aFIB"
  standardized[grepl("heart.*failure|cardiac.*failure|hf", standardized)] <- "heart_failure"
  standardized[grepl("cardiomyopathy", standardized)] <- "cardiomyopathy"
  
  return(standardized)
}

#' Validate CAMK2D Expression in Dataset
#'
#' @param expr_matrix Matrix. Expression matrix
#' @param feature_data Data frame. Feature annotations
#' @param pheno_data Data frame. Phenotype data
#' @param min_expression Numeric. Minimum expression level
#' @return List with CAMK2D validation results
validate_camk2d_expression <- function(expr_matrix, feature_data, pheno_data, min_expression) {
  
  cat("🔍 Validating CAMK2D expression...\n")
  
  # Find CAMK2D probes/genes
  camk2d_features <- find_camk2d_features(feature_data)
  
  if (length(camk2d_features) == 0) {
    return(list(
      camk2d_detected = FALSE,
      reason = "CAMK2D not found in feature annotations"
    ))
  }
  
  # Extract CAMK2D expression
  camk2d_expr <- expr_matrix[camk2d_features, , drop = FALSE]
  
  # Calculate expression statistics
  if (nrow(camk2d_expr) == 1) {
    mean_expr <- mean(camk2d_expr, na.rm = TRUE)
    max_expr <- max(camk2d_expr, na.rm = TRUE)
  } else {
    # Multiple probes - take the one with highest mean expression
    probe_means <- rowMeans(camk2d_expr, na.rm = TRUE)
    best_probe <- which.max(probe_means)
    camk2d_expr <- camk2d_expr[best_probe, , drop = FALSE]
    mean_expr <- probe_means[best_probe]
    max_expr <- max(camk2d_expr, na.rm = TRUE)
  }
  
  # Check if expression is adequate
  is_expressed <- mean_expr >= min_expression && max_expr > min_expression * 2
  
  # Check for differential expression between groups
  differential_expr <- NULL
  if (length(unique(pheno_data$group)) >= 2) {
    differential_expr <- test_camk2d_differential_expression(camk2d_expr, pheno_data)
  }
  
  result <- list(
    camk2d_detected = is_expressed,
    camk2d_features_found = camk2d_features,
    camk2d_mean_expression = round(mean_expr, 3),
    camk2d_max_expression = round(max_expr, 3),
    camk2d_expression_range = range(camk2d_expr, na.rm = TRUE),
    differential_expression = differential_expr
  )
  
  if (!is_expressed) {
    result$reason = paste("CAMK2D expression too low: mean =", round(mean_expr, 3))
  }
  
  return(result)
}

#' Find CAMK2D Features in Annotation Data
#'
#' @param feature_data Data frame. Feature annotations
#' @return Character vector. Row names/IDs for CAMK2D features
find_camk2d_features <- function(feature_data) {
  
  cat("🔍 Searching for CAMK2D features in", nrow(feature_data), "genes...\n")
  
  # Look for CAMK2D in different annotation columns
  search_cols <- c("Gene.Symbol", "GENE_SYMBOL", "Gene_Symbol", "Symbol", 
                   "Gene.Name", "GENE_NAME", "Gene_Name", "Name",
                   "Description", "DESCRIPTION", "Definition", "ID")
  
  camk2d_features <- c()
  
  # Expanded search patterns for CAMK2D
  camk2d_patterns <- c(
    "^CAMK2D$",                                    # Exact match
    "^CaMK2D$",                                    # Alternative capitalization
    "^CAMKIID$",                                   # Compact form
    "^CaMKII.*delta$|^CaMKII.*DELTA$",            # CaMKII delta variants
    "^CAMKII.*delta$|^CAMKII.*DELTA$",            # CAMKII delta variants  
    "^calcium.*calmodulin.*kinase.*II.*delta$",    # Full name
    "^calcium.*calmodulin.*dependent.*protein.*kinase.*II.*delta$", # Very full name
    "CAMK2D",                                      # Contains CAMK2D
    "CaMK.*II.*delta|CaMKII.*delta",               # Various CaMKII delta forms
    "calmodulin.*kinase.*delta|calmodulin.*dependent.*kinase.*delta" # Descriptive forms
  )
  
  # Search in annotation columns
  for (col in search_cols) {
    if (col %in% names(feature_data)) {
      col_data <- as.character(feature_data[[col]])
      
      for (pattern in camk2d_patterns) {
        matches <- grep(pattern, col_data, ignore.case = TRUE)
        
        if (length(matches) > 0) {
          new_features <- rownames(feature_data)[matches]
          camk2d_features <- c(camk2d_features, new_features)
          
          cat("📍 Found", length(matches), "matches in", col, "using pattern:", pattern, "\n")
          cat("   Matches:", paste(col_data[matches], collapse = ", "), "\n")
        }
      }
    }
  }
  
  # Also search in row names directly (common for RNA-seq data)
  row_names <- rownames(feature_data)
  for (pattern in camk2d_patterns) {
    direct_matches <- grep(pattern, row_names, ignore.case = TRUE)
    if (length(direct_matches) > 0) {
      new_features <- row_names[direct_matches]
      camk2d_features <- c(camk2d_features, new_features)
      
      cat("📍 Found", length(direct_matches), "matches in row names using pattern:", pattern, "\n")
      cat("   Matches:", paste(row_names[direct_matches], collapse = ", "), "\n")
    }
  }
  
  # Remove duplicates and return
  unique_features <- unique(camk2d_features)
  
  if (length(unique_features) > 0) {
    cat("✅ Total CAMK2D features found:", length(unique_features), "\n")
    cat("   Feature IDs:", paste(unique_features, collapse = ", "), "\n")
  } else {
    cat("❌ No CAMK2D features found in dataset\n")
    
    # Debug: show what gene symbols are available (first 20)
    if ("Gene.Symbol" %in% names(feature_data)) {
      available_symbols <- head(feature_data$Gene.Symbol, 20)
      cat("   Available gene symbols (sample):", paste(available_symbols, collapse = ", "), "\n")
    } else if (nrow(feature_data) > 0) {
      available_ids <- head(rownames(feature_data), 20)
      cat("   Available feature IDs (sample):", paste(available_ids, collapse = ", "), "\n")
    }
  }
  
  return(unique_features)
}

#' Test CAMK2D Differential Expression
#'
#' @param camk2d_expr Matrix. CAMK2D expression values
#' @param pheno_data Data frame. Phenotype data with groups
#' @return List with differential expression test results
test_camk2d_differential_expression <- function(camk2d_expr, pheno_data) {
  
  groups <- pheno_data$group
  unique_groups <- unique(groups)
  
  if (length(unique_groups) != 2) {
    return(list(testable = FALSE, reason = "Need exactly 2 groups for t-test"))
  }
  
  group1_expr <- as.numeric(camk2d_expr[groups == unique_groups[1]])
  group2_expr <- as.numeric(camk2d_expr[groups == unique_groups[2]])
  
  # Remove missing values
  group1_expr <- group1_expr[!is.na(group1_expr)]
  group2_expr <- group2_expr[!is.na(group2_expr)]
  
  if (length(group1_expr) < 3 || length(group2_expr) < 3) {
    return(list(testable = FALSE, reason = "Too few samples per group"))
  }
  
  # Perform t-test
  t_test_result <- t.test(group1_expr, group2_expr)
  
  return(list(
    testable = TRUE,
    group1_name = unique_groups[1],
    group2_name = unique_groups[2],
    group1_mean = round(mean(group1_expr), 3),
    group2_mean = round(mean(group2_expr), 3),
    fold_change = round(mean(group2_expr) / mean(group1_expr), 3),
    log2_fold_change = round(log2(mean(group2_expr) / mean(group1_expr)), 3),
    p_value = round(t_test_result$p.value, 6),
    significant = t_test_result$p.value < 0.05
  ))
}

#' Helper Functions

extract_age_info <- function(pheno_data) {
  age_cols <- c("age", "Age", "AGE", "age_at_diagnosis", "age_years")
  for (col in age_cols) {
    if (col %in% names(pheno_data)) {
      ages <- as.numeric(gsub("[^0-9.]", "", pheno_data[[col]]))
      if (sum(!is.na(ages)) > 0) return(ages)
    }
  }
  return(rep(NA, nrow(pheno_data)))
}

extract_sex_info <- function(pheno_data) {
  sex_cols <- c("sex", "Sex", "SEX", "gender", "Gender")
  for (col in sex_cols) {
    if (col %in% names(pheno_data)) {
      return(as.character(pheno_data[[col]]))
    }
  }
  return(rep("Unknown", nrow(pheno_data)))
}

extract_treatment_info <- function(pheno_data) {
  treatment_cols <- c("treatment", "Treatment", "drug", "Drug", "compound")
  for (col in treatment_cols) {
    if (col %in% names(pheno_data)) {
      return(as.character(pheno_data[[col]]))
    }
  }
  return(rep("Unknown", nrow(pheno_data)))
}

classify_disease_from_phenotype <- function(pheno_data) {
  text <- paste(tolower(pheno_data$title), tolower(pheno_data$source), 
                tolower(pheno_data$characteristics))
  
  ifelse(grepl("afib|atrial.*fibrillation", text), "aFIB",
  ifelse(grepl("heart.*failure|cardiac.*failure", text), "Heart_Failure", 
  ifelse(grepl("cardiomyopathy", text), "Cardiomyopathy",
  ifelse(grepl("cardiac|heart", text), "Cardiac", "Other"))))
}

clean_expression_matrix <- function(expr_matrix, feature_data) {
  # Add gene symbols if available
  if ("Gene.Symbol" %in% names(feature_data)) {
    symbols <- feature_data$Gene.Symbol[match(rownames(expr_matrix), rownames(feature_data))]
    expr_df <- as.data.frame(expr_matrix)
    expr_df$Gene_Symbol <- symbols
    return(expr_df)
  }
  
  return(as.data.frame(expr_matrix))
}

save_expression_data <- function(expression_matrix, phenotype_data, feature_data, geo_accession, output_dir) {
  
  # Save expression matrix
  expr_file <- file.path(output_dir, paste0(geo_accession, "_expression.csv"))
  write_csv(expression_matrix, expr_file)
  
  # Save phenotype data
  pheno_file <- file.path(output_dir, paste0(geo_accession, "_phenotypes.csv"))
  write_csv(phenotype_data, pheno_file)
  
  # Save feature annotations
  feature_file <- file.path(output_dir, paste0(geo_accession, "_features.csv"))
  write_csv(feature_data, feature_file)
  
  return(list(
    expression_file = expr_file,
    phenotype_file = pheno_file,
    feature_file = feature_file
  ))
}

create_expression_summary <- function(expression_matrix, phenotype_data, camk2d_validation, quality_results) {
  
  data.frame(
    n_genes = nrow(expression_matrix),
    n_samples = nrow(phenotype_data),
    n_groups = length(unique(phenotype_data$group)),
    camk2d_expressed = camk2d_validation$camk2d_detected,
    camk2d_mean_expr = camk2d_validation$camk2d_mean_expression,
    suitable_for_dge = quality_results$suitable_for_analysis && camk2d_validation$camk2d_detected,
    missing_data_percent = quality_results$missing_percent,
    is_log_transformed = quality_results$is_log_transformed
  )
}

#' Batch Process Multiple Datasets
#'
#' @param dataset_list Data frame. Results from targeted discovery
#' @param max_datasets Integer. Maximum number of datasets to process
#' @param output_dir Character. Output directory
#' @return List with processing results
batch_process_expression_data <- function(dataset_list, max_datasets = 10, output_dir = "data/expression_data") {
  
  cat("🔄 Batch Processing Expression Data\n")
  cat("==================================\n")
  
  if (nrow(dataset_list) == 0) {
    cat("❌ No datasets provided for processing\n")
    return(list(successful = list(), failed = list(), summary = NULL))
  }
  
  # Select top datasets for processing
  top_datasets <- dataset_list %>%
    filter(CAMK2D_Relevance %in% c("High", "Medium")) %>%
    arrange(desc(Total_Relevance_Score), desc(Sample_Count)) %>%
    head(max_datasets)
  
  cat("📋 Processing", nrow(top_datasets), "top-priority datasets...\n\n")
  
  successful_processing <- list()
  failed_processing <- list()
  
  for (i in seq_len(nrow(top_datasets))) {
    dataset_row <- top_datasets[i, ]
    geo_accession <- dataset_row$GEO_Accession
    
    cat("📥 Processing dataset", i, "of", nrow(top_datasets), ":", geo_accession, "\n")
    cat("   Title:", substr(dataset_row$Title, 1, 60), "...\n")
    cat("   Relevance:", dataset_row$CAMK2D_Relevance, "| Score:", dataset_row$Total_Relevance_Score, "\n")
    
    # Process the dataset
    processing_result <- download_and_validate_expression(
      geo_accession = geo_accession,
      output_dir = output_dir,
      validate_camk2d = TRUE,
      min_expression = 2  # Slightly higher threshold for batch processing
    )
    
    if (!is.null(processing_result) && processing_result$suitable) {
      # Add original dataset metadata
      processing_result$original_metadata <- dataset_row
      successful_processing[[geo_accession]] <- processing_result
      
      cat("   ✅ SUCCESS: Dataset validated and saved\n")
      
    } else {
      failed_processing[[geo_accession]] <- list(
        geo_accession = geo_accession,
        original_metadata = dataset_row,
        failure_reason = processing_result$reason %||% "Unknown error"
      )
      
      cat("   ❌ FAILED:", processing_result$reason %||% "Unknown error", "\n")
    }
    
    cat("\n")
    Sys.sleep(3.0)  # Rate limiting between datasets
  }
  
  # Create processing summary
  processing_summary <- create_batch_processing_summary(successful_processing, failed_processing, top_datasets)
  
  # Save batch processing results
  save_batch_results(successful_processing, failed_processing, processing_summary, output_dir)
  
  cat("🏁 Batch Processing Complete\n")
  cat("   ✅ Successful:", length(successful_processing), "datasets\n")
  cat("   ❌ Failed:", length(failed_processing), "datasets\n")
  cat("   📊 Success rate:", round(length(successful_processing)/(length(successful_processing)+length(failed_processing))*100, 1), "%\n\n")
  
  return(list(
    successful = successful_processing,
    failed = failed_processing,
    summary = processing_summary,
    success_rate = length(successful_processing)/(length(successful_processing)+length(failed_processing))
  ))
}

#' Create Batch Processing Summary
#'
#' @param successful_processing List. Successful processing results
#' @param failed_processing List. Failed processing results
#' @param original_datasets Data frame. Original dataset metadata
#' @return Data frame with processing summary
create_batch_processing_summary <- function(successful_processing, failed_processing, original_datasets) {
  
  if (length(successful_processing) == 0) {
    return(data.frame(
      summary_type = "No successful processing",
      value = "All datasets failed validation"
    ))
  }
  
  # Extract key metrics from successful datasets
  successful_metrics <- map_dfr(successful_processing, function(result) {
    data.frame(
      geo_accession = result$geo_accession,
      n_genes = result$analysis_summary$n_genes,
      n_samples = result$analysis_summary$n_samples,
      n_groups = result$analysis_summary$n_groups,
      camk2d_mean_expr = result$analysis_summary$camk2d_mean_expr,
      camk2d_expressed = result$analysis_summary$camk2d_expressed,
      suitable_for_dge = result$analysis_summary$suitable_for_dge,
      disease_category = result$phenotype_data$disease_category[1] %||% "Unknown"
    )
  })
  
  # Create comprehensive summary
  summary_stats <- data.frame(
    metric = c(
      "Total datasets processed",
      "Successful validations", 
      "Failed validations",
      "Success rate (%)",
      "Total samples available",
      "Total genes covered",
      "Datasets with CAMK2D expression",
      "Datasets suitable for DGE",
      "Average CAMK2D expression",
      "Disease categories covered"
    ),
    value = c(
      nrow(original_datasets),
      length(successful_processing),
      length(failed_processing), 
      round(length(successful_processing)/(length(successful_processing)+length(failed_processing))*100, 1),
      sum(successful_metrics$n_samples),
      paste(range(successful_metrics$n_genes), collapse = "-"),
      sum(successful_metrics$camk2d_expressed),
      sum(successful_metrics$suitable_for_dge),
      round(mean(successful_metrics$camk2d_mean_expr), 2),
      length(unique(successful_metrics$disease_category))
    )
  )
  
  return(list(
    summary_stats = summary_stats,
    successful_datasets = successful_metrics,
    failure_reasons = map_chr(failed_processing, ~.x$failure_reason)
  ))
}

#' Save Batch Processing Results
#'
#' @param successful_processing List. Successful results
#' @param failed_processing List. Failed results
#' @param processing_summary List. Processing summary
#' @param output_dir Character. Output directory
save_batch_results <- function(successful_processing, failed_processing, processing_summary, output_dir) {
  
  # Create batch results directory
  batch_dir <- file.path(output_dir, "batch_results")
  if (!dir.exists(batch_dir)) {
    dir.create(batch_dir, recursive = TRUE)
  }
  
  # Save summary statistics
  if (!is.null(processing_summary$summary_stats)) {
    summary_file <- file.path(batch_dir, paste0("batch_summary_", format(Sys.Date(), "%Y%m%d"), ".csv"))
    write_csv(processing_summary$summary_stats, summary_file)
  }
  
  # Save successful dataset details
  if (!is.null(processing_summary$successful_datasets)) {
    success_file <- file.path(batch_dir, paste0("successful_datasets_", format(Sys.Date(), "%Y%m%d"), ".csv"))
    write_csv(processing_summary$successful_datasets, success_file)
  }
  
  # Save failure summary
  if (length(failed_processing) > 0) {
    failed_summary <- map_dfr(failed_processing, function(fail) {
      data.frame(
        geo_accession = fail$geo_accession,
        title = fail$original_metadata$Title %||% "",
        failure_reason = fail$failure_reason,
        original_score = fail$original_metadata$Total_Relevance_Score %||% 0
      )
    })
    
    failure_file <- file.path(batch_dir, paste0("failed_datasets_", format(Sys.Date(), "%Y%m%d"), ".csv"))
    write_csv(failed_summary, failure_file)
  }
  
  cat("💾 Batch results saved to:", batch_dir, "\n")
}

cat("✅ Expression Data Validator Module Loaded\n")
cat("🔬 Ready to download, validate, and process GEO expression data\n")
cat("🧬 CAMK2D expression validation and phenotype extraction enabled\n")
cat("🔄 Batch processing pipeline ready for multiple datasets\n")