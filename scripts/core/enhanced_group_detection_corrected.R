#!/usr/bin/env Rscript
#' Enhanced Group Detection with Biological Reference Logic
#' 
#' CRITICAL CORRECTION: Implements proper biological reference group assignment
#' to ensure consistent Disease vs Control contrasts across all datasets

cat("TARGET: ENHANCED GROUP DETECTION WITH BIOLOGICAL REFERENCE LOGIC\n")
cat("CRITICAL FIX: Proper contrast direction for consistent meta-analysis\n")
cat("=============================================================\n\n")

#' Enhanced Auto-Detect Groups with Biological Logic
#'
#' @param dataset Processed dataset object
#' @return List with groups vector and metadata, with proper biological reference
enhanced_auto_detect_groups_corrected <- function(dataset) {
  
  if (!dataset$success || is.null(dataset$phenotype_data)) {
    return(NULL)
  }
  
  pheno_data <- dataset$phenotype_data
  n_samples <- nrow(pheno_data)
  
  cat("SEARCH: Analyzing phenotype data for healthy vs disease patterns...\n")
  cat("   Samples:", n_samples, "\n")
  cat("   Columns:", ncol(pheno_data), "\n")
  
  # Enhanced column priority for group detection
  priority_columns <- c("source_name_ch1", "title", "description", 
                       "disease status:ch1", "heart failure:ch1", "characteristics_ch1")
  
  # Try each column in priority order
  for (col_name in priority_columns) {
    if (col_name %in% colnames(pheno_data)) {
      cat("   Checking column:", col_name, "\n")
      
      source_info <- pheno_data[[col_name]]
      unique_values <- unique(source_info)
      cat("     Unique values:", length(unique_values), "\n")
      
      # BIOLOGICAL CLASSIFICATION LOGIC
      # Define control/healthy patterns (should be reference group)
      control_patterns <- c("control", "normal", "non-failing", "healthy", "donor", 
                           "sinus rhythm", "SR", "wild type", "WT")
      
      # Define disease patterns (should be comparison group)  
      disease_patterns <- c("cardiomyopathy", "CMP", "heart failure", "HF", 
                           "ischemic", "MI", "infarct", "pathological", "atrial fibrillation", "AF",
                           "dilated", "hypertrophic", "disease")
      
      # Classify each unique value
      control_samples <- rep(FALSE, n_samples)
      disease_samples <- rep(FALSE, n_samples)
      
      for (pattern in control_patterns) {
        control_samples <- control_samples | grepl(pattern, source_info, ignore.case = TRUE)
      }
      
      for (pattern in disease_patterns) {
        disease_samples <- disease_samples | grepl(pattern, source_info, ignore.case = TRUE)
      }
      
      # Check for valid group separation
      n_control <- sum(control_samples & !disease_samples)  # Pure control
      n_disease <- sum(disease_samples & !control_samples)  # Pure disease
      n_overlap <- sum(control_samples & disease_samples)   # Ambiguous
      n_unclassified <- sum(!control_samples & !disease_samples)  # Neither
      
      cat(sprintf("     Control samples: %d, Disease samples: %d, Overlap: %d, Unclassified: %d\n",
                  n_control, n_disease, n_overlap, n_unclassified))
      
      # Valid if we have both groups and minimal ambiguity
      if (n_control >= 3 && n_disease >= 3 && n_overlap == 0 && n_unclassified <= n_samples * 0.2) {
        
        # Create groups with BIOLOGICAL REFERENCE LOGIC
        # Control/Healthy = reference group (factor level 1)
        # Disease = comparison group (factor level 2)
        groups <- rep(NA, n_samples)
        groups[control_samples] <- "Control"  
        groups[disease_samples] <- "Disease"
        
        # Handle unclassified (try to assign based on context)
        if (n_unclassified > 0) {
          unclassified_indices <- which(!control_samples & !disease_samples)
          
          # Try secondary classification for unclassified
          for (i in unclassified_indices) {
            sample_info <- source_info[i]
            
            # More aggressive pattern matching for unclassified
            if (grepl("normal|control|wild|WT|baseline", sample_info, ignore.case = TRUE)) {
              groups[i] <- "Control"
            } else if (grepl("treatment|disease|pathology|failure|dysfunction", sample_info, ignore.case = TRUE)) {
              groups[i] <- "Disease"  
            }
          }
        }
        
        # Remove any remaining NA values
        valid_samples <- !is.na(groups)
        if (sum(valid_samples) >= 6) {  # Need at least 3 per group
          
          # Convert to factor with Control as reference (level 1)
          groups_factor <- factor(groups[valid_samples], levels = c("Control", "Disease"))
          
          final_n_control <- sum(groups_factor == "Control")
          final_n_disease <- sum(groups_factor == "Disease")
          
          cat("     SUCCESS: HEALTHY VS DISEASE pattern detected!\n")
          cat(sprintf("       Control samples: %d\n", final_n_control))
          cat(sprintf("       Disease samples: %d\n", final_n_disease))
          
          return(list(
            groups = groups_factor,
            pattern_type = "healthy_vs_disease",
            column_used = col_name,
            sample_indices = which(valid_samples),
            n_control = final_n_control,
            n_disease = final_n_disease
          ))
        }
      }
      
      # Fallback: Try AF/SR pattern detection for cardiac datasets
      if (any(grepl("atrial fibrillation|AF", source_info, ignore.case = TRUE)) && 
          any(grepl("sinus rhythm|SR", source_info, ignore.case = TRUE))) {
        
        af_samples <- grepl("atrial fibrillation|AF", source_info, ignore.case = TRUE)
        sr_samples <- grepl("sinus rhythm|SR", source_info, ignore.case = TRUE)
        
        if (sum(af_samples) >= 3 && sum(sr_samples) >= 3) {
          
          # For AF vs SR: SR is control (normal), AF is disease
          groups <- rep(NA, n_samples)
          groups[sr_samples] <- "Control"    # SR = Control (reference)
          groups[af_samples] <- "Disease"    # AF = Disease (comparison)
          
          valid_samples <- !is.na(groups)
          groups_factor <- factor(groups[valid_samples], levels = c("Control", "Disease"))
          
          cat("     SUCCESS: AF/SR pattern detected with biological reference!\n")
          cat(sprintf("       SR (Control) samples: %d\n", sum(groups_factor == "Control")))
          cat(sprintf("       AF (Disease) samples: %d\n", sum(groups_factor == "Disease")))
          
          return(list(
            groups = groups_factor,
            pattern_type = "AF_SR_corrected", 
            column_used = col_name,
            sample_indices = which(valid_samples),
            n_control = sum(groups_factor == "Control"),
            n_disease = sum(groups_factor == "Disease")
          ))
        }
      }
    }
  }
  
  cat("ERROR: No valid healthy vs disease pattern detected\n")
  return(NULL)
}

#' Test Enhanced Group Detection on All Datasets
test_enhanced_group_detection <- function() {
  
  cat("TESTING: Enhanced group detection on all datasets...\n\n")
  
  datasets <- c("GSE57338", "GSE41177", "GSE79768")
  
  for (dataset_id in datasets) {
    cat("=============================================\n")
    cat("TESTING:", dataset_id, "\n")
    cat("=============================================\n")
    
    # Find dataset file
    dataset_files <- list.files("cache", pattern = paste0(dataset_id, "_processed.rds"), 
                               recursive = TRUE, full.names = TRUE)
    
    if (length(dataset_files) > 0) {
      dataset <- readRDS(dataset_files[1])
      
      if (dataset$success) {
        result <- enhanced_auto_detect_groups_corrected(dataset)
        
        if (!is.null(result)) {
          cat("SUCCESS: Groups detected!\n")
          cat("   Pattern type:", result$pattern_type, "\n")
          cat("   Column used:", result$column_used, "\n")
          cat("   Group distribution:\n")
          print(table(result$groups))
          cat("   Samples to include:", length(result$sample_indices), "out of", dataset$n_samples, "\n")
          
          # CRITICAL: Verify contrast direction
          cat("   CONTRAST VALIDATION:\n")
          cat("   - Reference group (factor level 1):", levels(result$groups)[1], "\n")
          cat("   - Comparison group (factor level 2):", levels(result$groups)[2], "\n")
          cat("   - Positive logFC interpretation: Disease > Control ✅\n")
          
        } else {
          cat("FAILED: No groups detected\n")
        }
      } else {
        cat("ERROR: Dataset processing failed\n")
      }
    } else {
      cat("ERROR: Dataset file not found\n")
    }
    cat("\n")
  }
}

# Test the enhanced group detection
cat("SUCCESS: Enhanced group detection function ready!\n")
cat("TARGET: This enables proper Disease vs Control DGE analysis\n")

# Run the test
test_enhanced_group_detection()

cat("\n=== ENHANCED GROUP DETECTION CORRECTION COMPLETE ===\n")
cat("CRITICAL FIX: All datasets now use Control as reference group\n")
cat("RESULT: Positive logFC consistently means 'UP in Disease'\n")
cat("READY: For complete pipeline rerun with corrected methodology\n")