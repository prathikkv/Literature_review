#!/usr/bin/env Rscript
#' Enhanced Group Detection for Healthy vs Disease Analysis
#' 
#' Adds proper healthy vs disease group detection to the analysis functions

# Load required functions
source("functions/analysis.R")

cat("🎯 ENHANCED HEALTHY VS DISEASE GROUP DETECTION\n")
cat("==============================================\n\n")

# Enhanced auto_detect_groups function with healthy vs disease focus
enhanced_auto_detect_groups <- function(dataset) {
  
  pheno_data <- dataset$phenotype_data
  
  # Safety check for empty phenotype data
  if (is.null(pheno_data) || nrow(pheno_data) == 0) {
    return(NULL)
  }
  
  cat("🔍 Analyzing phenotype data for healthy vs disease patterns...\n")
  cat("   Samples:", nrow(pheno_data), "\n")
  cat("   Columns:", ncol(pheno_data), "\n")
  
  # Priority order: Look for healthy vs disease indicators first
  healthy_disease_columns <- c("source_name_ch1", "disease status:ch1", "heart failure:ch1", 
                               "disease", "condition", "group", "tissue", "treatment", 
                               "disease.state", "disease_state", "phenotype", "title")
  
  for (col in healthy_disease_columns) {
    if (col %in% colnames(pheno_data)) {
      cat("   Checking column:", col, "\n")
      
      values <- pheno_data[[col]]
      if (is.null(values) || length(values) == 0) next
      
      # Convert to character for pattern matching
      values_char <- as.character(values)
      unique_vals <- unique(values_char)
      
      cat("     Unique values:", length(unique_vals), "\n")
      
      # Healthy vs Disease Pattern Detection
      if (col == "source_name_ch1") {
        # Special handling for GSE57338 pattern
        healthy_pattern <- grepl("non-failing|control|normal", values_char, ignore.case = TRUE)
        disease_pattern <- grepl("dilated|ischemic|failure|cardiomyopathy", values_char, ignore.case = TRUE)
        
        if (any(healthy_pattern) && any(disease_pattern)) {
          groups <- ifelse(healthy_pattern, "Healthy", 
                          ifelse(disease_pattern, "Disease", "Other"))
          
          # Remove "Other" samples if present
          valid_samples <- groups %in% c("Healthy", "Disease")
          if (sum(valid_samples) >= 10 && sum(groups == "Healthy") >= 3 && sum(groups == "Disease") >= 3) {
            
            groups_final <- as.factor(groups[valid_samples])
            
            cat("     ✅ HEALTHY VS DISEASE pattern detected!\n")
            cat("       Healthy samples:", sum(groups == "Healthy"), "\n")
            cat("       Disease samples:", sum(groups == "Disease"), "\n")
            
            return(list(
              groups = groups_final,
              column = col,
              pattern_type = "healthy_vs_disease",
              sample_indices = which(valid_samples),
              group_labels = c("Healthy", "Disease")
            ))
          }
        }
      }
      
      # General healthy/control vs disease pattern detection
      healthy_keywords <- c("control", "normal", "healthy", "non-failing", "non_failing", "donor", "control")
      disease_keywords <- c("disease", "patient", "failure", "AF", "fibrillation", "ischemic", "dilated", "HF", "cardiomyopathy")
      
      # Count samples matching each pattern
      healthy_matches <- sapply(healthy_keywords, function(kw) 
        sum(grepl(kw, values_char, ignore.case = TRUE)))
      disease_matches <- sapply(disease_keywords, function(kw) 
        sum(grepl(kw, values_char, ignore.case = TRUE)))
      
      total_healthy <- sum(healthy_matches > 0)
      total_disease <- sum(disease_matches > 0)
      
      if (total_healthy > 0 && total_disease > 0) {
        # Create groups based on keyword matching
        is_healthy <- sapply(values_char, function(val) {
          any(sapply(healthy_keywords, function(kw) grepl(kw, val, ignore.case = TRUE)))
        })
        
        is_disease <- sapply(values_char, function(val) {
          any(sapply(disease_keywords, function(kw) grepl(kw, val, ignore.case = TRUE)))
        })
        
        # Create group labels
        groups <- rep("Other", length(values_char))
        groups[is_healthy] <- "Healthy"
        groups[is_disease] <- "Disease" 
        
        # Check if we have good group balance
        group_counts <- table(groups)
        if ("Healthy" %in% names(group_counts) && "Disease" %in% names(group_counts) &&
            group_counts["Healthy"] >= 3 && group_counts["Disease"] >= 3) {
          
          # Filter to just healthy and disease samples
          valid_samples <- groups %in% c("Healthy", "Disease")
          
          if (sum(valid_samples) >= 10) {
            groups_final <- as.factor(groups[valid_samples])
            
            cat("     ✅ General healthy vs disease pattern detected!\n")
            cat("       Healthy samples:", sum(groups == "Healthy"), "\n") 
            cat("       Disease samples:", sum(groups == "Disease"), "\n")
            
            return(list(
              groups = groups_final,
              column = col,
              pattern_type = "healthy_vs_disease_general",
              sample_indices = which(valid_samples),
              group_labels = c("Healthy", "Disease")
            ))
          }
        }
      }
      
      # Traditional group detection (fallback)
      groups <- as.factor(values)
      group_counts <- table(groups)
      if (length(group_counts) >= 2 && min(group_counts) >= 3 && max(group_counts) <= nrow(pheno_data) * 0.8) {
        cat("     ✅ Traditional groups detected:", names(group_counts), "\n")
        return(list(groups = groups, column = col, pattern_type = "traditional"))
      }
    }
  }
  
  # Fallback: AF/SR detection for atrial fibrillation datasets
  sample_info <- if (!is.null(pheno_data$title)) pheno_data$title else rownames(pheno_data)
  
  if (!is.null(sample_info) && length(sample_info) > 0) {
    if (any(grepl("AF|fibrillation", sample_info, ignore.case = TRUE)) && 
        any(grepl("SR|sinus|control", sample_info, ignore.case = TRUE))) {
      
      groups <- ifelse(grepl("AF|fibrillation", sample_info, ignore.case = TRUE), "AF", "SR")
      if (length(unique(groups)) >= 2) {
        cat("   ✅ AF/SR pattern detected (fallback)\n")
        return(list(groups = as.factor(groups), column = "auto_detected", pattern_type = "AF_SR"))
      }
    }
  }
  
  cat("   ❌ No suitable groups detected\n")
  return(NULL)
}

# Test the enhanced function on GSE57338
cat("🧪 Testing enhanced group detection on GSE57338...\n")
gse57338 <- readRDS("cache/comprehensive/GSE57338_processed.rds")
detected_groups <- enhanced_auto_detect_groups(gse57338)

if (!is.null(detected_groups)) {
  cat("✅ SUCCESS! Groups detected:\n")
  cat("   Pattern type:", detected_groups$pattern_type, "\n")
  cat("   Column used:", detected_groups$column, "\n") 
  cat("   Group distribution:\n")
  print(table(detected_groups$groups))
  
  if (!is.null(detected_groups$sample_indices)) {
    cat("   Samples to include:", length(detected_groups$sample_indices), "out of", nrow(gse57338$phenotype_data), "\n")
  }
} else {
  cat("❌ No groups detected\n")
}

cat("\n✅ Enhanced group detection function ready!\n")
cat("🎯 This enables proper healthy vs disease DGE analysis\n")