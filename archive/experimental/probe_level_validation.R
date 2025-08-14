#!/usr/bin/env Rscript
#' Probe-Level Validation Analysis
#' 
#' Deep dive into GSE57338 probe mapping and sample selection methodology

library(GEOquery)

cat("=== PROBE-LEVEL VALIDATION ANALYSIS ===\n\n")

# Try to download original GSE57338 data for comparison
cat("PHASE 1: ORIGINAL DATA DOWNLOAD\n")
cat("===============================\n")

tryCatch({
  cat("Downloading GSE57338 series matrix...\n")
  gset <- getGEO("GSE57338", GSEMatrix = TRUE, AnnotGPL = TRUE)
  
  if (length(gset) > 0) {
    gset <- gset[[1]]
    
    # Get raw expression data
    raw_expr <- exprs(gset)
    raw_pheno <- pData(gset)
    feature_info <- fData(gset)
    
    cat("Raw data dimensions:\n")
    cat("- Expression matrix:", nrow(raw_expr), "x", ncol(raw_expr), "\n")
    cat("- Feature annotations:", nrow(feature_info), "columns:", ncol(feature_info), "\n")
    cat("- Phenotype data:", nrow(raw_pheno), "columns:", ncol(raw_pheno), "\n\n")
    
    # Check platform information
    cat("Platform Information:\n")
    cat("- Platform ID:", gset@annotation, "\n")
    if ("platform_id" %in% colnames(raw_pheno)) {
      cat("- Platform details:", unique(raw_pheno$platform_id), "\n")
    }
    cat("\n")
    
    # Analyze CAMK2D probes
    cat("CAMK2D PROBE ANALYSIS:\n")
    cat("======================\n")
    
    # Look for CAMK2D in feature annotations
    gene_symbol_cols <- c("Gene.symbol", "Gene_Symbol", "GENE_SYMBOL", "gene_symbol", "Symbol", "symbol")
    camk2d_probes <- NULL
    gene_col_used <- NULL
    
    for (col in gene_symbol_cols) {
      if (col %in% colnames(feature_info)) {
        gene_symbols <- feature_info[[col]]
        camk2d_indices <- which(grepl("CAMK2D", gene_symbols, ignore.case = TRUE))
        if (length(camk2d_indices) > 0) {
          camk2d_probes <- rownames(feature_info)[camk2d_indices]
          gene_col_used <- col
          break
        }
      }
    }
    
    if (!is.null(camk2d_probes)) {
      cat("CAMK2D probes found:", length(camk2d_probes), "\n")
      cat("Gene symbol column used:", gene_col_used, "\n")
      cat("Probe IDs:", paste(camk2d_probes, collapse = ", "), "\n")
      
      # Get probe details
      for (probe in camk2d_probes) {
        probe_info <- feature_info[probe, ]
        cat(sprintf("  Probe %s: %s\n", probe, probe_info[[gene_col_used]]))
        
        # Additional annotation columns
        other_cols <- setdiff(colnames(feature_info), c("ID", gene_col_used))
        for (col in other_cols[1:min(3, length(other_cols))]) {
          if (!is.na(probe_info[[col]]) && probe_info[[col]] != "") {
            cat(sprintf("    %s: %s\n", col, probe_info[[col]]))
          }
        }
      }
      
      # Analyze CAMK2D expression across samples
      cat("\nCAMK2D PROBE EXPRESSION ANALYSIS:\n")
      
      # Sample group identification
      source_col <- NULL
      for (col in c("source_name_ch1", "title", "description")) {
        if (col %in% colnames(raw_pheno)) {
          source_col <- col
          break
        }
      }
      
      if (!is.null(source_col)) {
        source_info <- raw_pheno[[source_col]]
        cat("Sample source column:", source_col, "\n")
        
        # Identify sample groups
        dcm_samples <- grepl("dilated|CMP", source_info, ignore.case = TRUE)
        nf_samples <- grepl("non-failing", source_info, ignore.case = TRUE)
        
        cat("Sample groups identified:\n")
        cat("- DCM samples:", sum(dcm_samples), "\n")
        cat("- Non-failing samples:", sum(nf_samples), "\n")
        
        if (sum(dcm_samples) > 0 && sum(nf_samples) > 0) {
          # Analyze each CAMK2D probe
          for (probe in camk2d_probes) {
            probe_expr <- raw_expr[probe, ]
            
            dcm_expr <- probe_expr[dcm_samples]
            nf_expr <- probe_expr[nf_samples]
            
            # Calculate statistics
            dcm_mean <- mean(dcm_expr, na.rm = TRUE)
            nf_mean <- mean(nf_expr, na.rm = TRUE)
            fold_change <- dcm_mean / nf_mean
            log2_fc <- log2(fold_change)
            
            # T-test
            ttest_result <- t.test(dcm_expr, nf_expr)
            
            cat(sprintf("\\nProbe %s Results:\\n", probe))
            cat(sprintf("  DCM mean: %.4f\\n", dcm_mean))
            cat(sprintf("  NF mean: %.4f\\n", nf_mean))
            cat(sprintf("  Fold change: %.4f\\n", fold_change))
            cat(sprintf("  Log2 FC: %.4f\\n", log2_fc))
            cat(sprintf("  P-value: %.6f\\n", ttest_result$p.value))
            
            # Compare to publication
            cat("  Comparison to publication:\\n")
            cat(sprintf("    Publication FC: 2.1111, Our FC: %.4f\\n", fold_change))
            cat(sprintf("    Publication logFC: 1.078, Our logFC: %.4f\\n", log2_fc))
            cat(sprintf("    Publication p: 0.000521, Our p: %.6f\\n", ttest_result$p.value))
            
            # Check if this matches publication better
            if (abs(log2_fc - 1.078) < abs(0.0405 - 1.078)) {
              cat("    ✅ THIS PROBE CLOSER TO PUBLICATION THAN OUR PIPELINE\\n")
            } else {
              cat("    ❌ Still different from publication\\n")
            }
          }
        } else {
          cat("ERROR: Could not identify DCM and NF sample groups\\n")
        }
      } else {
        cat("ERROR: Could not find sample source information\\n")
      }
      
    } else {
      cat("ERROR: No CAMK2D probes found in feature annotations\\n")
      cat("Available gene symbol columns:\\n")
      print(gene_symbol_cols[gene_symbol_cols %in% colnames(feature_info)])
      
      if (ncol(feature_info) > 0) {
        cat("\\nFirst few feature annotation columns:\\n")
        print(head(colnames(feature_info), 10))
        
        # Try searching for CAMK2D in all text columns
        text_cols <- sapply(feature_info, is.character)
        for (col in names(text_cols)[text_cols]) {
          camk2d_matches <- grepl("CAMK2D", feature_info[[col]], ignore.case = TRUE)
          if (any(camk2d_matches)) {
            cat(sprintf("Found CAMK2D in column '%s': %d matches\\n", col, sum(camk2d_matches)))
          }
        }
      }
    }
    
  } else {
    cat("ERROR: Could not download GSE57338 data\\n")
  }
  
}, error = function(e) {
  cat("ERROR downloading GSE57338 data:", e$message, "\\n")
  
  # Alternative: Analyze our processed data more thoroughly
  cat("\\nFALLBACK: ANALYZING OUR PROCESSED DATA\\n")
  cat("======================================\\n")
  
  dataset <- readRDS("cache/comprehensive/GSE57338_processed.rds")
  
  if ("feature_data" %in% names(dataset) && !is.null(dataset$feature_data)) {
    feature_data <- dataset$feature_data
    cat("Feature data available in processed dataset\\n")
    cat("Columns:", paste(colnames(feature_data), collapse = ", "), "\\n")
    
    # Look for CAMK2D probe information
    if ("Gene.symbol" %in% colnames(feature_data)) {
      camk2d_features <- feature_data[feature_data$Gene.symbol == "CAMK2D", ]
      if (nrow(camk2d_features) > 0) {
        cat("CAMK2D probe details:\\n")
        print(camk2d_features)
      }
    }
  } else {
    cat("No detailed feature data available in processed dataset\\n")
    cat("Cannot perform probe-level validation\\n")
  }
})

cat("\\n=== PROBE-LEVEL VALIDATION COMPLETE ===\\n")
cat("Next: Review normalization methodology differences\\n")