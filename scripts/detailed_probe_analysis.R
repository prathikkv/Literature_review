#!/usr/bin/env Rscript
#' Detailed CAMK2D Probe Analysis - Raw vs Processed Comparison

library(GEOquery)

cat("=== DETAILED CAMK2D PROBE ANALYSIS ===\n\n")

# Download and analyze original GSE57338 data
gset <- getGEO("GSE57338", GSEMatrix = TRUE, AnnotGPL = TRUE)
gset <- gset[[1]]

raw_expr <- exprs(gset)
raw_pheno <- pData(gset)
feature_info <- fData(gset)

cat("Raw GSE57338 Data Analysis:\n")
cat("- Expression matrix:", nrow(raw_expr), "x", ncol(raw_expr), "\n")
cat("- Platform:", gset@annotation, "\n\n")

# Find CAMK2D probe(s)
camk2d_probes <- which(grepl("CAMK2D", feature_info$`Gene symbol`, ignore.case = TRUE))
cat("CAMK2D probes found:", length(camk2d_probes), "\n")

if (length(camk2d_probes) > 0) {
  for (i in camk2d_probes) {
    probe_id <- rownames(feature_info)[i]
    gene_symbol <- feature_info$`Gene symbol`[i]
    gene_title <- feature_info$`Gene title`[i]
    
    cat(sprintf("Probe: %s, Symbol: %s\n", probe_id, gene_symbol))
    cat(sprintf("Title: %s\n", gene_title))
  }
  
  # Use the first CAMK2D probe for analysis
  camk2d_probe_id <- rownames(feature_info)[camk2d_probes[1]]
  cat(sprintf("\nUsing probe %s for analysis\n\n", camk2d_probe_id))
  
  # Extract CAMK2D expression
  camk2d_raw_expr <- raw_expr[camk2d_probe_id, ]
  
  # Identify sample groups from phenotype data
  source_names <- raw_pheno$source_name_ch1
  
  # Sample group identification
  dcm_samples <- grepl("dilated CMP", source_names, ignore.case = TRUE)
  nf_samples <- grepl("non-failing", source_names, ignore.case = TRUE)
  ischemic_samples <- grepl("ischemic", source_names, ignore.case = TRUE)
  
  cat("Sample Groups in Raw Data:\n")
  cat("- DCM samples:", sum(dcm_samples), "\n")
  cat("- Non-failing samples:", sum(nf_samples), "\n") 
  cat("- Ischemic samples:", sum(ischemic_samples), "\n")
  cat("- Total samples:", length(source_names), "\n\n")
  
  # Extract expression for DCM vs NF comparison
  dcm_expr <- camk2d_raw_expr[dcm_samples]
  nf_expr <- camk2d_raw_expr[nf_samples]
  
  cat("RAW DATA CAMK2D ANALYSIS:\n")
  cat("=========================\n")
  cat("DCM expression stats:\n")
  cat("- Count:", length(dcm_expr), "\n")
  cat("- Mean:", round(mean(dcm_expr), 4), "\n")
  cat("- Range:", round(range(dcm_expr), 4), "\n")
  cat("- SD:", round(sd(dcm_expr), 4), "\n\n")
  
  cat("NF expression stats:\n")
  cat("- Count:", length(nf_expr), "\n")
  cat("- Mean:", round(mean(nf_expr), 4), "\n")
  cat("- Range:", round(range(nf_expr), 4), "\n") 
  cat("- SD:", round(sd(nf_expr), 4), "\n\n")
  
  # Calculate fold change and statistics
  raw_dcm_mean <- mean(dcm_expr)
  raw_nf_mean <- mean(nf_expr)
  raw_fc <- raw_dcm_mean / raw_nf_mean
  raw_log2fc <- log2(raw_fc)
  
  # Statistical test
  raw_ttest <- t.test(dcm_expr, nf_expr)
  
  cat("RAW DATA RESULTS:\n")
  cat("- DCM mean:", round(raw_dcm_mean, 4), "\n")
  cat("- NF mean:", round(raw_nf_mean, 4), "\n")
  cat("- Fold change:", round(raw_fc, 4), "\n")
  cat("- Log2 FC:", round(raw_log2fc, 4), "\n")
  cat("- P-value:", format(raw_ttest$p.value, scientific = TRUE), "\n\n")
  
  # Compare with publication
  cat("COMPARISON WITH PUBLICATION:\n")
  cat("============================\n")
  cat("Publication (DCM vs NF):\n")
  cat("- FPKM DCM: 59.7739\n")
  cat("- FPKM NF: 28.3137\n")
  cat("- FC: 2.1111\n") 
  cat("- Log2FC: 1.078\n")
  cat("- P-value: 0.000521\n\n")
  
  cat("Our raw data analysis:\n")
  cat("- Raw DCM:", round(raw_dcm_mean, 4), "\n")
  cat("- Raw NF:", round(raw_nf_mean, 4), "\n")
  cat("- FC:", round(raw_fc, 4), "\n")
  cat("- Log2FC:", round(raw_log2fc, 4), "\n")
  cat("- P-value:", format(raw_ttest$p.value, scientific = TRUE), "\n\n")
  
  # Scale comparison
  cat("SCALE ANALYSIS:\n")
  cat("===============\n")
  pub_scale_ratio <- 59.7739 / raw_dcm_mean
  cat("Publication FPKM / Our raw expression =", round(pub_scale_ratio, 4), "\n")
  cat("This suggests", ifelse(pub_scale_ratio > 1, "FPKM values are higher", "raw values are higher"), "scale\n\n")
  
  # Load our processed data for comparison
  cat("COMPARISON WITH OUR PROCESSED PIPELINE:\n")
  cat("=======================================\n")
  
  dataset <- readRDS("cache/comprehensive/GSE57338_processed.rds")
  processed_expr <- dataset$expression_matrix
  processed_pheno <- dataset$phenotype_data
  
  # Check if CAMK2D exists in processed data
  if ("CAMK2D" %in% rownames(processed_expr)) {
    processed_camk2d <- processed_expr["CAMK2D", ]
    
    # Same sample identification
    proc_source_names <- processed_pheno$source_name_ch1
    proc_dcm_samples <- grepl("dilated CMP", proc_source_names)
    proc_nf_samples <- grepl("non-failing", proc_source_names)
    
    proc_dcm_expr <- processed_camk2d[proc_dcm_samples]
    proc_nf_expr <- processed_camk2d[proc_nf_samples]
    
    proc_dcm_mean <- mean(proc_dcm_expr)
    proc_nf_mean <- mean(proc_nf_expr)
    proc_log2fc <- proc_dcm_mean - proc_nf_mean  # Already log2 scale
    
    cat("Processed data (log2 scale):\n")
    cat("- DCM mean:", round(proc_dcm_mean, 4), "\n")
    cat("- NF mean:", round(proc_nf_mean, 4), "\n")
    cat("- Log2FC:", round(proc_log2fc, 4), "\n")
    
    # Convert to linear for comparison
    proc_dcm_linear <- 2^proc_dcm_mean
    proc_nf_linear <- 2^proc_nf_mean
    proc_linear_fc <- proc_dcm_linear / proc_nf_linear
    
    cat("\nProcessed data (converted to linear):\n")
    cat("- DCM linear:", round(proc_dcm_linear, 4), "\n")
    cat("- NF linear:", round(proc_nf_linear, 4), "\n")
    cat("- Linear FC:", round(proc_linear_fc, 4), "\n\n")
    
    # COMPREHENSIVE COMPARISON
    cat("COMPREHENSIVE COMPARISON:\n")
    cat("=========================\n")
    cat("Data Source           | DCM Mean  | NF Mean   | FC      | Log2FC  | P-value\n")
    cat("---------------------|-----------|-----------|---------|---------|----------\n")
    cat(sprintf("Publication (FPKM)   | %8.4f | %8.4f | %7.4f | %7.4f | %.6f\n", 
               59.7739, 28.3137, 2.1111, 1.078, 0.000521))
    cat(sprintf("Raw probe data       | %8.4f | %8.4f | %7.4f | %7.4f | %.6f\n",
               raw_dcm_mean, raw_nf_mean, raw_fc, raw_log2fc, raw_ttest$p.value))
    cat(sprintf("Our pipeline (log2)  | %8.4f | %8.4f | %7.4f | %7.4f | %.6f\n",
               proc_dcm_mean, proc_nf_mean, proc_linear_fc, proc_log2fc, 0.127508))
    cat(sprintf("Our pipeline (linear)| %8.4f | %8.4f | %7.4f | %7.4f | %.6f\n",
               proc_dcm_linear, proc_nf_linear, proc_linear_fc, proc_log2fc, 0.127508))
    
    cat("\nKEY FINDINGS:\n")
    if (abs(raw_log2fc - 1.078) < abs(proc_log2fc - 1.078)) {
      cat("✅ RAW DATA CLOSER TO PUBLICATION than our processed pipeline\n")
      cat("❗ Our preprocessing pipeline may be dampening the signal\n")
    } else {
      cat("❌ Both raw and processed data differ significantly from publication\n")
    }
    
    if (raw_ttest$p.value < 0.05) {
      cat("✅ RAW DATA SHOWS SIGNIFICANT RESULT\n")
    } else {
      cat("❌ Raw data also non-significant (p =", format(raw_ttest$p.value, scientific=TRUE), ")\n")
    }
    
    # Direction consistency
    all_positive <- all(c(raw_log2fc, proc_log2fc, 1.078) > 0)
    if (all_positive) {
      cat("✅ ALL ANALYSES SHOW UPREGULATION (direction consistent)\n")
    }
    
  } else {
    cat("ERROR: CAMK2D not found in processed data\n")
  }
  
} else {
  cat("ERROR: No CAMK2D probes found\n")
}

cat("\n=== DETAILED PROBE ANALYSIS COMPLETE ===\n")
cat("Summary: Investigating whether raw data better matches publication\n")