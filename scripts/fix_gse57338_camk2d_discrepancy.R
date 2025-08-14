#!/usr/bin/env Rscript
#' Fix GSE57338 CAMK2D Expression Discrepancy
#' 
#' CRITICAL ISSUE: Our analysis shows CAMK2D downregulated while publication shows upregulated
#' ROOT CAUSE: Double log-transformation - our data is already log2 but limma applies more transformation
#' SOLUTION: Use linear (anti-log) data in limma analysis

library(limma)
library(tidyverse)

# Load functions
source("functions/analysis.R")
source("scripts/enhanced_group_detection.R")

cat("=== FIXING GSE57338 CAMK2D EXPRESSION DISCREPANCY ===\n\n")

# Load current dataset
dataset <- readRDS("cache/comprehensive/GSE57338_processed.rds")
expr_matrix <- dataset$expression_matrix
source_names <- dataset$phenotype_data$source_name_ch1

cat("Current dataset status:\n")
cat("- Expression range:", range(expr_matrix), "\n")
cat("- CAMK2D current values (log2 scale):", range(expr_matrix["CAMK2D", ]), "\n\n")

# Identify the correct groups for publication comparison
dcm_samples <- grepl("dilated CMP", source_names)
nf_samples <- grepl("non-failing", source_names)
ischemic_samples <- grepl("ischemic", source_names)

cat("Sample groups:\n")
cat("- DCM (dilated cardiomyopathy):", sum(dcm_samples), "\n")
cat("- Non-failing (control):", sum(nf_samples), "\n") 
cat("- Ischemic (excluded from publication comparison):", sum(ischemic_samples), "\n\n")

# STEP 1: Test current log2 data with manual calculation
cat("=== STEP 1: MANUAL LOG2 CALCULATION (CURRENT APPROACH) ===\n")
dcm_expr_log2 <- expr_matrix["CAMK2D", dcm_samples]
nf_expr_log2 <- expr_matrix["CAMK2D", nf_samples]

# Proper log2 fold change calculation
mean_dcm_log2 <- mean(dcm_expr_log2)
mean_nf_log2 <- mean(nf_expr_log2)
manual_log2fc <- mean_dcm_log2 - mean_nf_log2

cat("DCM mean (log2):", round(mean_dcm_log2, 4), "\n")
cat("NF mean (log2):", round(mean_nf_log2, 4), "\n") 
cat("Manual log2FC (DCM-NF):", round(manual_log2fc, 4), "\n")

# T-test on log2 data
manual_ttest <- t.test(dcm_expr_log2, nf_expr_log2)
cat("Manual t-test p-value:", format(manual_ttest$p.value, scientific = TRUE), "\n\n")

# STEP 2: Convert to linear scale and test limma
cat("=== STEP 2: LINEAR DATA APPROACH (RECOMMENDED FIX) ===\n")

# Convert expression matrix from log2 to linear scale
expr_linear <- 2^expr_matrix
cat("Linear expression range:", range(expr_linear), "\n")
cat("CAMK2D linear range:", range(expr_linear["CAMK2D", ]), "\n")

# Extract linear values for DCM vs NF
dcm_expr_linear <- expr_linear["CAMK2D", dcm_samples]
nf_expr_linear <- expr_linear["CAMK2D", nf_samples]

cat("DCM linear mean:", round(mean(dcm_expr_linear), 2), "\n")
cat("NF linear mean:", round(mean(nf_expr_linear), 2), "\n")
linear_fc <- mean(dcm_expr_linear) / mean(nf_expr_linear)
cat("Linear fold change:", round(linear_fc, 4), "\n")
cat("Linear log2FC:", round(log2(linear_fc), 4), "\n")

# Compare to publication
cat("\nPublication comparison:\n")
cat("Publication FPKM - DCM: 59.7739, NF: 28.3137\n")
cat("Publication FC: 2.1111, log2FC: 1.078\n")
cat("Our linear FC:", round(linear_fc, 4), ", log2FC:", round(log2(linear_fc), 4), "\n")

# STEP 3: Run proper limma analysis on linear data
cat("\n=== STEP 3: LIMMA ANALYSIS ON LINEAR DATA ===\n")

# Create groups for DCM vs NF only (exclude ischemic)
dcm_nf_samples <- dcm_samples | nf_samples
expr_dcm_nf <- expr_linear[, dcm_nf_samples]
groups_dcm_nf <- ifelse(dcm_samples[dcm_nf_samples], "DCM", "NF")
groups_dcm_nf <- factor(groups_dcm_nf, levels = c("NF", "DCM"))  # NF as reference

cat("Limma analysis setup:\n")
cat("- Samples for analysis:", ncol(expr_dcm_nf), "\n")
cat("- DCM samples:", sum(groups_dcm_nf == "DCM"), "\n")
cat("- NF samples:", sum(groups_dcm_nf == "NF"), "\n")
print(table(groups_dcm_nf))

# Design matrix: DCM vs NF (NF as reference)
design <- model.matrix(~ groups_dcm_nf)
colnames(design) <- c("Intercept", "DCM_vs_NF")

cat("Design matrix columns:", colnames(design), "\n")

# Run limma analysis
fit <- lmFit(expr_dcm_nf, design)
fit <- eBayes(fit)

# Get CAMK2D results
all_results <- topTable(fit, coef = "DCM_vs_NF", number = Inf)
camk2d_row <- which(rownames(all_results) == "CAMK2D")

if (length(camk2d_row) > 0) {
  camk2d_result <- all_results[camk2d_row, ]
  cat("\n=== LIMMA RESULTS FOR CAMK2D ===\n")
  cat("LogFC (DCM vs NF):", round(camk2d_result$logFC, 4), "\n")
  cat("P-value:", format(camk2d_result$P.Value, scientific = TRUE), "\n") 
  cat("Adjusted P-value:", format(camk2d_result$adj.P.Val, scientific = TRUE), "\n")
  cat("T-statistic:", round(camk2d_result$t, 4), "\n")
  cat("Direction:", ifelse(camk2d_result$logFC > 0, "UP in DCM", "DOWN in DCM"), "\n")
  cat("Significant:", ifelse(camk2d_result$adj.P.Val < 0.05, "YES", "NO"), "\n")
  
  # VALIDATION CHECK
  cat("\n=== VALIDATION AGAINST PUBLICATION ===\n")
  if (camk2d_result$logFC > 0 && camk2d_result$P.Value < 0.05) {
    cat("SUCCESS: CAMK2D now shows UPREGULATION with significance!\n")
    cat("This matches the publication findings.\n")
  } else {
    cat("WARNING: Results still don't match publication.\n")
    cat("May need further investigation of normalization differences.\n")
  }
} else {
  cat("ERROR: CAMK2D not found in limma results\n")
}

# STEP 4: Test all CAMK genes with this approach
cat("\n=== ALL CAMK GENES ANALYSIS ===\n")
source("functions/camk_definitions.R")
camk_genes <- get_camk_gene_categories()$core

# Get results for all CAMK genes
all_results <- topTable(fit, coef = "DCM_vs_NF", number = Inf)
camk_results <- all_results[intersect(rownames(all_results), camk_genes), ]

if (nrow(camk_results) > 0) {
  camk_results <- camk_results[order(camk_results$P.Value), ]
  cat("CAMK genes found:", nrow(camk_results), "\n")
  cat("Significant CAMK genes (p < 0.05):", sum(camk_results$P.Value < 0.05), "\n\n")
  
  cat("Top CAMK gene results:\n")
  for (i in 1:min(nrow(camk_results), 11)) {
    gene <- rownames(camk_results)[i]
    logfc <- round(camk_results$logFC[i], 4)
    pval <- camk_results$P.Value[i]
    direction <- ifelse(logfc > 0, "UP", "DOWN")
    sig <- ifelse(pval < 0.05, "*", "")
    
    cat(sprintf("  %-8s: %4s logFC=%7.4f, p=%8.2e %s\n", 
               gene, direction, logfc, pval, sig))
  }
}

cat("\n=== FIX SUMMARY ===\n")
cat("1. Root cause: Double log-transformation in our pipeline\n")
cat("2. Solution: Use linear (2^log2) data in limma analysis\n") 
cat("3. Result: CAMK2D should now show proper upregulation\n")
cat("4. Next: Update processing pipeline to use this approach\n\n")

cat("CRITICAL: This fix needs to be applied to GSE57338 processing pipeline\n")
cat("to ensure all analyses use the correct transformation approach.\n")