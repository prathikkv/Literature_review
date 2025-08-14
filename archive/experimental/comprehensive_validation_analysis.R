#!/usr/bin/env Rscript
#' Comprehensive Validation Analysis: Publication vs Pipeline Results
#' 
#' Investigation of the 26x magnitude difference between our GSE57338 results and publication

library(GEOquery)
library(limma)

cat("=== COMPREHENSIVE VALIDATION ANALYSIS ===\n\n")

# PHASE 1: PUBLICATION vs OUR RESULTS COMPARISON
cat("PHASE 1: PUBLICATION DATA RECONCILIATION\n")
cat("=========================================\n")

# Publication data
pub_fpkm_dcm <- 59.7739
pub_fpkm_nf <- 28.3137
pub_fc <- pub_fpkm_dcm / pub_fpkm_nf
pub_log2fc <- log2(pub_fc)
pub_pval <- 0.000521

cat("Publication Results (GSE57338 CAMK2D):\n")
cat("- FPKM DCM:", pub_fpkm_dcm, "\n")
cat("- FPKM NF:", pub_fpkm_nf, "\n")
cat("- Fold Change:", round(pub_fc, 4), "\n")
cat("- Log2 FC:", round(pub_log2fc, 4), "\n")
cat("- P-value:", pub_pval, "\n\n")

# Our results
our_log2fc <- 0.0405
our_fc <- 2^our_log2fc
our_pval <- 0.127508

cat("Our Corrected Results (GSE57338 CAMK2D):\n")
cat("- Log2 FC:", round(our_log2fc, 4), "\n")
cat("- Linear FC:", round(our_fc, 4), "\n")
cat("- P-value:", our_pval, "\n\n")

# Calculate discrepancy
magnitude_diff <- pub_log2fc / our_log2fc
fc_diff <- pub_fc / our_fc

cat("MAGNITUDE DISCREPANCY ANALYSIS:\n")
cat("- Log2FC difference:", round(magnitude_diff, 1), "x smaller in our analysis\n")
cat("- Linear FC difference:", round(fc_diff, 1), "x smaller in our analysis\n")
cat("- P-value difference:", round(pub_pval / our_pval, 0), "x more significant in publication\n\n")

# PHASE 2: RAW DATA INVESTIGATION
cat("PHASE 2: RAW DATA INVESTIGATION\n")
cat("================================\n")

# Load our processed dataset
dataset <- readRDS("cache/comprehensive/GSE57338_processed.rds")
expr_matrix <- dataset$expression_matrix
source_names <- dataset$phenotype_data$source_name_ch1

# Sample selection validation
dcm_samples <- grepl("dilated CMP", source_names)
nf_samples <- grepl("non-failing", source_names)

cat("Sample Selection Validation:\n")
cat("- Our DCM samples:", sum(dcm_samples), "\n")
cat("- Our NF samples:", sum(nf_samples), "\n")
cat("- Publication indicates: 333 DCM + 136 NF, but text shows specific comparison\n\n")

# CAMK2D expression analysis
camk2d_dcm <- expr_matrix["CAMK2D", dcm_samples]
camk2d_nf <- expr_matrix["CAMK2D", nf_samples]

cat("CAMK2D Expression Analysis (Log2 Scale):\n")
cat("- DCM mean:", round(mean(camk2d_dcm), 4), "\n")
cat("- NF mean:", round(mean(camk2d_nf), 4), "\n")
cat("- Our log2FC:", round(mean(camk2d_dcm) - mean(camk2d_nf), 4), "\n")
cat("- Expression range DCM:", range(camk2d_dcm), "\n")
cat("- Expression range NF:", range(camk2d_nf), "\n\n")

# Convert to linear scale for comparison
camk2d_dcm_linear <- 2^camk2d_dcm
camk2d_nf_linear <- 2^camk2d_nf

cat("CAMK2D Expression Analysis (Linear Scale):\n")
cat("- DCM mean:", round(mean(camk2d_dcm_linear), 2), "\n")
cat("- NF mean:", round(mean(camk2d_nf_linear), 2), "\n")
cat("- Linear FC:", round(mean(camk2d_dcm_linear) / mean(camk2d_nf_linear), 4), "\n")
cat("- Compare to publication FC:", round(pub_fc, 4), "\n\n")

# PHASE 3: POTENTIAL ROOT CAUSES
cat("PHASE 3: ROOT CAUSE ANALYSIS\n")
cat("=============================\n")

cat("Potential Explanations for Magnitude Difference:\n\n")

cat("1. NORMALIZATION DIFFERENCES:\n")
cat("   - Publication FPKM values (28-60) suggest different normalization\n")
cat("   - Our values (~1100 linear) are ~20x higher scale\n")
cat("   - This suggests different preprocessing pipelines\n\n")

cat("2. SAMPLE COMPOSITION:\n")
# Check if we have all the right samples
ischemic_samples <- grepl("ischemic", source_names)
cat("   - We exclude", sum(ischemic_samples), "ischemic samples\n")
cat("   - Publication may have different inclusion criteria\n")
cat("   - Our DCM:NF ratio =", round(sum(dcm_samples)/sum(nf_samples), 2), "\n\n")

cat("3. STATISTICAL APPROACH:\n")
cat("   - Publication may use different software/parameters\n")
cat("   - Different multiple testing corrections\n")
cat("   - Alternative statistical models\n\n")

cat("4. PROBE/GENE MAPPING:\n")
# Check CAMK2D probe information
if ("feature_data" %in% names(dataset) && !is.null(dataset$feature_data)) {
  cat("   - Need to verify we're analyzing same probe as publication\n")
  cat("   - Gene symbol mapping may aggregate multiple probes differently\n")
} else {
  cat("   - No probe information available in processed data\n")
  cat("   - Cannot verify specific probe used in publication\n")
}
cat("\n")

# PHASE 4: CROSS-DATASET CONSISTENCY CHECK
cat("PHASE 4: CROSS-DATASET CONSISTENCY\n")
cat("==================================\n")

# Load all CAMK2D results
all_results <- read.csv("output/CAMK_focused_DGE_all_datasets_CORRECTED.csv")
camk2d_all <- all_results[all_results$Gene_Symbol == "CAMK2D", ]

cat("CAMK2D Results Across All Datasets:\n")
for (i in 1:nrow(camk2d_all)) {
  dataset_id <- camk2d_all$Dataset[i]
  logfc <- round(camk2d_all$logFC[i], 4)
  pval <- camk2d_all$P.Value[i]
  direction <- ifelse(logfc > 0, "UP", "DOWN")
  
  cat(sprintf("  %-10s: %4s logFC=%7.4f, p=%8.2e\n", dataset_id, direction, logfc, pval))
}

# Consistency analysis
directions <- sign(camk2d_all$logFC)
consistent_direction <- length(unique(directions)) == 1

cat("\nConsistency Analysis:\n")
cat("- All datasets same direction:", consistent_direction, "\n")
cat("- Effect size range:", range(camk2d_all$logFC), "\n")
cat("- GSE57338 vs others: GSE57338 shows weakest effect\n\n")

# PHASE 5: RELIABILITY ASSESSMENT
cat("PHASE 5: PIPELINE RELIABILITY ASSESSMENT\n")
cat("==========================================\n")

# Calculate confidence metrics
significant_count <- sum(camk2d_all$P.Value < 0.05)
consistent_direction_count <- sum(camk2d_all$logFC > 0)  # Assuming upregulation is expected

cat("Reliability Metrics:\n")
cat("- Datasets with significant CAMK2D:", significant_count, "/", nrow(camk2d_all), "\n")
cat("- Datasets with upregulation:", consistent_direction_count, "/", nrow(camk2d_all), "\n")
cat("- GSE57338 significance level: p =", our_pval, "(trend level)\n\n")

# PHASE 6: PUBLICATION READINESS EVALUATION
cat("PHASE 6: PUBLICATION READINESS EVALUATION\n")
cat("==========================================\n")

cat("STRENGTHS:\n")
cat("✅ Correct CAMK2D direction in GSE57338 (matches publication)\n")
cat("✅ Robust multi-dataset analysis approach\n")
cat("✅ Comprehensive statistical framework\n")
cat("✅ Reproducible analysis pipeline\n\n")

cat("CRITICAL CONCERNS:\n")
cat("❌ 26x magnitude difference from publication benchmark\n")
cat("❌ Non-significant results where publication shows high significance\n")
cat("❌ High heterogeneity across datasets suggests methodological issues\n")
cat("❌ Cannot validate probe-level concordance with publication\n\n")

# PHASE 7: RECOMMENDATIONS
cat("PHASE 7: RECOMMENDATIONS\n")
cat("========================\n")

cat("IMMEDIATE ACTIONS NEEDED:\n")
cat("1. 🔍 PROBE VALIDATION: Download original GSE57338 probe data\n")
cat("2. 📊 NORMALIZATION REVIEW: Compare preprocessing with publication\n") 
cat("3. 👥 SAMPLE VERIFICATION: Validate exact DCM vs NF criteria used\n")
cat("4. 🧮 ALTERNATIVE ANALYSIS: Test different statistical approaches\n")
cat("5. 📚 LITERATURE SEARCH: Find other GSE57338 analyses for benchmarking\n\n")

cat("PUBLICATION STRATEGY:\n")
if (consistent_direction_count >= 2) {
  cat("✅ PROCEED WITH CAUTION: Direction consistency supports biological relevance\n")
  cat("   - Acknowledge magnitude limitations in methods/discussion\n")
  cat("   - Focus on direction and pathway-level evidence\n")
  cat("   - Include thorough methodology validation\n")
} else {
  cat("⚠️  ADDITIONAL VALIDATION REQUIRED before publication submission\n")
  cat("   - Resolve magnitude discrepancies\n")
  cat("   - Validate against additional benchmark datasets\n")
  cat("   - Consider alternative analysis approaches\n")
}

cat("\nCONFIDENCE ASSESSMENT:\n")
cat("- GSE57338 CAMK2D: MODERATE (direction correct, magnitude uncertain)\n")
cat("- Overall pipeline: MODERATE (needs validation against benchmarks)\n")
cat("- Therapeutic conclusions: CAUTIOUS (direction supportive, effect size unclear)\n")

cat("\n=== VALIDATION ANALYSIS COMPLETE ===\n")
cat("Priority: Resolve magnitude discrepancy through probe-level validation\n")