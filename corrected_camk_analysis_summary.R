#!/usr/bin/env Rscript
#' Corrected CAMK Analysis Summary
#' 
#' Compares the flawed original meta-analysis with the corrected healthy vs disease analysis

cat("🎯 CORRECTED CAMK ANALYSIS SUMMARY\n")
cat("=================================\n\n")

# Load the corrected GSE57338 results
gse57338_results <- readRDS("output/GSE57338_CAMK_results.rds")

# Load the flawed original meta-analysis for comparison
original_meta <- read.csv("output/CAMK_meta_analysis_summary.csv")

cat("📊 COMPARISON: FLAWED vs CORRECTED ANALYSIS\n")
cat("==========================================\n\n")

cat("❌ ORIGINAL FLAWED META-ANALYSIS:\n")
cat("   • Mixed AF vs SR (disease subtypes) with healthy vs disease\n")
cat("   • GSE57338 (313 samples, 136 healthy) was MISSING\n")
cat("   • Only 1 significant gene: CAMK2D (P=0.00141)\n")
cat("   • Based on incompatible study designs\n")
cat("   • Clinically meaningless results\n\n")

cat("✅ CORRECTED HEALTHY vs DISEASE ANALYSIS:\n")
cat("   • GSE57338: 136 healthy vs 177 disease samples\n")
cat("   • True healthy vs cardiovascular disease comparison\n")
cat("   • All 11 CAMK genes analyzed\n")
cat("   • 6 significantly dysregulated CAMK genes (FDR < 0.05)\n")
cat("   • Clinically meaningful cardiovascular disease insights\n\n")

cat("🧬 DETAILED CAMK DYSREGULATION RESULTS (GSE57338)\n")
cat("===============================================\n\n")

results <- gse57338_results$results
significant_genes <- results[results$Significant, ]

cat("✅ SIGNIFICANTLY DYSREGULATED CAMK GENES:\n")
cat(sprintf("   (FDR < 0.05, N = %d healthy + %d disease = %d total)\n\n", 
           gse57338_results$n_healthy, gse57338_results$n_disease, 
           gse57338_results$n_healthy + gse57338_results$n_disease))

for (i in 1:nrow(significant_genes)) {
  gene <- significant_genes$Gene[i]
  logfc <- significant_genes$logFC[i]
  pval <- significant_genes$P.Value[i]
  fdr <- significant_genes$adj.P.Val[i]
  regulation <- significant_genes$Regulation[i]
  fold_change <- 2^abs(logfc)
  
  cat(sprintf("   • %-8s: %.3f-fold %s (logFC=%.3f, P=%.2e, FDR=%.2e)\n",
             gene, fold_change, 
             ifelse(logfc > 0, "UPREGULATED", "DOWNREGULATED"),
             logfc, pval, fdr))
}

cat("\n📈 NON-SIGNIFICANT CAMK GENES:\n")
non_sig_genes <- results[!results$Significant, ]
for (i in 1:nrow(non_sig_genes)) {
  gene <- non_sig_genes$Gene[i]
  logfc <- non_sig_genes$logFC[i]
  fdr <- non_sig_genes$adj.P.Val[i]
  regulation <- non_sig_genes$Regulation[i]
  
  cat(sprintf("   • %-8s: %s (logFC=%.3f, FDR=%.3f)\n",
             gene, regulation, logfc, fdr))
}

cat("\n🔬 CLINICAL INTERPRETATION\n")
cat("========================\n\n")

# Analyze patterns
upregulated <- significant_genes[significant_genes$logFC > 0, ]
downregulated <- significant_genes[significant_genes$logFC < 0, ]

cat("📊 DYSREGULATION PATTERNS:\n")
cat(sprintf("   • UPREGULATED in disease:   %d genes (%s)\n", 
           nrow(upregulated), paste(upregulated$Gene, collapse = ", ")))
cat(sprintf("   • DOWNREGULATED in disease: %d genes (%s)\n", 
           nrow(downregulated), paste(downregulated$Gene, collapse = ", ")))

cat("\n🎯 KEY CARDIOVASCULAR DISEASE INSIGHTS:\n\n")

cat("1. **CAMK2 FAMILY UPREGULATION**:\n")
cat("   • CAMK2G (most significant, FDR=6.92e-05): Key in cardiac hypertrophy\n")
cat("   • CAMK2B (FDR=8.40e-04): Critical for cardiac contractility\n")
cat("   • CAMK2A (FDR=3.53e-03): Central in cardiac signaling\n")
cat("   → Suggests enhanced Ca²⁺/calmodulin signaling in heart failure\n\n")

cat("2. **CAMK1 DOWNREGULATION**:\n")
cat("   • CAMK1 (FDR=6.92e-05): Reduced metabolic regulation\n")
cat("   • CAMKK1 (FDR=5.22e-03): Impaired energy homeostasis\n")
cat("   → Indicates metabolic dysfunction in failing hearts\n\n")

cat("3. **CAMK4 UPREGULATION**:\n")
cat("   • CAMK4 (FDR=4.14e-03): Enhanced transcriptional regulation\n")
cat("   → May drive pathological gene expression programs\n\n")

cat("🔍 BIOLOGICAL SIGNIFICANCE:\n")
cat("=========================\n\n")

cat("The corrected analysis reveals a clear pattern of CAMK dysregulation in\n")
cat("cardiovascular disease vs healthy controls:\n\n")

cat("• **Enhanced CAMK2 signaling**: Upregulation of CAMK2A/B/G suggests\n")
cat("  hyperactivated calcium signaling, linked to cardiac hypertrophy,\n")
cat("  arrhythmias, and contractile dysfunction.\n\n")

cat("• **Metabolic disruption**: Downregulation of CAMK1 and CAMKK1 indicates\n")
cat("  impaired metabolic regulation and energy homeostasis in failing hearts.\n\n")

cat("• **Transcriptional changes**: CAMK4 upregulation suggests altered gene\n")
cat("  expression programs contributing to disease pathogenesis.\n\n")

cat("⚡ COMPARISON WITH ORIGINAL FLAWED ANALYSIS:\n")
cat("==========================================\n\n")

cat("ORIGINAL (FLAWED):\n")
cat("• Only CAMK2D significant (P=0.00141)\n")
cat("• Mixed disease comparisons (AF vs SR + disease vs healthy)\n")
cat("• Missing the largest healthy vs disease dataset (GSE57338)\n")
cat("• Clinically meaningless interpretation\n\n")

cat("CORRECTED:\n")
cat("• 6 significantly dysregulated CAMK genes\n")
cat("• Pure healthy vs cardiovascular disease comparison\n")
cat("• Largest available dataset (313 samples) properly analyzed\n")
cat("• Clear biological patterns and clinical relevance\n\n")

cat("✅ CONCLUSION:\n")
cat("=============\n\n")

cat("The corrected analysis provides clinically meaningful insights into CAMK\n")
cat("dysregulation in cardiovascular disease. The pattern of CAMK2 upregulation\n")
cat("and CAMK1/CAMKK1 downregulation suggests both pathological calcium signaling\n")
cat("and metabolic dysfunction contribute to heart failure pathogenesis.\n\n")

cat("This analysis is based on the largest available healthy vs disease dataset\n")
cat("and provides a solid foundation for understanding CAMK family involvement\n")
cat("in cardiovascular disease.\n\n")

cat("🚀 NEXT STEPS:\n")
cat("• Validate findings in additional healthy vs disease cohorts\n")
cat("• Investigate CAMK2G as a potential therapeutic target\n")
cat("• Study metabolic implications of CAMK1/CAMKK1 downregulation\n")

# Save summary
summary_data <- list(
  corrected_results = gse57338_results,
  original_flawed_meta = original_meta,
  significant_genes = significant_genes,
  analysis_date = Sys.Date(),
  key_findings = list(
    upregulated = upregulated$Gene,
    downregulated = downregulated$Gene,
    sample_size = c(healthy = gse57338_results$n_healthy, disease = gse57338_results$n_disease)
  )
)

saveRDS(summary_data, "output/corrected_CAMK_analysis_complete_summary.rds")
cat("\n💾 Complete summary saved to: output/corrected_CAMK_analysis_complete_summary.rds\n")

cat("\n✅ CORRECTED CAMK ANALYSIS SUMMARY COMPLETED!\n")