#!/usr/bin/env Rscript
#' Comprehensive Dataset Inventory and Classification
#' 
#' Documents which datasets have healthy controls vs disease subtypes vs unclear groups

cat("📋 COMPREHENSIVE DATASET INVENTORY\n")
cat("=================================\n\n")

# Create comprehensive inventory of all datasets and their comparison types
datasets_inventory <- data.frame(
  Dataset_ID = c("GSE57338", "GSE14975", "GSE31821", "GSE41177", "GSE79768", "GSE120895"),
  Total_Samples = c(313, 10, 6, 38, 26, 55),
  Comparison_Type = c(
    "Healthy vs Disease",
    "Healthy vs Disease", 
    "AF vs SR",
    "AF vs SR", 
    "AF vs SR",
    "No clear groups"
  ),
  Group_1 = c("136 Healthy", "5 Healthy", "4 AF", "32 AF", "14 AF", "Unknown"),
  Group_2 = c("177 Disease", "5 Disease", "2 SR", "6 SR", "12 SR", "Unknown"),
  Disease_Types = c(
    "Ischemic + Dilated CMP",
    "Heart failure", 
    "Atrial fibrillation",
    "Atrial fibrillation",
    "Atrial fibrillation", 
    "Mixed cardiovascular"
  ),
  Clinical_Relevance = c(
    "HIGH - True disease dysregulation",
    "MODERATE - Small sample size",
    "MODERATE - Disease subtype comparison", 
    "MODERATE - Disease subtype comparison",
    "MODERATE - Disease subtype comparison",
    "LOW - Unclear groups"
  ),
  Statistical_Power = c(
    "EXCELLENT - Large, well-powered",
    "LOW - Only 10 samples total",
    "LOW - Only 6 samples total",
    "MODERATE - 38 samples",
    "MODERATE - 26 samples", 
    "UNKNOWN - Groups unclear"
  ),
  CAMK_Analysis_Suitable = c(
    "YES - Primary dataset",
    "NO - Too small",
    "SECONDARY - Different question", 
    "SECONDARY - Different question",
    "SECONDARY - Different question",
    "NO - No clear groups"
  ),
  Notes = c(
    "Gold standard: 136 healthy vs 177 disease, all CAMK genes detected",
    "Healthy vs disease but insufficient power for reliable results",
    "Compares two disease states, not healthy vs disease",
    "Compares two disease states, not healthy vs disease", 
    "Compares two disease states, not healthy vs disease",
    "Group structure unclear, excluded from analysis"
  )
)

cat("📊 DATASET CLASSIFICATION SUMMARY\n")
cat("================================\n\n")

# Print formatted inventory
for (i in 1:nrow(datasets_inventory)) {
  dataset <- datasets_inventory[i, ]
  
  cat("📋 Dataset:", dataset$Dataset_ID, "\n")
  cat("   ", paste(rep("=", nchar(dataset$Dataset_ID) + 9), collapse = ""), "\n")
  cat("   📊 Total samples:", dataset$Total_Samples, "\n")
  cat("   🔬 Comparison type:", dataset$Comparison_Type, "\n")
  cat("   👥 Groups:", dataset$Group_1, "vs", dataset$Group_2, "\n")
  cat("   🏥 Disease context:", dataset$Disease_Types, "\n")
  cat("   📈 Clinical relevance:", dataset$Clinical_Relevance, "\n")
  cat("   📊 Statistical power:", dataset$Statistical_Power, "\n")
  cat("   🧬 CAMK analysis:", dataset$CAMK_Analysis_Suitable, "\n")
  cat("   💡 Notes:", dataset$Notes, "\n\n")
}

# Categorize datasets by analysis type
healthy_vs_disease <- datasets_inventory[datasets_inventory$Comparison_Type == "Healthy vs Disease", ]
disease_subtypes <- datasets_inventory[datasets_inventory$Comparison_Type == "AF vs SR", ]
excluded <- datasets_inventory[!datasets_inventory$Comparison_Type %in% c("Healthy vs Disease", "AF vs SR"), ]

cat("🎯 ANALYSIS CATEGORIES\n")
cat("=====================\n\n")

cat("✅ PRIMARY ANALYSIS: Healthy vs Disease\n")
cat("---------------------------------------\n")
if (nrow(healthy_vs_disease) > 0) {
  for (i in 1:nrow(healthy_vs_disease)) {
    ds <- healthy_vs_disease[i, ]
    power_symbol <- ifelse(ds$Statistical_Power == "EXCELLENT - Large, well-powered", "🟢", 
                          ifelse(grepl("LOW", ds$Statistical_Power), "🔴", "🟡"))
    cat(sprintf("   %s %-12s: %3d samples (%s) %s\n", 
               power_symbol, ds$Dataset_ID, ds$Total_Samples, ds$Comparison_Type, 
               ifelse(ds$CAMK_Analysis_Suitable == "YES - Primary dataset", "⭐ PRIMARY", "")))
  }
} else {
  cat("   ❌ No healthy vs disease datasets found\n")
}

cat("\n🟡 SECONDARY ANALYSIS: Disease Subtypes\n")
cat("--------------------------------------\n")
if (nrow(disease_subtypes) > 0) {
  for (i in 1:nrow(disease_subtypes)) {
    ds <- disease_subtypes[i, ]
    cat(sprintf("   🟡 %-12s: %3d samples (%s)\n", 
               ds$Dataset_ID, ds$Total_Samples, ds$Comparison_Type))
  }
  cat("   💡 These compare disease states, not healthy vs disease\n")
  cat("   💡 Useful for understanding disease heterogeneity\n")
} else {
  cat("   No disease subtype datasets found\n")
}

cat("\n❌ EXCLUDED FROM ANALYSIS\n")
cat("-------------------------\n")
if (nrow(excluded) > 0) {
  for (i in 1:nrow(excluded)) {
    ds <- excluded[i, ]
    cat(sprintf("   ❌ %-12s: %s\n", ds$Dataset_ID, ds$Notes))
  }
} else {
  cat("   No excluded datasets\n")
}

cat("\n🔍 KEY INSIGHTS\n")
cat("==============\n\n")

cat("1. **ONLY ONE ROBUST HEALTHY vs DISEASE DATASET**:\n")
cat("   • GSE57338 (313 samples) is the ONLY dataset with adequate sample size\n")
cat("   • GSE14975 too small (10 samples) for reliable analysis\n")
cat("   • This explains why GSE57338 alone provides better results than mixed meta-analysis\n\n")

cat("2. **MOST DATASETS COMPARE DISEASE SUBTYPES**:\n")
cat("   • GSE31821, GSE41177, GSE79768 all compare AF vs SR\n")
cat("   • These are disease-vs-disease comparisons, not healthy-vs-disease\n")
cat("   • Mixing these with healthy-vs-disease creates meaningless meta-analysis\n\n")

cat("3. **WHY THE ORIGINAL META-ANALYSIS WAS FLAWED**:\n")
cat("   • Mixed incompatible comparison types (healthy-vs-disease + AF-vs-SR)\n")
cat("   • Missing the largest healthy-vs-disease dataset (GSE57338)\n")
cat("   • Created clinically meaningless 'average' across different questions\n\n")

cat("4. **WHY GSE57338 ANALYSIS IS SUPERIOR**:\n")
cat("   • Pure healthy vs cardiovascular disease comparison\n")
cat("   • Large sample size (313 samples)\n")
cat("   • All 11 CAMK genes detected and analyzed\n")
cat("   • Clinically interpretable results\n\n")

# Calculate summary statistics
total_samples <- sum(datasets_inventory$Total_Samples, na.rm = TRUE)
healthy_disease_samples <- sum(healthy_vs_disease$Total_Samples, na.rm = TRUE)
disease_subtype_samples <- sum(disease_subtypes$Total_Samples, na.rm = TRUE)

cat("📊 SAMPLE SIZE SUMMARY\n")
cat("=====================\n")
cat(sprintf("   Total samples across all datasets: %d\n", total_samples))
cat(sprintf("   Healthy vs disease samples: %d (%.1f%%)\n", 
           healthy_disease_samples, 100 * healthy_disease_samples / total_samples))
cat(sprintf("   Disease subtype samples: %d (%.1f%%)\n", 
           disease_subtype_samples, 100 * disease_subtype_samples / total_samples))
cat(sprintf("   GSE57338 alone: %d samples (%.1f%% of all data)\n", 
           313, 100 * 313 / total_samples))

cat("\n💡 CONCLUSION\n")
cat("=============\n\n")
cat("The user's question reveals an important misconception. NOT all datasets\n")
cat("have healthy controls. In fact:\n\n")
cat("• Only 1 dataset (GSE57338) provides robust healthy vs disease comparison\n")
cat("• 3 datasets compare disease subtypes (AF vs SR)\n") 
cat("• 1 dataset is too small for meaningful analysis\n")
cat("• 1 dataset has unclear group structure\n\n")
cat("This is why the corrected analysis focusing on GSE57338 alone provides\n")
cat("more meaningful and clinically relevant results than the original flawed\n") 
cat("meta-analysis that mixed incompatible comparison types.\n")

# Save inventory
write.csv(datasets_inventory, "output/comprehensive_dataset_inventory.csv", row.names = FALSE)
saveRDS(datasets_inventory, "output/comprehensive_dataset_inventory.rds")

cat("\n💾 Dataset inventory saved to:\n")
cat("   • output/comprehensive_dataset_inventory.csv\n") 
cat("   • output/comprehensive_dataset_inventory.rds\n")

cat("\n✅ Comprehensive dataset inventory completed!\n")