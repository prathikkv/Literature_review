#!/usr/bin/env Rscript
#' Update Final DGE Results with Corrected GSE57338 Data

# Load results
original_results <- read.csv("output/CAMK_focused_DGE_all_datasets_FINAL.csv")
corrected_gse57338 <- read.csv("output/GSE57338_corrected_CAMK_results.csv")

cat("=== UPDATING FINAL DGE RESULTS ===\n")
cat("Original entries:", nrow(original_results), "\n")
cat("GSE57338 entries to replace:", sum(original_results$Dataset == "GSE57338"), "\n")
cat("New GSE57338 entries:", nrow(corrected_gse57338), "\n")

# Remove original GSE57338 results
non_gse57338 <- original_results[original_results$Dataset != "GSE57338", ]

# Ensure column compatibility
corrected_subset <- corrected_gse57338[, colnames(original_results)]

# Combine results
updated_results <- rbind(non_gse57338, corrected_subset)

# Save updated results
write.csv(updated_results, "output/CAMK_focused_DGE_all_datasets_CORRECTED.csv", row.names = FALSE)

cat("Updated entries:", nrow(updated_results), "\n")

# Validate CAMK2D correction
camk2d_original <- original_results[original_results$Gene_Symbol == "CAMK2D" & original_results$Dataset == "GSE57338", ]
camk2d_corrected <- updated_results[updated_results$Gene_Symbol == "CAMK2D" & updated_results$Dataset == "GSE57338", ]

cat("\n=== CAMK2D VALIDATION ===\n")
if (nrow(camk2d_original) > 0 && nrow(camk2d_corrected) > 0) {
  cat("Original CAMK2D:", round(camk2d_original$logFC, 4), "(", camk2d_original$Regulation, ")\n")
  cat("Corrected CAMK2D:", round(camk2d_corrected$logFC, 4), "(", camk2d_corrected$Regulation, ")\n")
  
  if (camk2d_original$logFC < 0 && camk2d_corrected$logFC > 0) {
    cat("✅ SUCCESS: CAMK2D direction corrected (DOWN → UP)\n")
  }
}

cat("\nUpdated file saved: output/CAMK_focused_DGE_all_datasets_CORRECTED.csv\n")