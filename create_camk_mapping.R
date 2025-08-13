#!/usr/bin/env Rscript
#' Create CAMK Gene Mapping for Entrez IDs
#' 
#' Maps CAMK gene symbols to their corresponding Entrez IDs to identify them in existing results

library(org.Hs.eg.db)
library(AnnotationDbi)

cat("🧬 CREATING CAMK GENE MAPPING\n")
cat("=============================\n\n")

# Define CAMK family genes
camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMKK1", "CAMKK2", "CAMK1", "CAMK1D", "CAMK1G", "CAMK4", "CAMKV")

cat("🎯 CAMK genes to map:", length(camk_genes), "\n")
cat("   ", paste(camk_genes, collapse = ", "), "\n\n")

# Get Entrez IDs for CAMK genes
tryCatch({
  cat("🔍 Mapping gene symbols to Entrez IDs...\n")
  
  camk_entrez_mapping <- mapIds(org.Hs.eg.db,
                               keys = camk_genes,
                               column = "ENTREZID", 
                               keytype = "SYMBOL",
                               multiVals = "first")
  
  # Create comprehensive mapping table
  mapping_table <- data.frame(
    Gene_Symbol = camk_genes,
    Entrez_ID = camk_entrez_mapping,
    Found = !is.na(camk_entrez_mapping),
    stringsAsFactors = FALSE
  )
  
  cat("✅ Mapping completed!\n\n")
  
  # Display results
  cat("📋 CAMK Gene to Entrez ID Mapping:\n")
  cat("==================================\n")
  
  for (i in 1:nrow(mapping_table)) {
    gene <- mapping_table$Gene_Symbol[i]
    entrez <- mapping_table$Entrez_ID[i]
    status <- if (mapping_table$Found[i]) "✅" else "❌"
    
    cat(sprintf("%-8s: %s %s\n", gene, entrez, status))
  }
  
  # Summary
  found_count <- sum(mapping_table$Found)
  cat(sprintf("\n📊 Summary: %d/%d CAMK genes mapped to Entrez IDs\n", found_count, length(camk_genes)))
  
  # Save mapping table
  write.csv(mapping_table, "output/CAMK_entrez_mapping.csv", row.names = FALSE)
  saveRDS(mapping_table, "output/CAMK_entrez_mapping.rds")
  
  cat("💾 Mapping saved to:\n")
  cat("   • output/CAMK_entrez_mapping.csv\n")
  cat("   • output/CAMK_entrez_mapping.rds\n\n")
  
  # Now search for CAMK genes in the existing DGE results
  cat("🔎 Searching for CAMK genes in existing DGE results...\n")
  
  if (file.exists("output/GSE57338_healthy_vs_disease_DGE.csv")) {
    dge_results <- read.csv("output/GSE57338_healthy_vs_disease_DGE.csv", row.names = 1)
    
    cat("   DGE results loaded:", nrow(dge_results), "genes\n")
    
    # Get list of Entrez IDs to search for
    search_entrez_ids <- mapping_table$Entrez_ID[mapping_table$Found]
    
    # Search in row names (which should be Entrez IDs)
    camk_matches <- rownames(dge_results)[rownames(dge_results) %in% search_entrez_ids]
    
    if (length(camk_matches) > 0) {
      cat("   ✅ Found", length(camk_matches), "CAMK genes in DGE results!\n\n")
      
      # Extract CAMK results
      camk_dge_results <- dge_results[camk_matches, , drop = FALSE]
      
      # Add gene symbols
      camk_dge_results$Entrez_ID <- rownames(camk_dge_results)
      camk_dge_results$Gene_Symbol <- sapply(camk_dge_results$Entrez_ID, function(entrez_id) {
        idx <- which(mapping_table$Entrez_ID == entrez_id)
        if (length(idx) > 0) mapping_table$Gene_Symbol[idx[1]] else "Unknown"
      })
      
      # Reorder columns to put gene info first
      col_order <- c("Gene_Symbol", "Entrez_ID", setdiff(colnames(camk_dge_results), c("Gene_Symbol", "Entrez_ID")))
      camk_dge_results <- camk_dge_results[, col_order]
      
      # Sort by significance
      camk_dge_results <- camk_dge_results[order(camk_dge_results$adj.P.Val), ]
      
      cat("🎯 CAMK FAMILY DGE RESULTS (Healthy vs Disease):\n")
      cat("===============================================\n")
      
      for (i in 1:nrow(camk_dge_results)) {
        gene <- camk_dge_results$Gene_Symbol[i]
        entrez <- camk_dge_results$Entrez_ID[i]
        logfc <- round(camk_dge_results$logFC[i], 3)
        pval <- camk_dge_results$P.Value[i]
        adj_pval <- camk_dge_results$adj.P.Val[i]
        
        # Determine significance and direction
        sig_level <- if (adj_pval < 0.001) "***" else if (adj_pval < 0.01) "**" else if (adj_pval < 0.05) "*" else ""
        direction <- if (logfc > 0) "↑" else "↓"
        reg_status <- if (logfc > 0) "UP in Disease" else "DOWN in Disease"
        
        cat(sprintf("%-8s (%s): %s logFC=%6.3f, P=%8.2e, FDR=%8.2e %s [%s]\n", 
                    gene, entrez, direction, logfc, pval, adj_pval, sig_level, reg_status))
      }
      
      # Count significant CAMK genes
      sig_camk <- sum(camk_dge_results$adj.P.Val < 0.05)
      cat(sprintf("\n🎯 Significant CAMK genes (FDR < 0.05): %d/%d\n", sig_camk, nrow(camk_dge_results)))
      
      # Save CAMK-specific results
      write.csv(camk_dge_results, "output/CAMK_family_healthy_vs_disease_DGE.csv", row.names = FALSE)
      cat("💾 CAMK results saved to: output/CAMK_family_healthy_vs_disease_DGE.csv\n")
      
    } else {
      cat("   ❌ No CAMK genes found in DGE results\n")
      cat("   This might indicate a gene ID format mismatch\n")
      
      # Show examples of what's in the results
      cat("   Examples of gene IDs in results:", paste(head(rownames(dge_results), 5), collapse = ", "), "\n")
    }
  } else {
    cat("   ❌ DGE results file not found\n")
  }
  
}, error = function(e) {
  cat("❌ Error in mapping:", e$message, "\n")
})

cat("\n✅ CAMK mapping completed!\n")