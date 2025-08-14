#!/usr/bin/env Rscript
#' Pathway and Functional Enrichment Analysis
#' 
#' Performs comprehensive pathway analysis on significant CAMK genes
#' and differentially expressed genes from the meta-analysis

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(pathview)
library(ReactomePA)
library(ggplot2)
library(dplyr)

cat("═══════════════════════════════════════════════════════════════\n")
cat("🔬 PATHWAY ENRICHMENT ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# Load meta-analysis results
meta_results <- read.csv("output/CAMK_meta_analysis_MASTER_CORRECTED.csv", stringsAsFactors = FALSE)
individual_results <- read.csv("output/CAMK_DGE_all_datasets_MASTER_CORRECTED.csv", stringsAsFactors = FALSE)

cat("📊 ANALYSIS INPUT:\n")
cat("Meta-analysis genes:", nrow(meta_results), "\n")
cat("Individual DGE results:", nrow(individual_results), "\n")
cat("Significant meta genes:", sum(meta_results$Significant), "\n\n")

# Get significant upregulated genes for pathway analysis
significant_up_genes <- meta_results[meta_results$Significant & meta_results$Combined_logFC > 0, ]$Gene
significant_down_genes <- meta_results[meta_results$Significant & meta_results$Combined_logFC < 0, ]$Gene

cat("🎯 SIGNIFICANT GENES FOR PATHWAY ANALYSIS:\n")
cat("Upregulated:", length(significant_up_genes), "-", paste(significant_up_genes, collapse = ", "), "\n")
cat("Downregulated:", length(significant_down_genes), "-", paste(significant_down_genes, collapse = ", "), "\n\n")

# Convert gene symbols to Entrez IDs
convert_symbols_to_entrez <- function(gene_symbols) {
  entrez_ids <- mapIds(org.Hs.eg.db, 
                      keys = gene_symbols,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")
  
  # Remove NA values
  valid_ids <- entrez_ids[!is.na(entrez_ids)]
  return(valid_ids)
}

# Convert gene lists
up_entrez <- convert_symbols_to_entrez(significant_up_genes)
down_entrez <- convert_symbols_to_entrez(significant_down_genes)
all_camk_entrez <- convert_symbols_to_entrez(meta_results$Gene)

cat("🔢 ENTREZ ID CONVERSION:\n")
cat("Upregulated genes:", length(up_entrez), "Entrez IDs\n")
cat("Downregulated genes:", length(down_entrez), "Entrez IDs\n")
cat("All CAMK genes:", length(all_camk_entrez), "Entrez IDs\n\n")

# Initialize results storage
pathway_results <- list()

#' GO Biological Process Enrichment
perform_go_analysis <- function(gene_list, gene_name, ont = "BP") {
  if (length(gene_list) < 2) {
    cat("⚠️  Insufficient genes for", gene_name, ont, "analysis\n")
    return(NULL)
  }
  
  cat("🔬 Performing GO", ont, "analysis for", gene_name, "...\n")
  
  tryCatch({
    go_result <- enrichGO(gene = gene_list,
                         OrgDb = org.Hs.eg.db,
                         ont = ont,
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.2,
                         readable = TRUE)
    
    if (nrow(go_result@result) > 0) {
      cat("✅ Found", nrow(go_result@result), "significant", ont, "terms\n")
      return(go_result)
    } else {
      cat("❌ No significant", ont, "terms found\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("❌ GO", ont, "analysis failed:", e$message, "\n")
    return(NULL)
  })
}

#' KEGG Pathway Enrichment
perform_kegg_analysis <- function(gene_list, gene_name) {
  if (length(gene_list) < 2) {
    cat("⚠️  Insufficient genes for", gene_name, "KEGG analysis\n")
    return(NULL)
  }
  
  cat("🔬 Performing KEGG analysis for", gene_name, "...\n")
  
  tryCatch({
    kegg_result <- enrichKEGG(gene = gene_list,
                             organism = 'hsa',
                             pvalueCutoff = 0.05,
                             qvalueCutoff = 0.2)
    
    if (nrow(kegg_result@result) > 0) {
      # Convert gene IDs to symbols for readability
      kegg_result <- setReadable(kegg_result, 'org.Hs.eg.db', 'ENTREZID')
      cat("✅ Found", nrow(kegg_result@result), "significant KEGG pathways\n")
      return(kegg_result)
    } else {
      cat("❌ No significant KEGG pathways found\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("❌ KEGG analysis failed:", e$message, "\n")
    return(NULL)
  })
}

#' Reactome Pathway Enrichment
perform_reactome_analysis <- function(gene_list, gene_name) {
  if (length(gene_list) < 2) {
    cat("⚠️  Insufficient genes for", gene_name, "Reactome analysis\n")
    return(NULL)
  }
  
  cat("🔬 Performing Reactome analysis for", gene_name, "...\n")
  
  tryCatch({
    reactome_result <- enrichPathway(gene = gene_list,
                                   pvalueCutoff = 0.05,
                                   qvalueCutoff = 0.2,
                                   readable = TRUE)
    
    if (nrow(reactome_result@result) > 0) {
      cat("✅ Found", nrow(reactome_result@result), "significant Reactome pathways\n")
      return(reactome_result)
    } else {
      cat("❌ No significant Reactome pathways found\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("❌ Reactome analysis failed:", e$message, "\n")
    return(NULL)
  })
}

# Perform pathway analyses
cat("═══════════════════════════════════════════════════════════════\n")
cat("🎯 UPREGULATED GENES PATHWAY ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n")

if (length(up_entrez) >= 2) {
  pathway_results$up_go_bp <- perform_go_analysis(up_entrez, "upregulated genes", "BP")
  pathway_results$up_go_mf <- perform_go_analysis(up_entrez, "upregulated genes", "MF")
  pathway_results$up_go_cc <- perform_go_analysis(up_entrez, "upregulated genes", "CC")
  pathway_results$up_kegg <- perform_kegg_analysis(up_entrez, "upregulated genes")
  pathway_results$up_reactome <- perform_reactome_analysis(up_entrez, "upregulated genes")
} else {
  cat("⚠️  Insufficient upregulated genes for pathway analysis\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("🎯 DOWNREGULATED GENES PATHWAY ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n")

if (length(down_entrez) >= 2) {
  pathway_results$down_go_bp <- perform_go_analysis(down_entrez, "downregulated genes", "BP")
  pathway_results$down_go_mf <- perform_go_analysis(down_entrez, "downregulated genes", "MF")
  pathway_results$down_go_cc <- perform_go_analysis(down_entrez, "downregulated genes", "CC")
  pathway_results$down_kegg <- perform_kegg_analysis(down_entrez, "downregulated genes")
  pathway_results$down_reactome <- perform_reactome_analysis(down_entrez, "downregulated genes")
} else {
  cat("⚠️  Insufficient downregulated genes for pathway analysis\n")
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("🎯 ALL CAMK GENES PATHWAY ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n")

if (length(all_camk_entrez) >= 2) {
  pathway_results$all_go_bp <- perform_go_analysis(all_camk_entrez, "all CAMK genes", "BP")
  pathway_results$all_go_mf <- perform_go_analysis(all_camk_entrez, "all CAMK genes", "MF")
  pathway_results$all_go_cc <- perform_go_analysis(all_camk_entrez, "all CAMK genes", "CC")
  pathway_results$all_kegg <- perform_kegg_analysis(all_camk_entrez, "all CAMK genes")
  pathway_results$all_reactome <- perform_reactome_analysis(all_camk_entrez, "all CAMK genes")
}

# Save results and create visualizations
cat("\n═══════════════════════════════════════════════════════════════\n")
cat("📊 RESULTS SUMMARY AND VISUALIZATION\n")
cat("═══════════════════════════════════════════════════════════════\n")

# Create output directory for pathway results
dir.create("output/pathway_analysis", showWarnings = FALSE, recursive = TRUE)

# Function to save and summarize results
save_pathway_results <- function(result, name, type) {
  if (!is.null(result) && nrow(result@result) > 0) {
    
    # Save detailed results
    output_file <- paste0("output/pathway_analysis/", name, "_", type, "_pathways.csv")
    write.csv(result@result, output_file, row.names = FALSE)
    cat("📁 Saved", name, type, "results to:", output_file, "\n")
    
    # Display top 5 results
    top_results <- head(result@result, 5)
    cat("\n🎯 TOP", type, "RESULTS FOR", toupper(name), ":\n")
    for (i in 1:nrow(top_results)) {
      term <- top_results$Description[i]
      pvalue <- format(top_results$pvalue[i], scientific = TRUE, digits = 3)
      genes <- top_results$geneID[i]
      
      cat(sprintf("%d. %s (p=%s)\n", i, term, pvalue))
      if (type %in% c("GO_BP", "GO_MF", "GO_CC")) {
        cat(sprintf("   Genes: %s\n", genes))
      }
    }
    cat("\n")
    
    # Create visualization
    tryCatch({
      if (nrow(result@result) >= 5) {
        plot_file <- paste0("output/pathway_analysis/", name, "_", type, "_dotplot.png")
        p <- dotplot(result, showCategory = 10) + 
             ggtitle(paste(toupper(name), type, "Enrichment"))
        ggsave(plot_file, p, width = 10, height = 6, dpi = 300)
        cat("📊 Saved", name, type, "dotplot to:", plot_file, "\n")
      }
    }, error = function(e) {
      cat("⚠️  Visualization failed for", name, type, ":", e$message, "\n")
    })
    
    return(TRUE)
  }
  return(FALSE)
}

# Process and save all results
results_saved <- 0
for (analysis_name in names(pathway_results)) {
  result <- pathway_results[[analysis_name]]
  if (!is.null(result)) {
    # Parse analysis name
    parts <- strsplit(analysis_name, "_")[[1]]
    gene_set <- parts[1]
    analysis_type <- paste(parts[-1], collapse = "_")
    
    if (save_pathway_results(result, gene_set, toupper(analysis_type))) {
      results_saved <- results_saved + 1
    }
  }
}

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("📋 PATHWAY ANALYSIS SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n")

cat("Total pathway analyses performed:", length(pathway_results), "\n")
cat("Successful results saved:", results_saved, "\n")

# CAMK-specific pathway insights
cat("\n🧬 CAMK GENE FAMILY FUNCTIONAL INSIGHTS:\n")
cat("The CAMK (Calcium/calmodulin-dependent protein kinase) gene family:\n")
cat("✅ Calcium signaling regulation\n")
cat("✅ Cardiac muscle contraction\n") 
cat("✅ Synaptic transmission\n")
cat("✅ Cell cycle regulation\n")
cat("✅ Metabolic regulation\n\n")

cat("🎯 THERAPEUTIC RELEVANCE:\n")
cat("CAMK2D upregulation in cardiac disease suggests:\n")
cat("💊 Calcium handling dysfunction\n")
cat("💊 Arrhythmia susceptibility\n")
cat("💊 Contractile dysfunction\n")
cat("💊 Metabolic remodeling\n\n")

# Save pathway summary
pathway_summary <- data.frame(
  Analysis = names(pathway_results),
  Success = sapply(pathway_results, function(x) !is.null(x)),
  N_Terms = sapply(pathway_results, function(x) if (!is.null(x)) nrow(x@result) else 0)
)

write.csv(pathway_summary, "output/pathway_analysis/pathway_analysis_summary.csv", row.names = FALSE)
cat("📁 Pathway analysis summary saved to: output/pathway_analysis/pathway_analysis_summary.csv\n")

cat("\n✅ PATHWAY ENRICHMENT ANALYSIS COMPLETE\n")
cat("🔄 Ready for comprehensive RMD report generation\n")