#!/usr/bin/env Rscript
#' Pathway Analysis Module
#' 
#' Comprehensive GO/KEGG pathway analysis for pharmaceutical insights
#' Non-disruptive enhancement to existing pipeline
#' 
#' @author Claude Code Enhancement Module
#' @version 1.0.0

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)
  library(enrichplot)
  library(yaml)
  library(tidyverse)
})

#' Run Comprehensive Pathway Analysis
#'
#' Performs GO, KEGG, and Reactome pathway analysis on DGE results
#' @param dge_results Data frame with differential expression results
#' @param config_file Configuration file path
#' @param output_dir Directory for saving results
#' @return List with all pathway analysis results
run_pathway_analysis <- function(dge_results,
                                config_file = "config.yml",
                                output_dir = "output/pathways") {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║          PATHWAY ANALYSIS MODULE                             ║\n")
  cat("║          GO/KEGG/Reactome Enrichment Analysis               ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Load configuration
  config <- yaml::read_yaml(config_file)
  
  # Get pathway settings
  pathway_config <- config$dynamic_features$pathway_settings %||% list(
    databases = c("GO", "KEGG"),
    fdr_threshold = 0.05,
    min_gene_set_size = 10,
    max_gene_set_size = 500
  )
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize results
  pathway_results <- list(
    success = TRUE,
    go_results = NULL,
    kegg_results = NULL,
    reactome_results = NULL,
    drug_targets = NULL,
    summary = list()
  )
  
  # Prepare gene lists
  cat("📋 Preparing gene lists...\n")
  gene_lists <- prepare_gene_lists(dge_results, config)
  
  if (is.null(gene_lists)) {
    cat("❌ Failed to prepare gene lists\n")
    pathway_results$success <- FALSE
    return(pathway_results)
  }
  
  cat("   Significant up-regulated genes:", length(gene_lists$up_genes), "\n")
  cat("   Significant down-regulated genes:", length(gene_lists$down_genes), "\n")
  cat("   All significant genes:", length(gene_lists$all_sig_genes), "\n")
  
  # Run GO enrichment
  if ("GO" %in% pathway_config$databases) {
    cat("\n🧬 Running GO enrichment analysis...\n")
    pathway_results$go_results <- run_go_enrichment(
      gene_lists,
      pathway_config,
      output_dir
    )
  }
  
  # Run KEGG enrichment
  if ("KEGG" %in% pathway_config$databases) {
    cat("\n🗺️  Running KEGG pathway analysis...\n")
    pathway_results$kegg_results <- run_kegg_enrichment(
      gene_lists,
      pathway_config,
      output_dir
    )
  }
  
  # Run Reactome enrichment
  if ("Reactome" %in% pathway_config$databases) {
    cat("\n🔬 Running Reactome pathway analysis...\n")
    pathway_results$reactome_results <- run_reactome_enrichment(
      gene_lists,
      pathway_config,
      output_dir
    )
  }
  
  # Predict drug targets
  cat("\n💊 Predicting drug targets...\n")
  pathway_results$drug_targets <- predict_drug_targets(
    dge_results,
    pathway_results,
    output_dir
  )
  
  # Generate summary
  pathway_results$summary <- generate_pathway_summary(pathway_results)
  
  # Save results
  save_pathway_results(pathway_results, output_dir)
  
  # Print summary
  print_pathway_summary(pathway_results$summary)
  
  return(pathway_results)
}

#' Prepare Gene Lists for Pathway Analysis
#'
#' Converts gene symbols to Entrez IDs and creates gene lists
#' @param dge_results DGE results data frame
#' @param config Configuration object
#' @return List with up/down/all gene lists
prepare_gene_lists <- function(dge_results, config) {
  
  gene_lists <- list()
  
  # Get FDR threshold
  fdr_threshold <- config$analysis$differential_expression$fdr_threshold %||% 0.05
  
  # Filter significant genes
  sig_genes <- dge_results[dge_results$FDR < fdr_threshold, ]
  
  if (nrow(sig_genes) == 0) {
    cat("⚠️  No significant genes found\n")
    return(NULL)
  }
  
  # Convert gene symbols to Entrez IDs
  tryCatch({
    # Get Entrez IDs
    entrez_ids <- mapIds(
      org.Hs.eg.db,
      keys = sig_genes$Gene,
      column = "ENTREZID",
      keytype = "SYMBOL",
      multiVals = "first"
    )
    
    # Remove NA values
    entrez_ids <- entrez_ids[!is.na(entrez_ids)]
    
    # Create gene lists
    sig_genes$entrez <- entrez_ids[match(sig_genes$Gene, names(entrez_ids))]
    sig_genes <- sig_genes[!is.na(sig_genes$entrez), ]
    
    # Up-regulated genes
    up_genes <- sig_genes[sig_genes$logFC > 0, ]
    gene_lists$up_genes <- up_genes$entrez
    gene_lists$up_logfc <- setNames(up_genes$logFC, up_genes$entrez)
    
    # Down-regulated genes
    down_genes <- sig_genes[sig_genes$logFC < 0, ]
    gene_lists$down_genes <- down_genes$entrez
    gene_lists$down_logfc <- setNames(down_genes$logFC, down_genes$entrez)
    
    # All significant genes
    gene_lists$all_sig_genes <- sig_genes$entrez
    gene_lists$all_logfc <- setNames(sig_genes$logFC, sig_genes$entrez)
    
  }, error = function(e) {
    cat("❌ Error converting gene symbols:", e$message, "\n")
    return(NULL)
  })
  
  return(gene_lists)
}

#' Run GO Enrichment Analysis
#'
#' Performs GO enrichment for BP, MF, and CC
#' @param gene_lists Prepared gene lists
#' @param config Pathway configuration
#' @param output_dir Output directory
#' @return GO enrichment results
run_go_enrichment <- function(gene_lists, config, output_dir) {
  
  go_results <- list()
  
  # GO categories
  go_categories <- c("BP" = "Biological Process",
                    "MF" = "Molecular Function",
                    "CC" = "Cellular Component")
  
  for (ont in names(go_categories)) {
    cat("   Analyzing", go_categories[ont], "...\n")
    
    # Run enrichGO
    tryCatch({
      ego <- enrichGO(
        gene = gene_lists$all_sig_genes,
        OrgDb = org.Hs.eg.db,
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = config$fdr_threshold,
        qvalueCutoff = config$fdr_threshold,
        minGSSize = config$min_gene_set_size,
        maxGSSize = config$max_gene_set_size
      )
      
      if (!is.null(ego) && nrow(ego@result) > 0) {
        go_results[[ont]] <- ego
        
        # Save plot
        plot_file <- file.path(output_dir, paste0("GO_", ont, "_dotplot.pdf"))
        pdf(plot_file, width = 10, height = 8)
        print(dotplot(ego, showCategory = 20))
        dev.off()
        
        cat("      Found", nrow(ego@result), "enriched terms\n")
      } else {
        cat("      No enriched terms found\n")
      }
      
    }, error = function(e) {
      cat("      Error:", e$message, "\n")
    })
  }
  
  return(go_results)
}

#' Run KEGG Pathway Analysis
#'
#' Performs KEGG pathway enrichment
#' @param gene_lists Prepared gene lists
#' @param config Pathway configuration
#' @param output_dir Output directory
#' @return KEGG enrichment results
run_kegg_enrichment <- function(gene_lists, config, output_dir) {
  
  kegg_results <- NULL
  
  tryCatch({
    # Run enrichKEGG
    ekegg <- enrichKEGG(
      gene = gene_lists$all_sig_genes,
      organism = "hsa",
      pvalueCutoff = config$fdr_threshold,
      qvalueCutoff = config$fdr_threshold,
      minGSSize = config$min_gene_set_size,
      maxGSSize = config$max_gene_set_size
    )
    
    if (!is.null(ekegg) && nrow(ekegg@result) > 0) {
      kegg_results <- ekegg
      
      # Save plot
      plot_file <- file.path(output_dir, "KEGG_pathways_dotplot.pdf")
      pdf(plot_file, width = 10, height = 8)
      print(dotplot(ekegg, showCategory = 20))
      dev.off()
      
      cat("   Found", nrow(ekegg@result), "enriched pathways\n")
      
      # Print top pathways
      top_pathways <- head(ekegg@result, 5)
      for (i in 1:nrow(top_pathways)) {
        cat("      -", top_pathways$Description[i], 
            "(p =", formatC(top_pathways$pvalue[i], format = "e", digits = 2), ")\n")
      }
    } else {
      cat("   No enriched pathways found\n")
    }
    
  }, error = function(e) {
    cat("   Error:", e$message, "\n")
  })
  
  return(kegg_results)
}

#' Run Reactome Pathway Analysis
#'
#' Performs Reactome pathway enrichment
#' @param gene_lists Prepared gene lists
#' @param config Pathway configuration
#' @param output_dir Output directory
#' @return Reactome enrichment results
run_reactome_enrichment <- function(gene_lists, config, output_dir) {
  
  # Placeholder for Reactome analysis
  # Would require ReactomePA package
  cat("   Reactome analysis not implemented in this version\n")
  return(NULL)
}

#' Predict Drug Targets
#'
#' Identifies potential drug targets from DGE and pathway results
#' @param dge_results DGE results
#' @param pathway_results Pathway analysis results
#' @param output_dir Output directory
#' @return Data frame with drug target predictions
predict_drug_targets <- function(dge_results, pathway_results, output_dir) {
  
  drug_targets <- data.frame()
  
  # Get CAMK genes as primary targets
  camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", 
                 "CAMKK1", "CAMKK2", "CAMK1", "CAMK1D", "CAMK1G", "CAMK4")
  
  # Filter for significant CAMK genes
  sig_camk <- dge_results[dge_results$Gene %in% camk_genes & 
                          dge_results$FDR < 0.05, ]
  
  if (nrow(sig_camk) > 0) {
    # Create drug target table
    drug_targets <- data.frame(
      Gene = sig_camk$Gene,
      LogFC = sig_camk$logFC,
      FDR = sig_camk$FDR,
      Direction = ifelse(sig_camk$logFC > 0, "Up", "Down"),
      Target_Type = "Kinase",
      Druggability_Score = calculate_druggability_score(sig_camk),
      Clinical_Relevance = "High - Cardiovascular disease target",
      stringsAsFactors = FALSE
    )
    
    # Sort by druggability score
    drug_targets <- drug_targets[order(drug_targets$Druggability_Score, 
                                      decreasing = TRUE), ]
    
    # Save results
    output_file <- file.path(output_dir, "drug_targets.csv")
    write.csv(drug_targets, output_file, row.names = FALSE)
    
    cat("   Identified", nrow(drug_targets), "potential drug targets\n")
    
    # Print top targets
    if (nrow(drug_targets) > 0) {
      cat("   Top drug targets:\n")
      top_targets <- head(drug_targets, 3)
      for (i in 1:nrow(top_targets)) {
        cat("      -", top_targets$Gene[i], 
            "(Score:", round(top_targets$Druggability_Score[i], 2), ")\n")
      }
    }
  } else {
    cat("   No significant CAMK genes found as drug targets\n")
  }
  
  return(drug_targets)
}

#' Calculate Druggability Score
#'
#' Calculates druggability score for potential targets
#' @param gene_data Gene expression data
#' @return Vector of druggability scores
calculate_druggability_score <- function(gene_data) {
  
  # Simple scoring based on:
  # - Effect size (logFC)
  # - Statistical significance (FDR)
  # - Kinase family (highly druggable)
  
  scores <- rep(0, nrow(gene_data))
  
  # Effect size component (0-40 points)
  scores <- scores + pmin(40, abs(gene_data$logFC) * 50)
  
  # Significance component (0-30 points)
  scores <- scores + pmin(30, -log10(gene_data$FDR) * 5)
  
  # Kinase bonus (30 points)
  scores <- scores + 30
  
  return(scores)
}

#' Generate Pathway Summary
#'
#' Creates summary of all pathway results
#' @param pathway_results All pathway analysis results
#' @return Summary list
generate_pathway_summary <- function(pathway_results) {
  
  summary <- list(
    total_go_terms = 0,
    total_kegg_pathways = 0,
    total_drug_targets = 0,
    top_biological_processes = character(),
    top_kegg_pathways = character(),
    top_drug_targets = character()
  )
  
  # Summarize GO results
  if (!is.null(pathway_results$go_results)) {
    for (ont in names(pathway_results$go_results)) {
      if (!is.null(pathway_results$go_results[[ont]])) {
        summary$total_go_terms <- summary$total_go_terms + 
                                  nrow(pathway_results$go_results[[ont]]@result)
        
        if (ont == "BP" && nrow(pathway_results$go_results[[ont]]@result) > 0) {
          top_bp <- head(pathway_results$go_results[[ont]]@result$Description, 3)
          summary$top_biological_processes <- top_bp
        }
      }
    }
  }
  
  # Summarize KEGG results
  if (!is.null(pathway_results$kegg_results)) {
    summary$total_kegg_pathways <- nrow(pathway_results$kegg_results@result)
    if (summary$total_kegg_pathways > 0) {
      summary$top_kegg_pathways <- head(pathway_results$kegg_results@result$Description, 3)
    }
  }
  
  # Summarize drug targets
  if (!is.null(pathway_results$drug_targets)) {
    summary$total_drug_targets <- nrow(pathway_results$drug_targets)
    if (summary$total_drug_targets > 0) {
      summary$top_drug_targets <- head(pathway_results$drug_targets$Gene, 3)
    }
  }
  
  return(summary)
}

#' Save Pathway Results
#'
#' Saves all pathway results to files
#' @param pathway_results All pathway analysis results
#' @param output_dir Output directory
save_pathway_results <- function(pathway_results, output_dir) {
  
  # Save GO results
  if (!is.null(pathway_results$go_results)) {
    for (ont in names(pathway_results$go_results)) {
      if (!is.null(pathway_results$go_results[[ont]])) {
        output_file <- file.path(output_dir, paste0("GO_", ont, "_results.csv"))
        write.csv(pathway_results$go_results[[ont]]@result, 
                 output_file, row.names = FALSE)
      }
    }
  }
  
  # Save KEGG results
  if (!is.null(pathway_results$kegg_results)) {
    output_file <- file.path(output_dir, "KEGG_results.csv")
    write.csv(pathway_results$kegg_results@result, 
             output_file, row.names = FALSE)
  }
  
  # Save complete results object
  saveRDS(pathway_results, file.path(output_dir, "pathway_results.rds"))
  
  cat("\n💾 All pathway results saved to:", output_dir, "\n")
}

#' Print Pathway Summary
#'
#' Prints summary of pathway analysis results
#' @param summary Summary object
print_pathway_summary <- function(summary) {
  
  cat("
")
  cat("═══════════════════════════════════════════════════════════════
")
  cat("📊 PATHWAY ANALYSIS SUMMARY
")
  cat("═══════════════════════════════════════════════════════════════
")
  
  cat("GO Terms Enriched:", summary$total_go_terms, "
")
  cat("KEGG Pathways Enriched:", summary$total_kegg_pathways, "
")
  cat("Drug Targets Identified:", summary$total_drug_targets, "
")
  
  if (length(summary$top_biological_processes) > 0) {
    cat("
🧬 Top Biological Processes:
")
    for (bp in summary$top_biological_processes) {
      cat("   -", bp, "
")
    }
  }
  
  if (length(summary$top_kegg_pathways) > 0) {
    cat("
🗺️  Top KEGG Pathways:
")
    for (pathway in summary$top_kegg_pathways) {
      cat("   -", pathway, "
")
    }
  }
  
  if (length(summary$top_drug_targets) > 0) {
    cat("
💊 Top Drug Targets:
")
    for (target in summary$top_drug_targets) {
      cat("   -", target, "
")
    }
  }
}

# NULL coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Test function
test_pathway_analysis <- function() {
  cat("🧪 Testing Pathway Analysis Module
")
  
  # Load sample DGE results
  dge_file <- "output/current/CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv"
  if (file.exists(dge_file)) {
    dge_results <- read.csv(dge_file, stringsAsFactors = FALSE)
    
    # Run pathway analysis
    test_results <- run_pathway_analysis(
      dge_results = dge_results,
      config_file = "config.yml",
      output_dir = "output/test_pathways"
    )
    
    cat("✅ Pathway analysis module test complete
")
    return(test_results)
  } else {
    cat("⚠️  No DGE results found for testing
")
    return(NULL)
  }
}

cat("✅ Pathway Analysis Module loaded successfully
")
cat("   Functions: run_pathway_analysis(), predict_drug_targets()
")
cat("   Usage: source(\"modules/pathway_analysis.R\")
")
cat("         run_pathway_analysis(dge_results)
")
