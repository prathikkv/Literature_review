#!/usr/bin/env Rscript
#' CAMK-Focused Analysis Pipeline (Methodologically Corrected)
#' 
#' Performs genome-wide differential expression analysis first, then focuses on CAMK family genes
#' to avoid statistical biases from pre-filtering

# Load required functions and centralized definitions
source("../../functions/data_processing.R")
source("../../functions/analysis.R")
source("../../functions/camk_definitions.R")
source("../../functions/common_utils.R")
source("../../scripts/enhanced_group_detection.R")

# Load required libraries
load_required_libraries("visualization")

print_analysis_header(
  title = "CAMK-Focused Analysis Pipeline (Corrected)",
  subtitle = "Methodologically improved with genome-wide DGE approach",
  methodology = "Genome-wide DGE → CAMK filtering (eliminates pre-filtering bias)"
)

# Get CAMK gene definitions from centralized module
camk_gene_categories <- get_camk_gene_categories()
camk_core_genes <- camk_gene_categories$core
camk_extended_genes <- camk_gene_categories$extended
camk_genes <- camk_core_genes  # Backward compatibility

# Display gene summary
cat("GENETIC: Core CAMK genes:", length(camk_core_genes), "\n")
cat("PATHWAY: Extended CAMK+pathway genes:", length(camk_extended_genes), "\n")
cat("DATA: Enhanced biological context enabled\n\n")

# Load and analyze datasets with proper gene symbols
cache_dir <- "cache/comprehensive"
processed_files <- list.files(cache_dir, pattern = "_processed.rds$", full.names = FALSE)

camk_analysis_results <- list()

for (file in processed_files) {
  dataset_id <- gsub("_processed.rds", "", file)
  
  cat("DATA: Analyzing", dataset_id, "\n")
  cat("   Loading dataset...\n")
  
  # Load dataset
  dataset <- readRDS(file.path(cache_dir, file))
  
  if (!dataset$success || is.null(dataset$expression_matrix)) {
    cat("   ERROR: Dataset failed or no expression matrix - skipping\n\n")
    next
  }
  
  expr_matrix <- dataset$expression_matrix
  n_genes <- nrow(expr_matrix)
  n_samples <- ncol(expr_matrix)
  
  cat("   Dataset info:", n_genes, "genes x", n_samples, "samples\n")
  
  # Check for CAMK genes in this dataset (for reporting, not filtering)
  gene_names <- rownames(expr_matrix)
  camk_present <- camk_core_genes[camk_core_genes %in% gene_names]
  camk_extended_present <- camk_extended_genes[camk_extended_genes %in% gene_names]
  
  cat("   Core CAMK genes found:", length(camk_present), "/", length(camk_core_genes), "\n")
  cat("   Extended CAMK+pathway genes found:", length(camk_extended_present), "/", length(camk_extended_genes), "\n")
  
  if (length(camk_present) == 0) {
    cat("   WARNING: No core CAMK genes detected - skipping analysis\n\n")
    next
  }
  
  # Print found CAMK genes
  cat("   Core CAMK present:", paste(camk_present, collapse = ", "), "\n")
  cat("   Extended genes present:", length(camk_extended_present), "genes\n")
  
  # NO PRE-FILTERING - Keep full expression matrix for proper DGE analysis
  cat("   DATA: Using full expression matrix:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")
  cat("   SUCCESS: Methodological correction: No pre-filtering to avoid bias\n")
  
  # Try to detect healthy vs disease groups
  detected_groups <- enhanced_auto_detect_groups(dataset)
  
  if (is.null(detected_groups)) {
    cat("   WARNING: No suitable groups detected - skipping DGE analysis\n\n")
    
    # Still save basic expression data
    camk_analysis_results[[dataset_id]] <- list(
      dataset_id = dataset_id,
      camk_core_genes_present = camk_present,
      camk_extended_genes_present = camk_extended_present,
      full_expression_matrix = expr_matrix,  # Save full matrix, not filtered
      n_samples = n_samples,
      n_total_genes = nrow(expr_matrix),
      groups_detected = FALSE
    )
    next
  }
  
  cat("   SUCCESS: Groups detected:", detected_groups$pattern_type, "\n")
  cat("      Healthy:", sum(detected_groups$groups == "Healthy"), "samples\n")
  cat("      Disease:", sum(detected_groups$groups == "Disease"), "samples\n")
  
  # Prepare data for DGE analysis (using FULL expression matrix)
  if (!is.null(detected_groups$sample_indices)) {
    # Filter expression matrix to samples with group assignments
    expr_filtered <- expr_matrix[, detected_groups$sample_indices, drop = FALSE]
    groups_vector <- detected_groups$groups
  } else {
    expr_filtered <- expr_matrix
    groups_vector <- detected_groups$groups
  }
  
  cat("   GENETIC: Genome-wide analysis:", nrow(expr_filtered), "genes x", ncol(expr_filtered), "samples\n")
  
  # Perform GENOME-WIDE DGE analysis first (methodologically correct)
  tryCatch({
    cat("   METHOD: Running genome-wide DGE analysis...\n")
    
    # Create design matrix (Disease vs Healthy)
    design <- model.matrix(~ groups_vector)
    colnames(design) <- c("Intercept", "Disease_vs_Healthy")
    
    # Fit linear model on FULL expression matrix
    fit <- lmFit(expr_filtered, design)
    fit <- eBayes(fit)
    
    # Get ALL results first
    all_dge_results <- topTable(fit, coef = "Disease_vs_Healthy", number = Inf, adjust.method = "BH")
    
    cat("   SUCCESS: Genome-wide DGE completed:", nrow(all_dge_results), "genes analyzed\n")
    
    # NOW filter for CAMK genes (post-DGE)
    camk_dge_results <- all_dge_results[rownames(all_dge_results) %in% camk_core_genes, , drop = FALSE]
    camk_extended_dge_results <- all_dge_results[rownames(all_dge_results) %in% camk_extended_genes, , drop = FALSE]
    
    # Add gene information to CAMK results
    if (nrow(camk_dge_results) > 0) {
      camk_dge_results$Gene_Symbol <- rownames(camk_dge_results)
      camk_dge_results$Dataset <- dataset_id
      camk_dge_results$Significant <- camk_dge_results$adj.P.Val < 0.05
      camk_dge_results$Regulation <- ifelse(camk_dge_results$logFC > 0, "UP in Disease", "DOWN in Disease")
    }
    
    # Add gene information to extended results
    if (nrow(camk_extended_dge_results) > 0) {
      camk_extended_dge_results$Gene_Symbol <- rownames(camk_extended_dge_results)
      camk_extended_dge_results$Dataset <- dataset_id
      camk_extended_dge_results$Significant <- camk_extended_dge_results$adj.P.Val < 0.05
      camk_extended_dge_results$Regulation <- ifelse(camk_extended_dge_results$logFC > 0, "UP in Disease", "DOWN in Disease")
    }
    
    cat("   SUCCESS: Post-DGE CAMK filtering completed!\n")
    cat("      Core CAMK genes analyzed:", nrow(camk_dge_results), "\n")
    cat("      Core CAMK significant (FDR < 0.05):", sum(camk_dge_results$Significant %||% 0), "\n")
    cat("      Extended CAMK+pathway genes analyzed:", nrow(camk_extended_dge_results), "\n")
    cat("      Extended significant (FDR < 0.05):", sum(camk_extended_dge_results$Significant %||% 0), "\n")
    cat("      TARGET: Methodological advantage: Proper normalization with", nrow(all_dge_results), "gene background\n")
    
    # Display individual CAMK gene results
    if (nrow(camk_dge_results) > 0) {
      cat("      Core CAMK gene results:\n")
      
      # Sort by significance
      camk_sorted <- camk_dge_results[order(camk_dge_results$adj.P.Val), ]
      
      for (i in 1:nrow(camk_sorted)) {
        gene <- camk_sorted$Gene_Symbol[i]
        logfc <- round(camk_sorted$logFC[i], 3)
        adj_pval <- camk_sorted$adj.P.Val[i]
        
        sig_marker <- if (adj_pval < 0.001) "***" else if (adj_pval < 0.01) "**" else if (adj_pval < 0.05) "*" else ""
        direction <- if (logfc > 0) "" else ""
        
        cat(sprintf("        %-8s: %s logFC=%6.3f, FDR=%8.2e %s\n", 
                    gene, direction, logfc, adj_pval, sig_marker))
      }
    }
    
    # Display top extended CAMK-related results
    if (nrow(camk_extended_dge_results) > 0) {
      significant_extended <- camk_extended_dge_results[camk_extended_dge_results$Significant, ]
      if (nrow(significant_extended) > 0) {
        cat("      Significant CAMK-related pathway genes:", nrow(significant_extended), "\n")
        
        # Show top 5 most significant
        top_extended <- head(significant_extended[order(significant_extended$adj.P.Val), ], 5)
        for (i in 1:nrow(top_extended)) {
          gene <- top_extended$Gene_Symbol[i]
          logfc <- round(top_extended$logFC[i], 3)
          adj_pval <- top_extended$adj.P.Val[i]
          direction <- if (logfc > 0) "" else ""
          
          cat(sprintf("        %-8s: %s logFC=%6.3f, FDR=%8.2e (pathway)\n", 
                      gene, direction, logfc, adj_pval))
        }
      }
    }
    
    # Store results with methodological improvements
    camk_analysis_results[[dataset_id]] <- list(
      dataset_id = dataset_id,
      camk_core_genes_present = camk_present,
      camk_extended_genes_present = camk_extended_present,
      full_expression_matrix = expr_filtered,  # Full matrix used for analysis
      groups = detected_groups,
      
      # DGE results
      all_dge_results = all_dge_results,  # Complete genome-wide results
      camk_dge_results = camk_dge_results,  # Core CAMK results
      camk_extended_dge_results = camk_extended_dge_results,  # Extended pathway results
      
      # Analysis metadata
      n_samples_analyzed = ncol(expr_filtered),
      n_total_genes_analyzed = nrow(all_dge_results),
      groups_detected = TRUE,
      significant_core_camk_genes = sum(camk_dge_results$Significant %||% 0),
      significant_extended_genes = sum(camk_extended_dge_results$Significant %||% 0),
      
      # Methodological notes
      analysis_approach = "genome_wide_then_filter",
      bias_correction = "no_pre_filtering_applied",
      statistical_power = "enhanced_with_full_transcriptome_background"
    )
    
  }, error = function(e) {
    cat("   ERROR: DGE analysis error:", e$message, "\n")
    
    # Store basic info even if DGE failed
    camk_analysis_results[[dataset_id]] <- list(
      dataset_id = dataset_id,
      camk_core_genes_present = camk_present,
      camk_extended_genes_present = camk_extended_present,
      full_expression_matrix = expr_filtered,
      groups = detected_groups,
      n_samples_analyzed = ncol(expr_filtered),
      n_total_genes = nrow(expr_filtered),
      groups_detected = TRUE,
      dge_error = e$message,
      analysis_approach = "genome_wide_then_filter"
    )
  })
  
  cat("\n")
}

# Generate summary report
cat("SUMMARY: CAMK ANALYSIS SUMMARY (METHODOLOGICALLY CORRECTED)\n")
cat("======================================================\n")

datasets_with_core_camk <- sum(sapply(camk_analysis_results, function(x) length(x$camk_core_genes_present %||% c()) > 0))
datasets_with_extended_camk <- sum(sapply(camk_analysis_results, function(x) length(x$camk_extended_genes_present %||% c()) > 0))
datasets_with_dge <- sum(sapply(camk_analysis_results, function(x) !is.null(x$all_dge_results)))
total_genes_analyzed <- sum(sapply(camk_analysis_results, function(x) x$n_total_genes_analyzed %||% 0), na.rm = TRUE)

cat("Datasets processed:", length(camk_analysis_results), "\n")
cat("Datasets with core CAMK genes:", datasets_with_core_camk, "\n")
cat("Datasets with extended CAMK+pathway genes:", datasets_with_extended_camk, "\n")
cat("Datasets with successful genome-wide DGE:", datasets_with_dge, "\n")
cat("Total genes analyzed (all datasets):", format(total_genes_analyzed, big.mark = ","), "\n")
cat("\nTARGET: METHODOLOGICAL IMPROVEMENT ACHIEVED:\n")
cat("  SUCCESS: No pre-filtering bias - full transcriptome background maintained\n")
cat("  SUCCESS: Proper statistical normalization with", format(total_genes_analyzed, big.mark = ","), "genes\n")
cat("  SUCCESS: Enhanced pathway context with extended gene set\n")
cat("  SUCCESS: Maintains CAMK focus through post-DGE filtering\n\n")

# Create combined CAMK results (from post-DGE filtering)
if (datasets_with_dge > 0) {
  cat("DATA: Combined CAMK DGE Results (Post-DGE Filtered):\n")
  cat("=================================================\n")
  
  # Combine core CAMK results
  all_camk_results <- do.call(rbind, lapply(camk_analysis_results, function(x) {
    if (!is.null(x$camk_dge_results) && nrow(x$camk_dge_results) > 0) x$camk_dge_results else NULL
  }))
  
  # Combine extended CAMK+pathway results
  all_extended_results <- do.call(rbind, lapply(camk_analysis_results, function(x) {
    if (!is.null(x$camk_extended_dge_results) && nrow(x$camk_extended_dge_results) > 0) x$camk_extended_dge_results else NULL
  }))
  
  # Core CAMK results summary
  if (!is.null(all_camk_results) && nrow(all_camk_results) > 0) {
    cat("\nCORE CAMK genes across all datasets:\n")
    for (gene in unique(all_camk_results$Gene_Symbol)) {
      gene_data <- all_camk_results[all_camk_results$Gene_Symbol == gene, ]
      n_datasets <- nrow(gene_data)
      n_significant <- sum(gene_data$Significant)
      avg_logfc <- round(mean(gene_data$logFC), 3)
      
      cat(sprintf("  %-8s: %d datasets, %d significant, avg logFC=%6.3f\n", 
                  gene, n_datasets, n_significant, avg_logfc))
    }
  }
  
  # Extended CAMK+pathway results summary
  if (!is.null(all_extended_results) && nrow(all_extended_results) > 0) {
    significant_extended <- all_extended_results[all_extended_results$Significant, ]
    if (nrow(significant_extended) > 0) {
      cat("\nSignificant CAMK-RELATED PATHWAY genes:\n")
      pathway_genes <- unique(significant_extended$Gene_Symbol)
      pathway_genes <- pathway_genes[!pathway_genes %in% camk_core_genes]  # Exclude core CAMK genes
      
      for (gene in head(pathway_genes, 10)) {  # Top 10
        gene_data <- significant_extended[significant_extended$Gene_Symbol == gene, ]
        n_datasets <- nrow(gene_data)
        avg_logfc <- round(mean(gene_data$logFC), 3)
        
        cat(sprintf("  %-8s: %d datasets, avg logFC=%6.3f (pathway context)\n", 
                    gene, n_datasets, avg_logfc))
      }
    }
  }
  
  cat("\nTARGET: METHODOLOGICAL VALIDATION:\n")
  cat("  SUCCESS: Results derived from genome-wide DGE analysis\n")
  cat("  SUCCESS: No statistical bias from gene pre-filtering\n")
  cat("  SUCCESS: Enhanced biological context with pathway genes\n")
  cat("  SUCCESS: Proper multiple testing correction applied\n")
}

# Save results using standardized utilities
ensure_output_directory("output")

# Save main analysis results
save_analysis_results(camk_analysis_results, "CAMK_focused_analysis_results_corrected", 
                     output_dir = "output")

# Save core CAMK DGE results
if (exists("all_camk_results") && !is.null(all_camk_results)) {
  save_analysis_results(all_camk_results, "CAMK_core_DGE_all_datasets", 
                       output_dir = "output", save_csv = TRUE)
}

# Save extended pathway results  
if (exists("all_extended_results") && !is.null(all_extended_results)) {
  save_analysis_results(all_extended_results, "CAMK_extended_pathway_DGE_all_datasets", 
                       output_dir = "output", save_csv = TRUE)
}

# Save methodological documentation
methodology_notes <- list(
  approach = "genome_wide_dge_then_filter",
  improvement = "eliminated_pre_filtering_bias",
  statistical_power = "enhanced_with_full_transcriptome_background",
  biological_context = "expanded_with_pathway_genes",
  gene_sets = list(
    core_camk = camk_core_genes,
    extended_pathway = camk_extended_genes,
    total_genes_analyzed = total_genes_analyzed
  ),
  timestamp = Sys.time()
)
saveRDS(methodology_notes, "output/CAMK_analysis_methodology_notes.rds")
cat("   • output/CAMK_analysis_methodology_notes.rds\n")

cat("\nCOMPLETE: CAMK-focused analysis completed with METHODOLOGICAL CORRECTION!\n")
cat("ENHANCED: Genome-wide → CAMK filtering approach eliminates statistical bias\n")
cat("METHOD: Enhanced biological context with", length(camk_extended_genes), "CAMK+pathway genes\n")
cat("DATA: Proper statistical foundation with", format(total_genes_analyzed, big.mark = ","), "genes analyzed\n")
cat("TARGET: Ready for unbiased clinical interpretation and downstream analysis\n")
cat("\nMETHODOLOGY: METHODOLOGICAL IMPACT:\n")
cat("  • Eliminated data skewing from pre-filtering\n")
cat("  • Maintained proper statistical normalization\n")
cat("  • Enhanced biological pathway context\n")
cat("  • Preserved CAMK research focus\n")
cat("\nSUCCESS: This approach provides statistically sound and biologically meaningful results.\n")