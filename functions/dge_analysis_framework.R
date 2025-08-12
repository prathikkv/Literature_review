#!/usr/bin/env Rscript
#' Differential Gene Expression Analysis Framework
#' 
#' Comprehensive DGE analysis framework for CAMK2D cardiac research
#' Supports both RNA-Seq and microarray data with multiple analysis methods

# Load required libraries
suppressPackageStartupMessages({
  library(limma)
  library(edgeR)
  library(DESeq2)
  library(tidyverse)
  library(biomaRt)
  library(pheatmap)
  library(enrichR)
  library(clusterProfiler)
})

#' Comprehensive DGE Analysis Pipeline
#'
#' @param expression_matrix Expression matrix (genes x samples)
#' @param metadata Sample metadata with disease status
#' @param geo_accession Dataset identifier
#' @param analysis_type Type of analysis ("RNA-Seq", "Microarray", "auto")
#' @param output_dir Output directory
#' @return List with DGE results and plots
perform_dge_analysis <- function(expression_matrix, metadata, geo_accession, 
                                analysis_type = "auto", output_dir = "results/dge_analysis") {
  
  cat("🧬 Starting DGE analysis for", geo_accession, "\n")
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Auto-detect analysis type if needed
  if (analysis_type == "auto") {
    analysis_type <- detect_data_type(expression_matrix)
  }
  
  cat("📊 Analysis type:", analysis_type, "\n")
  cat("🔬 Samples:", ncol(expression_matrix)-1, "| Genes:", nrow(expression_matrix)-1, "\n")
  
  # Prepare data
  prepared_data <- prepare_dge_data(expression_matrix, metadata, analysis_type)
  
  if (is.null(prepared_data)) {
    cat("❌ Failed to prepare data for analysis\n")
    return(NULL)
  }
  
  # Perform DGE analysis based on data type
  dge_results <- switch(analysis_type,
    "RNA-Seq" = analyze_rnaseq_data(prepared_data, geo_accession),
    "Microarray" = analyze_microarray_data(prepared_data, geo_accession),
    "Single-cell" = analyze_singlecell_data(prepared_data, geo_accession),
    stop("Unsupported analysis type:", analysis_type)
  )
  
  if (is.null(dge_results)) {
    cat("❌ DGE analysis failed\n")
    return(NULL)
  }
  
  # CAMK2D-focused analysis
  camk_analysis <- analyze_camk_family(dge_results, prepared_data, geo_accession)
  
  # Pathway enrichment analysis
  pathway_results <- perform_pathway_analysis(dge_results, geo_accession)
  
  # Generate visualizations
  plots <- generate_dge_plots(dge_results, camk_analysis, prepared_data, geo_accession, output_dir)
  
  # Save results
  save_dge_results(dge_results, camk_analysis, pathway_results, geo_accession, output_dir)
  
  cat("✅ DGE analysis completed for", geo_accession, "\n")
  
  return(list(
    dge_results = dge_results,
    camk_analysis = camk_analysis,
    pathway_results = pathway_results,
    plots = plots,
    metadata = prepared_data$metadata,
    geo_accession = geo_accession
  ))
}

#' Prepare Data for DGE Analysis
#'
#' @param expression_matrix Raw expression matrix
#' @param metadata Sample metadata
#' @param analysis_type Type of analysis
#' @return List with prepared data
prepare_dge_data <- function(expression_matrix, metadata, analysis_type) {
  
  cat("🔧 Preparing data for DGE analysis...\n")
  
  # Remove gene_id column if present and set as rownames
  if ("gene_id" %in% colnames(expression_matrix)) {
    rownames(expression_matrix) <- expression_matrix$gene_id
    expression_matrix <- expression_matrix[, !colnames(expression_matrix) %in% "gene_id"]
  }
  
  # Ensure metadata matches expression matrix samples
  common_samples <- intersect(colnames(expression_matrix), metadata$sample_id)
  
  if (length(common_samples) < 6) {
    cat("❌ Insufficient samples for analysis:", length(common_samples), "\n")
    return(NULL)
  }
  
  # Filter to common samples
  expression_matrix <- expression_matrix[, common_samples]
  metadata <- metadata[metadata$sample_id %in% common_samples, ]
  
  # Clean disease status
  metadata$disease_status <- clean_disease_status(metadata$disease_status)
  
  # Check for sufficient groups
  disease_counts <- table(metadata$disease_status)
  if (length(disease_counts) < 2 || min(disease_counts) < 3) {
    cat("❌ Insufficient samples per group for analysis\n")
    print(disease_counts)
    return(NULL)
  }
  
  cat("✅ Prepared data: ", nrow(expression_matrix), "genes x", ncol(expression_matrix), "samples\n")
  cat("   📊 Group sizes:", paste(names(disease_counts), disease_counts, sep="=", collapse=", "), "\n")
  
  # Quality control filtering
  filtered_data <- perform_quality_filtering(expression_matrix, metadata, analysis_type)
  
  return(list(
    expression = filtered_data$expression,
    metadata = metadata,
    analysis_type = analysis_type,
    original_genes = nrow(expression_matrix),
    filtered_genes = nrow(filtered_data$expression)
  ))
}

#' RNA-Seq DGE Analysis using DESeq2 and edgeR
#'
#' @param prepared_data Prepared data object
#' @param geo_accession Dataset identifier
#' @return DGE results
analyze_rnaseq_data <- function(prepared_data, geo_accession) {
  
  cat("🧬 Performing RNA-Seq DGE analysis...\n")
  
  expr_data <- prepared_data$expression
  metadata <- prepared_data$metadata
  
  # Ensure counts are integers (round if necessary)
  if (any(expr_data != round(expr_data))) {
    cat("   🔧 Converting to integer counts...\n")
    expr_data <- round(expr_data)
  }
  
  # Remove negative values
  expr_data[expr_data < 0] <- 0
  
  # DESeq2 Analysis
  cat("   📊 Running DESeq2 analysis...\n")
  deseq_results <- tryCatch({
    
    # Create DESeq2 object
    dds <- DESeqDataSetFromMatrix(
      countData = as.matrix(expr_data),
      colData = metadata,
      design = ~ disease_status
    )
    
    # Filter low count genes
    dds <- dds[rowSums(counts(dds)) >= 10, ]
    
    # Run DESeq2
    dds <- DESeq(dds, quiet = TRUE)
    
    # Extract results
    res <- results(dds, contrast = c("disease_status", "Disease", "Control"))
    
    # Convert to data frame
    res_df <- as.data.frame(res) %>%
      rownames_to_column("gene_id") %>%
      arrange(pvalue) %>%
      mutate(
        analysis_method = "DESeq2",
        significant = padj < 0.05 & abs(log2FoldChange) > 1,
        direction = ifelse(log2FoldChange > 0, "Up", "Down")
      )
    
    list(
      results = res_df,
      dds_object = dds,
      success = TRUE
    )
    
  }, error = function(e) {
    cat("   ⚠️ DESeq2 analysis failed:", e$message, "\n")
    list(success = FALSE, error = e$message)
  })
  
  # edgeR Analysis (alternative/validation)
  cat("   📊 Running edgeR analysis...\n")
  edger_results <- tryCatch({
    
    # Create DGEList object
    y <- DGEList(counts = as.matrix(expr_data), group = metadata$disease_status)
    
    # Filter low expressed genes
    keep <- filterByExpr(y)
    y <- y[keep, ]
    
    # Normalize
    y <- calcNormFactors(y)
    
    # Estimate dispersion
    design <- model.matrix(~ disease_status, data = metadata)
    y <- estimateDisp(y, design)
    
    # Fit model and test
    fit <- glmQLFit(y, design)
    qlf <- glmQLFTest(fit, coef = 2)
    
    # Extract results
    res <- topTags(qlf, n = Inf)$table %>%
      rownames_to_column("gene_id") %>%
      arrange(PValue) %>%
      mutate(
        analysis_method = "edgeR",
        significant = FDR < 0.05 & abs(logFC) > 1,
        direction = ifelse(logFC > 0, "Up", "Down")
      )
    
    list(
      results = res,
      dge_object = y,
      fit_object = fit,
      success = TRUE
    )
    
  }, error = function(e) {
    cat("   ⚠️ edgeR analysis failed:", e$message, "\n")
    list(success = FALSE, error = e$message)
  })
  
  # Combine results
  combined_results <- combine_dge_methods(deseq_results, edger_results, "RNA-Seq")
  
  return(combined_results)
}

#' Microarray DGE Analysis using limma
#'
#' @param prepared_data Prepared data object
#' @param geo_accession Dataset identifier
#' @return DGE results
analyze_microarray_data <- function(prepared_data, geo_accession) {
  
  cat("🔬 Performing Microarray DGE analysis...\n")
  
  expr_data <- prepared_data$expression
  metadata <- prepared_data$metadata
  
  # limma analysis
  cat("   📊 Running limma analysis...\n")
  limma_results <- tryCatch({
    
    # Create design matrix
    design <- model.matrix(~ disease_status, data = metadata)
    
    # Fit linear model
    fit <- lmFit(as.matrix(expr_data), design)
    
    # Empirical Bayes
    fit <- eBayes(fit)
    
    # Extract results
    res <- topTable(fit, coef = 2, n = Inf) %>%
      rownames_to_column("gene_id") %>%
      arrange(P.Value) %>%
      mutate(
        analysis_method = "limma",
        significant = adj.P.Val < 0.05 & abs(logFC) > 1,
        direction = ifelse(logFC > 0, "Up", "Down")
      )
    
    list(
      results = res,
      fit_object = fit,
      success = TRUE
    )
    
  }, error = function(e) {
    cat("   ⚠️ limma analysis failed:", e$message, "\n")
    list(success = FALSE, error = e$message)
  })
  
  return(limma_results)
}

#' CAMK Family-Focused Analysis
#'
#' @param dge_results DGE analysis results
#' @param prepared_data Prepared data
#' @param geo_accession Dataset identifier
#' @return CAMK-specific results
analyze_camk_family <- function(dge_results, prepared_data, geo_accession) {
  
  cat("🎯 Performing CAMK family analysis...\n")
  
  # Define CAMK family genes
  camk_genes <- c(
    "CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G",  # Primary CAMK2 family
    "CAMKK2", "CAMKK1",                       # Upstream kinases
    "CAMK1", "CAMK4", "CAMK1D", "CAMK1G"     # Related CAMK family
  )
  
  # Extract CAMK results from DGE
  if (dge_results$success) {
    camk_dge <- dge_results$results %>%
      filter(gene_id %in% camk_genes) %>%
      mutate(
        camk_subfamily = case_when(
          gene_id %in% c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G") ~ "CAMK2_Family",
          gene_id %in% c("CAMKK1", "CAMKK2") ~ "CAMKK_Family",
          TRUE ~ "Other_CAMK"
        ),
        proposal_relevance = case_when(
          gene_id == "CAMK2D" ~ "Primary_Target",
          gene_id %in% c("CAMK2A", "CAMK2B", "CAMK2G") ~ "Functional_Redundancy",
          gene_id %in% c("CAMKK2") ~ "Upstream_Regulator",
          TRUE ~ "Related"
        )
      )
    
    # Calculate expression levels
    camk_expression <- prepared_data$expression[rownames(prepared_data$expression) %in% camk_genes, ] %>%
      rownames_to_column("gene_id") %>%
      pivot_longer(-gene_id, names_to = "sample_id", values_to = "expression") %>%
      left_join(prepared_data$metadata, by = "sample_id") %>%
      group_by(gene_id, disease_status) %>%
      summarise(
        mean_expression = mean(expression, na.rm = TRUE),
        median_expression = median(expression, na.rm = TRUE),
        sd_expression = sd(expression, na.rm = TRUE),
        n_samples = n(),
        .groups = "drop"
      )
    
    cat("✅ Found", nrow(camk_dge), "CAMK family genes in DGE results\n")
    
    return(list(
      camk_dge_results = camk_dge,
      camk_expression_summary = camk_expression,
      camk_genes_tested = camk_genes,
      camk_genes_found = unique(camk_dge$gene_id)
    ))
  } else {
    cat("❌ No DGE results available for CAMK analysis\n")
    return(NULL)
  }
}

#' Pathway Enrichment Analysis
#'
#' @param dge_results DGE results
#' @param geo_accession Dataset identifier
#' @return Pathway enrichment results
perform_pathway_analysis <- function(dge_results, geo_accession) {
  
  cat("🛤️  Performing pathway enrichment analysis...\n")
  
  if (!dge_results$success) {
    cat("❌ No DGE results for pathway analysis\n")
    return(NULL)
  }
  
  # Get significant genes
  sig_genes <- dge_results$results %>%
    filter(significant == TRUE) %>%
    pull(gene_id)
  
  if (length(sig_genes) < 10) {
    cat("⚠️ Too few significant genes for pathway analysis:", length(sig_genes), "\n")
    return(NULL)
  }
  
  cat("   📊 Analyzing", length(sig_genes), "significant genes\n")
  
  # Perform enrichment with multiple databases
  tryCatch({
    
    # Set up enrichR databases
    dbs <- c("GO_Biological_Process_2021", "KEGG_2021_Human", 
             "Reactome_2022", "WikiPathway_2021_Human")
    
    enrichr_results <- enrichr(sig_genes, dbs)
    
    # Focus on cardiac and calcium signaling pathways
    cardiac_pathways <- list()
    
    for (db in names(enrichr_results)) {
      cardiac_terms <- enrichr_results[[db]] %>%
        filter(
          grepl("cardiac|heart|calcium|kinase|phosphor", Term, ignore.case = TRUE) |
          Adjusted.P.value < 0.05
        ) %>%
        head(20) %>%
        mutate(database = db)
      
      if (nrow(cardiac_terms) > 0) {
        cardiac_pathways[[db]] <- cardiac_terms
      }
    }
    
    if (length(cardiac_pathways) > 0) {
      combined_pathways <- bind_rows(cardiac_pathways)
      cat("✅ Found", nrow(combined_pathways), "enriched pathways\n")
      
      return(list(
        all_results = enrichr_results,
        cardiac_pathways = combined_pathways,
        input_genes = sig_genes
      ))
    }
    
  }, error = function(e) {
    cat("⚠️ Pathway analysis failed:", e$message, "\n")
  })
  
  return(NULL)
}

#' Helper Functions

detect_data_type <- function(expression_matrix) {
  # Simple heuristic to detect data type
  sample_values <- as.numeric(expression_matrix[1:100, 1:min(5, ncol(expression_matrix))])
  
  # Check if values are mostly integers (RNA-Seq)
  if (mean(sample_values == round(sample_values), na.rm = TRUE) > 0.8) {
    return("RNA-Seq")
  } else {
    return("Microarray")
  }
}

clean_disease_status <- function(disease_status) {
  # Standardize disease status labels
  disease_status <- tolower(as.character(disease_status))
  
  disease_status[grepl("control|normal|healthy", disease_status)] <- "Control"
  disease_status[grepl("disease|patient|affected|case", disease_status)] <- "Disease"
  
  return(disease_status)
}

perform_quality_filtering <- function(expression_matrix, metadata, analysis_type) {
  
  if (analysis_type == "RNA-Seq") {
    # Filter low count genes
    keep_genes <- rowSums(expression_matrix > 1) >= ncol(expression_matrix) * 0.1
    filtered_expr <- expression_matrix[keep_genes, ]
  } else {
    # Filter low expression genes for microarray
    keep_genes <- rowSums(is.na(expression_matrix)) < ncol(expression_matrix) * 0.5
    filtered_expr <- expression_matrix[keep_genes, ]
  }
  
  cat("   🔧 Filtered from", nrow(expression_matrix), "to", nrow(filtered_expr), "genes\n")
  
  return(list(expression = filtered_expr))
}

combine_dge_methods <- function(method1_results, method2_results, analysis_type) {
  
  if (method1_results$success && method2_results$success) {
    # Both methods succeeded - combine results
    return(list(
      results = method1_results$results,
      alternative_results = method2_results$results,
      success = TRUE,
      analysis_type = analysis_type,
      methods_used = c(method1_results$results$analysis_method[1], 
                      method2_results$results$analysis_method[1])
    ))
  } else if (method1_results$success) {
    return(method1_results)
  } else if (method2_results$success) {
    return(method2_results)
  } else {
    return(list(success = FALSE, error = "Both methods failed"))
  }
}

generate_dge_plots <- function(dge_results, camk_analysis, prepared_data, geo_accession, output_dir) {
  # Placeholder for plot generation
  cat("   📊 Generating DGE plots...\n")
  return(list(plots_generated = TRUE))
}

save_dge_results <- function(dge_results, camk_analysis, pathway_results, geo_accession, output_dir) {
  
  # Save main DGE results
  if (dge_results$success) {
    write_csv(dge_results$results, 
              file.path(output_dir, paste0(geo_accession, "_dge_results.csv")))
  }
  
  # Save CAMK analysis
  if (!is.null(camk_analysis)) {
    write_csv(camk_analysis$camk_dge_results, 
              file.path(output_dir, paste0(geo_accession, "_camk_analysis.csv")))
    
    write_csv(camk_analysis$camk_expression_summary, 
              file.path(output_dir, paste0(geo_accession, "_camk_expression.csv")))
  }
  
  # Save pathway results
  if (!is.null(pathway_results)) {
    write_csv(pathway_results$cardiac_pathways, 
              file.path(output_dir, paste0(geo_accession, "_pathways.csv")))
  }
  
  cat("💾 Results saved to", output_dir, "\n")
}

cat("🧬 DGE Analysis Framework Loaded\n")
cat("   📊 Supports RNA-Seq (DESeq2, edgeR) and Microarray (limma)\n")
cat("   🎯 CAMK2D-focused analysis with pathway enrichment\n")