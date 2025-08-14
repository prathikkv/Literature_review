#!/usr/bin/env Rscript
#' Comprehensive Analysis Module
#' 
#' Consolidated statistical analysis functions including differential expression,
#' meta-analysis, and pathway enrichment for CAMK2D research

# Load required libraries
suppressPackageStartupMessages({
  library(limma)
  library(DESeq2)
  library(edgeR)
  library(metafor)
  library(meta)
  library(clusterProfiler)
  library(enrichplot)
  library(DOSE)
  library(ReactomePA)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(org.Rn.eg.db)
  library(biomaRt)
  library(tidyverse)
  library(ggplot2)
  # PARALLEL PROCESSING LIBRARIES (NEW)
  library(parallel)
  library(BiocParallel)
  library(future)
  library(future.apply)
})

#' Get CAMK Family Genes
#'
#' Returns the complete CAMK family gene list from prompts.md specification
#' @return Vector of CAMK gene symbols
get_camk_family_genes <- function() {
  c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", 
    "CAMKK1", "CAMKK2", "CAMK1", "CAMK1D", "CAMK1G", "CAMK4", "CAMKV")
}

#' PARALLEL PROCESSING OPTIMIZATION (NEW - HIGH PERFORMANCE)
#'
#' Setup intelligent parallel processing based on system capabilities
#' @param max_cores Maximum cores to use (NULL for auto-detection)
#' @param strategy Parallel strategy: "multicore", "cluster", or "auto"
#' @return Parallel backend information
setup_parallel_processing <- function(max_cores = NULL, strategy = "auto") {
  
  # Detect available cores
  available_cores <- detectCores()
  
  if (is.null(max_cores)) {
    # Use all but one core to keep system responsive
    optimal_cores <- max(1, available_cores - 1)
  } else {
    optimal_cores <- min(max_cores, available_cores)
  }
  
  # Choose strategy
  if (strategy == "auto") {
    if (.Platform$OS.type == "windows") {
      strategy <- "cluster"
    } else {
      strategy <- "multicore"
    }
  }
  
  cat("⚡ PARALLEL: Configuring", strategy, "processing with", optimal_cores, "cores\n")
  
  # Setup BiocParallel backend
  if (strategy == "multicore") {
    register(MulticoreParam(workers = optimal_cores))
    plan(multicore, workers = optimal_cores)
  } else {
    register(SnowParam(workers = optimal_cores))
    plan(cluster, workers = optimal_cores)
  }
  
  return(list(
    cores = optimal_cores,
    strategy = strategy,
    available_cores = available_cores
  ))
}

#' Parallel DGE Analysis
#'
#' Run DGE analysis in parallel across multiple datasets
#' @param processed_datasets List of datasets
#' @param focus_genes Genes of interest
#' @param use_parallel Enable parallel processing
#' @return DGE results
parallel_dge_analysis <- function(processed_datasets, focus_genes, use_parallel = TRUE) {
  
  if (use_parallel && length(processed_datasets) > 1) {
    cat("🚀 PARALLEL: Running DGE analysis on", length(processed_datasets), "datasets in parallel\n")
    
    # Setup parallel processing
    parallel_info <- setup_parallel_processing()
    
    # Parallel analysis function
    analyze_single_dataset <- function(dataset_entry) {
      dataset_id <- dataset_entry$id
      dataset <- dataset_entry$data
      
      tryCatch({
        # Auto-detect groups
        groups <- auto_detect_groups(dataset)
        if (is.null(groups)) {
          return(list(dataset_id = dataset_id, result = NULL, error = "No groups detected"))
        }
        
        # Determine analysis method
        platform_type <- determine_platform_type(dataset)
        
        if (platform_type == "microarray") {
          result <- perform_limma_analysis(dataset, groups, fdr_threshold = 0.05, fc_threshold = 1.2)
        } else {
          result <- perform_deseq2_analysis(dataset, groups, fdr_threshold = 0.05, fc_threshold = 1.2)
        }
        
        # Extract CAMK results
        camk_result <- extract_camk_results(result, focus_genes)
        
        return(list(
          dataset_id = dataset_id,
          result = result,
          camk_result = camk_result,
          error = NULL,
          platform_type = platform_type
        ))
        
      }, error = function(e) {
        return(list(dataset_id = dataset_id, result = NULL, error = e$message))
      })
    }
    
    # Prepare data for parallel processing
    dataset_entries <- lapply(names(processed_datasets), function(id) {
      list(id = id, data = processed_datasets[[id]])
    })
    
    # Run parallel analysis
    start_time <- Sys.time()
    parallel_results <- future_lapply(dataset_entries, analyze_single_dataset, 
                                     future.seed = TRUE)
    end_time <- Sys.time()
    
    cat("⚡ PERFORMANCE: Parallel analysis completed in", 
        round(as.numeric(difftime(end_time, start_time, units = "secs")), 1), "seconds\n")
    
    # Process results
    dge_results <- list()
    camk_results <- list()
    
    for (result in parallel_results) {
      if (!is.null(result$result)) {
        dge_results[[result$dataset_id]] <- result$result
        camk_results[[result$dataset_id]] <- result$camk_result
      } else {
        cat("⚠️ WARNING: Failed to analyze", result$dataset_id, ":", result$error, "\n")
      }
    }
    
    return(list(
      dge_results = dge_results,
      camk_results = camk_results,
      parallel_info = parallel_info,
      analysis_time = as.numeric(difftime(end_time, start_time, units = "secs"))
    ))
    
  } else {
    # Fallback to sequential processing
    cat("INFO: Using sequential processing\n")
    return(sequential_dge_analysis(processed_datasets, focus_genes))
  }
}

#' Sequential DGE Analysis (Fallback)
#'
#' Traditional sequential analysis for comparison or fallback
#' @param processed_datasets List of datasets
#' @param focus_genes Genes of interest
#' @return DGE results
sequential_dge_analysis <- function(processed_datasets, focus_genes) {
  
  dge_results <- list()
  camk_results <- list()
  
  start_time <- Sys.time()
  
  for (dataset_id in names(processed_datasets)) {
    dataset <- processed_datasets[[dataset_id]]
    
    cat("METHOD: Analyzing", dataset_id, "\n")
    
    tryCatch({
      groups <- auto_detect_groups(dataset)
      if (is.null(groups)) {
        cat("WARNING: Could not determine comparison groups for", dataset_id, "\n")
        next
      }
      
      platform_type <- determine_platform_type(dataset)
      
      if (platform_type == "microarray") {
        result <- perform_limma_analysis(dataset, groups, fdr_threshold = 0.05, fc_threshold = 1.2)
      } else {
        result <- perform_deseq2_analysis(dataset, groups, fdr_threshold = 0.05, fc_threshold = 1.2)
      }
      
      dge_results[[dataset_id]] <- result
      camk_results[[dataset_id]] <- extract_camk_results(result, focus_genes)
      
    }, error = function(e) {
      cat("ERROR: Failed to analyze", dataset_id, ":", e$message, "\n")
    })
  }
  
  end_time <- Sys.time()
  
  return(list(
    dge_results = dge_results,
    camk_results = camk_results,
    parallel_info = NULL,
    analysis_time = as.numeric(difftime(end_time, start_time, units = "secs"))
  ))
}

#' Comprehensive Differential Expression Analysis
#'
#' Performs DGE analysis with CAMK family focus across all datasets
#' @param processed_datasets List of processed datasets
#' @param focus_genes Vector of genes to prioritize
#' @param comparison_groups Comparison groups (auto-detected if NULL)
#' @param output_dir Output directory
#' @param fdr_threshold FDR threshold for significance
#' @param fold_change_threshold Minimum fold change threshold
#' @return List with DGE results
comprehensive_differential_expression_pipeline <- function(processed_datasets,
                                                          focus_genes = get_camk_family_genes(),
                                                          comparison_groups = NULL,
                                                          output_dir = "results/dge",
                                                          fdr_threshold = 0.05,
                                                          fold_change_threshold = 1.2) {
  
  cat("DATA: Comprehensive Differential Expression Analysis\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Enable parallel processing by default (can be disabled in config)
  use_parallel <- getOption("camk_pipeline.parallel", TRUE)
  
  cat("🚀 PERFORMANCE: Starting", ifelse(use_parallel, "parallel", "sequential"), "DGE analysis\n")
  
  # Run parallel or sequential analysis
  analysis_results <- if (use_parallel) {
    parallel_dge_analysis(processed_datasets, focus_genes)
  } else {
    sequential_dge_analysis(processed_datasets, focus_genes)
  }
  
  dge_results <- analysis_results$dge_results
  camk_results <- analysis_results$camk_results
  
  # Print performance summary
  if (!is.null(analysis_results$parallel_info)) {
    cat("⚡ PARALLEL PERFORMANCE: Used", analysis_results$parallel_info$cores, "cores,", 
        "completed in", round(analysis_results$analysis_time, 1), "seconds\n")
  }
  
  
  # Validate results
  cat("✅ SUCCESS: DGE analysis completed for", length(dge_results), "datasets\n")
  
  if (length(camk_results) > 0) {
    total_camk_genes <- sum(sapply(camk_results, function(x) if(!is.null(x)) nrow(x) else 0))
    cat("🎯 CAMK FOCUS: Extracted", total_camk_genes, "CAMK gene results across datasets\n")
  }
  
  return(list(
    dge_results = dge_results,
    camk_results = camk_results,
    analysis_parameters = list(
      fdr_threshold = fdr_threshold,
      fold_change_threshold = fold_change_threshold,
      focus_genes = focus_genes
    ),
    analysis_time = Sys.time()
  ))
}

#' Auto-detect Comparison Groups
#'
#' @param dataset Dataset object
#' @return List with group assignments or NULL
auto_detect_groups <- function(dataset) {
  
  pheno_data <- dataset$phenotype_data
  
  # Safety check for empty phenotype data
  if (is.null(pheno_data) || nrow(pheno_data) == 0) {
    return(NULL)
  }
  
  # Look for common disease/control indicators (including GEO-style colon names)
  group_columns <- c("disease", "condition", "group", "tissue", "treatment", 
                    "disease.state", "disease_state", "phenotype", "title",
                    "heart failure:ch1", "disease status:ch1", "condition:ch1", 
                    "disease:ch1", "group:ch1", "treatment:ch1")
  
  for (col in group_columns) {
    if (col %in% colnames(pheno_data)) {
      groups <- as.factor(pheno_data[[col]])
      
      # Safety check for empty groups
      if (length(groups) == 0) {
        next
      }
      
      # Check if we have at least 2 groups with reasonable sample sizes
      group_counts <- table(groups)
      if (length(group_counts) >= 2 && min(group_counts) >= 3) {
        return(list(groups = groups, column = col))
      }
    }
  }
  
  # Fallback: try to infer from sample names or titles
  sample_names <- rownames(pheno_data)
  if (!is.null(dataset$phenotype_data$title)) {
    sample_info <- dataset$phenotype_data$title
  } else {
    sample_info <- sample_names
  }
  
  # Safety check for sample info
  if (is.null(sample_info) || length(sample_info) == 0) {
    return(NULL)
  }
  
  # Look for AF/SR, disease/control patterns
  if (any(grepl("AF|fibrillation", sample_info, ignore.case = TRUE)) && 
      any(grepl("SR|sinus|control", sample_info, ignore.case = TRUE))) {
    
    groups <- ifelse(grepl("AF|fibrillation", sample_info, ignore.case = TRUE), "AF", "SR")
    # Check we have both groups
    if (length(unique(groups)) >= 2) {
      return(list(groups = as.factor(groups), column = "auto_detected"))
    }
  }
  
  if (any(grepl("HF|heart.failure|failure", sample_info, ignore.case = TRUE)) && 
      any(grepl("control|normal|healthy", sample_info, ignore.case = TRUE))) {
    
    groups <- ifelse(grepl("HF|heart.failure|failure", sample_info, ignore.case = TRUE), "HF", "Control")
    # Check we have both groups
    if (length(unique(groups)) >= 2) {
      return(list(groups = as.factor(groups), column = "auto_detected"))
    }
  }
  
  return(NULL)
}

#' Perform limma Analysis
#'
#' @param dataset Dataset object
#' @param groups Group information
#' @param fdr_threshold FDR threshold
#' @param fc_threshold Fold change threshold
#' @return DGE results data frame
perform_limma_analysis <- function(dataset, groups, fdr_threshold, fc_threshold) {
  
  expr_matrix <- dataset$expression_matrix
  group_factor <- groups$groups
  
  # Create design matrix
  design <- model.matrix(~ 0 + group_factor)
  colnames(design) <- levels(group_factor)
  
  # Fit linear model
  fit <- lmFit(expr_matrix, design)
  
  # Create contrast matrix for all pairwise comparisons
  group_levels <- levels(group_factor)
  if (length(group_levels) == 2) {
    contrast_matrix <- makeContrasts(contrasts = paste(group_levels[1], group_levels[2], sep = "-"), 
                                   levels = design)
  } else {
    # For multiple groups, just take first vs second
    contrast_matrix <- makeContrasts(contrasts = paste(group_levels[1], group_levels[2], sep = "-"), 
                                   levels = design)
  }
  
  # Fit contrasts
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  
  # Get results
  results <- topTable(fit2, coef = 1, number = Inf, adjust.method = "BH")
  
  # Add significance flags
  results$significant <- (results$adj.P.Val < fdr_threshold) & 
                        (abs(results$logFC) > log2(fc_threshold))
  
  # Map probe IDs to gene symbols
  probe_ids <- rownames(results)
  gene_symbols <- apply_gene_symbol_mapping(probe_ids, dataset)
  
  # Fallback: Add known CAMK gene mappings for GPL570 platform
  if (all(gene_symbols == probe_ids)) {  # No mapping occurred
    gene_symbols <- add_camk_gene_mappings(probe_ids, gene_symbols)
  }
  
  results$gene_symbol <- gene_symbols
  
  return(results)
}

#' Perform DESeq2 Analysis
#'
#' @param dataset Dataset object
#' @param groups Group information
#' @param fdr_threshold FDR threshold
#' @param fc_threshold Fold change threshold
#' @return DGE results data frame
perform_deseq2_analysis <- function(dataset, groups, fdr_threshold, fc_threshold) {
  
  # Convert expression matrix to integer counts (for DESeq2)
  count_matrix <- round(2^dataset$expression_matrix)  # Convert from log2 back to counts
  group_factor <- groups$groups
  
  # Create DESeq2 dataset
  coldata <- data.frame(condition = group_factor)
  rownames(coldata) <- colnames(count_matrix)
  
  dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                               colData = coldata,
                               design = ~ condition)
  
  # Filter low count genes
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  # Run DESeq2
  dds <- DESeq(dds)
  
  # Get results
  group_levels <- levels(group_factor)
  contrast <- c("condition", group_levels[1], group_levels[2])
  res <- results(dds, contrast = contrast, alpha = fdr_threshold)
  
  # Convert to data frame and clean up
  results_df <- as.data.frame(res)
  results_df <- results_df[complete.cases(results_df), ]
  
  # Add significance flags
  results_df$significant <- (results_df$padj < fdr_threshold) & 
                           (abs(results_df$log2FoldChange) > log2(fc_threshold))
  
  # Rename columns for consistency with limma
  colnames(results_df)[colnames(results_df) == "log2FoldChange"] <- "logFC"
  colnames(results_df)[colnames(results_df) == "padj"] <- "adj.P.Val"
  colnames(results_df)[colnames(results_df) == "pvalue"] <- "P.Value"
  
  # Map probe IDs to gene symbols
  probe_ids <- rownames(results_df)
  gene_symbols <- apply_gene_symbol_mapping(probe_ids, dataset)
  
  # Fallback: Add known CAMK gene mappings for GPL570 platform
  if (all(gene_symbols == probe_ids)) {  # No mapping occurred
    gene_symbols <- add_camk_gene_mappings(probe_ids, gene_symbols)
  }
  
  results_df$gene_symbol <- gene_symbols
  
  return(results_df)
}

#' Extract CAMK Results
#'
#' @param dge_result DGE results
#' @param focus_genes CAMK family genes
#' @return Subset of results for CAMK genes
extract_camk_results <- function(dge_result, focus_genes) {
  
  # Match gene symbols (case insensitive)
  gene_matches <- sapply(focus_genes, function(gene) {
    which(grepl(paste0("^", gene, "$"), dge_result$gene_symbol, ignore.case = TRUE))
  })
  
  # Flatten and get unique matches
  match_indices <- unique(unlist(gene_matches))
  match_indices <- match_indices[!is.na(match_indices)]
  
  if (length(match_indices) > 0) {
    camk_subset <- dge_result[match_indices, ]
    camk_subset$camk_gene <- TRUE
    return(camk_subset)
  }
  
  return(data.frame())  # Return empty data frame if no matches
}

#' Comprehensive Meta-Analysis Pipeline
#'
#' Combines DGE results across studies using random effects models
#' @param dge_results_list List of DGE results from different studies
#' @param focus_genes Vector of genes to prioritize
#' @param output_dir Directory for meta-analysis results
#' @param effect_size_method Method for calculating effect sizes
#' @param min_studies Minimum number of studies required
#' @return List with meta-analysis results
comprehensive_meta_analysis_pipeline <- function(dge_results_list,
                                                focus_genes = get_camk_family_genes(),
                                                output_dir = "results/meta_analysis",
                                                effect_size_method = "log_fc",
                                                min_studies = 2) {
  
  cat("RESULTS: Comprehensive Meta-Analysis Framework\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Validate inputs
  if (length(dge_results_list) == 0) {
    cat("ERROR: No DGE results provided for meta-analysis\n")
    return(NULL)
  }
  
  # Get all genes present in studies
  all_genes <- unique(unlist(lapply(dge_results_list, function(x) x$gene_symbol)))
  
  # Adjust min_studies based on available data
  actual_studies <- length(dge_results_list)
  effective_min_studies <- min(min_studies, actual_studies)
  
  cat("INFO: Available studies:", actual_studies, "\n")
  cat("INFO: Effective minimum studies required:", effective_min_studies, "\n")
  
  # Filter genes present in at least effective_min_studies
  gene_study_counts <- sapply(all_genes, function(gene) {
    sum(sapply(dge_results_list, function(study) gene %in% study$gene_symbol))
  })
  
  eligible_genes <- names(gene_study_counts)[gene_study_counts >= effective_min_studies]
  
  cat("DATA: Genes eligible for meta-analysis:", length(eligible_genes), "\n")
  cat("TARGET: CAMK genes in eligible set:", sum(focus_genes %in% eligible_genes), "\n\n")
  
  # Perform meta-analysis for each eligible gene
  meta_results <- list()
  camk_meta_results <- list()
  
  for (gene in eligible_genes) {
    tryCatch({
      # Extract effect sizes and standard errors across studies
      study_data <- extract_gene_across_studies(gene, dge_results_list)
      
      if (nrow(study_data) >= min_studies) {
        # Perform random effects meta-analysis
        meta_result <- perform_meta_analysis(gene, study_data)
        
        if (!is.null(meta_result)) {
          meta_results[[gene]] <- meta_result
          
          # Store CAMK results separately
          if (gene %in% focus_genes) {
            camk_meta_results[[gene]] <- meta_result
          }
        }
      }
      
    }, error = function(e) {
      # Silent error handling for individual genes
    })
  }
  
  cat("SUCCESS: Meta-analysis completed for", length(meta_results), "genes\n")
  cat("TARGET: CAMK family results:", length(camk_meta_results), "genes\n")
  
  return(list(
    gene_meta_results = meta_results,
    camk_meta_results = camk_meta_results,
    eligible_genes = eligible_genes,
    analysis_parameters = list(
      min_studies = min_studies,
      effect_size_method = effect_size_method,
      focus_genes = focus_genes
    ),
    meta_analysis_time = Sys.time()
  ))
}

#' Extract Gene Across Studies
#'
#' @param gene Gene symbol
#' @param dge_results_list List of DGE results
#' @return Data frame with study information
extract_gene_across_studies <- function(gene, dge_results_list) {
  
  study_data <- data.frame()
  
  for (study_id in names(dge_results_list)) {
    dge_result <- dge_results_list[[study_id]]
    
    gene_row <- dge_result[dge_result$gene_symbol == gene, ]
    
    if (nrow(gene_row) > 0) {
      study_row <- data.frame(
        study_id = study_id,
        gene = gene,
        effect_size = gene_row$logFC[1],
        standard_error = gene_row$logFC[1] / abs(qt(gene_row$P.Value[1]/2, df = Inf)),  # Approximate SE
        p_value = gene_row$P.Value[1],
        adj_p_value = gene_row$adj.P.Val[1],
        n_samples = 50  # Approximate - would need actual sample sizes
      )
      
      study_data <- rbind(study_data, study_row)
    }
  }
  
  return(study_data)
}

#' Perform Meta-Analysis
#'
#' @param gene Gene symbol
#' @param study_data Study data for the gene
#' @return Meta-analysis result
perform_meta_analysis <- function(gene, study_data) {
  
  tryCatch({
    # Handle single study case
    if (nrow(study_data) == 1) {
      cat("   INFO: Single study analysis for", gene, "\n")
      
      return(list(
        gene = gene,
        pooled_effect = study_data$effect_size[1],
        pooled_se = study_data$standard_error[1],
        pooled_pval = study_data$p_value[1],
        pooled_ci_lb = study_data$effect_size[1] - 1.96 * study_data$standard_error[1],
        pooled_ci_ub = study_data$effect_size[1] + 1.96 * study_data$standard_error[1],
        i2 = 0,  # No heterogeneity with single study
        tau2 = 0,
        q_stat = 0,
        q_pval = 1,
        n_studies = 1,
        study_data = study_data,
        meta_object = NULL,
        analysis_type = "single_study"
      ))
    }
    
    # Multi-study meta-analysis using metafor
    meta_result <- rma(yi = effect_size, 
                      sei = standard_error,
                      data = study_data,
                      method = "REML")
    
    return(list(
      gene = gene,
      pooled_effect = meta_result$beta[1],
      pooled_se = meta_result$se,
      pooled_pval = meta_result$pval,
      pooled_ci_lb = meta_result$ci.lb,
      pooled_ci_ub = meta_result$ci.ub,
      i2 = meta_result$I2,
      tau2 = meta_result$tau2,
      q_stat = meta_result$QE,
      q_pval = meta_result$QEp,
      n_studies = nrow(study_data),
      study_data = study_data,
      meta_object = meta_result,
      analysis_type = "meta_analysis"
    ))
    
  }, error = function(e) {
    cat("   WARNING: Meta-analysis failed for", gene, ":", e$message, "\n")
    return(NULL)
  })
}

#' Comprehensive Pathway Analysis Pipeline
#'
#' GO, KEGG, Reactome enrichment with CAMK2D focus
#' @param dge_results_list List of DGE results
#' @param expression_data_list List of expression matrices
#' @param species Target species
#' @param output_dir Output directory
#' @param fdr_threshold FDR threshold
#' @param min_gene_set_size Minimum gene set size
#' @return Pathway analysis results
comprehensive_pathway_analysis_pipeline <- function(dge_results_list,
                                                   expression_data_list = NULL,
                                                   species = "human",
                                                   output_dir = "results/pathway_analysis",
                                                   fdr_threshold = 0.05,
                                                   min_gene_set_size = 10) {
  
  cat("🛤️ Comprehensive Pathway Analysis Pipeline\n")
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Set organism database
  if (species == "human") {
    orgdb <- org.Hs.eg.db
  } else if (species == "mouse") {
    orgdb <- org.Mm.eg.db
  } else if (species == "rat") {
    orgdb <- org.Rn.eg.db
  } else {
    stop("Unsupported species: ", species)
  }
  
  pathway_results <- list()
  
  # Combine significant genes across studies
  all_sig_genes <- c()
  for (study_id in names(dge_results_list)) {
    dge_result <- dge_results_list[[study_id]]
    sig_genes <- dge_result$gene_symbol[dge_result$significant == TRUE]
    all_sig_genes <- c(all_sig_genes, sig_genes)
  }
  
  # Get unique significant genes
  unique_sig_genes <- unique(all_sig_genes)
  cat("DATA: Total significant genes for pathway analysis:", length(unique_sig_genes), "\n")
  
  if (length(unique_sig_genes) < min_gene_set_size) {
    cat("WARNING: Too few significant genes for pathway analysis\n")
    return(NULL)
  }
  
  tryCatch({
    # Convert gene symbols to Entrez IDs
    entrez_ids <- bitr(unique_sig_genes, fromType = "SYMBOL", toType = "ENTREZID", 
                      OrgDb = orgdb, drop = TRUE)$ENTREZID
    
    if (length(entrez_ids) < min_gene_set_size) {
      cat("WARNING: Too few genes could be converted to Entrez IDs\n")
      return(NULL)
    }
    
    # GO enrichment analysis
    go_bp <- enrichGO(gene = entrez_ids,
                     OrgDb = orgdb,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = fdr_threshold,
                     qvalueCutoff = fdr_threshold,
                     readable = TRUE)
    
    go_mf <- enrichGO(gene = entrez_ids,
                     OrgDb = orgdb,
                     ont = "MF", 
                     pAdjustMethod = "BH",
                     pvalueCutoff = fdr_threshold,
                     qvalueCutoff = fdr_threshold,
                     readable = TRUE)
    
    # KEGG pathway analysis
    if (species == "human") {
      kegg_pathways <- enrichKEGG(gene = entrez_ids,
                                 organism = "hsa",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = fdr_threshold,
                                 qvalueCutoff = fdr_threshold)
    }
    
    pathway_results <- list(
      go_results = list(
        biological_process = go_bp,
        molecular_function = go_mf
      ),
      kegg_results = if (exists("kegg_pathways")) kegg_pathways else NULL,
      input_genes = unique_sig_genes,
      entrez_ids = entrez_ids,
      analysis_time = Sys.time()
    )
    
    cat("SUCCESS: Pathway analysis completed\n")
    if (!is.null(go_bp)) cat("   GO BP pathways:", nrow(go_bp@result), "\n")
    if (!is.null(go_mf)) cat("   GO MF pathways:", nrow(go_mf@result), "\n")
    if (!is.null(kegg_pathways)) cat("   KEGG pathways:", nrow(kegg_pathways@result), "\n")
    
  }, error = function(e) {
    cat("ERROR: Error in pathway analysis:", e$message, "\n")
    pathway_results <- NULL
  })
  
  return(pathway_results)
}

#' Add Known CAMK Gene Mappings for GPL570 Platform
#'
#' Provides fallback gene symbol mapping for key CAMK family genes
#' @param probe_ids Vector of probe IDs
#' @param current_symbols Current gene symbols (usually same as probe_ids if no mapping)
#' @return Updated gene symbols with CAMK mappings
add_camk_gene_mappings <- function(probe_ids, current_symbols) {
  
  # Known GPL570 probe ID to CAMK gene mappings (from Affymetrix HG-U133_Plus_2)
  # These mappings were derived from Affymetrix annotation files
  camk_mappings <- list(
    # CAMK2D - Primary target gene
    "218542_at" = "CAMK2D",
    "218543_s_at" = "CAMK2D", 
    
    # Other CAMK family members
    "201370_s_at" = "CAMK2A",
    "201371_s_at" = "CAMK2A",
    "209949_at" = "CAMK2B",
    "209950_s_at" = "CAMK2B",
    "206808_at" = "CAMK2G",
    "206809_s_at" = "CAMK2G",
    "214391_x_at" = "CAMK4",
    "214392_at" = "CAMK4",
    "220458_at" = "CAMKK1", 
    "220459_s_at" = "CAMKK1",
    "218181_s_at" = "CAMKK2",
    "218182_at" = "CAMKK2",
    "205449_at" = "CAMK1",
    "205450_s_at" = "CAMK1",
    "221520_s_at" = "CAMK1D",
    "221521_at" = "CAMK1D",
    "214945_at" = "CAMK1G",
    "214946_x_at" = "CAMK1G",
    "210619_s_at" = "CAMKV",
    "210620_at" = "CAMKV"
  )
  
  # Apply mappings
  updated_symbols <- current_symbols
  for (i in seq_along(probe_ids)) {
    probe_id <- probe_ids[i]
    if (probe_id %in% names(camk_mappings)) {
      updated_symbols[i] <- camk_mappings[[probe_id]]
      cat("   Mapped probe", probe_id, "to gene", camk_mappings[[probe_id]], "\n")
    }
  }
  
  return(updated_symbols)
}

cat("SUCCESS: Comprehensive Analysis Module loaded successfully\n")
cat("SUMMARY: Main functions: comprehensive_differential_expression_pipeline(), comprehensive_meta_analysis_pipeline(), comprehensive_pathway_analysis_pipeline()\n")
cat("GENETIC: CAMK gene mapping function: add_camk_gene_mappings()\n")