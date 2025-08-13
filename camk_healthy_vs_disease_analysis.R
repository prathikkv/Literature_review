#!/usr/bin/env Rscript
#' CAMK Healthy vs Disease Analysis - Corrected Pipeline
#' 
#' Focuses exclusively on datasets with true healthy vs disease comparisons
#' Excludes AF vs SR and other disease subtype comparisons

cat("🎯 CAMK HEALTHY vs DISEASE ANALYSIS\n")
cat("==================================\n\n")

# Load required functions
source("functions/analysis.R")
source("enhanced_group_detection.R")

# Define CAMK genes
camk_genes <- c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", "CAMKK1", "CAMKK2", 
                "CAMK1", "CAMK1D", "CAMK1G", "CAMK4", "CAMKV")

cat("🧬 CAMK genes to analyze:", length(camk_genes), "\n")
cat("   ", paste(camk_genes, collapse = ", "), "\n\n")

# Initialize results storage
healthy_disease_results <- list()

# Dataset IDs to check for healthy vs disease patterns
datasets_to_check <- c("GSE57338", "GSE14975", "GSE120895")

cat("📊 ANALYZING DATASETS FOR HEALTHY VS DISEASE PATTERNS\n")
cat("====================================================\n\n")

for (dataset_id in datasets_to_check) {
  
  cat("📋 Processing:", dataset_id, "\n")
  cat("   ", paste(rep("=", nchar(dataset_id) + 12), collapse = ""), "\n")
  
  # Check if dataset exists
  cache_file <- paste0("cache/comprehensive/", dataset_id, "_processed.rds")
  
  if (!file.exists(cache_file)) {
    cat("   ❌ Dataset not found in cache\n\n")
    next
  }
  
  # Load dataset
  tryCatch({
    dataset <- readRDS(cache_file)
    cat("   ✅ Dataset loaded successfully\n")
    cat("   📊 Samples:", nrow(dataset$phenotype_data), "\n")
    cat("   🧬 Genes:", nrow(dataset$expression_data), "\n")
    
    # Apply enhanced group detection
    group_info <- enhanced_auto_detect_groups(dataset)
    
    if (is.null(group_info)) {
      cat("   ❌ No suitable groups detected\n\n")
      next
    }
    
    # Check if this is a healthy vs disease comparison
    if (!group_info$pattern_type %in% c("healthy_vs_disease", "healthy_vs_disease_general")) {
      cat("   ⚠️ Not a healthy vs disease comparison (", group_info$pattern_type, ")\n")
      cat("   🔄 Skipping for healthy vs disease analysis\n\n")
      next
    }
    
    cat("   ✅ HEALTHY vs DISEASE pattern confirmed!\n")
    cat("   📊 Groups:", paste(names(table(group_info$groups)), collapse = " vs "), "\n")
    cat("   👥 Sample distribution:\n")
    group_table <- table(group_info$groups)
    for (group_name in names(group_table)) {
      cat(sprintf("      %-12s: %3d samples\n", group_name, group_table[group_name]))
    }
    
    # Filter to samples with valid groups
    if (!is.null(group_info$sample_indices)) {
      valid_samples <- group_info$sample_indices
      filtered_expression <- dataset$expression_data[, valid_samples]
      filtered_phenotype <- dataset$phenotype_data[valid_samples, ]
      groups <- group_info$groups
    } else {
      filtered_expression <- dataset$expression_data
      filtered_phenotype <- dataset$phenotype_data
      groups <- group_info$groups
    }
    
    # Filter to CAMK genes only
    gene_symbols <- rownames(filtered_expression)
    camk_indices <- which(gene_symbols %in% camk_genes)
    
    if (length(camk_indices) == 0) {
      cat("   ❌ No CAMK genes found in expression data\n\n")
      next
    }
    
    camk_expression <- filtered_expression[camk_indices, ]
    camk_genes_present <- gene_symbols[camk_indices]
    
    cat("   🧬 CAMK genes found:", length(camk_genes_present), "/", length(camk_genes), "\n")
    cat("      ", paste(camk_genes_present, collapse = ", "), "\n")
    
    # Perform differential expression analysis
    cat("   🔬 Running differential expression analysis...\n")
    
    # Set up design matrix - Disease vs Healthy (so positive logFC = upregulated in disease)
    groups_factor <- factor(groups, levels = c("Healthy", "Disease"))
    design_matrix <- model.matrix(~ groups_factor)
    colnames(design_matrix) <- c("Intercept", "Disease_vs_Healthy")
    
    # Fit linear model
    fit <- lmFit(camk_expression, design_matrix)
    fit <- eBayes(fit)
    
    # Extract results for Disease vs Healthy comparison
    results <- topTable(fit, coef = "Disease_vs_Healthy", number = Inf, sort.by = "P")
    
    # Add gene symbols and format results
    results$Gene <- rownames(results)
    results$Dataset <- dataset_id
    results$N_Healthy <- sum(groups == "Healthy")
    results$N_Disease <- sum(groups == "Disease")
    results$Total_N <- length(groups)
    
    # Determine regulation direction
    results$Regulation <- ifelse(results$logFC > 0, "UP in Disease", "DOWN in Disease")
    results$Significant <- results$adj.P.Val < 0.05
    
    cat("   📈 Results summary:\n")
    cat("      Significant genes (FDR < 0.05):", sum(results$Significant), "\n")
    if (sum(results$Significant) > 0) {
      sig_results <- results[results$Significant, ]
      cat("      Significant CAMK genes:\n")
      for (i in 1:nrow(sig_results)) {
        cat(sprintf("        • %-8s: logFC=%.3f, P=%.2e (%s)\n", 
                   sig_results$Gene[i], sig_results$logFC[i], 
                   sig_results$adj.P.Val[i], sig_results$Regulation[i]))
      }
    }
    
    # Store results
    healthy_disease_results[[dataset_id]] <- list(
      results = results,
      groups = group_info,
      camk_genes_present = camk_genes_present,
      n_healthy = sum(groups == "Healthy"),
      n_disease = sum(groups == "Disease"),
      total_samples = length(groups)
    )
    
    cat("   ✅ Analysis completed successfully\n\n")
    
  }, error = function(e) {
    cat("   ❌ Error processing dataset:", e$message, "\n\n")
  })
}

# Summary of healthy vs disease datasets
cat("📊 HEALTHY vs DISEASE ANALYSIS SUMMARY\n")
cat("=====================================\n\n")

if (length(healthy_disease_results) == 0) {
  cat("❌ No datasets with healthy vs disease comparisons analyzed successfully\n")
} else {
  cat("✅ Datasets with healthy vs disease comparisons:\n\n")
  
  for (dataset_id in names(healthy_disease_results)) {
    result_info <- healthy_disease_results[[dataset_id]]
    cat(sprintf("📋 %-12s: %3d healthy + %3d disease = %3d total samples\n", 
               dataset_id, result_info$n_healthy, result_info$n_disease, result_info$total_samples))
    cat(sprintf("   🧬 CAMK genes: %d/11 available\n", length(result_info$camk_genes_present)))
    
    sig_genes <- sum(result_info$results$Significant)
    if (sig_genes > 0) {
      cat(sprintf("   📈 Significant: %d CAMK genes (FDR < 0.05)\n", sig_genes))
      sig_results <- result_info$results[result_info$results$Significant, ]
      for (i in 1:nrow(sig_results)) {
        cat(sprintf("      • %-8s: %.3f (%s)\n", 
                   sig_results$Gene[i], sig_results$logFC[i], sig_results$Regulation[i]))
      }
    } else {
      cat("   📈 No significant CAMK genes (FDR < 0.05)\n")
    }
    cat("\n")
  }
}

# Save results
output_file <- "output/CAMK_healthy_vs_disease_analysis_results.rds"
saveRDS(healthy_disease_results, output_file)
cat("💾 Results saved to:", output_file, "\n")

# Create meta-analysis if multiple datasets available
if (length(healthy_disease_results) > 1) {
  cat("\n🔬 META-ANALYSIS OF HEALTHY vs DISEASE STUDIES\n")
  cat("=============================================\n\n")
  
  # Prepare data for meta-analysis
  all_genes <- unique(unlist(lapply(healthy_disease_results, function(x) x$camk_genes_present)))
  
  meta_results <- data.frame()
  
  for (gene in all_genes) {
    gene_data <- data.frame()
    
    for (dataset_id in names(healthy_disease_results)) {
      result_info <- healthy_disease_results[[dataset_id]]
      
      if (gene %in% result_info$results$Gene) {
        gene_result <- result_info$results[result_info$results$Gene == gene, ]
        
        gene_data <- rbind(gene_data, data.frame(
          Dataset = dataset_id,
          logFC = gene_result$logFC,
          SE = sqrt(gene_result$t^2 / (gene_result$t^2 + result_info$total_samples - 2)) * abs(gene_result$logFC) / abs(gene_result$t),
          P_Value = gene_result$P.Value,
          N = result_info$total_samples
        ))
      }
    }
    
    if (nrow(gene_data) >= 2) {
      # Perform random-effects meta-analysis
      tryCatch({
        meta_model <- rma(yi = gene_data$logFC, sei = gene_data$SE, method = "REML")
        
        meta_results <- rbind(meta_results, data.frame(
          Gene = gene,
          N_Studies = nrow(gene_data),
          Combined_logFC = meta_model$beta[1],
          Combined_SE = meta_model$se,
          Combined_P_Value = meta_model$pval,
          CI_Lower = meta_model$ci.lb,
          CI_Upper = meta_model$ci.ub,
          Heterogeneity_I2 = meta_model$I2,
          Heterogeneity_P = ifelse(is.null(meta_model$QEp), NA, meta_model$QEp),
          Regulation = ifelse(meta_model$beta[1] > 0, "UP in Disease", "DOWN in Disease"),
          Significant = meta_model$pval < 0.05,
          Datasets = paste(gene_data$Dataset, collapse = ", ")
        ))
      }, error = function(e) {
        cat("   ⚠️ Meta-analysis failed for", gene, ":", e$message, "\n")
      })
    }
  }
  
  if (nrow(meta_results) > 0) {
    # Sort by p-value
    meta_results <- meta_results[order(meta_results$Combined_P_Value), ]
    
    cat("📊 Meta-analysis results (healthy vs disease only):\n")
    cat("=================================================\n\n")
    
    significant_genes <- meta_results[meta_results$Significant, ]
    if (nrow(significant_genes) > 0) {
      cat("✅ Significantly dysregulated CAMK genes:\n")
      for (i in 1:nrow(significant_genes)) {
        cat(sprintf("   • %-8s: logFC=%.3f, P=%.2e (%s, %d studies)\n",
                   significant_genes$Gene[i], significant_genes$Combined_logFC[i],
                   significant_genes$Combined_P_Value[i], significant_genes$Regulation[i],
                   significant_genes$N_Studies[i]))
      }
    } else {
      cat("❌ No significantly dysregulated CAMK genes in meta-analysis\n")
    }
    
    # Save meta-analysis results
    meta_file <- "output/CAMK_healthy_vs_disease_meta_analysis.csv"
    write.csv(meta_results, meta_file, row.names = FALSE)
    cat("\n💾 Meta-analysis results saved to:", meta_file, "\n")
    
    # Save complete results
    complete_meta_file <- "output/CAMK_healthy_vs_disease_meta_complete.rds"
    saveRDS(list(
      individual_results = healthy_disease_results,
      meta_analysis = meta_results
    ), complete_meta_file)
    cat("💾 Complete results saved to:", complete_meta_file, "\n")
  }
}

cat("\n✅ CAMK HEALTHY vs DISEASE ANALYSIS COMPLETED!\n")
cat("==============================================\n")
cat("🎯 This analysis focused exclusively on healthy vs disease comparisons\n")
cat("🔬 Meta-analysis combined only compatible study designs\n")
cat("📈 Results represent true disease dysregulation vs healthy controls\n")