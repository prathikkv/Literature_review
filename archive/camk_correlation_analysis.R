#!/usr/bin/env Rscript
#' CAMK Correlation and Expression Pattern Analysis
#' 
#' Analyzes inter-gene correlations and expression patterns within CAMK family

library(corrplot)
library(pheatmap)
library(ggplot2)
library(reshape2)

cat("🔗 CAMK CORRELATION & EXPRESSION PATTERN ANALYSIS\n")
cat("==================================================\n\n")

# Load CAMK-focused results
if (!file.exists("output/CAMK_focused_analysis_results.rds")) {
  stop("❌ CAMK analysis results not found. Please run camk_focused_analysis.R first.")
}

camk_results <- readRDS("output/CAMK_focused_analysis_results.rds")

cat("📊 Loaded CAMK analysis results for", length(camk_results), "datasets\n\n")

# Function to create CAMK correlation matrix
create_camk_correlation_matrix <- function(dataset_results, dataset_id) {
  
  if (is.null(dataset_results$camk_expression_matrix)) {
    return(NULL)
  }
  
  expr_matrix <- dataset_results$camk_expression_matrix
  
  if (nrow(expr_matrix) < 2 || ncol(expr_matrix) < 3) {
    cat("   ⚠️ Insufficient data for correlation analysis\n")
    return(NULL)
  }
  
  # Calculate correlation matrix
  cor_matrix <- cor(t(expr_matrix), use = "complete.obs")
  
  cat("   ✅ Correlation matrix created:", nrow(cor_matrix), "x", ncol(cor_matrix), "\n")
  
  return(cor_matrix)
}

# Function to create CAMK expression heatmap
create_camk_expression_heatmap <- function(dataset_results, dataset_id) {
  
  if (is.null(dataset_results$camk_expression_matrix)) {
    return(NULL)
  }
  
  expr_matrix <- dataset_results$camk_expression_matrix
  
  # Check if we have group information for better annotation
  sample_annotation <- NULL
  group_colors <- NULL
  
  if (!is.null(dataset_results$groups)) {
    # Handle different group data structures
    if (is.list(dataset_results$groups) && !is.null(dataset_results$groups$groups)) {
      groups <- dataset_results$groups$groups
    } else if (is.factor(dataset_results$groups) || is.character(dataset_results$groups)) {
      groups <- dataset_results$groups
    } else {
      groups <- NULL
    }
    
    if (!is.null(groups) && length(groups) == ncol(expr_matrix)) {
      # Create annotation for samples
      sample_annotation <- data.frame(
        Group = groups,
        row.names = colnames(expr_matrix)
      )
      
      # Create color mapping for groups
      group_colors <- list(Group = c("Healthy" = "#2E86C1", "Disease" = "#E74C3C", "AF" = "#E74C3C", "SR" = "#2E86C1"))
    }
  }
  
  # Create heatmap
  tryCatch({
    heatmap_file <- paste0("output/CAMK_expression_heatmap_", dataset_id, ".png")
    
    png(heatmap_file, width = 800, height = 600, res = 100)
    
    pheatmap(expr_matrix,
             main = paste("CAMK Family Expression -", dataset_id),
             cluster_rows = TRUE,
             cluster_cols = TRUE,
             scale = "row",
             annotation_col = sample_annotation,
             annotation_colors = group_colors,
             show_colnames = FALSE,
             fontsize_row = 10,
             color = colorRampPalette(c("#2E86C1", "white", "#E74C3C"))(100))
    
    dev.off()
    
    cat("   📊 Heatmap saved:", heatmap_file, "\n")
    
  }, error = function(e) {
    cat("   ❌ Heatmap creation error:", e$message, "\n")
  })
}

# Function to analyze CAMK expression patterns
analyze_camk_expression_patterns <- function(dataset_results, dataset_id) {
  
  if (is.null(dataset_results$camk_expression_matrix) || is.null(dataset_results$groups)) {
    return(NULL)
  }
  
  expr_matrix <- dataset_results$camk_expression_matrix
  
  # Handle different group data structures
  if (is.list(dataset_results$groups) && !is.null(dataset_results$groups$groups)) {
    groups <- dataset_results$groups$groups
  } else if (is.factor(dataset_results$groups) || is.character(dataset_results$groups)) {
    groups <- dataset_results$groups
  } else {
    return(NULL)
  }
  
  # Calculate mean expression by group for each CAMK gene
  group_means <- list()
  
  for (group in unique(groups)) {
    group_samples <- which(groups == group)
    if (length(group_samples) > 0) {
      group_means[[group]] <- apply(expr_matrix[, group_samples, drop = FALSE], 1, mean, na.rm = TRUE)
    }
  }
  
  if (length(group_means) < 2) {
    return(NULL)
  }
  
  # Create comparison data frame
  pattern_data <- data.frame(
    Gene = rownames(expr_matrix),
    stringsAsFactors = FALSE
  )
  
  for (group_name in names(group_means)) {
    pattern_data[[group_name]] <- group_means[[group_name]]
  }
  
  # Calculate fold changes (if we have two groups)
  if (length(group_means) == 2) {
    group_names <- names(group_means)
    pattern_data$log2_FoldChange <- log2(group_means[[group_names[1]]] / group_means[[group_names[2]]])
    pattern_data$FoldChange_Direction <- ifelse(pattern_data$log2_FoldChange > 0, 
                                               paste(group_names[1], "Higher"), 
                                               paste(group_names[2], "Higher"))
  }
  
  return(pattern_data)
}

# Process each dataset
correlation_matrices <- list()
expression_patterns <- list()

for (dataset_id in names(camk_results)) {
  cat("🔬 Processing", dataset_id, "\n")
  
  dataset_results <- camk_results[[dataset_id]]
  
  # Skip if no CAMK genes or no expression data
  if (length(dataset_results$camk_genes_present) == 0) {
    cat("   ⚠️ No CAMK genes present - skipping\n\n")
    next
  }
  
  cat("   CAMK genes present:", length(dataset_results$camk_genes_present), "\n")
  cat("   ", paste(dataset_results$camk_genes_present, collapse = ", "), "\n")
  
  # Create correlation matrix
  cor_matrix <- create_camk_correlation_matrix(dataset_results, dataset_id)
  if (!is.null(cor_matrix)) {
    correlation_matrices[[dataset_id]] <- cor_matrix
    
    # Save correlation plot
    tryCatch({
      png(paste0("output/CAMK_correlation_", dataset_id, ".png"), width = 600, height = 600, res = 100)
      corrplot(cor_matrix, 
               method = "color", 
               type = "upper",
               order = "hclust",
               title = paste("CAMK Gene Correlations -", dataset_id),
               mar = c(0,0,2,0),
               tl.cex = 0.8,
               tl.col = "black")
      dev.off()
      cat("   📊 Correlation plot saved: output/CAMK_correlation_", dataset_id, ".png\n")
    }, error = function(e) {
      cat("   ❌ Correlation plot error:", e$message, "\n")
    })
  }
  
  # Create expression heatmap
  create_camk_expression_heatmap(dataset_results, dataset_id)
  
  # Analyze expression patterns
  patterns <- analyze_camk_expression_patterns(dataset_results, dataset_id)
  if (!is.null(patterns)) {
    expression_patterns[[dataset_id]] <- patterns
    
    cat("   📈 Expression pattern analysis completed\n")
    
    # Print top patterns if fold changes available
    if ("log2_FoldChange" %in% colnames(patterns)) {
      cat("      Top CAMK expression changes:\n")
      patterns_sorted <- patterns[order(abs(patterns$log2_FoldChange), decreasing = TRUE), ]
      for (i in 1:min(5, nrow(patterns_sorted))) {
        gene <- patterns_sorted$Gene[i]
        fc <- round(patterns_sorted$log2_FoldChange[i], 2)
        direction <- patterns_sorted$FoldChange_Direction[i]
        cat(sprintf("        %-8s: log2FC = %6.2f [%s]\n", gene, fc, direction))
      }
    }
  }
  
  cat("\n")
}

# Create combined CAMK correlation analysis
cat("🔗 COMBINED CAMK CORRELATION ANALYSIS\n")
cat("====================================\n")

if (length(correlation_matrices) > 0) {
  # Calculate average correlation across datasets
  all_genes <- unique(unlist(lapply(correlation_matrices, rownames)))
  
  cat("CAMK genes analyzed across datasets:", length(all_genes), "\n")
  cat("Datasets with correlation data:", length(correlation_matrices), "\n\n")
  
  # Find most consistent correlations
  consistent_correlations <- list()
  
  for (i in 1:(length(all_genes)-1)) {
    for (j in (i+1):length(all_genes)) {
      gene1 <- all_genes[i]
      gene2 <- all_genes[j]
      
      # Get correlations from all datasets where both genes are present
      cors <- c()
      for (dataset_id in names(correlation_matrices)) {
        cor_mat <- correlation_matrices[[dataset_id]]
        if (gene1 %in% rownames(cor_mat) && gene2 %in% colnames(cor_mat)) {
          cors <- c(cors, cor_mat[gene1, gene2])
        }
      }
      
      if (length(cors) >= 2) {
        avg_cor <- mean(cors, na.rm = TRUE)
        sd_cor <- sd(cors, na.rm = TRUE)
        n_datasets <- length(cors)
        
        consistent_correlations[[paste(gene1, gene2, sep = "_vs_")]] <- list(
          gene1 = gene1,
          gene2 = gene2,
          mean_correlation = avg_cor,
          sd_correlation = sd_cor,
          n_datasets = n_datasets,
          correlations = cors
        )
      }
    }
  }
  
  # Sort by absolute mean correlation
  consistent_correlations <- consistent_correlations[order(sapply(consistent_correlations, function(x) abs(x$mean_correlation)), decreasing = TRUE)]
  
  cat("🎯 TOP CAMK GENE CORRELATIONS (across multiple datasets):\n")
  cat("========================================================\n")
  
  for (i in 1:min(10, length(consistent_correlations))) {
    corr_info <- consistent_correlations[[i]]
    
    cat(sprintf("%-8s vs %-8s: r=%.3f ± %.3f (%d datasets)\n",
                corr_info$gene1, corr_info$gene2, 
                corr_info$mean_correlation, corr_info$sd_correlation,
                corr_info$n_datasets))
  }
}

# Save all results
save_data <- list(
  correlation_matrices = correlation_matrices,
  expression_patterns = expression_patterns,
  consistent_correlations = if(exists("consistent_correlations")) consistent_correlations else NULL
)

saveRDS(save_data, "output/CAMK_correlation_expression_analysis.rds")

# Create summary report
cat(sprintf("\n📋 SUMMARY REPORT\n"))
cat("================\n")
cat(sprintf("Datasets analyzed: %d\n", length(camk_results)))
cat(sprintf("Correlation matrices created: %d\n", length(correlation_matrices)))
cat(sprintf("Expression pattern analyses: %d\n", length(expression_patterns)))

if (exists("consistent_correlations")) {
  cat(sprintf("Consistent correlations identified: %d\n", length(consistent_correlations)))
}

cat("\n💾 Files created:\n")
camk_files <- list.files("output", pattern = "CAMK_.*\\.(png|rds|csv)$", full.names = FALSE)
for (file in camk_files) {
  cat("   •", file, "\n")
}

cat("\n🎉 CAMK correlation and expression pattern analysis completed!\n")
cat("✨ Ready for clinical interpretation and publication\n")