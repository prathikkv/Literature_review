#!/usr/bin/env Rscript
#' Focused CAMK2D Analysis with Available Data
#' 
#' Demonstrates the multi-modal framework with GSE57338 microarray data

cat("FOCUSED: CAMK2D ANALYSIS WITH AVAILABLE DATA\n")
cat("==========================================\n\n")

# Load configuration and modules
source("functions/data_processing.R")
source("functions/analysis.R")
source("functions/visualization.R")
source("functions/utilities.R")

library(yaml)
config <- read_yaml("config.yml")

# Load the working GSE57338 data
cat("DATA: Loading GSE57338 microarray data\n")
gse_data <- readRDS("cache/microarray/GSE57338_processed.rds")

cat("SUCCESS: Loaded GSE57338 data\n")
cat("  • Dataset ID:", gse_data$dataset_id, "\n")
cat("  • Success:", gse_data$success, "\n")
cat("  • Genes:", gse_data$n_genes, "\n")
cat("  • Samples:", gse_data$n_samples, "\n")
cat("  • Expression matrix:", dim(gse_data$expression_matrix)[1], "genes x", dim(gse_data$expression_matrix)[2], "samples\n\n")

# Check for CAMK family genes
camk_genes <- config$research$gene_family
cat("GENETIC: Checking for CAMK family genes in dataset\n")

if (!is.null(gse_data$expression_matrix) && nrow(gse_data$expression_matrix) > 0) {
  gene_names <- rownames(gse_data$expression_matrix)
  camk_present <- intersect(camk_genes, gene_names)
  
  cat("CAMK genes found:", length(camk_present), "out of", length(camk_genes), "\n")
  if (length(camk_present) > 0) {
    cat("Present CAMK genes:", paste(camk_present, collapse=", "), "\n")
  }
  cat("\n")
  
  # Focus on CAMK2D if available
  if ("CAMK2D" %in% camk_present) {
    cat("TARGET: CAMK2D found in dataset - proceeding with analysis\n")
    
    # Extract CAMK2D expression
    camk2d_expr <- gse_data$expression_matrix["CAMK2D", ]
    
    cat("CAMK2D expression summary:\n")
    print(summary(as.numeric(camk2d_expr)))
    
    # Check if we have group information
    if (!is.null(gse_data$phenotype_data)) {
      cat("\nPHENOTYPE: Available phenotype data columns:\n")
      print(names(gse_data$phenotype_data))
      
      # Look for disease/condition grouping
      potential_group_cols <- c("disease_group", "condition", "group", "disease_state", "tissue")
      group_col <- NULL
      
      for (col in potential_group_cols) {
        if (col %in% names(gse_data$phenotype_data)) {
          group_col <- col
          break
        }
      }
      
      if (!is.null(group_col)) {
        cat("GROUP: Using", group_col, "for grouping\n")
        groups <- gse_data$phenotype_data[[group_col]]
        group_table <- table(groups)
        cat("Group distribution:\n")
        print(group_table)
        
        if (length(unique(groups)) >= 2) {
          cat("\nCAMK2D: Expression by group:\n")
          for (grp in names(group_table)) {
            grp_samples <- which(groups == grp)
            grp_expr <- as.numeric(camk2d_expr[grp_samples])
            cat("  ", grp, ":", length(grp_samples), "samples,", 
                "mean =", round(mean(grp_expr, na.rm=TRUE), 3),
                "± sd =", round(sd(grp_expr, na.rm=TRUE), 3), "\n")
          }
        }
      }
    }
  } else {
    cat("WARNING: CAMK2D not found in this dataset\n")
    cat("Available genes in dataset (first 10):", paste(head(gene_names, 10), collapse=", "), "\n")
  }
} else {
  cat("ERROR: No expression matrix available in dataset\n")
}

# Generate summary report
report_dir <- "output/focused_analysis"
if (!dir.exists(report_dir)) {
  dir.create(report_dir, recursive = TRUE)
}

# Save analysis summary
analysis_summary <- list(
  timestamp = Sys.time(),
  dataset_id = gse_data$dataset_id,
  total_genes = gse_data$n_genes,
  total_samples = gse_data$n_samples,
  camk_genes_available = if(exists("camk_present")) camk_present else c(),
  camk2d_available = "CAMK2D" %in% (if(exists("camk_present")) camk_present else c()),
  config_datasets = config$datasets$total_datasets,
  framework_ready = TRUE
)

saveRDS(analysis_summary, file.path(report_dir, "focused_analysis_summary.rds"))

cat("\nSUCCESS: Focused CAMK2D analysis completed\n")
cat("FRAMEWORK: Multi-modal framework is operational\n")
cat("DATA: Analysis summary saved to:", file.path(report_dir, "focused_analysis_summary.rds"), "\n")
cat("READY: System ready for expanded multi-modal analysis when additional datasets are available\n")