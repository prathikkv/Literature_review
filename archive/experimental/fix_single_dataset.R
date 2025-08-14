#!/usr/bin/env Rscript
#' Fix Single Dataset Gene Mapping
#' 
#' Process one dataset at a time to avoid timeout

# Load functions
source("functions/data_processing.R")
source("functions/camk_definitions.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  cat("Usage: Rscript fix_single_dataset.R <dataset_id>\n")
  quit(status = 1)
}

dataset_id <- args[1]
camk_genes <- get_camk_gene_categories()$core

cat("=== Fixing", dataset_id, "===\n")

# Find dataset file
dataset_file <- list.files("cache", pattern = paste0(dataset_id, "_processed.rds"), 
                         recursive = TRUE, full.names = TRUE)[1]

if (is.null(dataset_file) || !file.exists(dataset_file)) {
  cat("ERROR: Dataset file not found\n")
  quit(status = 1)
}

# Load dataset
dataset <- readRDS(dataset_file)

if (!dataset$success || is.null(dataset$expression_matrix)) {
  cat("ERROR: Dataset failed or no expression matrix\n")
  quit(status = 1)
}

expr_matrix <- dataset$expression_matrix
cat("Original:", nrow(expr_matrix), "genes x", ncol(expr_matrix), "samples\n")

# Check CAMK genes before mapping
camk_before <- intersect(rownames(expr_matrix), camk_genes)
cat("CAMK genes before:", length(camk_before), "/", length(camk_genes), "\n")

if (length(camk_before) == length(camk_genes)) {
  cat("All CAMK genes already present\n")
  quit(status = 0)
}

# Apply gene mapping
cat("Applying gene symbol mapping...\n")
mapped_matrix <- apply_gene_symbol_mapping(expr_matrix, platform = "GPL570")

cat("After mapping:", nrow(mapped_matrix), "genes x", ncol(mapped_matrix), "samples\n")

# Check CAMK genes after mapping
camk_after <- intersect(rownames(mapped_matrix), camk_genes)
cat("CAMK genes after:", length(camk_after), "/", length(camk_genes), "\n")

if (length(camk_after) > length(camk_before)) {
  # Update dataset
  dataset$expression_matrix <- mapped_matrix
  dataset$n_genes <- nrow(mapped_matrix)
  
  # Save updated dataset
  saveRDS(dataset, dataset_file)
  cat("SUCCESS:", dataset_id, "updated with", nrow(mapped_matrix), "genes and", length(camk_after), "CAMK genes\n")
} else {
  cat("WARNING: No improvement after mapping\n")
}

cat("=== Complete ===\n")