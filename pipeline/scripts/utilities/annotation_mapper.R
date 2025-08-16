#!/usr/bin/env Rscript
#' Annotation Mapper Utilities
#' 
#' Functions to map probe IDs to gene symbols using platform annotation files

#' Load GPL570 Annotation
#'
#' Loads and parses the GPL570 platform annotation file
#' @param annotation_file Path to GPL570.soft.gz file
#' @return Data frame with probe ID to gene symbol mapping
load_gpl570_annotation <- function(annotation_file = "cache/downloads/GPL570.soft.gz") {
  
  if (!file.exists(annotation_file)) {
    stop("GPL570 annotation file not found: ", annotation_file)
  }
  
  cat("📋 Loading GPL570 annotation from:", annotation_file, "\n")
  
  # Read the file and find the annotation table
  lines <- readLines(annotation_file)
  
  # Find the start of the annotation table
  table_start <- which(grepl("!platform_table_begin", lines))
  if (length(table_start) == 0) {
    stop("Could not find platform table in GPL570 file")
  }
  
  table_start <- table_start[1] + 1  # Skip the header line
  
  # Find the end of the annotation table
  table_end <- which(grepl("!platform_table_end", lines))
  if (length(table_end) == 0) {
    table_end <- length(lines)
  } else {
    table_end <- table_end[1] - 1
  }
  
  cat("📊 Parsing annotation table: lines", table_start, "to", table_end, "\n")
  
  # Extract the header
  header_line <- lines[table_start]
  headers <- strsplit(header_line, "\t")[[1]]
  
  # Find Gene Symbol column (column 11 based on our analysis)
  gene_symbol_col <- which(headers == "Gene Symbol")
  if (length(gene_symbol_col) == 0) {
    stop("Could not find 'Gene Symbol' column in annotation")
  }
  
  cat("📋 Gene Symbol column found at position:", gene_symbol_col, "\n")
  
  # Parse the annotation data
  annotation_data <- list()
  probe_count <- 0
  
  for (i in (table_start + 1):table_end) {
    fields <- strsplit(lines[i], "\t")[[1]]
    
    if (length(fields) >= gene_symbol_col) {
      probe_id <- fields[1]
      gene_symbols <- fields[gene_symbol_col]
      
      # Handle multiple gene symbols (separated by ///)
      if (gene_symbols != "" && !is.na(gene_symbols) && gene_symbols != "---") {
        symbols <- strsplit(gene_symbols, " /// ")[[1]]
        
        # Use the first non-empty symbol
        for (symbol in symbols) {
          symbol <- trimws(symbol)
          if (symbol != "" && symbol != "---") {
            annotation_data[[probe_id]] <- symbol
            probe_count <- probe_count + 1
            break
          }
        }
      }
    }
  }
  
  cat("✅ Loaded", probe_count, "probe-to-gene mappings\n")
  
  # Convert to data frame
  probe_to_gene <- data.frame(
    ProbeID = names(annotation_data),
    GeneSymbol = unlist(annotation_data),
    stringsAsFactors = FALSE
  )
  
  return(probe_to_gene)
}

#' Map Probe IDs to Gene Symbols
#'
#' Maps probe IDs in expression matrix to gene symbols using annotation
#' @param expr_matrix Expression matrix with probe IDs as rownames
#' @param annotation_df Annotation data frame from load_gpl570_annotation()
#' @param method Method for handling multiple probes per gene ("mean", "max", "first")
#' @return Expression matrix with gene symbols as rownames
map_probes_to_genes <- function(expr_matrix, annotation_df, method = "mean") {
  
  cat("🧬 PROBE TO GENE SYMBOL MAPPING\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  probe_ids <- rownames(expr_matrix)
  cat("📊 Input matrix:", nrow(expr_matrix), "probes x", ncol(expr_matrix), "samples\n")
  
  # Find probes that have gene symbol annotations
  mapped_probes <- intersect(probe_ids, annotation_df$ProbeID)
  unmapped_probes <- setdiff(probe_ids, annotation_df$ProbeID)
  
  cat("📋 Probe mapping summary:\n")
  cat("   ✅ Mapped probes:", length(mapped_probes), "\n")
  cat("   ❌ Unmapped probes:", length(unmapped_probes), "\n")
  cat("   📊 Mapping coverage:", round(length(mapped_probes)/length(probe_ids)*100, 1), "%\n")
  
  if (length(mapped_probes) == 0) {
    warning("No probes could be mapped to gene symbols")
    return(expr_matrix)
  }
  
  # Create mapping for mapped probes
  probe_to_gene_map <- annotation_df[annotation_df$ProbeID %in% mapped_probes, ]
  
  # Filter expression matrix to mapped probes only
  mapped_expr <- expr_matrix[mapped_probes, , drop = FALSE]
  
  # Add gene symbols
  gene_symbols <- probe_to_gene_map$GeneSymbol[match(rownames(mapped_expr), probe_to_gene_map$ProbeID)]
  
  # Handle multiple probes per gene
  unique_genes <- unique(gene_symbols)
  cat("📊 Gene consolidation:\n")
  cat("   🧬 Unique genes after mapping:", length(unique_genes), "\n")
  
  # Create gene-level expression matrix
  gene_expr_matrix <- matrix(
    nrow = length(unique_genes), 
    ncol = ncol(mapped_expr),
    dimnames = list(unique_genes, colnames(mapped_expr))
  )
  
  consolidation_methods <- character()
  
  for (gene in unique_genes) {
    probe_indices <- which(gene_symbols == gene)
    
    if (length(probe_indices) == 1) {
      # Single probe for this gene
      gene_expr_matrix[gene, ] <- mapped_expr[probe_indices, ]
      consolidation_methods <- c(consolidation_methods, "single")
    } else {
      # Multiple probes for this gene
      probe_data <- mapped_expr[probe_indices, , drop = FALSE]
      
      if (method == "mean") {
        gene_expr_matrix[gene, ] <- colMeans(probe_data, na.rm = TRUE)
      } else if (method == "max") {
        gene_expr_matrix[gene, ] <- apply(probe_data, 2, max, na.rm = TRUE)
      } else if (method == "first") {
        gene_expr_matrix[gene, ] <- probe_data[1, ]
      } else {
        stop("Unknown consolidation method: ", method)
      }
      
      consolidation_methods <- c(consolidation_methods, paste0("multi(", length(probe_indices), ")"))
    }
  }
  
  # Summary of consolidation
  consolidation_summary <- table(consolidation_methods)
  cat("   📊 Consolidation summary:\n")
  for (i in 1:length(consolidation_summary)) {
    cat("      ", names(consolidation_summary)[i], ":", consolidation_summary[i], "genes\n")
  }
  
  cat("   🔧 Consolidation method:", method, "\n")
  
  cat("✅ MAPPING COMPLETED\n")
  cat("   📊 Final matrix:", nrow(gene_expr_matrix), "genes x", ncol(gene_expr_matrix), "samples\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  return(gene_expr_matrix)
}

#' Get CAMK Gene Coverage
#'
#' Checks which CAMK genes are present after mapping
#' @param gene_expr_matrix Expression matrix with gene symbols
#' @param camk_genes Vector of CAMK gene symbols to check
#' @return List with coverage information
check_camk_gene_coverage <- function(gene_expr_matrix, camk_genes) {
  
  available_genes <- rownames(gene_expr_matrix)
  
  present_camk <- intersect(available_genes, camk_genes)
  missing_camk <- setdiff(camk_genes, available_genes)
  
  coverage_pct <- length(present_camk) / length(camk_genes) * 100
  
  result <- list(
    present_genes = present_camk,
    missing_genes = missing_camk,
    coverage_count = length(present_camk),
    total_camk_genes = length(camk_genes),
    coverage_percentage = coverage_pct
  )
  
  cat("🧬 CAMK GENE COVERAGE ANALYSIS:\n")
  cat("   ✅ Present CAMK genes:", length(present_camk), "/", length(camk_genes), 
      "(", round(coverage_pct, 1), "%)\n")
  
  if (length(present_camk) > 0) {
    cat("   🧬 Present:", paste(present_camk, collapse = ", "), "\n")
  }
  
  if (length(missing_camk) > 0) {
    cat("   ❌ Missing:", paste(missing_camk, collapse = ", "), "\n")
  }
  
  return(result)
}

cat("✅ ANNOTATION MAPPER utilities loaded successfully\n")
cat("📋 Functions: load_gpl570_annotation(), map_probes_to_genes(), check_camk_gene_coverage()\n\n")