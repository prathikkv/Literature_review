#!/usr/bin/env Rscript
#' CAMK Gene Definitions - Central Configuration
#' 
#' Centralized definitions for all CAMK family genes and related pathway genes
#' to eliminate redundancy across analysis scripts

#' Get Core CAMK Family Genes
#' 
#' Returns the 11 core CAMK family genes as specified in original analysis
#' @return Character vector of core CAMK gene symbols
get_camk_core_genes <- function() {
  c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G", 
    "CAMKK1", "CAMKK2", "CAMK1", "CAMK1D", "CAMK1G", "CAMK4", "CAMKV")
}

#' Get CAMK Pathway Context Genes
#' 
#' Returns genes involved in CAMK-related pathways to provide biological context
#' @return Character vector of pathway gene symbols
get_camk_pathway_genes <- function() {
  c(
    # Calcium signaling pathway
    "CALM1", "CALM2", "CALM3", "CALR", "CACNA1C", "CACNA1D", 
    "RYR2", "PLN", "SERCA2A", "ATP2A2",
    
    # Cardiac contraction pathway
    "MYH6", "MYH7", "TNNT2", "TNNI3", "TPM1", "ACTC1",
    
    # CAMK substrates/targets
    "HDAC4", "HDAC5", "MEF2A", "MEF2C", "CREB1", "FOXO1"
  )
}

#' Get Extended CAMK Gene Set
#' 
#' Returns combined core CAMK + pathway genes for enhanced biological context
#' @return Character vector of all genes (core + pathway)
get_camk_extended_genes <- function() {
  c(get_camk_core_genes(), get_camk_pathway_genes())
}

#' Get CAMK Gene Categories
#' 
#' Returns categorized gene lists for specialized analyses
#' @return Named list with gene categories
get_camk_gene_categories <- function() {
  list(
    core = get_camk_core_genes(),
    pathway = get_camk_pathway_genes(),
    extended = get_camk_extended_genes(),
    
    # Specialized subcategories
    camk2_family = c("CAMK2D", "CAMK2A", "CAMK2B", "CAMK2G"),
    calcium_signaling = c("CALM1", "CALM2", "CALM3", "CALR", "CACNA1C", "CACNA1D", "RYR2", "PLN", "ATP2A2"),
    cardiac_contraction = c("MYH6", "MYH7", "TNNT2", "TNNI3", "TPM1", "ACTC1"),
    transcription_factors = c("HDAC4", "HDAC5", "MEF2A", "MEF2C", "CREB1", "FOXO1")
  )
}

#' Validate CAMK Gene Coverage
#' 
#' Check what percentage of CAMK genes are present in a dataset
#' @param available_genes Character vector of available gene symbols
#' @param gene_set Character vector of target genes to check (default: core CAMK)
#' @return List with coverage statistics
validate_camk_coverage <- function(available_genes, gene_set = get_camk_core_genes()) {
  
  genes_found <- gene_set[gene_set %in% available_genes]
  genes_missing <- gene_set[!gene_set %in% available_genes]
  
  list(
    total_target_genes = length(gene_set),
    genes_found = genes_found,
    genes_missing = genes_missing,
    coverage_count = length(genes_found),
    coverage_percentage = round(length(genes_found) / length(gene_set) * 100, 1),
    is_adequate = length(genes_found) >= (length(gene_set) * 0.5)  # At least 50% coverage
  )
}

#' Print CAMK Gene Summary
#' 
#' Utility function to display CAMK gene information
#' @param show_details Logical, whether to show detailed breakdown
print_camk_summary <- function(show_details = FALSE) {
  core_genes <- get_camk_core_genes()
  pathway_genes <- get_camk_pathway_genes()
  
  cat("GENETIC: CAMK Gene Definitions Summary\n")
  cat("===============================\n")
  cat("Core CAMK family genes:", length(core_genes), "\n")
  cat("Pathway context genes:", length(pathway_genes), "\n") 
  cat("Total extended set:", length(core_genes) + length(pathway_genes), "\n\n")
  
  if (show_details) {
    cat("SUMMARY: Core CAMK Genes:\n")
    cat("   ", paste(core_genes, collapse = ", "), "\n\n")
    
    categories <- get_camk_gene_categories()
    cat("SUMMARY: Gene Categories:\n")
    for (category in names(categories)) {
      if (!category %in% c("core", "pathway", "extended")) {
        cat("   ", category, ":", length(categories[[category]]), "genes\n")
      }
    }
  }
}

# Backward compatibility - create variables for existing scripts
camk_core_genes <- get_camk_core_genes()
camk_pathway_genes <- get_camk_pathway_genes() 
camk_extended_genes <- get_camk_extended_genes()

# Legacy variable names for compatibility
camk_genes <- camk_core_genes  # Most common usage
camk_related_genes <- camk_pathway_genes  # Alternative naming