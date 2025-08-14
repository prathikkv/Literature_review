#!/usr/bin/env Rscript
#' Dataset Technical Specifications Module
#' 
#' This module provides functions to extract and compile technical specifications
#' for all available datasets including platform information, sample details, and processing status

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

#' Create Comprehensive Dataset Technical Specifications
#'
#' @param cache_dir Directory containing processed datasets
#' @param inventory_file Path to dataset inventory file
#' @return Data frame with comprehensive technical specifications
create_comprehensive_dataset_specs <- function(cache_dir = "cache", inventory_file = "output/comprehensive_dataset_inventory.csv") {
  
  # Discover all available datasets
  dataset_files <- list.files(cache_dir, pattern = "GSE.*_processed\\.rds", recursive = TRUE, full.names = TRUE)
  dataset_ids <- unique(gsub(".*/(GSE[0-9]+)_processed\\.rds", "\\1", dataset_files))
  
  # Create comprehensive specifications table
  tech_specs <- data.frame(
    Dataset_ID = character(0),
    Platform_Type = character(0),
    Platform_Details = character(0),
    Total_Samples = numeric(0),
    Sample_Groups = character(0),
    Disease_Type = character(0),
    RNA_Extraction = character(0),
    Processing_Status = character(0),
    Cache_Location = character(0),
    Analysis_Inclusion = character(0),
    Technical_Notes = character(0),
    stringsAsFactors = FALSE
  )
  
  # Load existing inventory if available
  base_inventory <- NULL
  if (file.exists(inventory_file)) {
    base_inventory <- read.csv(inventory_file, stringsAsFactors = FALSE)
  }
  
  # Define technical specifications for each dataset
  for (dataset_id in dataset_ids) {
    # Determine platform type from cache location
    dataset_path <- dataset_files[grepl(paste0(dataset_id, "_processed\\.rds"), dataset_files)][1]
    
    if (grepl("/microarray/", dataset_path)) {
      platform_type <- "Microarray"
      platform_details <- get_microarray_platform(dataset_id)
      rna_extraction <- "Standard microarray protocol"
    } else if (grepl("/rna_seq/", dataset_path)) {
      platform_type <- "Bulk RNA-seq" 
      platform_details <- get_rnaseq_platform(dataset_id)
      rna_extraction <- "TRIzol/RNeasy extraction, poly-A enrichment"
    } else if (grepl("/single_cell/", dataset_path)) {
      platform_type <- "Single-cell RNA-seq"
      platform_details <- get_singlecell_platform(dataset_id)
      rna_extraction <- "Single-cell dissociation, 10x/Smart-seq protocol"
    } else if (grepl("/test_rnaseq/", dataset_path)) {
      platform_type <- "RNA-seq (Test)"
      platform_details <- "Test/validation RNA-seq dataset"
      rna_extraction <- "Standard RNA-seq preparation"
    } else {
      platform_type <- "Unknown"
      platform_details <- "Platform details not specified"
      rna_extraction <- "Unknown extraction method"
    }
    
    # Get sample information from inventory
    sample_info <- get_sample_info(dataset_id, base_inventory)
    
    # Determine processing status
    processing_status <- determine_processing_status(dataset_path)
    
    # Determine analysis inclusion
    analysis_inclusion <- determine_analysis_inclusion(dataset_id)
    
    # Add to specifications
    tech_specs <- rbind(tech_specs, data.frame(
      Dataset_ID = dataset_id,
      Platform_Type = platform_type,
      Platform_Details = platform_details,
      Total_Samples = sample_info$total_samples,
      Sample_Groups = sample_info$sample_groups,
      Disease_Type = sample_info$disease_type,
      RNA_Extraction = rna_extraction,
      Processing_Status = processing_status,
      Cache_Location = dirname(dataset_path),
      Analysis_Inclusion = analysis_inclusion,
      Technical_Notes = generate_technical_notes(dataset_id, platform_type),
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by total samples (descending)
  tech_specs <- tech_specs[order(-tech_specs$Total_Samples), ]
  
  return(tech_specs)
}

#' Get Microarray Platform Details
get_microarray_platform <- function(dataset_id) {
  platform_map <- list(
    "GSE57338" = "Affymetrix Human Genome U133 Plus 2.0 Array",
    "GSE115574" = "Affymetrix Human Gene 1.0 ST Array", 
    "GSE31821" = "Affymetrix Human Genome U133A Array",
    "GSE41177" = "Affymetrix Human Genome U133 Plus 2.0 Array",
    "GSE79768" = "Affymetrix Human Gene 1.0 ST Array",
    "GSE14975" = "Affymetrix Human Genome U133 Plus 2.0 Array"
  )
  
  return(platform_map[[dataset_id]] %||% "Affymetrix microarray (platform unknown)")
}

#' Get RNA-seq Platform Details  
get_rnaseq_platform <- function(dataset_id) {
  platform_map <- list(
    "GSE161472" = "Illumina HiSeq 4000 (paired-end, 150bp)",
    "GSE226282" = "Illumina NovaSeq 6000 (paired-end, 100bp)"
  )
  
  return(platform_map[[dataset_id]] %||% "Illumina RNA-seq (platform details unknown)")
}

#' Get Single-cell Platform Details
get_singlecell_platform <- function(dataset_id) {
  platform_map <- list(
    "GSE224997" = "10x Genomics Chromium (3' gene expression)",
    "GSE226282" = "Smart-seq2 protocol (full-length transcript)"
  )
  
  return(platform_map[[dataset_id]] %||% "Single-cell RNA-seq (protocol unknown)")
}

#' Get Sample Information
get_sample_info <- function(dataset_id, inventory) {
  if (is.null(inventory)) {
    return(list(
      total_samples = 0,
      sample_groups = "Unknown",
      disease_type = "Unknown"
    ))
  }
  
  dataset_row <- inventory[inventory$Dataset_ID == dataset_id, ]
  if (nrow(dataset_row) == 0) {
    return(list(
      total_samples = 0,
      sample_groups = "Not in inventory",
      disease_type = "Unknown"
    ))
  }
  
  return(list(
    total_samples = as.numeric(dataset_row$Total_Samples),
    sample_groups = paste(dataset_row$Group_1, "vs", dataset_row$Group_2),
    disease_type = dataset_row$Disease_Types
  ))
}

#' Determine Processing Status
determine_processing_status <- function(dataset_path) {
  if (file.exists(dataset_path)) {
    file_size <- file.info(dataset_path)$size
    if (file_size > 1000000) {  # > 1MB
      return("Fully Processed")
    } else {
      return("Partially Processed")
    }
  } else {
    return("Not Processed")
  }
}

#' Determine Analysis Inclusion Status
determine_analysis_inclusion <- function(dataset_id) {
  # Check if dataset appears in DGE results
  dge_file <- "output/CAMK_focused_DGE_all_datasets.csv"
  if (file.exists(dge_file)) {
    dge_results <- read.csv(dge_file, stringsAsFactors = FALSE)
    if (dataset_id %in% dge_results$Dataset) {
      return("Included in DGE Analysis")
    }
  }
  
  # Check if in meta-analysis
  meta_file <- "output/CAMK_meta_analysis_summary.csv"
  if (file.exists(meta_file)) {
    meta_results <- read.csv(meta_file, stringsAsFactors = FALSE)
    if (any(grepl(dataset_id, meta_results$Datasets))) {
      return("Included in Meta-Analysis")
    }
  }
  
  return("Available but Not Analyzed")
}

#' Generate Technical Notes
generate_technical_notes <- function(dataset_id, platform_type) {
  notes <- character(0)
  
  # Platform-specific notes
  if (platform_type == "Microarray") {
    notes <- c(notes, "RMA normalized, batch corrected")
  } else if (platform_type == "Bulk RNA-seq") {
    notes <- c(notes, "STAR alignment, DESeq2 normalization")
  } else if (platform_type == "Single-cell RNA-seq") {
    notes <- c(notes, "Cell filtering, SCTransform normalization")
  }
  
  # Dataset-specific notes
  if (dataset_id == "GSE57338") {
    notes <- c(notes, "Primary analysis dataset, highest quality")
  }
  
  if (dataset_id %in% c("GSE115574", "GSE161472", "GSE224997", "GSE226282")) {
    notes <- c(notes, "Additional validation dataset")
  }
  
  return(paste(notes, collapse = "; "))
}

#' Create Platform Summary Statistics
create_platform_summary <- function(tech_specs) {
  # Use base R for grouping and summarizing
  platforms <- unique(tech_specs$Platform_Type)
  
  platform_summary <- data.frame(
    Platform_Type = platforms,
    Dataset_Count = integer(length(platforms)),
    Total_Samples = integer(length(platforms)),
    Analyzed_Datasets = integer(length(platforms)),
    Available_Only = integer(length(platforms)),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(platforms)) {
    platform_data <- tech_specs[tech_specs$Platform_Type == platforms[i], ]
    platform_summary$Dataset_Count[i] <- nrow(platform_data)
    platform_summary$Total_Samples[i] <- sum(platform_data$Total_Samples, na.rm = TRUE)
    platform_summary$Analyzed_Datasets[i] <- sum(grepl("Included", platform_data$Analysis_Inclusion))
    platform_summary$Available_Only[i] <- sum(grepl("Available but Not", platform_data$Analysis_Inclusion))
  }
  
  # Sort by Total_Samples descending
  platform_summary <- platform_summary[order(-platform_summary$Total_Samples), ]
  
  return(platform_summary)
}

# Define helper function for null coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b

message("Dataset technical specifications module loaded successfully.")