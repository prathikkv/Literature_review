#!/usr/bin/env Rscript
#' Comprehensive Multi-Dataset CAMK Analysis Investigation
#' 
#' Systematic investigation and classification of 15 datasets for CAMK family analysis
#' Integration framework for multi-technology validation (microarray + RNA-seq + scRNA-seq)

cat("SEARCH: COMPREHENSIVE DATASET INVESTIGATION FOR CAMK ANALYSIS\n")
cat("=====================================================\n\n")

# Load required libraries
required_packages <- c(
  "GEOquery", "tidyverse", "limma", "affy", "oligo", "pd.hg.u133.plus.2",
  "annotate", "hgu133plus2.db", "biomaRt", "org.Hs.eg.db", 
  "ReactomePA", "clusterProfiler", "WGCNA", "igraph", "corrplot",
  "pheatmap", "VennDiagram", "UpSetR", "ComplexHeatmap", "circlize"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    if (pkg %in% c("ReactomePA", "clusterProfiler", "org.Hs.eg.db", "annotate", 
                   "hgu133plus2.db", "pd.hg.u133.plus.2", "affy", "oligo", "GEOquery")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg, ask = FALSE)
    } else {
      install.packages(pkg, dependencies = TRUE)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

# Load utility functions
source("../../functions/data_processing.R")
source("../../functions/analysis.R") 
source("../../functions/utilities.R")

# =============================================================================
# DATASET CLASSIFICATION AND INVESTIGATION FRAMEWORK
# =============================================================================

#' Comprehensive Dataset Investigation Function
#'
#' @param dataset_ids Vector of GEO dataset IDs to investigate
#' @param output_dir Output directory for results
#' @return List containing dataset classifications and metadata
investigate_datasets_comprehensive <- function(dataset_ids, output_dir = "results/dataset_investigation") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  cat("SEARCH: Starting comprehensive investigation of", length(dataset_ids), "datasets\n\n")
  
  # Initialize results structure
  investigation_results <- list(
    dataset_metadata = list(),
    classification = data.frame(
      Dataset_ID = character(),
      Status = character(),
      Platform = character(), 
      Samples = numeric(),
      Data_Type = character(),
      Species = character(),
      Disease_Context = character(),
      Comparison_Type = character(),
      Scientific_Value = character(),
      CAMK_Relevance = character(),
      Integration_Tier = character(),
      stringsAsFactors = FALSE
    ),
    investigation_summary = list()
  )
  
  # Known high-value datasets (already integrated)
  known_datasets <- list(
    "GSE57338" = list(
      status = "INTEGRATED", platform = "GPL11532", samples = 313,
      data_type = "Microarray", species = "Human", 
      disease_context = "Heart Failure", comparison = "Healthy vs HF",
      value = "EXTREMELY HIGH", camk_relevance = "PRIMARY DISCOVERY",
      tier = "TIER 1"
    ),
    "GSE148507" = list(
      status = "INTEGRATED", platform = "10x Genomics", samples = 386,
      data_type = "Single-cell RNA-seq", species = "Human",
      disease_context = "Atrial Fibrillation", comparison = "AF vs SR", 
      value = "EXTREMELY HIGH", camk_relevance = "CELL-TYPE RESOLUTION",
      tier = "TIER 2"
    ),
    "GSE148506" = list(
      status = "INTEGRATED", platform = "SmartSeq2", samples = 384,
      data_type = "Single-cell RNA-seq", species = "Human",
      disease_context = "Atrial Fibrillation", comparison = "AF vs Healthy",
      value = "EXTREMELY HIGH", camk_relevance = "FIBROBLAST-SPECIFIC", 
      tier = "TIER 2"
    ),
    "GSE299292" = list(
      status = "INTEGRATED", platform = "GPL1261", samples = 16,
      data_type = "Microarray", species = "Mouse",
      disease_context = "Myocardial Infarction", comparison = "Control vs MI",
      value = "MODERATE-HIGH", camk_relevance = "TEMPORAL DYNAMICS",
      tier = "TIER 3"
    ),
    "GSE297444" = list(
      status = "INTEGRATED", platform = "GPL1355", samples = 26, 
      data_type = "Microarray", species = "Rat",
      disease_context = "Pressure Overload", comparison = "Sham vs TAC",
      value = "HIGH", camk_relevance = "CROSS-SPECIES VALIDATION",
      tier = "TIER 3"
    )
  )
  
  # Investigate each dataset
  for (dataset_id in dataset_ids) {
    cat("DATA: Investigating", dataset_id, "\n")
    cat(paste(rep("-", 40), collapse = ""), "\n")
    
    if (dataset_id %in% names(known_datasets)) {
      # Use known information for already integrated datasets
      known_info <- known_datasets[[dataset_id]]
      investigation_results$classification <- rbind(
        investigation_results$classification,
        data.frame(
          Dataset_ID = dataset_id,
          Status = known_info$status,
          Platform = known_info$platform,
          Samples = known_info$samples,
          Data_Type = known_info$data_type,
          Species = known_info$species,
          Disease_Context = known_info$disease_context,
          Comparison_Type = known_info$comparison,
          Scientific_Value = known_info$value,
          CAMK_Relevance = known_info$camk_relevance,
          Integration_Tier = known_info$tier,
          stringsAsFactors = FALSE
        )
      )
      
      cat("   SUCCESS: Status:", known_info$status, "\n")
      cat("   SUMMARY: Platform:", known_info$platform, "\n")
      cat("   DATA: Samples:", known_info$samples, "\n")
      cat("   GENETIC: Data Type:", known_info$data_type, "\n")
      cat("   TARGET: Tier:", known_info$tier, "\n\n")
      
    } else {
      # Investigate new datasets
      dataset_info <- investigate_individual_dataset(dataset_id)
      
      if (!is.null(dataset_info)) {
        investigation_results$classification <- rbind(
          investigation_results$classification,
          data.frame(
            Dataset_ID = dataset_id,
            Status = dataset_info$status,
            Platform = dataset_info$platform,
            Samples = dataset_info$samples,
            Data_Type = dataset_info$data_type,
            Species = dataset_info$species,
            Disease_Context = dataset_info$disease_context,
            Comparison_Type = dataset_info$comparison_type,
            Scientific_Value = dataset_info$scientific_value,
            CAMK_Relevance = dataset_info$camk_relevance,
            Integration_Tier = dataset_info$tier,
            stringsAsFactors = FALSE
          )
        )
        
        investigation_results$dataset_metadata[[dataset_id]] <- dataset_info
        
        cat("   SUMMARY: Status:", dataset_info$status, "\n")
        cat("   GENETIC: Platform:", dataset_info$platform, "\n")
        cat("   DATA: Samples:", dataset_info$samples, "\n")
        cat("   TARGET: Scientific Value:", dataset_info$scientific_value, "\n\n")
        
      } else {
        # Handle failed investigation
        investigation_results$classification <- rbind(
          investigation_results$classification,
          data.frame(
            Dataset_ID = dataset_id,
            Status = "INVESTIGATION_FAILED",
            Platform = "Unknown",
            Samples = 0,
            Data_Type = "Unknown",
            Species = "Unknown", 
            Disease_Context = "Unknown",
            Comparison_Type = "Unknown",
            Scientific_Value = "UNKNOWN",
            CAMK_Relevance = "UNKNOWN",
            Integration_Tier = "UNKNOWN",
            stringsAsFactors = FALSE
          )
        )
        
        cat("   ERROR: Investigation failed for", dataset_id, "\n\n")
      }
    }
  }
  
  # Generate investigation summary
  investigation_results$investigation_summary <- generate_investigation_summary(investigation_results$classification)
  
  # Save results
  saveRDS(investigation_results, file.path(output_dir, "comprehensive_dataset_investigation.rds"))
  write.csv(investigation_results$classification, 
            file.path(output_dir, "dataset_classification_table.csv"),
            row.names = FALSE)
  
  # Print summary
  print_investigation_summary(investigation_results$investigation_summary)
  
  return(investigation_results)
}

#' Individual Dataset Investigation Function
#'
#' @param dataset_id GEO dataset ID to investigate
#' @return List containing dataset metadata and classification
investigate_individual_dataset <- function(dataset_id) {
  
  dataset_info <- tryCatch({
    
    # Get GEO metadata
    cat("   SEARCH: Fetching GEO metadata...\n")
    gset <- getGEO(dataset_id, GSEMatrix = FALSE, getGPL = FALSE)
    
    # Extract basic information
    title <- Meta(gset)$title %||% "No title available"
    summary <- Meta(gset)$summary %||% "No summary available"  
    organism <- Meta(gset)$taxon %||% "Unknown"
    platform_id <- Meta(gset)$platform_id %||% "Unknown"
    submission_date <- Meta(gset)$submission_date %||% "Unknown"
    
    # Get sample information
    sample_names <- names(GSMList(gset))
    n_samples <- length(sample_names)
    
    # Determine data type based on platform and metadata
    data_type <- classify_data_type(platform_id, title, summary)
    
    # Classify disease context
    disease_context <- classify_disease_context(title, summary)
    
    # Assess CAMK relevance
    camk_relevance <- assess_camk_relevance(title, summary, disease_context)
    
    # Determine scientific value
    scientific_value <- assess_scientific_value(n_samples, data_type, disease_context, camk_relevance)
    
    # Assign tier
    tier <- assign_integration_tier(scientific_value, camk_relevance, data_type, n_samples)
    
    # Classify comparison type
    comparison_type <- classify_comparison_type(title, summary, sample_names)
    
    # Map species
    species <- map_organism_to_species(organism)
    
    list(
      status = "INVESTIGATED",
      platform = platform_id,
      samples = n_samples,
      data_type = data_type,
      species = species,
      disease_context = disease_context,
      comparison_type = comparison_type,
      scientific_value = scientific_value,
      camk_relevance = camk_relevance,
      tier = tier,
      title = title,
      summary = summary,
      submission_date = submission_date,
      sample_names = sample_names
    )
    
  }, error = function(e) {
    cat("   WARNING: Error investigating", dataset_id, ":", e$message, "\n")
    return(NULL)
  })
  
  return(dataset_info)
}

#' Classify Data Type based on Platform and Metadata
classify_data_type <- function(platform_id, title, summary) {
  
  # Single-cell indicators
  sc_indicators <- c("single.cell", "single cell", "scRNA", "10x", "drop.seq", 
                     "smartseq", "smart.seq", "cell.seq", "inDrop", "MARS-seq")
  
  # RNA-seq indicators  
  rnaseq_indicators <- c("RNA.seq", "RNA seq", "transcriptome", "next.generation", 
                         "illumina", "sequencing", "RNAseq", "rna-seq")
  
  # Combine title and summary for searching
  text_search <- paste(tolower(title), tolower(summary))
  
  # Check for single-cell
  if (any(sapply(sc_indicators, function(x) grepl(x, text_search, ignore.case = TRUE)))) {
    return("Single-cell RNA-seq")
  }
  
  # Check for RNA-seq
  if (any(sapply(rnaseq_indicators, function(x) grepl(x, text_search, ignore.case = TRUE)))) {
    return("RNA-seq")
  }
  
  # Check platform ID for microarray
  if (grepl("GPL", platform_id)) {
    return("Microarray")
  }
  
  return("Unknown")
}

#' Classify Disease Context
classify_disease_context <- function(title, summary) {
  
  text_search <- paste(tolower(title), tolower(summary))
  
  # Disease context mapping
  disease_contexts <- list(
    "Heart Failure" = c("heart failure", "cardiac failure", "cardiomyopathy", "dcm", "hcm"),
    "Atrial Fibrillation" = c("atrial fibrillation", "afib", "af", "atrial", "fibrillation"), 
    "Myocardial Infarction" = c("myocardial infarction", "heart attack", "mi", "ischemia", "ischemic"),
    "Coronary Artery Disease" = c("coronary", "atherosclerosis", "cad", "coronary artery"),
    "Cardiac Hypertrophy" = c("hypertrophy", "cardiac hypertrophy", "ventricular hypertrophy"),
    "Arrhythmia" = c("arrhythmia", "rhythm", "conduction", "electrical"),
    "Cardiovascular" = c("cardiovascular", "cardiac", "heart", "cardio")
  )
  
  for (disease in names(disease_contexts)) {
    if (any(sapply(disease_contexts[[disease]], function(x) grepl(x, text_search)))) {
      return(disease)
    }
  }
  
  return("Unknown")
}

#' Assess CAMK Relevance
assess_camk_relevance <- function(title, summary, disease_context) {
  
  text_search <- paste(tolower(title), tolower(summary))
  
  # Direct CAMK mentions
  camk_indicators <- c("camk", "calcium.calmodulin", "ca2+.calmodulin", "camkii", "camk2")
  
  if (any(sapply(camk_indicators, function(x) grepl(x, text_search)))) {
    return("DIRECT CAMK STUDY")
  }
  
  # High relevance contexts
  high_relevance_contexts <- c("Heart Failure", "Atrial Fibrillation", "Cardiac Hypertrophy", "Arrhythmia")
  
  if (disease_context %in% high_relevance_contexts) {
    return("HIGH RELEVANCE")
  }
  
  # Moderate relevance
  if (disease_context %in% c("Myocardial Infarction", "Coronary Artery Disease", "Cardiovascular")) {
    return("MODERATE RELEVANCE")
  }
  
  return("LOW RELEVANCE")
}

#' Assess Scientific Value
assess_scientific_value <- function(n_samples, data_type, disease_context, camk_relevance) {
  
  score <- 0
  
  # Sample size scoring
  if (n_samples >= 300) score <- score + 3
  else if (n_samples >= 100) score <- score + 2  
  else if (n_samples >= 50) score <- score + 1
  
  # Data type scoring
  if (data_type == "Single-cell RNA-seq") score <- score + 3
  else if (data_type == "RNA-seq") score <- score + 2
  else if (data_type == "Microarray") score <- score + 1
  
  # Disease context scoring
  if (disease_context %in% c("Heart Failure", "Atrial Fibrillation")) score <- score + 2
  else if (disease_context != "Unknown") score <- score + 1
  
  # CAMK relevance scoring
  if (camk_relevance == "DIRECT CAMK STUDY") score <- score + 3
  else if (camk_relevance == "HIGH RELEVANCE") score <- score + 2
  else if (camk_relevance == "MODERATE RELEVANCE") score <- score + 1
  
  # Convert score to category
  if (score >= 8) return("EXTREMELY HIGH")
  else if (score >= 6) return("HIGH")
  else if (score >= 4) return("MODERATE")
  else if (score >= 2) return("LOW")
  else return("VERY LOW")
}

#' Assign Integration Tier
assign_integration_tier <- function(scientific_value, camk_relevance, data_type, n_samples) {
  
  if (scientific_value == "EXTREMELY HIGH" && camk_relevance %in% c("DIRECT CAMK STUDY", "HIGH RELEVANCE")) {
    if (data_type == "Single-cell RNA-seq") return("TIER 2")
    else return("TIER 1")
  }
  
  if (scientific_value %in% c("HIGH", "EXTREMELY HIGH")) {
    if (data_type == "Single-cell RNA-seq") return("TIER 2")
    else return("TIER 1")
  }
  
  if (scientific_value == "MODERATE") return("TIER 3")
  
  return("TIER 4")
}

#' Additional helper functions
map_organism_to_species <- function(organism) {
  species_map <- list(
    "Homo sapiens" = "Human",
    "Mus musculus" = "Mouse", 
    "Rattus norvegicus" = "Rat"
  )
  
  return(species_map[[organism]] %||% organism %||% "Unknown")
}

classify_comparison_type <- function(title, summary, sample_names) {
  
  text_search <- paste(tolower(title), tolower(summary), paste(sample_names, collapse = " "))
  
  # Common comparison types
  if (grepl("healthy.*vs.*disease|control.*vs.*disease|normal.*vs.*disease", text_search)) {
    return("Healthy vs Disease")
  }
  
  if (grepl("before.*vs.*after|pre.*vs.*post|baseline.*vs.*treatment", text_search)) {
    return("Before vs After Treatment")
  }
  
  if (grepl("time.*course|temporal|longitudinal|0.*day.*1.*day", text_search)) {
    return("Time Course")
  }
  
  return("Unknown")
}

#' Generate Investigation Summary
generate_investigation_summary <- function(classification_df) {
  
  summary_stats <- list(
    total_datasets = nrow(classification_df),
    by_status = table(classification_df$Status),
    by_tier = table(classification_df$Integration_Tier),
    by_data_type = table(classification_df$Data_Type),
    by_species = table(classification_df$Species),
    by_scientific_value = table(classification_df$Scientific_Value),
    total_samples = sum(classification_df$Samples, na.rm = TRUE),
    high_value_datasets = sum(classification_df$Scientific_Value %in% c("HIGH", "EXTREMELY HIGH"), na.rm = TRUE)
  )
  
  return(summary_stats)
}

#' Print Investigation Summary
print_investigation_summary <- function(summary_stats) {
  
  cat("\nTARGET: COMPREHENSIVE DATASET INVESTIGATION SUMMARY\n")
  cat(paste(rep("=", 55), collapse = ""), "\n\n")
  
  cat("DATA: **Overall Statistics:**\n")
  cat("   • Total datasets investigated:", summary_stats$total_datasets, "\n")
  cat("   • Total samples across all datasets:", summary_stats$total_samples, "\n")
  cat("   • High-value datasets identified:", summary_stats$high_value_datasets, "\n\n")
  
  cat("TARGET: **Integration Tier Distribution:**\n")
  for (tier in names(summary_stats$by_tier)) {
    cat("   •", tier, ":", summary_stats$by_tier[tier], "datasets\n")
  }
  cat("\n")
  
  cat("GENETIC: **Data Type Distribution:**\n")
  for (dtype in names(summary_stats$by_data_type)) {
    cat("   •", dtype, ":", summary_stats$by_data_type[dtype], "datasets\n")
  }
  cat("\n")
  
  cat("METHOD: **Species Distribution:**\n")
  for (species in names(summary_stats$by_species)) {
    cat("   •", species, ":", summary_stats$by_species[species], "datasets\n")
  }
  cat("\n")
  
  cat("STAR: **Scientific Value Assessment:**\n")
  for (value in names(summary_stats$by_scientific_value)) {
    cat("   •", value, ":", summary_stats$by_scientific_value[value], "datasets\n")
  }
  cat("\n")
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

# Define the 15 datasets for investigation
target_datasets <- c(
  # New datasets to investigate
  "GSE181114", "GSE59876", "GSE57345", "GSE95143", "GSE670", 
  "GSE49937", "GSE216264", "GSE42579", "GSE197049", "GSE31619",
  "GSE29744",  # Note: Check if this differs from GSE299292
  
  # Already integrated datasets (for completeness)
  "GSE57338", "GSE148507", "GSE148506", "GSE299292"
)

cat("TARGET: **TARGET DATASETS FOR CAMK ANALYSIS:**\n")
for (i in seq_along(target_datasets)) {
  cat("  ", i, ".", target_datasets[i], "\n")
}
cat("\n")

# Run comprehensive investigation
cat("LAUNCH: Starting comprehensive dataset investigation...\n\n")

investigation_results <- investigate_datasets_comprehensive(
  dataset_ids = target_datasets,
  output_dir = "results/comprehensive_dataset_investigation"
)

cat("\nSUCCESS: Comprehensive dataset investigation completed!\n")
cat("SAVED: Results saved to: results/comprehensive_dataset_investigation/\n")
cat("DATA: Classification table: dataset_classification_table.csv\n")
cat("SAVED: Full results: comprehensive_dataset_investigation.rds\n\n")

cat("TARGET: Ready for multi-dataset CAMK analysis integration!\n")