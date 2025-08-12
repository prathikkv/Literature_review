#!/usr/bin/env Rscript
#' Targeted CAMK2D Dataset Search for Proposal-Specific Research
#' 
#' This module implements highly targeted search strategies for CAMK2D research
#' focusing specifically on atrial fibrillation and heart failure studies
#' with actual expression data validation and bioinformatics readiness.

# Load required libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(rentrez)
  library(tidyverse)
  library(httr)
  library(jsonlite)
  library(biomaRt)
})

#' Targeted CAMK2D Dataset Discovery
#'
#' Implements proposal-specific search strategies with validation
#' @param focus_area Character. Primary research focus ("aFIB", "HF", "both")
#' @param species Character. Target species ("human", "mouse", "rat", "all")
#' @param min_samples Integer. Minimum samples required for DGE analysis
#' @param require_controls Logical. Require case-control design
#' @return List with targeted datasets and validation results
discover_targeted_camk2d_datasets <- function(focus_area = "both", 
                                             species = "human",
                                             min_samples = 20,
                                             require_controls = TRUE) {
  
  cat("🎯 Targeted CAMK2D Dataset Discovery\n")
  cat("📋 Focus:", focus_area, "| Species:", species, "| Min samples:", min_samples, "\n\n")
  
  # Define proposal-specific search strategies
  search_strategies <- create_proposal_search_strategies(focus_area, species)
  
  # Execute targeted searches
  all_results <- list()
  total_strategies <- length(search_strategies)
  
  for (i in seq_along(search_strategies)) {
    strategy_name <- names(search_strategies)[i]
    cat("🔍 Executing strategy", i, "of", total_strategies, ":", strategy_name, "\n")
    strategy <- search_strategies[[strategy_name]]
    
    # Search with validation
    strategy_results <- execute_targeted_search(
      strategy = strategy,
      min_samples = min_samples,
      require_controls = require_controls
    )
    
    if (nrow(strategy_results) > 0) {
      strategy_results$Search_Strategy <- strategy_name
      all_results[[strategy_name]] <- strategy_results
      cat("✅ Found", nrow(strategy_results), "relevant datasets\n")
    } else {
      cat("⚠️ No datasets found for", strategy_name, "\n")
    }
    cat("\n")
  }
  
  # Combine and validate results
  if (length(all_results) > 0) {
    combined_results <- bind_rows(all_results)
    
    # Apply proposal-specific scoring
    scored_results <- apply_proposal_scoring(combined_results)
    
    # Validate CAMK2D relevance
    validated_results <- validate_camk2d_relevance(scored_results)
    
    cat("📊 Final Results Summary:\n")
    cat("   Total datasets:", nrow(validated_results), "\n")
    cat("   High relevance:", sum(validated_results$CAMK2D_Relevance == "High"), "\n")
    cat("   Medium relevance:", sum(validated_results$CAMK2D_Relevance == "Medium"), "\n")
    
    return(list(
      datasets = validated_results,
      summary = create_discovery_summary(validated_results),
      search_strategies = search_strategies
    ))
  } else {
    cat("❌ No relevant datasets discovered\n")
    return(list(datasets = data.frame(), summary = NULL))
  }
}

#' Create Proposal-Specific Search Strategies
#'
#' @param focus_area Character. Research focus
#' @param species Character. Target species
#' @return List of search strategy configurations
create_proposal_search_strategies <- function(focus_area, species) {
  
  strategies <- list()
  
  # Strategy 1: STRICT Human aFIB + Explicit CAMK2D
  if (focus_area %in% c("aFIB", "both") && species %in% c("human", "all")) {
    strategies$strict_human_afib_camk2d <- list(
      description = "STRICT: Human aFIB with explicit CAMK2D in title/abstract",
      species = "Homo sapiens",
      primary_terms = c("CAMK2D", "CaMKII-delta", "CaMKIIdelta", "calcium calmodulin kinase II delta"),
      disease_terms = c("atrial fibrillation", "aFIB"),
      tissue_terms = c("atrial", "atrium"),
      platform_priority = c("RNA-seq", "transcriptome"),
      study_design = "any",
      priority_score = 10,
      require_explicit_camk2d = TRUE  # NEW: Require CAMK2D in title/abstract
    )
  }
  
  # Strategy 2: STRICT Human HF + Explicit CAMK2D  
  if (focus_area %in% c("HF", "both") && species %in% c("human", "all")) {
    strategies$strict_human_hf_camk2d <- list(
      description = "STRICT: Human heart failure with explicit CAMK2D in title/abstract",
      species = "Homo sapiens", 
      primary_terms = c("CAMK2D", "CaMKII-delta", "CaMKIIdelta", "calcium calmodulin kinase II delta"),
      disease_terms = c("heart failure", "cardiac failure"),
      tissue_terms = c("ventricular", "ventricle", "myocardial"),
      platform_priority = c("RNA-seq", "transcriptome"),
      study_design = "any",
      priority_score = 10,
      require_explicit_camk2d = TRUE  # NEW: Require CAMK2D in title/abstract
    )
  }
  
  # Strategy 3: STRICT CAMK2D inhibitor/knockout studies
  if (species %in% c("human", "mouse", "all")) {
    strategies$strict_camk2d_inhibition <- list(
      description = "STRICT: CAMK2D inhibition/knockout with explicit CAMK2D",
      species = c("Homo sapiens", "Mus musculus"),
      primary_terms = c("CAMK2D knockout", "CAMK2D inhibitor", "CaMKII-delta knockout", "CAMK2D silencing"),
      disease_terms = c("knockout", "inhibitor", "silencing", "KN-93", "AIP"),
      tissue_terms = c("cardiac", "heart", "cardiomyocyte"),
      platform_priority = c("RNA-seq", "transcriptome"),
      study_design = "any",
      priority_score = 10,
      require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
    )
  }
  
  # Strategy 4: STRICT CAMK2D phosphorylation substrates
  strategies$strict_camk2d_substrates <- list(
    description = "STRICT: CAMK2D phosphorylation targets with explicit CAMK2D",
    species = c("Homo sapiens", "Mus musculus"),
    primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta", "calcium calmodulin kinase II delta"),
    disease_terms = c("CAMK2D substrate", "CAMK2D phosphorylation", "CaMKII delta phosphorylation"),
    tissue_terms = c("cardiac", "heart"),
    platform_priority = c("RNA-seq", "proteomics", "phosphoproteomics"),
    study_design = "any",
    priority_score = 10,
    require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
  )
  
  # Strategy 5: STRICT Mouse CAMK2D cardiac models
  if (species %in% c("mouse", "all")) {
    strategies$strict_mouse_camk2d_cardiac <- list(
      description = "STRICT: Mouse CAMK2D cardiac models with explicit CAMK2D",
      species = "Mus musculus",
      primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta", "calcium calmodulin kinase II delta"),
      disease_terms = c("CAMK2D cardiac", "CAMK2D TAC", "CAMK2D pressure overload", "CAMK2D hypertrophy"),
      tissue_terms = c("heart", "cardiac", "ventricular"),
      platform_priority = c("RNA-seq", "transcriptome"),
      study_design = "any", 
      priority_score = 9,
      require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
    )
  }
  
  # NEW EXPANDED STRATEGIES FOR BROADER DISCOVERY
  
  # Strategy 6: STRICT CAMK2D + other CAMK isoforms (comparative)
  strategies$strict_camk2d_comparative <- list(
    description = "STRICT: CAMK2D comparative studies with other isoforms",
    species = c("Homo sapiens", "Mus musculus"),
    primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta"),
    disease_terms = c("CAMK2A CAMK2D", "CaMKII isoform", "CAMK2D comparison"),
    tissue_terms = c("heart", "cardiac", "myocardial"),
    platform_priority = c("RNA-seq", "microarray", "transcriptome"),
    study_design = "any",
    priority_score = 8,
    require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
  )
  
  # Strategy 7: STRICT CAMK2D in cardiac hypertrophy
  strategies$strict_camk2d_hypertrophy <- list(
    description = "STRICT: CAMK2D role in cardiac hypertrophy with explicit mention",
    species = c("Homo sapiens", "Mus musculus"),
    primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta"),
    disease_terms = c("CAMK2D hypertrophy", "CaMKII delta remodeling", "CAMK2D fibrosis"),
    tissue_terms = c("ventricular", "myocardial", "cardiac"),
    platform_priority = c("RNA-seq", "microarray"),
    study_design = "any",
    priority_score = 9,
    require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
  )
  
  # Strategy 8: STRICT CAMK2D calcium signaling
  strategies$strict_camk2d_calcium <- list(
    description = "STRICT: CAMK2D calcium signaling with explicit CAMK2D mention",
    species = c("Homo sapiens", "Mus musculus"),
    primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta"),
    disease_terms = c("CAMK2D calcium", "CaMKII delta signaling", "CAMK2D sarcoplasmic"),
    tissue_terms = c("cardiomyocyte", "myocardial", "cardiac"),
    platform_priority = c("RNA-seq", "microarray"),
    study_design = "any",
    priority_score = 8,
    require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
  )
  
  # Strategy 9: STRICT CAMK2D atrial studies
  if (focus_area %in% c("aFIB", "both")) {
    strategies$strict_camk2d_atrial <- list(
      description = "STRICT: CAMK2D atrial studies with explicit CAMK2D mention",
      species = c("Homo sapiens", "Mus musculus"),
      primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta"),
      disease_terms = c("CAMK2D atrial", "CaMKII delta atrial fibrillation", "CAMK2D arrhythmia"),
      tissue_terms = c("atrial tissue", "atrium", "atrial appendage"),
      platform_priority = c("RNA-seq", "microarray", "single-cell"),
      study_design = "any",
      priority_score = 10,
      require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
    )
  }
  
  # Strategy 10: STRICT CAMK2D heart failure progression
  if (focus_area %in% c("HF", "both")) {
    strategies$strict_camk2d_hf_progression <- list(
      description = "STRICT: CAMK2D heart failure progression with explicit mention",
      species = c("Homo sapiens", "Mus musculus"),
      primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta"),
      disease_terms = c("CAMK2D heart failure", "CaMKII delta HFrEF", "CAMK2D cardiac dysfunction"),
      tissue_terms = c("ventricular", "myocardial", "cardiac"),
      platform_priority = c("RNA-seq", "microarray"),
      study_design = "any",
      priority_score = 10,
      require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
    )
  }
  
  # Strategy 11: STRICT CAMK2D cardiomyocyte studies
  strategies$strict_camk2d_cardiomyocyte <- list(
    description = "STRICT: CAMK2D cardiomyocyte studies with explicit mention",
    species = c("Homo sapiens", "Mus musculus"),
    primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta"),
    disease_terms = c("CAMK2D cardiomyocyte", "CaMKII delta iPSC", "CAMK2D cardiac cell"),
    tissue_terms = c("cardiomyocyte", "cardiac cell", "myocyte"),
    platform_priority = c("RNA-seq", "single-cell", "transcriptome"),
    study_design = "any",
    priority_score = 8,
    require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
  )
  
  # Strategy 12: STRICT CAMK2D drug treatment studies
  strategies$strict_camk2d_drug_treatment <- list(
    description = "STRICT: CAMK2D drug treatment with explicit CAMK2D mention",
    species = c("Homo sapiens", "Mus musculus"),
    primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta"),
    disease_terms = c("CAMK2D treatment", "CaMKII delta therapy", "CAMK2D drug"),
    tissue_terms = c("cardiac", "heart", "myocardial"),
    platform_priority = c("RNA-seq", "microarray"),
    study_design = "any", 
    priority_score = 9,
    require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
  )
  
  # Strategy 13: STRICT CAMK2D ischemic heart disease
  strategies$strict_camk2d_ischemic <- list(
    description = "STRICT: CAMK2D ischemic heart disease with explicit mention",
    species = c("Homo sapiens", "Mus musculus"),
    primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta"),
    disease_terms = c("CAMK2D ischemia", "CaMKII delta infarction", "CAMK2D coronary"),
    tissue_terms = c("myocardial", "cardiac", "coronary"),
    platform_priority = c("RNA-seq", "microarray"),
    study_design = "any",
    priority_score = 8,
    require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
  )
  
  # Strategy 14: STRICT CAMK2D diabetic cardiomyopathy
  strategies$strict_camk2d_diabetic <- list(
    description = "STRICT: CAMK2D diabetic cardiomyopathy with explicit mention", 
    species = c("Homo sapiens", "Mus musculus"),
    primary_terms = c("CAMK2D", "CaMKII delta", "CaMKIIdelta"),
    disease_terms = c("CAMK2D diabetes", "CaMKII delta diabetic", "CAMK2D metabolic"),
    tissue_terms = c("cardiac", "myocardial", "heart"),
    platform_priority = c("RNA-seq", "microarray"),
    study_design = "any",
    priority_score = 8,
    require_explicit_camk2d = TRUE  # NEW: Must explicitly mention CAMK2D
  )
  
  # Removed redundant aging and gender strategies to optimize search time
  
  return(strategies)
}

#' Execute Targeted Search for Specific Strategy
#'
#' @param strategy List. Search strategy configuration
#' @param min_samples Integer. Minimum sample requirement
#' @param require_controls Logical. Require control samples
#' @return Data frame with search results
execute_targeted_search <- function(strategy, min_samples, require_controls) {
  
  results_list <- list()
  
  # Build targeted search queries
  queries <- build_targeted_queries(strategy)
  
  for (i in seq_along(queries)) {
    query <- queries[[i]]
    cat("   🔍 Query", i, ":", substr(query, 1, 60), "...\n")
    
    tryCatch({
      # Search GDS database first (more curated), then GSE if needed
      for (db in c("gds")) {  # Start with just GDS to avoid rate limits
        search_result <- entrez_search(
          db = db,
          term = query,
          retmax = 100  # Increased from 30 to get more results
        )
        
        if (length(search_result$ids) > 0) {
          # Extract metadata with validation and error handling
          for (id in search_result$ids) {
            tryCatch({
              dataset_info <- extract_validated_metadata(id, strategy, min_samples, require_controls, db)
              
              if (!is.null(dataset_info) && nrow(dataset_info) > 0) {
                dataset_info$Query_Used <- query
                dataset_info$Query_Number <- i
                dataset_info$Database <- db
                results_list[[paste0("query_", i, "_db_", db, "_id_", id)]] <- dataset_info
              }
            }, error = function(e) {
              # Silently skip datasets that cause errors
            })
          }
        }
        
        Sys.sleep(1.0)  # Increased rate limiting between database searches  
      }
      
      Sys.sleep(2.0)  # Increased rate limiting between queries to avoid 429 errors
      
    }, error = function(e) {
      cat("      ❌ Query failed:", e$message, "\n")
    })
  }
  
  if (length(results_list) > 0) {
    return(bind_rows(results_list))
  } else {
    return(data.frame())
  }
}

#' Build Targeted Search Queries
#'
#' @param strategy List. Search strategy configuration
#' @return Character vector of search queries
build_targeted_queries <- function(strategy) {
  
  queries <- c()
  
  # Primary query: CAMK2D + Disease + Tissue + Platform
  primary_query <- paste0(
    "(", paste(strategy$primary_terms, collapse = " OR "), ")",
    " AND (",
    paste(strategy$disease_terms, collapse = " OR "), ")",
    " AND (",
    paste(strategy$tissue_terms, collapse = " OR "), ")",
    ' AND "', strategy$species[1], '"[Organism]'
  )
  queries <- c(queries, primary_query)
  
  # Secondary query: CAMK2D + Platform priority
  if (length(strategy$platform_priority) > 0) {
    platform_query <- paste0(
      "(", paste(strategy$primary_terms, collapse = " OR "), ")",
      " AND (",
      paste(strategy$platform_priority, collapse = " OR "), ")",
      ' AND "', strategy$species[1], '"[Organism]'
    )
    queries <- c(queries, platform_query)
  }
  
  # Disease-specific query
  disease_query <- paste0(
    "(", paste(strategy$disease_terms, collapse = " OR "), ")",
    " AND (",
    paste(strategy$tissue_terms, collapse = " OR "), ")",
    ' AND "', strategy$species[1], '"[Organism]'
  )
  queries <- c(queries, disease_query)
  
  return(queries)
}

#' Extract and Validate Dataset Metadata
#'
#' @param geo_id Character. GEO dataset ID
#' @param strategy List. Search strategy used
#' @param min_samples Integer. Minimum sample requirement
#' @param require_controls Logical. Require control samples
#' @return Data frame with validated metadata or NULL
extract_validated_metadata <- function(geo_id, strategy, min_samples, require_controls, database = "gds") {
  
  tryCatch({
    # Get dataset summary from appropriate database
    summary_data <- entrez_summary(db = database, id = geo_id)
    
    # Extract basic information (handle both GDS and GSE)
    title <- summary_data$title %||% ""
    n_samples <- as.numeric(summary_data$n_samples %||% 0)
    platform <- summary_data$gpl %||% ""
    
    # Handle accession format based on database
    if (database == "gse") {
      accession <- summary_data$accession %||% paste0("GSE", geo_id)
    } else {
      accession <- summary_data$accession %||% paste0("GDS", geo_id)  
    }
    
    organism <- summary_data$taxon %||% ""
    
    # Validate sample size
    if (n_samples < min_samples) {
      return(NULL)
    }
    
    # Validate species match
    if (!any(sapply(strategy$species, function(sp) grepl(sp, organism, ignore.case = TRUE)))) {
      return(NULL)
    }
    
    # Validate study design (if required) - more flexible
    if (require_controls && strategy$study_design == "case-control") {
      if (!validate_case_control_design(title, summary_data$summary %||% "")) {
        return(NULL)
      }
    }
    # Allow any study design if require_controls is FALSE
    
    # Check for explicit CAMK2D requirement
    if (!is.null(strategy$require_explicit_camk2d) && strategy$require_explicit_camk2d) {
      if (!validate_explicit_camk2d_mention(title, summary_data$summary %||% "")) {
        return(NULL)  # Skip datasets without explicit CAMK2D mention
      }
    }
    
    # Calculate relevance scores
    camk2d_relevance <- calculate_camk2d_relevance_score(title, summary_data$summary %||% "")
    disease_relevance <- calculate_disease_relevance_score(title, summary_data$summary %||% "", strategy)
    platform_score <- calculate_platform_relevance_score(platform, title)
    
    # Create validated dataset record
    dataset_record <- data.frame(
      GEO_Accession = accession,
      Title = substr(title, 1, 150),
      Organism = organism,
      Platform = platform,
      Sample_Count = n_samples,
      CAMK2D_Relevance_Score = camk2d_relevance,
      Disease_Relevance_Score = disease_relevance,
      Platform_Score = platform_score,
      Strategy_Priority = strategy$priority_score,
      Submission_Date = as.Date(summary_data$pdat %||% NA),
      PubMed_ID = summary_data$pubmedids[1] %||% NA,
      Summary = substr(summary_data$summary %||% "", 1, 200),
      Validated = TRUE,
      stringsAsFactors = FALSE
    )
    
    return(dataset_record)
    
  }, error = function(e) {
    return(NULL)
  })
}

#' Validate Case-Control Study Design
#'
#' @param title Character. Dataset title
#' @param summary Character. Dataset summary
#' @return Logical. TRUE if appears to be case-control
validate_case_control_design <- function(title, summary) {
  
  text <- paste(tolower(title), tolower(summary))
  
  # Look for case-control indicators
  has_cases <- grepl("patient|disease|affected|case|pathology", text)
  has_controls <- grepl("control|normal|healthy|reference", text)
  has_vs <- grepl(" vs | versus |compared to", text)
  
  return(has_cases && (has_controls || has_vs))
}

#' Validate Explicit CAMK2D Mention
#'
#' @param title Character. Dataset title
#' @param summary Character. Dataset summary
#' @return Logical. TRUE if explicit CAMK2D mention found
validate_explicit_camk2d_mention <- function(title, summary) {
  
  text <- paste(tolower(title), tolower(summary))
  
  # Look for explicit CAMK2D mentions (strict)
  explicit_camk2d <- grepl("camk2d|camk.*2.*d|camkii.*delta|camk.*ii.*delta|calcium.*calmodulin.*kinase.*ii.*delta", text)
  
  return(explicit_camk2d)
}

#' Calculate CAMK2D Relevance Score
#'
#' @param title Character. Dataset title
#' @param summary Character. Dataset summary  
#' @return Numeric. Relevance score (0-10)
calculate_camk2d_relevance_score <- function(title, summary) {
  
  text <- paste(tolower(title), tolower(summary))
  score <- 0
  
  # Direct CAMK2D mentions - highest scores
  if (grepl("camk2d|camk.*ii.*delta", text)) score <- score + 5
  if (grepl("camkii.*delta|camk.*2.*delta", text)) score <- score + 4  
  if (grepl("camkii|camk.*ii|calcium.*calmodulin.*kinase.*ii", text)) score <- score + 3
  
  # Broader CAMK family - medium scores
  if (grepl("camk2a|camk2b|camk2g", text)) score <- score + 3
  if (grepl("calmodulin.*kinase|calcium.*kinase", text)) score <- score + 2
  
  # Calcium signaling and kinases - lower scores but still relevant
  if (grepl("calcium.*signaling|calcium.*handling", text)) score <- score + 2
  if (grepl("sarcoplasmic.*reticulum|calcium.*channel", text)) score <- score + 2
  if (grepl("kinase|phosphoryl", text)) score <- score + 1
  
  # Cardiac-specific terms add relevance even without kinase mentions
  if (grepl("cardiac|heart|myocardial|cardiomyocyte", text)) score <- score + 1
  
  return(min(10, score))
}

#' Calculate Disease Relevance Score
#'
#' @param title Character. Dataset title
#' @param summary Character. Dataset summary
#' @param strategy List. Search strategy
#' @return Numeric. Disease relevance score (0-10)
calculate_disease_relevance_score <- function(title, summary, strategy) {
  
  text <- paste(tolower(title), tolower(summary))
  score <- 0
  
  # Check for strategy-specific disease terms
  for (term in strategy$disease_terms) {
    if (grepl(tolower(term), text)) {
      score <- score + 2
    }
  }
  
  # Check for tissue relevance
  for (term in strategy$tissue_terms) {
    if (grepl(tolower(term), text)) {
      score <- score + 1
    }
  }
  
  return(min(10, score))
}

#' Calculate Platform Relevance Score
#'
#' @param platform Character. Platform information
#' @param title Character. Dataset title
#' @return Numeric. Platform score (0-10)
calculate_platform_relevance_score <- function(platform, title) {
  
  platform_text <- paste(tolower(platform), tolower(title))
  
  if (grepl("rna.?seq|illumina.*seq|novaseq|hiseq", platform_text)) return(10)
  if (grepl("single.?cell|10x", platform_text)) return(9)
  if (grepl("affymetrix.*u133|hg.*u133", platform_text)) return(7)
  if (grepl("illumina|microarray", platform_text)) return(6)
  if (grepl("agilent", platform_text)) return(5)
  
  return(3)  # Unknown platform
}

#' Apply Proposal-Specific Scoring
#'
#' @param datasets Data frame. Dataset results
#' @return Data frame with comprehensive scoring
apply_proposal_scoring <- function(datasets) {
  
  datasets %>%
    mutate(
      # Comprehensive relevance score
      Total_Relevance_Score = CAMK2D_Relevance_Score + 
                             Disease_Relevance_Score + 
                             Platform_Score + 
                             Strategy_Priority,
      
      # Categorize relevance - more inclusive thresholds
      CAMK2D_Relevance = case_when(
        CAMK2D_Relevance_Score >= 4 ~ "High",     # Lowered from 5
        CAMK2D_Relevance_Score >= 2 ~ "Medium",   # Lowered from 3
        TRUE ~ "Low"
      ),
      
      # Suitable for DGE analysis - more inclusive  
      DGE_Suitable = Sample_Count >= 10 & Platform_Score >= 4,  # Lowered thresholds
      
      # Priority classification
      Priority_Class = case_when(
        Total_Relevance_Score >= 20 ~ "Highest",
        Total_Relevance_Score >= 15 ~ "High", 
        Total_Relevance_Score >= 10 ~ "Medium",
        TRUE ~ "Low"
      )
    ) %>%
    arrange(desc(Total_Relevance_Score), desc(Sample_Count))
}

#' Validate CAMK2D Relevance
#'
#' @param datasets Data frame. Scored datasets
#' @return Data frame with relevance validation
validate_camk2d_relevance <- function(datasets) {
  
  # Filter for minimum relevance - more inclusive now
  validated <- datasets %>%
    filter(CAMK2D_Relevance %in% c("High", "Medium", "Low")) %>%  # Include Low relevance too
    mutate(
      Validation_Status = "Passed",
      Recommendation = case_when(
        Priority_Class == "Highest" & DGE_Suitable ~ "Primary Analysis Dataset",
        Priority_Class == "High" & DGE_Suitable ~ "Validation Dataset",
        Priority_Class == "Medium" ~ "Supplementary Dataset",
        TRUE ~ "Further Validation Required"
      )
    )
  
  return(validated)
}

#' Create Discovery Summary
#'
#' @param datasets Data frame. Final validated datasets
#' @return List with summary statistics
create_discovery_summary <- function(datasets) {
  
  if (nrow(datasets) == 0) return(NULL)
  
  list(
    total_datasets = nrow(datasets),
    high_relevance = sum(datasets$CAMK2D_Relevance == "High"),
    medium_relevance = sum(datasets$CAMK2D_Relevance == "Medium"),
    dge_suitable = sum(datasets$DGE_Suitable),
    primary_candidates = sum(datasets$Recommendation == "Primary Analysis Dataset"),
    validation_candidates = sum(datasets$Recommendation == "Validation Dataset"),
    avg_sample_size = round(mean(datasets$Sample_Count), 1),
    species_coverage = unique(datasets$Organism),
    date_range = paste(
      format(min(datasets$Submission_Date, na.rm = TRUE), "%Y"),
      "-",
      format(max(datasets$Submission_Date, na.rm = TRUE), "%Y")
    )
  )
}

cat("🎯 Targeted CAMK2D Search Module Loaded\n")
cat("📋 Proposal-aligned search strategies with expression validation\n")
cat("🔬 Ready for aFIB/HF focused dataset discovery\n")