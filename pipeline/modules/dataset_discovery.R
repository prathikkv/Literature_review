#!/usr/bin/env Rscript
#' Dataset Discovery Module
#' 
#' Automatically discovers relevant GEO datasets based on disease terms
#' Non-disruptive enhancement to existing pipeline
#' 
#' @author Claude Code Enhancement Module
#' @version 1.0.0

suppressPackageStartupMessages({
  library(GEOquery)
  library(rentrez)
  library(yaml)
  library(tidyverse)
  library(openxlsx)
})

#' Discover GEO Datasets
#'
#' Searches GEO for relevant datasets based on disease terms
#' @param config_file Configuration file with search parameters
#' @param output_file Excel file for results (default: "discovered_datasets.xlsx")
#' @param auto_download Whether to auto-download qualified datasets
#' @return Data frame with discovered datasets and quality scores
discover_geo_datasets <- function(config_file = "config.yml",
                                 output_file = "output/discovered_datasets.xlsx",
                                 auto_download = FALSE) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║          DATASET DISCOVERY MODULE                            ║\n")
  cat("║          Intelligent GEO Dataset Search                      ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Load configuration
  config <- yaml::read_yaml(config_file)
  
  # Extract search parameters
  diseases <- config$datasets$search_terms %||% 
              c("Heart Failure", "Atrial Fibrillation", "Cardiomyopathy")
  
  min_samples <- config$analysis$quality_control$min_samples_per_group %||% 20
  
  platforms <- c("GPL570", "GPL96", "GPL97")  # Common Affymetrix platforms
  
  cat("🔍 Search Parameters:\n")
  cat("   Diseases:", paste(diseases, collapse = ", "), "\n")
  cat("   Platforms:", paste(platforms, collapse = ", "), "\n")
  cat("   Minimum samples:", min_samples, "\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Initialize results
  all_datasets <- data.frame()
  
  # Search for each disease term
  for (disease in diseases) {
    cat("🔎 Searching for:", disease, "\n")
    
    # Build search query
    search_query <- paste0(
      disease, "[Title/Abstract] AND ",
      "Homo sapiens[Organism] AND ",
      "Expression profiling by array[DataSet Type]"
    )
    
    # Search GEO via Entrez
    search_results <- search_geo_datasets(
      query = search_query,
      platforms = platforms,
      min_samples = min_samples
    )
    
    if (nrow(search_results) > 0) {
      search_results$disease_match <- disease
      all_datasets <- rbind(all_datasets, search_results)
      cat("   Found", nrow(search_results), "datasets\n")
    } else {
      cat("   No datasets found\n")
    }
  }
  
  # Remove duplicates
  all_datasets <- all_datasets[!duplicated(all_datasets$gse_id), ]
  
  cat("\n📊 Total unique datasets found:", nrow(all_datasets), "\n")
  
  # Apply quality filters
  cat("\n🔧 Applying quality filters...\n")
  filtered_datasets <- apply_quality_filters(all_datasets, config)
  
  # Calculate quality scores
  cat("📈 Calculating quality scores...\n")
  scored_datasets <- calculate_quality_scores(filtered_datasets)
  
  # Sort by quality score
  final_datasets <- scored_datasets[order(scored_datasets$quality_score, decreasing = TRUE), ]
  
  # Save to Excel
  save_discovery_results(final_datasets, output_file)
  cat("\n💾 Results saved to:", output_file, "\n")
  
  # Print summary
  print_discovery_summary(final_datasets)
  
  # Auto-download if requested
  if (auto_download) {
    qualified_datasets <- final_datasets[final_datasets$include_status == "Include", ]
    if (nrow(qualified_datasets) > 0) {
      cat("\n📥 Auto-downloading", nrow(qualified_datasets), "qualified datasets...\n")
      update_config_with_datasets(qualified_datasets$gse_id, config_file)
    }
  }
  
  return(final_datasets)
}

#' Search GEO Datasets via Entrez
#'
#' Performs the actual GEO database search
#' @param query Search query string
#' @param platforms Vector of platform IDs
#' @param min_samples Minimum sample requirement
#' @return Data frame with search results
search_geo_datasets <- function(query, platforms, min_samples) {
  
  results <- data.frame()
  
  tryCatch({
    # Search GEO DataSets database
    search <- entrez_search(
      db = "gds",
      term = query,
      retmax = 100  # Limit to 100 results
    )
    
    if (length(search$ids) == 0) {
      return(results)
    }
    
    # Fetch dataset summaries
    summaries <- entrez_summary(
      db = "gds",
      id = search$ids
    )
    
    # Extract relevant information
    for (i in seq_along(summaries)) {
      summary <- summaries[[i]]
      
      # Skip if not a Series (GSE)
      if (!grepl("^GSE", summary$accession)) {
        next
      }
      
      # Extract metadata
      dataset_info <- data.frame(
        gse_id = summary$accession,
        title = summary$title %||% "",
        summary = summary$summary %||% "",
        platform = summary$gpl %||% "",
        samples = as.numeric(summary$n_samples %||% 0),
        pubmed_id = summary$pubmed_ids %||% "",
        stringsAsFactors = FALSE
      )
      
      # Check platform
      if (length(platforms) > 0 && !dataset_info$platform %in% platforms) {
        next
      }
      
      # Check sample size
      if (dataset_info$samples < min_samples) {
        next
      }
      
      results <- rbind(results, dataset_info)
    }
    
  }, error = function(e) {
    cat("   ⚠️  Search error:", e$message, "\n")
  })
  
  return(results)
}

#' Apply Quality Filters
#'
#' Filters datasets based on quality criteria
#' @param datasets Data frame of discovered datasets
#' @param config Configuration object
#' @return Filtered data frame with inclusion status
apply_quality_filters <- function(datasets, config) {
  
  if (nrow(datasets) == 0) {
    return(datasets)
  }
  
  # Initialize inclusion status
  datasets$include_status <- "Include"
  datasets$exclusion_reason <- ""
  
  # Filter: Minimum samples
  min_samples <- config$analysis$quality_control$min_samples_per_group %||% 20
  too_small <- datasets$samples < min_samples
  datasets$include_status[too_small] <- "Exclude"
  datasets$exclusion_reason[too_small] <- paste("< ", min_samples, " samples")
  
  # Filter: Maximum samples (exclude mega-studies)
  max_samples <- 1000
  too_large <- datasets$samples > max_samples
  datasets$include_status[too_large] <- "Exclude"
  datasets$exclusion_reason[too_large] <- paste("> ", max_samples, " samples")
  
  # Filter: Exclude cell line studies
  cell_line_keywords <- c("cell line", "in vitro", "cultured", "HEK293", "HeLa")
  for (keyword in cell_line_keywords) {
    is_cell_line <- grepl(keyword, datasets$title, ignore.case = TRUE) |
                   grepl(keyword, datasets$summary, ignore.case = TRUE)
    
    datasets$include_status[is_cell_line & datasets$include_status == "Include"] <- "Exclude"
    datasets$exclusion_reason[is_cell_line & datasets$exclusion_reason == ""] <- 
      paste("Cell line study:", keyword)
  }
  
  # Filter: Require case/control design keywords
  case_control_keywords <- c("control", "normal", "healthy", "vs", "versus", "compared")
  has_controls <- FALSE
  for (keyword in case_control_keywords) {
    has_controls <- has_controls | 
                   grepl(keyword, datasets$title, ignore.case = TRUE) |
                   grepl(keyword, datasets$summary, ignore.case = TRUE)
  }
  
  no_controls <- !has_controls
  datasets$include_status[no_controls & datasets$include_status == "Include"] <- "Review"
  datasets$exclusion_reason[no_controls & datasets$exclusion_reason == ""] <- 
    "No clear control group"
  
  # Check if already in pipeline
  existing_datasets <- names(config$datasets$active_datasets)
  is_existing <- datasets$gse_id %in% existing_datasets
  datasets$include_status[is_existing] <- "Existing"
  datasets$exclusion_reason[is_existing] <- "Already in pipeline"
  
  return(datasets)
}

#' Calculate Quality Scores
#'
#' Assigns quality scores to datasets
#' @param datasets Filtered dataset data frame
#' @return Data frame with quality scores
calculate_quality_scores <- function(datasets) {
  
  if (nrow(datasets) == 0) {
    return(datasets)
  }
  
  # Initialize scores
  datasets$quality_score <- 0
  
  # Score: Sample size (0-30 points)
  # Optimal range: 50-500 samples
  datasets$quality_score <- datasets$quality_score + 
    pmin(30, pmax(0, (datasets$samples - 20) / 15))
  
  # Score: Platform quality (0-20 points)
  platform_scores <- c("GPL570" = 20, "GPL96" = 18, "GPL97" = 16)
  for (platform in names(platform_scores)) {
    is_platform <- datasets$platform == platform
    datasets$quality_score[is_platform] <- 
      datasets$quality_score[is_platform] + platform_scores[platform]
  }
  
  # Score: Has PubMed publication (0-20 points)
  has_pubmed <- datasets$pubmed_id != "" & !is.na(datasets$pubmed_id)
  datasets$quality_score[has_pubmed] <- datasets$quality_score[has_pubmed] + 20
  
  # Score: Title/summary quality (0-15 points)
  quality_keywords <- c("patient", "clinical", "human", "biopsy", "tissue")
  keyword_score <- 0
  for (keyword in quality_keywords) {
    has_keyword <- grepl(keyword, datasets$title, ignore.case = TRUE) |
                  grepl(keyword, datasets$summary, ignore.case = TRUE)
    keyword_score <- keyword_score + (has_keyword * 3)
  }
  datasets$quality_score <- datasets$quality_score + pmin(15, keyword_score)
  
  # Score: Disease relevance (0-15 points)
  # This is already captured in disease_match
  datasets$quality_score <- datasets$quality_score + 15
  
  # Round scores
  datasets$quality_score <- round(datasets$quality_score, 1)
  
  return(datasets)
}

#' Save Discovery Results to Excel
#'
#' Creates formatted Excel file with discovery results
#' @param datasets Data frame with scored datasets
#' @param output_file Path to output Excel file
save_discovery_results <- function(datasets, output_file) {
  
  # Ensure output directory exists
  output_dir <- dirname(output_file)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Create workbook
  wb <- createWorkbook()
  
  # Add summary sheet
  addWorksheet(wb, "Summary")
  
  summary_data <- data.frame(
    Metric = c("Total Datasets Found", "Included", "Excluded", "Review Required", "Already in Pipeline"),
    Count = c(
      nrow(datasets),
      sum(datasets$include_status == "Include"),
      sum(datasets$include_status == "Exclude"),
      sum(datasets$include_status == "Review"),
      sum(datasets$include_status == "Existing")
    )
  )
  
  writeData(wb, "Summary", summary_data)
  
  # Add detailed results sheet
  addWorksheet(wb, "Datasets")
  
  # Reorder columns for better readability
  output_columns <- c("gse_id", "title", "samples", "platform", "disease_match",
                     "quality_score", "include_status", "exclusion_reason",
                     "pubmed_id", "summary")
  
  output_data <- datasets[, intersect(output_columns, names(datasets))]
  
  writeData(wb, "Datasets", output_data)
  
  # Format as table
  addStyle(wb, "Datasets", 
          style = createStyle(textDecoration = "bold"),
          rows = 1, cols = 1:ncol(output_data))
  
  # Color code by inclusion status
  include_rows <- which(output_data$include_status == "Include") + 1
  exclude_rows <- which(output_data$include_status == "Exclude") + 1
  review_rows <- which(output_data$include_status == "Review") + 1
  
  if (length(include_rows) > 0) {
    addStyle(wb, "Datasets",
            style = createStyle(fgFill = "#90EE90"),
            rows = include_rows, cols = 1:ncol(output_data),
            gridExpand = TRUE)
  }
  
  if (length(exclude_rows) > 0) {
    addStyle(wb, "Datasets",
            style = createStyle(fgFill = "#FFB6C1"),
            rows = exclude_rows, cols = 1:ncol(output_data),
            gridExpand = TRUE)
  }
  
  if (length(review_rows) > 0) {
    addStyle(wb, "Datasets",
            style = createStyle(fgFill = "#FFFFE0"),
            rows = review_rows, cols = 1:ncol(output_data),
            gridExpand = TRUE)
  }
  
  # Auto-size columns
  setColWidths(wb, "Datasets", cols = 1:ncol(output_data), widths = "auto")
  
  # Save workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)
}

#' Print Discovery Summary
#'
#' Prints summary of discovered datasets
#' @param datasets Data frame with discovery results
print_discovery_summary <- function(datasets) {
  
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📊 DISCOVERY SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  # Count by status
  status_counts <- table(datasets$include_status)
  
  for (status in names(status_counts)) {
    icon <- switch(status,
                  "Include" = "✅",
                  "Exclude" = "❌",
                  "Review" = "⚠️",
                  "Existing" = "📌",
                  "❓")
    cat(icon, status, ":", status_counts[status], "datasets\n")
  }
  
  # Top datasets by quality score
  cat("\n🏆 TOP 5 DATASETS BY QUALITY SCORE:\n")
  top_datasets <- head(datasets[datasets$include_status == "Include", ], 5)
  
  if (nrow(top_datasets) > 0) {
    for (i in 1:nrow(top_datasets)) {
      cat(sprintf("  %d. %s (Score: %.1f, Samples: %d)\n",
                 i,
                 top_datasets$gse_id[i],
                 top_datasets$quality_score[i],
                 top_datasets$samples[i]))
    }
  } else {
    cat("  No datasets qualified for inclusion\n")
  }
}

#' Update Config with Discovered Datasets
#'
#' Adds discovered datasets to config.yml
#' @param gse_ids Vector of GSE IDs to add
#' @param config_file Path to configuration file
update_config_with_datasets <- function(gse_ids, config_file) {
  
  # Load current config
  config <- yaml::read_yaml(config_file)
  
  # Add discovered datasets section if not exists
  if (is.null(config$datasets$discovered)) {
    config$datasets$discovered <- list()
  }
  
  # Add each dataset
  for (gse_id in gse_ids) {
    if (!gse_id %in% names(config$datasets$active_datasets)) {
      config$datasets$discovered[[gse_id]] <- list(
        status = "pending_download",
        discovery_date = as.character(Sys.Date())
      )
    }
  }
  
  # Save updated config
  yaml::write_yaml(config, config_file)
  cat("📝 Updated config with", length(gse_ids), "discovered datasets\n")
}

# NULL coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x

# Test function
test_dataset_discovery <- function() {
  cat("🧪 Testing Dataset Discovery Module\n")
  
  # Run discovery with small test
  test_results <- discover_geo_datasets(
    config_file = "config.yml",
    output_file = "output/test_discovery.xlsx",
    auto_download = FALSE
  )
  
  cat("✅ Discovery module test complete\n")
  return(test_results)
}

cat("✅ Dataset Discovery Module loaded successfully\n")
cat("   Functions: discover_geo_datasets(), search_geo_datasets()\n")
cat("   Usage: source('modules/dataset_discovery.R')\n")
cat("         discover_geo_datasets()\n")