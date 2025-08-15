#!/usr/bin/env Rscript
#' Gene Family Discovery Module
#' 
#' Automatically discovers gene families from primary gene input
#' Non-disruptive enhancement to existing pipeline
#' 
#' @author Claude Code Enhancement Module
#' @version 1.0.0

# Load required libraries
suppressPackageStartupMessages({
  library(httr)
  library(jsonlite)
  library(xml2)
  library(tidyverse)
  library(yaml)
  library(rentrez)
  library(biomaRt)
})

#' Discover Gene Family
#'
#' Automatically discovers gene family members for a primary gene
#' @param primary_gene Primary gene symbol (e.g., "CAMK2D")
#' @param config_file Configuration file path
#' @param output_file Output file for gene family report
#' @param update_config Whether to update config.yml with discoveries
#' @return List with gene family members and evidence scores
discover_gene_family <- function(primary_gene,
                                config_file = "config.yml",
                                output_file = "output/gene_family_report.csv",
                                update_config = FALSE) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║           GENE FAMILY DISCOVERY MODULE                       ║\n")
  cat("║           Intelligent Family Member Identification           ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  cat("🧬 Primary Gene:", primary_gene, "\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Initialize results
  family_results <- list(
    primary_gene = primary_gene,
    family_members = character(),
    evidence_scores = list(),
    discovery_methods = list(),
    total_members = 0,
    timestamp = Sys.time()
  )
  
  # 1. NCBI Gene Database Search
  cat("🔍 METHOD 1: NCBI Gene Database Search\n")
  ncbi_family <- discover_via_ncbi(primary_gene)
  family_results$discovery_methods$ncbi <- ncbi_family
  
  # 2. UniProt Family Classification
  cat("\n🔍 METHOD 2: UniProt Protein Family Classification\n")
  uniprot_family <- discover_via_uniprot(primary_gene)
  family_results$discovery_methods$uniprot <- uniprot_family
  
  # 3. KEGG Pathway Co-membership
  cat("\n🔍 METHOD 3: KEGG Pathway Co-membership\n")
  kegg_family <- discover_via_kegg(primary_gene)
  family_results$discovery_methods$kegg <- kegg_family
  
  # 4. GO Term Similarity
  cat("\n🔍 METHOD 4: GO Term Functional Similarity\n")
  go_family <- discover_via_go(primary_gene)
  family_results$discovery_methods$go <- go_family
  
  # 5. Literature Co-occurrence (PubMed)
  cat("\n🔍 METHOD 5: PubMed Literature Co-occurrence\n")
  pubmed_family <- discover_via_pubmed(primary_gene)
  family_results$discovery_methods$pubmed <- pubmed_family
  
  # 6. BioMart Homology Search
  cat("\n🔍 METHOD 6: BioMart Homology and Paralog Search\n")
  biomart_family <- discover_via_biomart(primary_gene)
  family_results$discovery_methods$biomart <- biomart_family
  
  # Combine and score all discoveries
  cat("\n📊 COMBINING AND SCORING DISCOVERIES\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  all_genes <- combine_discoveries(family_results$discovery_methods)
  scored_genes <- calculate_evidence_scores(all_genes, primary_gene)
  
  # Filter and rank
  final_family <- scored_genes %>%
    filter(evidence_score >= 3) %>%  # Minimum evidence threshold
    arrange(desc(evidence_score))
  
  family_results$family_members <- final_family$gene
  family_results$evidence_scores <- final_family
  family_results$total_members <- nrow(final_family)
  
  # Save report
  if (!is.null(output_file)) {
    save_family_report(final_family, output_file, primary_gene)
  }
  
  # Update config if requested
  if (update_config && file.exists(config_file)) {
    update_config_with_family(config_file, final_family$gene, primary_gene)
  }
  
  # Print summary
  print_family_summary(family_results)
  
  return(family_results)
}

#' Discover via NCBI Gene Database
#'
#' Searches NCBI Gene for family members
#' @param gene Primary gene symbol
#' @return Vector of related genes
discover_via_ncbi <- function(gene) {
  tryCatch({
    # Search for gene family name
    search_query <- paste0(gene, "[Gene Name] AND Homo sapiens[Organism]")
    search_results <- entrez_search(db = "gene", term = search_query, retmax = 1)
    
    if (length(search_results$ids) > 0) {
      # Get gene summary
      gene_summary <- entrez_summary(db = "gene", id = search_results$ids[1])
      
      # Extract family information from summary
      gene_name <- gene_summary$name
      gene_desc <- gene_summary$description
      
      # Look for family patterns (e.g., CAMK family)
      family_pattern <- gsub("[0-9]+.*", "", gene)  # Remove numbers and suffixes
      
      # Search for all family members
      family_query <- paste0(family_pattern, "[Gene Name] AND Homo sapiens[Organism]")
      family_search <- entrez_search(db = "gene", term = family_query, retmax = 50)
      
      if (length(family_search$ids) > 0) {
        # Get all family member names
        family_summaries <- entrez_summary(db = "gene", id = family_search$ids)
        family_genes <- sapply(family_summaries, function(x) x$name)
        
        cat("  ✅ Found", length(family_genes), "potential family members\n")
        return(unique(family_genes))
      }
    }
    
    cat("  ⚠️  No family members found via NCBI\n")
    return(character())
    
  }, error = function(e) {
    cat("  ❌ NCBI search failed:", e$message, "\n")
    return(character())
  })
}

#' Discover via UniProt
#'
#' Searches UniProt for protein family members
#' @param gene Primary gene symbol
#' @return Vector of related genes
discover_via_uniprot <- function(gene) {
  tryCatch({
    # UniProt REST API query
    base_url <- "https://rest.uniprot.org/uniprotkb/search"
    query <- paste0("gene:", gene, " AND organism_id:9606")
    
    response <- GET(base_url, query = list(
      query = query,
      format = "json",
      size = 10
    ))
    
    if (status_code(response) == 200) {
      content <- content(response, "parsed")
      
      if (length(content$results) > 0) {
        # Extract family information
        protein_families <- character()
        
        for (result in content$results) {
          # Look for family annotations
          if (!is.null(result$proteinDescription)) {
            families <- result$proteinDescription$recommendedName$fullName
            protein_families <- c(protein_families, families)
          }
        }
        
        # Search for other proteins in same family
        if (length(protein_families) > 0) {
          family_pattern <- gsub(".*family.*", "", protein_families[1])
          
          # Get family members
          family_query <- paste0("family:", family_pattern, " AND organism_id:9606")
          family_response <- GET(base_url, query = list(
            query = family_query,
            format = "json",
            size = 50
          ))
          
          if (status_code(family_response) == 200) {
            family_content <- content(family_response, "parsed")
            family_genes <- character()
            
            for (result in family_content$results) {
              if (!is.null(result$genes)) {
                for (gene_info in result$genes) {
                  family_genes <- c(family_genes, gene_info$geneName$value)
                }
              }
            }
            
            cat("  ✅ Found", length(unique(family_genes)), "family members via UniProt\n")
            return(unique(family_genes))
          }
        }
      }
    }
    
    cat("  ⚠️  No family members found via UniProt\n")
    return(character())
    
  }, error = function(e) {
    cat("  ❌ UniProt search failed:", e$message, "\n")
    return(character())
  })
}

#' Discover via KEGG Pathways
#'
#' Finds genes in same KEGG pathways
#' @param gene Primary gene symbol
#' @return Vector of related genes
discover_via_kegg <- function(gene) {
  tryCatch({
    # KEGG API endpoints
    base_url <- "https://rest.kegg.jp"
    
    # Find KEGG ID for gene
    gene_search <- GET(paste0(base_url, "/find/hsa/", gene))
    
    if (status_code(gene_search) == 200) {
      gene_content <- content(gene_search, "text")
      
      if (nchar(gene_content) > 0) {
        # Extract KEGG gene ID
        gene_lines <- strsplit(gene_content, "\n")[[1]]
        kegg_id <- gsub("\t.*", "", gene_lines[1])
        
        # Get pathways for this gene
        pathway_response <- GET(paste0(base_url, "/link/pathway/", kegg_id))
        
        if (status_code(pathway_response) == 200) {
          pathway_content <- content(pathway_response, "text")
          pathway_lines <- strsplit(pathway_content, "\n")[[1]]
          pathways <- unique(gsub(".*\t", "", pathway_lines))
          
          # Get all genes in these pathways
          pathway_genes <- character()
          
          for (pathway in pathways[1:min(5, length(pathways))]) {  # Limit to top 5 pathways
            if (nchar(pathway) > 0) {
              genes_response <- GET(paste0(base_url, "/link/hsa/", pathway))
              
              if (status_code(genes_response) == 200) {
                genes_content <- content(genes_response, "text")
                genes_lines <- strsplit(genes_content, "\n")[[1]]
                
                for (line in genes_lines) {
                  if (nchar(line) > 0) {
                    gene_id <- gsub(".*\t", "", line)
                    # Convert KEGG ID to gene symbol
                    info_response <- GET(paste0(base_url, "/get/", gene_id))
                    if (status_code(info_response) == 200) {
                      info_content <- content(info_response, "text")
                      # Extract gene symbol from NAME field
                      name_match <- regmatches(info_content, 
                                              regexpr("NAME\\s+([A-Z0-9]+)", info_content))
                      if (length(name_match) > 0) {
                        symbol <- gsub("NAME\\s+", "", name_match[1])
                        pathway_genes <- c(pathway_genes, symbol)
                      }
                    }
                  }
                }
              }
            }
          }
          
          pathway_genes <- unique(pathway_genes)
          cat("  ✅ Found", length(pathway_genes), "co-pathway genes via KEGG\n")
          return(pathway_genes)
        }
      }
    }
    
    cat("  ⚠️  No pathway members found via KEGG\n")
    return(character())
    
  }, error = function(e) {
    cat("  ❌ KEGG search failed:", e$message, "\n")
    return(character())
  })
}

#' Discover via GO Terms
#'
#' Finds genes with similar GO terms
#' @param gene Primary gene symbol
#' @return Vector of related genes
discover_via_go <- function(gene) {
  tryCatch({
    # For GO analysis, we'll use biomaRt
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Get GO terms for primary gene
    go_terms <- getBM(
      attributes = c("go_id", "name_1006", "namespace_1003"),
      filters = "hgnc_symbol",
      values = gene,
      mart = ensembl
    )
    
    if (nrow(go_terms) > 0) {
      # Find genes with similar GO terms
      similar_genes <- character()
      
      # Focus on molecular function GO terms
      mf_terms <- go_terms %>%
        filter(namespace_1003 == "molecular_function") %>%
        pull(go_id)
      
      if (length(mf_terms) > 0) {
        # Get genes with same molecular functions
        for (go_term in mf_terms[1:min(5, length(mf_terms))]) {
          genes_with_term <- getBM(
            attributes = c("hgnc_symbol"),
            filters = "go",
            values = go_term,
            mart = ensembl
          )
          
          similar_genes <- c(similar_genes, genes_with_term$hgnc_symbol)
        }
      }
      
      similar_genes <- unique(similar_genes)
      cat("  ✅ Found", length(similar_genes), "functionally similar genes via GO\n")
      return(similar_genes)
    }
    
    cat("  ⚠️  No GO-similar genes found\n")
    return(character())
    
  }, error = function(e) {
    cat("  ❌ GO search failed:", e$message, "\n")
    return(character())
  })
}

#' Discover via PubMed Literature
#'
#' Finds genes frequently co-mentioned in literature
#' @param gene Primary gene symbol
#' @return Vector of related genes
discover_via_pubmed <- function(gene) {
  tryCatch({
    # Search for papers mentioning the gene
    search_query <- paste0(gene, "[Title/Abstract] AND gene family[Title/Abstract]")
    search_results <- entrez_search(db = "pubmed", term = search_query, retmax = 50)
    
    if (length(search_results$ids) > 0) {
      # Get abstracts
      abstracts <- entrez_fetch(db = "pubmed", id = search_results$ids, 
                               rettype = "abstract", retmode = "text")
      
      # Look for gene family patterns
      family_pattern <- gsub("[0-9]+.*", "", gene)
      
      # Common gene family member patterns
      patterns <- c(
        paste0(family_pattern, "[0-9]+[A-Z]*"),  # e.g., CAMK1, CAMK2A
        paste0(family_pattern, "[A-Z][0-9]*"),    # e.g., CAMKA1
        paste0(family_pattern, "-[A-Z]+")         # e.g., CAMK-II
      )
      
      mentioned_genes <- character()
      
      for (pattern in patterns) {
        matches <- gregexpr(pattern, abstracts, ignore.case = FALSE)
        extracted <- regmatches(abstracts, matches)
        mentioned_genes <- c(mentioned_genes, unlist(extracted))
      }
      
      # Clean and count mentions
      mentioned_genes <- unique(toupper(mentioned_genes))
      mentioned_genes <- mentioned_genes[mentioned_genes != toupper(gene)]
      
      cat("  ✅ Found", length(mentioned_genes), "co-mentioned genes in literature\n")
      return(mentioned_genes)
    }
    
    cat("  ⚠️  No literature co-mentions found\n")
    return(character())
    
  }, error = function(e) {
    cat("  ❌ PubMed search failed:", e$message, "\n")
    return(character())
  })
}

#' Discover via BioMart
#'
#' Finds paralogs and homologs
#' @param gene Primary gene symbol
#' @return Vector of related genes
discover_via_biomart <- function(gene) {
  tryCatch({
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    # Get paralogs
    paralogs <- getBM(
      attributes = c("hsapiens_paralog_associated_gene_name"),
      filters = "hgnc_symbol",
      values = gene,
      mart = ensembl
    )
    
    if (nrow(paralogs) > 0) {
      paralog_genes <- unique(paralogs$hsapiens_paralog_associated_gene_name)
      paralog_genes <- paralog_genes[nchar(paralog_genes) > 0]
      
      cat("  ✅ Found", length(paralog_genes), "paralogs via BioMart\n")
      return(paralog_genes)
    }
    
    cat("  ⚠️  No paralogs found via BioMart\n")
    return(character())
    
  }, error = function(e) {
    cat("  ❌ BioMart search failed:", e$message, "\n")
    return(character())
  })
}

#' Combine Discoveries from All Methods
#'
#' Combines gene lists from all discovery methods
#' @param discovery_methods List of discovery results
#' @return Data frame with genes and their discovery sources
combine_discoveries <- function(discovery_methods) {
  all_genes <- data.frame(
    gene = character(),
    sources = character(),
    source_count = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Collect all unique genes
  unique_genes <- character()
  
  for (method in names(discovery_methods)) {
    genes <- discovery_methods[[method]]
    unique_genes <- unique(c(unique_genes, genes))
  }
  
  # Count sources for each gene
  for (gene in unique_genes) {
    sources <- character()
    
    for (method in names(discovery_methods)) {
      if (gene %in% discovery_methods[[method]]) {
        sources <- c(sources, method)
      }
    }
    
    all_genes <- rbind(all_genes, data.frame(
      gene = gene,
      sources = paste(sources, collapse = ","),
      source_count = length(sources),
      stringsAsFactors = FALSE
    ))
  }
  
  return(all_genes)
}

#' Calculate Evidence Scores
#'
#' Calculates evidence scores for gene family members
#' @param gene_data Data frame with genes and sources
#' @param primary_gene Primary gene for context
#' @return Data frame with evidence scores
calculate_evidence_scores <- function(gene_data, primary_gene) {
  if (nrow(gene_data) == 0) {
    return(data.frame(
      gene = character(),
      evidence_score = numeric(),
      sources = character(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Weight different sources
  source_weights <- list(
    ncbi = 2.0,      # High confidence
    uniprot = 2.0,   # High confidence
    kegg = 1.5,      # Medium-high confidence
    go = 1.5,        # Medium-high confidence
    biomart = 2.0,   # High confidence (paralogs)
    pubmed = 1.0     # Medium confidence
  )
  
  gene_data$evidence_score <- 0
  
  for (i in 1:nrow(gene_data)) {
    sources <- strsplit(gene_data$sources[i], ",")[[1]]
    score <- 0
    
    for (source in sources) {
      score <- score + (source_weights[[source]] %||% 1.0)
    }
    
    # Bonus for multiple sources
    if (length(sources) > 1) {
      score <- score * (1 + 0.1 * (length(sources) - 1))
    }
    
    # Pattern matching bonus
    family_pattern <- gsub("[0-9]+.*", "", primary_gene)
    if (grepl(family_pattern, gene_data$gene[i], ignore.case = TRUE)) {
      score <- score * 1.2
    }
    
    gene_data$evidence_score[i] <- round(score, 2)
  }
  
  # Add evidence interpretation
  gene_data$evidence_level <- cut(
    gene_data$evidence_score,
    breaks = c(-Inf, 2, 4, 6, Inf),
    labels = c("Low", "Moderate", "High", "Very High")
  )
  
  return(gene_data %>% arrange(desc(evidence_score)))
}

#' Save Family Report
#'
#' Saves gene family discovery report
#' @param family_data Data frame with family members
#' @param output_file Output file path
#' @param primary_gene Primary gene name
save_family_report <- function(family_data, output_file, primary_gene) {
  # Ensure output directory exists
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Add metadata
  report_data <- family_data %>%
    mutate(
      primary_gene = primary_gene,
      discovery_date = Sys.Date()
    )
  
  # Save to CSV
  write.csv(report_data, output_file, row.names = FALSE)
  cat("\n📄 Gene family report saved to:", output_file, "\n")
}

#' Update Config with Family
#'
#' Updates configuration file with discovered gene family
#' @param config_file Configuration file path
#' @param family_genes Vector of family gene symbols
#' @param primary_gene Primary gene symbol
update_config_with_family <- function(config_file, family_genes, primary_gene) {
  tryCatch({
    config <- yaml::read_yaml(config_file)
    
    # Update or create gene family section
    if (is.null(config$gene_analysis)) {
      config$gene_analysis <- list()
    }
    
    config$gene_analysis$primary_gene <- primary_gene
    config$gene_analysis$gene_family <- unique(c(primary_gene, family_genes))
    config$gene_analysis$family_discovery_date <- as.character(Sys.Date())
    config$gene_analysis$family_size <- length(config$gene_analysis$gene_family)
    
    # Write updated config
    yaml::write_yaml(config, config_file)
    cat("✅ Configuration updated with", length(family_genes), "family members\n")
    
  }, error = function(e) {
    cat("⚠️  Could not update configuration:", e$message, "\n")
  })
}

#' Print Family Summary
#'
#' Prints summary of discovered gene family
#' @param results Gene family discovery results
print_family_summary <- function(results) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📊 GENE FAMILY DISCOVERY SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("🧬 Primary Gene:", results$primary_gene, "\n")
  cat("👥 Family Members Found:", results$total_members, "\n")
  
  if (results$total_members > 0) {
    # Show top members
    top_members <- head(results$evidence_scores, 10)
    
    cat("\n🏆 TOP FAMILY MEMBERS:\n")
    cat("───────────────────────────────────────────────────────────────\n")
    
    for (i in 1:nrow(top_members)) {
      cat(sprintf("%2d. %-10s (Score: %4.1f, Evidence: %s, Sources: %s)\n",
                  i,
                  top_members$gene[i],
                  top_members$evidence_score[i],
                  top_members$evidence_level[i],
                  top_members$sources[i]))
    }
    
    # Method summary
    cat("\n📈 DISCOVERY METHOD CONTRIBUTIONS:\n")
    cat("───────────────────────────────────────────────────────────────\n")
    
    for (method in names(results$discovery_methods)) {
      count <- length(results$discovery_methods[[method]])
      if (count > 0) {
        cat(sprintf("  %-15s: %3d genes\n", 
                   paste0(toupper(substr(method, 1, 1)), substr(method, 2, nchar(method))),
                   count))
      }
    }
  }
  
  cat("\n✅ Gene family discovery completed successfully!\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
}

# Module loading confirmation
cat("✅ Gene Family Discovery Module loaded successfully\n")
cat("   Function: discover_gene_family()\n")
cat("   Version: 1.0.0\n\n")