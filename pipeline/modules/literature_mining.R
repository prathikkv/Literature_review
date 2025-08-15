#!/usr/bin/env Rscript
#' Literature Mining and Hyperlink Module
#' 
#' Intelligent literature mining with dynamic hyperlinks
#' Non-disruptive enhancement to existing pipeline
#' 
#' @author Claude Code Enhancement Module
#' @version 1.0.0

# Load required libraries
suppressPackageStartupMessages({
  library(rentrez)
  library(httr)
  library(jsonlite)
  library(xml2)
  library(tidyverse)
  library(lubridate)
})

#' Mine PubMed Literature
#'
#' Mines PubMed for gene and disease combination literature
#' @param gene Gene symbol or list of genes
#' @param diseases Vector of disease terms
#' @param max_papers Maximum number of papers to retrieve (default 20)
#' @param output_file Optional output file for literature summary
#' @return List with literature mining results
mine_pubmed_literature <- function(gene, 
                                 diseases,
                                 max_papers = 20,
                                 output_file = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║           PUBMED LITERATURE MINING MODULE                    ║\n")
  cat("║           Intelligent Literature Analysis                    ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Handle gene input (single or multiple)
  if (is.character(gene)) {
    genes <- gene
  } else {
    genes <- as.character(gene)
  }
  
  cat("🧬 Genes:", paste(genes, collapse = ", "), "\n")
  cat("🏥 Diseases:", paste(diseases, collapse = ", "), "\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
  
  # Initialize results
  literature_results <- list(
    genes = genes,
    diseases = diseases,
    papers = list(),
    key_findings = list(),
    therapeutic_relevance = list(),
    mechanisms = list(),
    citations = list(),
    hyperlinks = list(),
    summary = character(),
    timestamp = Sys.time()
  )
  
  # Search for each gene-disease combination
  all_papers <- list()
  
  for (gene_symbol in genes) {
    for (disease in diseases) {
      cat("🔍 Searching:", gene_symbol, "+", disease, "\n")
      
      # Build search query
      search_query <- paste0(
        "(", gene_symbol, "[Title/Abstract] OR ", gene_symbol, "[MeSH Terms]) AND ",
        "(", disease, "[Title/Abstract] OR ", disease, "[MeSH Terms]) AND ",
        "Homo sapiens[Organism]"
      )
      
      # Search PubMed
      search_results <- tryCatch({
        entrez_search(db = "pubmed", 
                     term = search_query, 
                     retmax = max_papers,
                     sort = "relevance")
      }, error = function(e) {
        cat("  ⚠️  Search failed:", e$message, "\n")
        return(list(ids = character()))
      })
      
      if (length(search_results$ids) > 0) {
        cat("  ✅ Found", length(search_results$ids), "papers\n")
        
        # Fetch paper details
        papers <- fetch_paper_details(search_results$ids)
        
        # Add gene and disease tags
        for (i in seq_along(papers)) {
          papers[[i]]$gene <- gene_symbol
          papers[[i]]$disease <- disease
        }
        
        all_papers <- c(all_papers, papers)
      }
    }
  }
  
  cat("\n📚 Total papers found:", length(all_papers), "\n")
  
  if (length(all_papers) > 0) {
    # Process and analyze papers
    cat("\n📊 ANALYZING LITERATURE\n")
    cat("═══════════════════════════════════════════════════════════════\n")
    
    # Extract key information
    literature_results$papers <- all_papers
    literature_results$key_findings <- extract_key_findings(all_papers)
    literature_results$therapeutic_relevance <- extract_therapeutic_relevance(all_papers)
    literature_results$mechanisms <- extract_mechanisms(all_papers)
    
    # Generate citations with hyperlinks
    literature_results$citations <- generate_citations(all_papers)
    literature_results$hyperlinks <- generate_hyperlinks(all_papers)
    
    # Create summary
    literature_results$summary <- create_literature_summary(literature_results)
    
    # Save if output file specified
    if (!is.null(output_file)) {
      save_literature_report(literature_results, output_file)
    }
  }
  
  # Print summary
  print_literature_summary(literature_results)
  
  return(literature_results)
}

#' Extract Clinical Trials
#'
#' Searches and extracts clinical trials from ClinicalTrials.gov
#' @param gene Gene symbol or list
#' @param diseases Disease terms
#' @param status Trial status filter (e.g., "Recruiting", "Completed")
#' @param output_file Optional output file
#' @return List with clinical trials data
extract_clinical_trials <- function(gene,
                                  diseases,
                                  status = NULL,
                                  output_file = NULL) {
  
  cat("\n")
  cat("╔══════════════════════════════════════════════════════════════╗\n")
  cat("║           CLINICAL TRIALS EXTRACTION MODULE                  ║\n")
  cat("║           ClinicalTrials.gov Database Search                 ║\n")
  cat("╚══════════════════════════════════════════════════════════════╝\n")
  cat("\n")
  
  # Handle gene input
  if (is.character(gene)) {
    genes <- gene
  } else {
    genes <- as.character(gene)
  }
  
  # Initialize results
  trials_results <- list(
    genes = genes,
    diseases = diseases,
    trials = list(),
    active_trials = list(),
    completed_trials = list(),
    trial_phases = list(),
    interventions = list(),
    hyperlinks = list(),
    summary = character(),
    timestamp = Sys.time()
  )
  
  # ClinicalTrials.gov API endpoint
  base_url <- "https://clinicaltrials.gov/api/query/study_fields"
  
  all_trials <- list()
  
  for (gene_symbol in genes) {
    for (disease in diseases) {
      cat("🔍 Searching trials:", gene_symbol, "+", disease, "\n")
      
      # Build search expression
      search_expr <- paste0(gene_symbol, " AND ", disease)
      
      # API parameters
      params <- list(
        expr = search_expr,
        fields = "NCTId,BriefTitle,Condition,InterventionName,Phase,OverallStatus,StartDate,CompletionDate,EnrollmentCount",
        min_rnk = 1,
        max_rnk = 50,
        fmt = "json"
      )
      
      # Add status filter if specified
      if (!is.null(status)) {
        params$expr <- paste0(params$expr, " AND ", status)
      }
      
      # Make API request
      response <- tryCatch({
        GET(base_url, query = params)
      }, error = function(e) {
        cat("  ⚠️  API request failed:", e$message, "\n")
        return(NULL)
      })
      
      if (!is.null(response) && status_code(response) == 200) {
        content <- content(response, "parsed")
        
        if (!is.null(content$StudyFieldsResponse$StudyFields)) {
          trials <- content$StudyFieldsResponse$StudyFields
          
          cat("  ✅ Found", length(trials), "trials\n")
          
          # Process trials
          for (trial in trials) {
            trial_data <- process_trial_data(trial, gene_symbol, disease)
            all_trials <- append(all_trials, list(trial_data))
          }
        }
      }
    }
  }
  
  cat("\n🏥 Total trials found:", length(all_trials), "\n")
  
  if (length(all_trials) > 0) {
    # Organize trials
    trials_results$trials <- all_trials
    
    # Separate by status
    trials_results$active_trials <- all_trials[sapply(all_trials, function(x) 
      x$status %in% c("Recruiting", "Active, not recruiting", "Enrolling by invitation"))]
    
    trials_results$completed_trials <- all_trials[sapply(all_trials, function(x) 
      x$status == "Completed")]
    
    # Extract phases
    trials_results$trial_phases <- table(sapply(all_trials, function(x) x$phase))
    
    # Generate hyperlinks
    trials_results$hyperlinks <- generate_trial_hyperlinks(all_trials)
    
    # Create summary
    trials_results$summary <- create_trials_summary(trials_results)
    
    # Save if requested
    if (!is.null(output_file)) {
      save_trials_report(trials_results, output_file)
    }
  }
  
  # Print summary
  print_trials_summary(trials_results)
  
  return(trials_results)
}

#' Generate Hyperlinks
#'
#' Generates hyperlinks for various biological entities
#' @param text_content Text or data to add hyperlinks to
#' @param link_type Type of links to generate
#' @return Text with embedded hyperlinks or list of hyperlinks
generate_hyperlinks <- function(text_content, link_type = "all") {
  
  hyperlinks <- list()
  
  # Define link templates
  link_templates <- list(
    gene = "https://www.ncbi.nlm.nih.gov/gene/?term=",
    pubmed = "https://pubmed.ncbi.nlm.nih.gov/",
    kegg_pathway = "https://www.kegg.jp/kegg-bin/show_pathway?",
    reactome = "https://reactome.org/content/detail/",
    go_term = "http://amigo.geneontology.org/amigo/term/",
    clinicaltrials = "https://clinicaltrials.gov/ct2/show/",
    disease_ontology = "https://disease-ontology.org/?id=",
    drugbank = "https://go.drugbank.com/drugs/",
    uniprot = "https://www.uniprot.org/uniprot/"
  )
  
  # Process based on content type
  if (is.character(text_content)) {
    # Single text string - embed hyperlinks
    processed_text <- text_content
    
    if (link_type %in% c("all", "gene")) {
      # Find and link gene symbols (uppercase letters/numbers pattern)
      gene_pattern <- "\\b[A-Z][A-Z0-9]{2,}\\b"
      genes <- unique(unlist(regmatches(processed_text, 
                                       gregexpr(gene_pattern, processed_text))))
      
      for (gene in genes) {
        gene_url <- paste0(link_templates$gene, gene)
        hyperlinks[[paste0("gene_", gene)]] <- gene_url
      }
    }
    
    if (link_type %in% c("all", "pubmed")) {
      # Find and link PMIDs
      pmid_pattern <- "PMID:?\\s*([0-9]+)"
      pmids <- unique(unlist(regmatches(processed_text, 
                                       gregexpr(pmid_pattern, processed_text))))
      
      for (pmid in pmids) {
        pmid_num <- gsub(".*?([0-9]+).*", "\\1", pmid)
        pmid_url <- paste0(link_templates$pubmed, pmid_num)
        hyperlinks[[paste0("pmid_", pmid_num)]] <- pmid_url
      }
    }
    
  } else if (is.list(text_content)) {
    # List of items - generate appropriate links
    for (item in text_content) {
      if (!is.null(item$pmid)) {
        hyperlinks[[paste0("pmid_", item$pmid)]] <- 
          paste0(link_templates$pubmed, item$pmid)
      }
      
      if (!is.null(item$nct_id)) {
        hyperlinks[[paste0("trial_", item$nct_id)]] <- 
          paste0(link_templates$clinicaltrials, item$nct_id)
      }
      
      if (!is.null(item$gene)) {
        hyperlinks[[paste0("gene_", item$gene)]] <- 
          paste0(link_templates$gene, item$gene)
      }
    }
  }
  
  return(hyperlinks)
}

#' Fetch Paper Details
#'
#' Fetches detailed information for PubMed papers
#' @param pmids Vector of PubMed IDs
#' @return List of paper details
fetch_paper_details <- function(pmids) {
  papers <- list()
  
  # Fetch summaries
  summaries <- tryCatch({
    entrez_summary(db = "pubmed", id = pmids)
  }, error = function(e) {
    return(list())
  })
  
  # Process each summary
  for (i in seq_along(summaries)) {
    if (is.list(summaries[[i]])) {
      paper <- list(
        pmid = summaries[[i]]$uid,
        title = summaries[[i]]$title,
        authors = paste(summaries[[i]]$authors$name[1:min(3, length(summaries[[i]]$authors$name))], 
                       collapse = ", "),
        journal = summaries[[i]]$source,
        year = as.numeric(substr(summaries[[i]]$pubdate, 1, 4)),
        doi = summaries[[i]]$elocationid,
        abstract = NA  # Will fetch if needed
      )
      
      papers <- append(papers, list(paper))
    }
  }
  
  return(papers)
}

#' Extract Key Findings
#'
#' Extracts key findings from literature
#' @param papers List of paper data
#' @return List of key findings
extract_key_findings <- function(papers) {
  findings <- list()
  
  # Keywords to look for
  keywords <- c("significant", "associated", "upregulated", "downregulated", 
               "increased", "decreased", "correlated", "therapeutic", "target")
  
  for (paper in papers) {
    if (!is.null(paper$title)) {
      # Check title for key terms
      for (keyword in keywords) {
        if (grepl(keyword, paper$title, ignore.case = TRUE)) {
          findings <- append(findings, list(list(
            pmid = paper$pmid,
            finding = paper$title,
            gene = paper$gene,
            disease = paper$disease,
            year = paper$year
          )))
          break
        }
      }
    }
  }
  
  return(findings)
}

#' Extract Therapeutic Relevance
#'
#' Extracts therapeutic relevance from papers
#' @param papers List of papers
#' @return List of therapeutic findings
extract_therapeutic_relevance <- function(papers) {
  therapeutic <- list()
  
  # Therapeutic keywords
  therapy_keywords <- c("therapeutic", "treatment", "drug", "inhibitor", 
                       "agonist", "antagonist", "clinical", "therapy")
  
  for (paper in papers) {
    if (!is.null(paper$title)) {
      for (keyword in therapy_keywords) {
        if (grepl(keyword, paper$title, ignore.case = TRUE)) {
          therapeutic <- append(therapeutic, list(list(
            pmid = paper$pmid,
            relevance = paper$title,
            gene = paper$gene,
            disease = paper$disease,
            keyword = keyword
          )))
          break
        }
      }
    }
  }
  
  return(therapeutic)
}

#' Extract Mechanisms
#'
#' Extracts mechanism information
#' @param papers List of papers
#' @return List of mechanisms
extract_mechanisms <- function(papers) {
  mechanisms <- list()
  
  # Mechanism keywords
  mechanism_keywords <- c("mechanism", "pathway", "signaling", "regulation", 
                         "activation", "inhibition", "phosphorylation")
  
  for (paper in papers) {
    if (!is.null(paper$title)) {
      for (keyword in mechanism_keywords) {
        if (grepl(keyword, paper$title, ignore.case = TRUE)) {
          mechanisms <- append(mechanisms, list(list(
            pmid = paper$pmid,
            mechanism = paper$title,
            gene = paper$gene,
            disease = paper$disease,
            keyword = keyword
          )))
          break
        }
      }
    }
  }
  
  return(mechanisms)
}

#' Generate Citations
#'
#' Generates formatted citations
#' @param papers List of papers
#' @return List of formatted citations
generate_citations <- function(papers) {
  citations <- list()
  
  for (paper in papers) {
    citation <- paste0(
      paper$authors, " (",
      paper$year, "). ",
      paper$title, " ",
      paper$journal, ". ",
      "PMID: ", paper$pmid
    )
    
    citations <- append(citations, citation)
  }
  
  return(citations)
}

#' Process Trial Data
#'
#' Processes raw trial data from API
#' @param trial Raw trial data
#' @param gene Associated gene
#' @param disease Associated disease
#' @return Processed trial data
process_trial_data <- function(trial, gene, disease) {
  trial_data <- list(
    nct_id = trial$NCTId[[1]],
    title = trial$BriefTitle[[1]],
    condition = paste(trial$Condition, collapse = "; "),
    intervention = paste(trial$InterventionName, collapse = "; "),
    phase = trial$Phase[[1]],
    status = trial$OverallStatus[[1]],
    start_date = trial$StartDate[[1]],
    completion_date = trial$CompletionDate[[1]],
    enrollment = trial$EnrollmentCount[[1]],
    gene = gene,
    disease = disease
  )
  
  return(trial_data)
}

#' Generate Trial Hyperlinks
#'
#' Generates hyperlinks for clinical trials
#' @param trials List of trials
#' @return List of hyperlinks
generate_trial_hyperlinks <- function(trials) {
  hyperlinks <- list()
  
  for (trial in trials) {
    url <- paste0("https://clinicaltrials.gov/ct2/show/", trial$nct_id)
    hyperlinks[[trial$nct_id]] <- url
  }
  
  return(hyperlinks)
}

#' Create Literature Summary
#'
#' Creates summary of literature findings
#' @param results Literature results
#' @return Summary text
create_literature_summary <- function(results) {
  summary <- character()
  
  summary <- c(summary, paste("Analyzed", length(results$papers), "papers for",
                              paste(results$genes, collapse = ", "),
                              "in", paste(results$diseases, collapse = ", ")))
  
  if (length(results$key_findings) > 0) {
    summary <- c(summary, paste("\nKey findings:", length(results$key_findings), "significant associations"))
  }
  
  if (length(results$therapeutic_relevance) > 0) {
    summary <- c(summary, paste("Therapeutic relevance:", length(results$therapeutic_relevance), "papers"))
  }
  
  if (length(results$mechanisms) > 0) {
    summary <- c(summary, paste("Mechanisms described:", length(results$mechanisms), "papers"))
  }
  
  return(paste(summary, collapse = "\n"))
}

#' Create Trials Summary
#'
#' Creates summary of clinical trials
#' @param results Trials results
#' @return Summary text
create_trials_summary <- function(results) {
  summary <- character()
  
  summary <- c(summary, paste("Found", length(results$trials), "clinical trials"))
  
  if (length(results$active_trials) > 0) {
    summary <- c(summary, paste("Active trials:", length(results$active_trials)))
  }
  
  if (length(results$completed_trials) > 0) {
    summary <- c(summary, paste("Completed trials:", length(results$completed_trials)))
  }
  
  # Phase distribution
  if (length(results$trial_phases) > 0) {
    phase_text <- paste(names(results$trial_phases), ":", results$trial_phases, collapse = ", ")
    summary <- c(summary, paste("Phase distribution:", phase_text))
  }
  
  return(paste(summary, collapse = "\n"))
}

#' Save Literature Report
#'
#' Saves literature mining results
#' @param results Literature results
#' @param output_file Output file path
save_literature_report <- function(results, output_file) {
  # Ensure output directory exists
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Create data frame for saving
  papers_df <- do.call(rbind, lapply(results$papers, as.data.frame))
  
  # Save to CSV
  write.csv(papers_df, output_file, row.names = FALSE)
  cat("\n📄 Literature report saved to:", output_file, "\n")
}

#' Save Trials Report
#'
#' Saves clinical trials results
#' @param results Trials results
#' @param output_file Output file path
save_trials_report <- function(results, output_file) {
  # Ensure output directory exists
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Create data frame
  trials_df <- do.call(rbind, lapply(results$trials, as.data.frame))
  
  # Save to CSV
  write.csv(trials_df, output_file, row.names = FALSE)
  cat("\n📄 Trials report saved to:", output_file, "\n")
}

#' Print Literature Summary
#'
#' Prints literature mining summary
#' @param results Literature results
print_literature_summary <- function(results) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("📚 LITERATURE MINING SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  cat("📊 Papers analyzed:", length(results$papers), "\n")
  cat("🔑 Key findings:", length(results$key_findings), "\n")
  cat("💊 Therapeutic relevance:", length(results$therapeutic_relevance), "\n")
  cat("⚙️  Mechanisms identified:", length(results$mechanisms), "\n")
  cat("🔗 Hyperlinks generated:", length(results$hyperlinks), "\n")
  
  if (length(results$key_findings) > 0) {
    cat("\n🎯 TOP FINDINGS:\n")
    cat("───────────────────────────────────────────────────────────────\n")
    
    for (i in 1:min(3, length(results$key_findings))) {
      finding <- results$key_findings[[i]]
      cat(sprintf("%d. %s (%s)\n   Gene: %s | Disease: %s | Year: %d | PMID: %s\n",
                 i,
                 substr(finding$finding, 1, 60),
                 if(nchar(finding$finding) > 60) "..." else "",
                 finding$gene,
                 finding$disease,
                 finding$year,
                 finding$pmid))
    }
  }
  
  cat("\n✅ Literature mining completed successfully!\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
}

#' Print Trials Summary
#'
#' Prints clinical trials summary
#' @param results Trials results
print_trials_summary <- function(results) {
  cat("\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  cat("🏥 CLINICAL TRIALS SUMMARY\n")
  cat("═══════════════════════════════════════════════════════════════\n")
  
  cat("📊 Total trials:", length(results$trials), "\n")
  cat("✅ Active trials:", length(results$active_trials), "\n")
  cat("✓ Completed trials:", length(results$completed_trials), "\n")
  
  if (length(results$trial_phases) > 0) {
    cat("\n📈 PHASE DISTRIBUTION:\n")
    for (phase in names(results$trial_phases)) {
      cat(sprintf("  %s: %d trials\n", phase, results$trial_phases[[phase]]))
    }
  }
  
  if (length(results$active_trials) > 0) {
    cat("\n🎯 ACTIVE TRIALS:\n")
    cat("───────────────────────────────────────────────────────────────\n")
    
    for (i in 1:min(3, length(results$active_trials))) {
      trial <- results$active_trials[[i]]
      cat(sprintf("%d. %s\n   NCT: %s | Phase: %s | Status: %s\n",
                 i,
                 substr(trial$title, 1, 60),
                 trial$nct_id,
                 trial$phase,
                 trial$status))
    }
  }
  
  cat("\n✅ Clinical trials extraction completed!\n")
  cat("═══════════════════════════════════════════════════════════════\n\n")
}

# Module loading confirmation
cat("✅ Literature Mining Module loaded successfully\n")
cat("   Functions: mine_pubmed_literature(), extract_clinical_trials(), generate_hyperlinks()\n")
cat("   Version: 1.0.0\n\n")