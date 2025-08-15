#!/usr/bin/env Rscript
#' Setup Script for CAMK2D Analysis Pipeline
#' 
#' Installs all required packages and configures the analysis environment

cat("FIX: CAMK2D ANALYSIS PIPELINE SETUP\n")
cat("==================================\n")
cat("TIME: Started:", Sys.time(), "\n\n")

# Set options for installation
options(repos = c(CRAN = "https://cran.r-project.org"))
options(timeout = 600)

# Function to install packages with error handling
install_package_safely <- function(package_name, source = "CRAN") {
  cat("PACKAGE: Installing", package_name, "from", source, "...")
  
  tryCatch({
    if (source == "CRAN") {
      if (!requireNamespace(package_name, quietly = TRUE)) {
        install.packages(package_name, dependencies = TRUE)
      }
    } else if (source == "Bioconductor") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      if (!requireNamespace(package_name, quietly = TRUE)) {
        BiocManager::install(package_name, dependencies = TRUE)
      }
    }
    cat(" SUCCESS:\n")
    return(TRUE)
  }, error = function(e) {
    cat(" ERROR: Error:", e$message, "\n")
    return(FALSE)
  })
}

# Core CRAN packages
cran_packages <- c(
  "tidyverse", "dplyr", "ggplot2", "readr", "stringr",
  "httr", "jsonlite", "yaml", "openxlsx", "writexl",
  "data.table", "R.utils", "parallel", "doParallel",
  "plotly", "DT", "knitr", "rmarkdown", "gridExtra",
  "RColorBrewer", "pheatmap", "corrplot", "VennDiagram",
  "metafor", "meta", "forestplot"
)

# Bioconductor packages  
bioc_packages <- c(
  "GEOquery", "ArrayExpress", "Biobase", "BiocGenerics",
  "limma", "DESeq2", "edgeR", "biomaRt", "AnnotationDbi",
  "clusterProfiler", "enrichplot", "DOSE", "ReactomePA",
  "org.Hs.eg.db", "org.Mm.eg.db", "org.Rn.eg.db",
  "ComplexHeatmap", "SummarizedExperiment", "GenomicRanges",
  "WGCNA", "KEGGREST", "recount3", "rhdf5", "sva", "ChemmineR"
)

cat("FIX: Installing CRAN packages...\n")
cran_success <- sum(sapply(cran_packages, function(pkg) install_package_safely(pkg, "CRAN")))

cat("\nGENETIC: Installing Bioconductor packages...\n")
bioc_success <- sum(sapply(bioc_packages, function(pkg) install_package_safely(pkg, "Bioconductor")))

# Test core package loading
cat("\nTEST: Testing core package loading...\n")
core_test_packages <- c("GEOquery", "limma", "tidyverse", "metafor", "clusterProfiler")
working_packages <- sum(sapply(core_test_packages, function(pkg) requireNamespace(pkg, quietly = TRUE)))

cat("\nDATA: SETUP SUMMARY\n")
cat("================\n")
cat("CRAN packages installed:", cran_success, "/", length(cran_packages), "\n")
cat("Bioconductor packages installed:", bioc_success, "/", length(bioc_packages), "\n") 
cat("Core packages working:", working_packages, "/", length(core_test_packages), "\n")

if (working_packages == length(core_test_packages)) {
  cat("\nCELEBRATE: SUCCESS: Setup completed successfully!\n")
  cat("SUCCESS: Pipeline is ready for execution\n")
  cat("LAUNCH: Next step: Rscript run_pipeline.R\n")
} else {
  cat("\nWARNING: WARNING: Some packages failed to install\n")
  cat("FIX: Please check error messages above\n")
}

cat("\nTIME: Setup completed:", Sys.time(), "\n")