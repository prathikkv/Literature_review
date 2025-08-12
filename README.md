# CAMK2D Research Platform

A production-ready bioinformatics platform for CAMK2D (Calcium/calmodulin-dependent protein kinase II delta) research, specifically designed for atrial fibrillation and heart failure studies.

## Overview

This platform provides a complete pipeline for CAMK2D research:

- **Literature Analysis**: Automated mining of CAMK2D research papers with relevance scoring
- **Dataset Discovery**: Targeted search of GEO database for CAMK2D-relevant expression datasets  
- **Expression Validation**: Download and validation of expression data with CAMK2D detection
- **Analysis-Ready Output**: Processed datasets ready for downstream bioinformatics analysis

## Quick Start

### Simple Execution
```r
# Run the complete platform with production settings
source("run_platform.R")
```

### Command Line
```bash
Rscript run_platform.R
```

## Platform Components

### Core Pipeline
- `01_literature_processing.Rmd` - Literature analysis and paper discovery
- `02_cross_species_discovery.Rmd` - GEO dataset discovery with CAMK2D focus
- `03_integrated_discovery_validation.Rmd` - Expression data validation and integration
- `run_platform.R` - Main execution script with production configuration

### Functions
- `functions/targeted_camk2d_search.R` - CAMK2D-specific search strategies
- `functions/expression_data_validator.R` - Expression data download and validation
- `functions/dge_analysis_framework.R` - Differential expression analysis

### Configuration
- `config.yml` - Production configuration file
- `run_integrated_camk2d_platform.R` - Legacy execution script (use `run_platform.R`)

## Configuration

Edit `config.yml` to customize the platform:

```yaml
# Research Focus
research:
  focus_area: "both"              # "aFIB", "HF", "both"
  species: "human"

# Processing Limits  
datasets:
  max_datasets: 12               # Max datasets to process
  min_samples: 10                # Min samples per dataset
  min_camk2d_expression: 2       # Min CAMK2D expression level

# Expression Settings
expression:
  validation_enabled: true       # Download expression data
  cache_downloads: true         # Cache for performance
  max_retries: 3               # Retry failed downloads
```

## Output Structure

```
output/
├── CAMK2D_Literature_Results_YYYYMMDD.csv      # Literature analysis
├── CAMK2D_Discovery_Results_YYYYMMDD.csv       # Dataset discovery
├── Literature_Analysis_Report.html             # Literature report
├── Dataset_Discovery_Report.html               # Discovery report
├── Integrated_Validation_Report.html           # Validation report
└── expression_data/                            # Expression matrices
    ├── GSE12345_expression.csv
    ├── GSE12345_phenotypes.csv
    └── GSE12345_features.csv
```

## System Requirements

### Software
- R (>= 4.0.0)
- RStudio (recommended)

### R Packages
```r
# Bioconductor packages
BiocManager::install(c("GEOquery", "Biobase", "limma"))

# CRAN packages
install.packages(c("tidyverse", "rmarkdown", "yaml", 
                   "rentrez", "httr", "plotly", "openxlsx",
                   "DT", "knitr", "kableExtra"))
```

### System
- Internet connection for database queries
- ~1GB free disk space for cache and downloads

## Features

### Enhanced Performance
- **Caching**: Downloaded data is cached to avoid re-downloads
- **Retry Logic**: Automatic retries for failed network requests
- **Progress Tracking**: Detailed logging of all operations

### Quality Control
- CAMK2D expression validation in all datasets
- Sample size and data quality filtering
- Case-control design validation
- Phenotype data extraction and standardization

### Production Ready
- Comprehensive error handling
- Configurable parameters
- Clean folder structure
- Detailed documentation

## Troubleshooting

### Common Issues

**Network Timeouts**
- Increase timeout in `config.yml`: `timeout_seconds: 1200`
- Check internet connection and firewall settings

**Memory Issues**
- Reduce `max_datasets` in configuration
- Clear cache: `rm -rf cache/`

**Missing Dependencies**
- Install required packages (see System Requirements)
- Update Bioconductor: `BiocManager::install()`

### Getting Help

1. Check the logs for detailed error messages
2. Verify configuration settings in `config.yml`
3. Test network connectivity to NCBI GEO
4. Ensure all required packages are installed

## Advanced Usage

### Custom Analysis
```r
# Run specific components
rmarkdown::render("01_literature_processing.Rmd")
rmarkdown::render("02_cross_species_discovery.Rmd")

# Custom parameters
rmarkdown::render("03_integrated_discovery_validation.Rmd",
  params = list(
    focus_area = "aFIB",
    max_datasets = 5,
    min_samples = 20
  )
)
```

### Batch Processing
```r
# Process multiple configurations
configs <- list(
  list(focus_area = "aFIB", max_datasets = 5),
  list(focus_area = "HF", max_datasets = 5)
)

for (config in configs) {
  # Run with specific config
}
```

## Platform Status: Production Ready

✅ **Literature Analysis**: Functional  
✅ **Dataset Discovery**: Functional  
✅ **Expression Validation**: Functional  
✅ **Quality Control**: Implemented  
✅ **Error Handling**: Robust  
✅ **Documentation**: Complete  

---

**Ready for CAMK2D research with validated expression data and comprehensive analysis capabilities.**