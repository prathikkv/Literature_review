# CAMK2D Platform Usage Guide

## Quick Start

### Option 1: Simple Execution (Recommended)
```r
# Run with production settings
source("run_platform.R")
```

### Option 2: Command Line
```bash
Rscript run_platform.R
```

### Option 3: Legacy Script
```r
# Use the original script (still functional)
source("run_integrated_camk2d_platform.R")
```

## Configuration

Edit `config.yml` to customize behavior:

```yaml
# Basic settings most users need to modify
research:
  focus_area: "both"        # "aFIB", "HF", or "both"
  
datasets:
  max_datasets: 12          # Reduce for faster testing
  min_samples: 10           # Increase for higher quality
  
expression:
  validation_enabled: true  # Set false to skip downloads
```

## Step-by-Step Execution

If you want to run individual components:

```r
# Step 1: Literature Analysis
rmarkdown::render("01_literature_processing.Rmd")

# Step 2: Dataset Discovery  
rmarkdown::render("02_cross_species_discovery.Rmd")

# Step 3: Expression Validation
rmarkdown::render("03_integrated_discovery_validation.Rmd")
```

## Troubleshooting

### Platform Won't Start
- Check R packages: `BiocManager::install(c("GEOquery", "Biobase", "tidyverse"))`
- Verify config.yml exists and is valid YAML
- Check working directory: `getwd()`

### Downloads Failing
- Check internet connection
- Increase timeout in config.yml: `timeout_seconds: 1200`
- Clear cache: `unlink("cache", recursive = TRUE)`

### Memory Issues
- Reduce max_datasets to 5 or fewer
- Close other R sessions
- Restart R session: `.rs.restartR()` in RStudio

### Slow Performance
- Enable caching: `cache_downloads: true` in config.yml
- Use cached results from previous runs
- Run with fewer datasets initially

## Expected Runtime

- **Literature Analysis**: 2-5 minutes
- **Dataset Discovery**: 5-15 minutes (depending on network)
- **Expression Validation**: 10-60 minutes (depending on datasets)
- **Total**: 15-80 minutes for complete pipeline

## Output Files

Check these locations for results:

```bash
# Key result files
output/CAMK2D_Literature_Results_*.csv
output/CAMK2D_Discovery_Results_*.csv
output/expression_data/GSE*_expression.csv

# Reports (open in browser)
output/Literature_Analysis_Report.html
output/Dataset_Discovery_Report.html
output/Integrated_Discovery_Validation_Report.html
```

## Common Workflows

### Testing New Configuration
```r
# 1. Edit config.yml
# 2. Test with small number of datasets first
config.yml -> datasets.max_datasets: 3

# 3. Run platform
source("run_platform.R")

# 4. Check outputs
list.files("output", pattern = "*.csv")
```

### Production Run
```r
# 1. Use default config.yml settings
# 2. Ensure adequate disk space (~1GB)
# 3. Run full platform
source("run_platform.R")

# 4. Wait for completion (30-60 minutes)
# 5. Check all outputs generated
```

### Reprocessing Data
```r
# Option 1: Use existing discovery results
# (Faster - skips dataset discovery)
source("run_platform.R")

# Option 2: Fresh discovery
# Clear outputs first
unlink("output", recursive = TRUE)
source("run_platform.R")
```

## Getting Help

1. **Read error messages carefully** - they usually indicate the specific issue
2. **Check log output** - shows progress and where failures occur  
3. **Verify prerequisites** - R packages, internet connection, disk space
4. **Start small** - test with max_datasets: 3 before full runs
5. **Check existing results** - platform uses cached data when available

## Advanced Usage

### Custom Parameters
```r
# Run with specific focus
rmarkdown::render("03_integrated_discovery_validation.Rmd", 
  params = list(
    focus_area = "aFIB",
    max_datasets = 5,
    min_camk2d_expression = 1.5
  )
)
```

### Programmatic Access
```r
# Load results for analysis
discovery_results <- read.csv("output/CAMK2D_Discovery_Results_20250812.csv")
literature_results <- read.csv("output/CAMK2D_Literature_Results_20250812.csv")

# Access expression data
expr_files <- list.files("output/expression_data", pattern = "_expression.csv", full.names = TRUE)
first_dataset <- read.csv(expr_files[1])
```

---

**For detailed information, see README.md**