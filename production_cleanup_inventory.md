# Production Cleanup File Inventory

Generated: 2025-08-14
Purpose: Complete catalog of all files for production-level cleanup

## 🎯 CORE PRODUCTION FILES (MUST KEEP)

### Main Pipeline Scripts
- `/scripts/core/comprehensive_6_dataset_pipeline.R` - **MAIN ANALYSIS PIPELINE** 
  - Status: ⚠️ HAS PATH BUG (line 12: wrong path to enhanced_group_detection_corrected.R)
  - Dependencies: functions/camk_definitions.R, functions/analysis.R, enhanced_group_detection_corrected.R
  - Purpose: Processes all 6 datasets, generates DGE results

- `/scripts/core/enhanced_group_detection_corrected.R` - **GROUP DETECTION MODULE**
  - Status: ✅ REQUIRED by main pipeline
  - Dependencies: None identified
  - Purpose: Detects disease/control groups in datasets

- `/scripts/core/fixed_meta_analysis.R` - **META-ANALYSIS MODULE**
  - Status: ✅ REQUIRED
  - Dependencies: metafor, tidyverse
  - Purpose: Performs fixed-effects meta-analysis

- `/scripts/core/comprehensive_meta_analysis.R` - **ADVANCED META-ANALYSIS**
  - Status: ❓ UNCLEAR if used in current pipeline
  - Dependencies: metafor, tidyverse
  - Purpose: Alternative meta-analysis approach

### Utility Scripts
- `/scripts/utilities/setup.R` - **PACKAGE INSTALLER**
  - Status: ✅ ESSENTIAL for environment setup
  - Dependencies: None (installs all dependencies)
  - Purpose: Installs all required CRAN and Bioconductor packages

- `/scripts/utilities/run_pipeline.R` - **PIPELINE RUNNER**
  - Status: ❓ UNCLEAR if currently functional
  - Dependencies: yaml, all functions/*.R files
  - Purpose: Orchestrates complete pipeline execution

- `/scripts/utilities/validate_structure.R` - **VALIDATION UTILITY**
  - Status: ❓ UNCLEAR if used
  - Dependencies: Unknown
  - Purpose: Unknown validation functions

### Core Functions (ALL REQUIRED)
- `/functions/camk_definitions.R` - **GENE DEFINITIONS**
  - Status: ✅ REFERENCED by main pipeline
  - Purpose: Defines CAMK gene families and categories

- `/functions/analysis.R` - **CORE ANALYSIS FUNCTIONS**
  - Status: ✅ REFERENCED by main pipeline  
  - Purpose: Core differential expression analysis functions

- `/functions/data_processing.R` - **DATA PROCESSING**
  - Status: ✅ LIKELY REQUIRED
  - Purpose: Dataset download and preprocessing functions

- `/functions/multimodal_download.R` - **DOWNLOAD MODULE**
  - Status: ✅ LIKELY REQUIRED
  - Purpose: GEO/ArrayExpress data download

- `/functions/common_utils.R` - **UTILITY FUNCTIONS**
  - Status: ❓ UNCLEAR dependency status
  - Purpose: Common utility functions

- `/functions/dataset_technical_specs.R` - **DATASET SPECS**
  - Status: ❓ UNCLEAR dependency status
  - Purpose: Technical specifications for datasets

- `/functions/load_pipeline_results.R` - **RESULT LOADER**
  - Status: ❓ UNCLEAR dependency status
  - Purpose: Pipeline result loading functions

- `/functions/ml_biomarker_prediction.R` - **ML ANALYSIS**
  - Status: ❓ UNCLEAR if used in current pipeline
  - Purpose: Machine learning biomarker prediction

- `/functions/pathway_analysis.R` - **PATHWAY ANALYSIS**
  - Status: ❓ UNCLEAR if used in current pipeline
  - Purpose: Pathway enrichment analysis

- `/functions/survival_analysis.R` - **SURVIVAL ANALYSIS**
  - Status: ❓ UNCLEAR if used in current pipeline
  - Purpose: Survival analysis functions

- `/functions/utilities.R` - **GENERAL UTILITIES**
  - Status: ❓ UNCLEAR dependency status
  - Purpose: General utility functions

- `/functions/visualization.R` - **VISUALIZATION**
  - Status: ❓ UNCLEAR if used in current pipeline
  - Purpose: Plotting and visualization functions

### Essential Data Files
- `/cache/microarray/GSE57338_processed.rds` - **PROCESSED DATASET**
  - Status: ✅ REQUIRED (313 samples, heart failure)
  - Size: Large binary file
  - Purpose: Primary heart failure dataset

- `/cache/microarray/GSE41177_processed.rds` - **PROCESSED DATASET**
  - Status: ✅ REQUIRED (38 samples, atrial fibrillation)
  - Size: Large binary file  
  - Purpose: AF dataset with sample imbalance

- `/cache/microarray/GSE79768_processed.rds` - **PROCESSED DATASET**
  - Status: ✅ REQUIRED (26 samples, atrial fibrillation)
  - Size: Large binary file
  - Purpose: AF dataset for meta-analysis

- `/cache/microarray/GSE115574_processed.rds` - **PROCESSED DATASET**
  - Status: ✅ REQUIRED (59 samples, atrial fibrillation)
  - Size: Large binary file
  - Purpose: Largest AF dataset

### Current Results (REQUIRED BY RMD REPORT)
- `/output/current/CAMK_meta_analysis_FINAL.csv` - **FINAL META-ANALYSIS**
  - Status: ✅ REQUIRED by RMD report (line 78)
  - Purpose: Meta-analysis results for all CAMK genes

- `/output/current/CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv` - **COMPREHENSIVE DGE**
  - Status: ✅ REQUIRED by RMD report (line 81)
  - Purpose: Individual dataset results for all CAMK genes

- `/output/current/dataset_processing_summary_6_datasets.csv` - **DATASET SUMMARY**
  - Status: ✅ REQUIRED by RMD report (line 84)
  - Purpose: Dataset processing summary and statistics

- `/output/current/methodology_comparison_analysis.csv` - **METHODOLOGY ANALYSIS**
  - Status: ✅ REQUIRED by RMD report (line 252)
  - Purpose: Methodology validation results

- `/output/current/comprehensive_dataset_inventory.csv` - **DATASET INVENTORY**
  - Status: ❓ UNCLEAR if used by RMD
  - Purpose: Complete dataset inventory

- `/output/current/detailed_methodology_analysis.md` - **METHODOLOGY DOCUMENTATION**
  - Status: ❓ UNCLEAR if used by RMD
  - Purpose: Detailed methodology analysis

### Main Report
- `/reports/CAMK_Analysis_Professional_Report.Rmd` - **MAIN REPORT**
  - Status: ✅ PRIMARY OUTPUT DOCUMENT
  - Dependencies: All files in /output/current/
  - Purpose: Generates final HTML report

### Reference Data
- `/data/comprehensive_camk2d_literature.csv` - **LITERATURE DATA**
  - Status: ❓ UNCLEAR if used
  - Purpose: CAMK2D literature references

- `/data/verified_camk2d_literature.csv` - **VERIFIED LITERATURE**
  - Status: ❓ UNCLEAR if used
  - Purpose: Verified literature references

- `/data/ortholog_mappings/` - **ORTHOLOG DATA**
  - Status: ❓ UNCLEAR if used
  - Purpose: Cross-species gene mappings

### Configuration Files
- `/claude.md` - **CLAUDE INSTRUCTIONS**
  - Status: ✅ KEEP (user instructions)
  - Purpose: Project instructions and rules

- `/config.yml` - **CONFIGURATION**
  - Status: ❓ UNCLEAR if used by pipeline
  - Purpose: Pipeline configuration

- `/renv.lock` - **R ENVIRONMENT**
  - Status: ✅ KEEP (package management)
  - Purpose: R package version control

## 🗑️ CLEANUP CANDIDATES (MOVE TO ARCHIVE)

### Archive Directory (ENTIRE FOLDER)
- `/archive/` - **HISTORICAL DEVELOPMENT CODE**
  - Status: 🗑️ MOVE TO ARCHIVE_STORAGE
  - Contains: 35+ experimental scripts, old analysis attempts
  - Size: Large (multiple subdirectories)
  - Purpose: Historical development, not needed for production

### Old Results (MOVE TO ARCHIVE)
- `/output/archive/` - **OLD RESULTS**
  - Status: 🗑️ MOVE TO ARCHIVE_STORAGE
  - Contains: Legacy CSV files, old visualizations
  - Purpose: Historical results, replaced by current/

### Legacy Reports (MOVE TO ARCHIVE)  
- `/reports/CAMK_Enhanced_Professional_Report.html` - **OLD REPORT**
  - Status: 🗑️ MOVE TO ARCHIVE_STORAGE
  - Purpose: Superseded by Professional_Report

- `/reports/CAMK_Professional_Analysis_Report.html` - **OLD REPORT**
  - Status: 🗑️ MOVE TO ARCHIVE_STORAGE  
  - Purpose: Superseded by Professional_Report

- `/reports/CAMK_Scientifically_Validated_Report.html` - **OLD REPORT**
  - Status: 🗑️ MOVE TO ARCHIVE_STORAGE
  - Purpose: Superseded by Professional_Report

- `/reports/archive/` - **ARCHIVED REPORTS**
  - Status: 🗑️ MOVE TO ARCHIVE_STORAGE
  - Purpose: Already archived reports

### Unused Cache Files
- `/cache/comprehensive/` - **ALTERNATIVE CACHE**
  - Status: 🗑️ LIKELY UNUSED (only 1 file)
  - Purpose: Alternative processing approach

- `/cache/intelligent/` - **INTELLIGENT CACHE**  
  - Status: 🗑️ LIKELY UNUSED (batch correction experiments)
  - Purpose: Experimental batch correction

- `/cache/GPL11532.soft.gz` - **PLATFORM ANNOTATION**
  - Status: ❓ UNCLEAR if needed
  - Purpose: Microarray platform annotation

### Duplicate Output Files
Multiple CSV files in `/output/` root that appear to be duplicates of files in `/output/current/`:
- `CAMK_DGE_all_6_datasets_COMPREHENSIVE.csv` - **DUPLICATE**
- `CAMK_meta_analysis_FINAL.csv` - **DUPLICATE**
- Various other analysis outputs that may be superseded

### Documentation (KEEP ESSENTIAL, ARCHIVE DETAILED)
- `/docs/` - **PROJECT DOCUMENTATION**
  - Status: ✅/🗑️ MIXED (keep essential, archive detailed development docs)
  - Keep: README.md, PROJECT_STRUCTURE.md, PRODUCTION_READY.md
  - Archive: Detailed development summaries

### Notebooks
- `/notebooks/CAMK2D_Quick_Tutorial.Rmd` - **TUTORIAL**
  - Status: ❓ UNCLEAR if needed for production
  - Purpose: User tutorial

## 🔍 TESTING REQUIRED

Files marked with ❓ need dependency testing to determine if they're actually used by the current pipeline.

**Next Steps:**
1. Install packages via setup.R
2. Test pipeline end-to-end with current files
3. Monitor which files are actually loaded/used
4. Make final keep/archive decisions based on actual usage