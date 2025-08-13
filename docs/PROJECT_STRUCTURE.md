# CAMK Analysis Project Structure

## Overview
This project has been reorganized for better maintainability and clarity. All files are now logically grouped by their purpose.

## Directory Structure

```
├── analysis/                    # Main analysis scripts (organized by purpose)
│   ├── core/                   # Core CAMK analysis scripts
│   │   ├── camk_focused_analysis.R
│   │   ├── camk_healthy_vs_disease_analysis.R
│   │   └── camk2d_independent_analysis.R
│   ├── exploratory/           # Investigation and exploration scripts
│   │   ├── comprehensive_dataset_investigation.R
│   │   ├── comprehensive_dataset_inventory.R
│   │   ├── camk_family_interconnection_analysis.R
│   │   ├── multi_technology_camk_analysis.R
│   │   ├── corrected_camk_analysis_summary.R
│   │   └── enhanced_pipeline_integration.R
│   ├── visualization/         # Plotting and visualization scripts
│   │   ├── create_final_visualizations.R
│   │   ├── comprehensive_camk_analysis_report.R
│   │   ├── final_methodology_comparison_report.R
│   │   └── methodology_correction_clinical_interpretation.R
│   └── validation/           # Verification and testing scripts
│       ├── verify_dataset_comparisons.R
│       └── validate.R
├── scripts/                   # Utility and setup scripts
│   ├── setup.R
│   ├── run_pipeline.R
│   ├── enhanced_group_detection.R
│   └── create_camk_mapping.R
├── docs/                     # All documentation files
│   ├── PROJECT_STRUCTURE.md (this file)
│   ├── README.md
│   ├── METHODOLOGY_CORRECTION_SUMMARY.md
│   ├── REFACTORING_SUMMARY.md
│   └── [other .md files]
├── notebooks/                # R Markdown files and generated content
│   ├── CAMK2D_Analysis_Documentation.Rmd
│   ├── CAMK2D_Quick_Tutorial.Rmd
│   ├── CAMK2D_Analysis_Documentation_cache/
│   └── CAMK2D_Analysis_Documentation_files/
├── functions/                # Reusable function modules
│   ├── camk_definitions.R    # Centralized CAMK gene definitions
│   ├── common_utils.R        # Shared utility functions
│   ├── data_processing.R     # Data processing functions
│   ├── analysis.R            # Analysis functions
│   ├── visualization.R       # Visualization functions
│   └── utilities.R           # General utilities
├── archive/                  # Archived/backup files
│   ├── README.md
│   ├── README_archived_scripts.md
│   ├── CAMK2D_Analysis_Documentation_backup_20250813_153343.Rmd
│   ├── camk_correlation_analysis.R
│   └── camk_meta_analysis.R
├── cache/                    # Cached data files
├── data/                     # Input data files
├── output/                   # Analysis outputs
├── results/                  # Final results
├── config.yml               # Project configuration
├── renv.lock                # R package versions
└── CAM.Rproj                # RStudio project file
```

## Usage Guidelines

### Running Core Analyses
```r
# Change to project root directory first
setwd("/path/to/project/root")

# Run core CAMK analysis
source("analysis/core/camk_focused_analysis.R")

# Run healthy vs disease specific analysis
source("analysis/core/camk_healthy_vs_disease_analysis.R")

# Run CAMK2D independent analysis
source("analysis/core/camk2d_independent_analysis.R")
```

### Running the Pipeline
```r
# From project root
source("scripts/run_pipeline.R")
```

### Path References
All scripts have been updated to use relative paths from their new locations:
- Scripts in `analysis/core/`, `analysis/exploratory/`, `analysis/visualization/`, `analysis/validation/` reference functions with `../../functions/`
- Scripts in `scripts/` reference functions with `../functions/`
- All inter-script references have been updated accordingly

## Key Benefits of New Structure

1. **Clear Separation of Concerns**: Each directory has a specific purpose
2. **Easy Navigation**: Find functionality quickly by category
3. **Maintainable**: Changes to functions affect all dependent scripts consistently
4. **Professional**: Industry-standard project organization
5. **Scalable**: Easy to add new analyses in appropriate categories

## Migration Notes

- All `source()` paths have been updated to work from new locations
- Functions are centralized in `functions/` directory
- Documentation is consolidated in `docs/`
- Archives are unified in single `archive/` directory
- Core functionality remains unchanged - only organization improved

## File Dependencies

### Core Analysis Dependencies
- All core analysis scripts depend on: `functions/camk_definitions.R`, `functions/common_utils.R`
- Path detection uses: `scripts/enhanced_group_detection.R`

### Function Module Dependencies
- `functions/camk_definitions.R`: Standalone (no dependencies)
- `functions/common_utils.R`: Standalone (no dependencies)  
- Other function modules may have interdependencies

## Validation

After reorganization:
- ✅ All core functionality preserved
- ✅ All file paths updated correctly
- ✅ Function dependencies maintained
- ✅ Analysis outputs unchanged
- ✅ Documentation consolidated