# Archived Redundant Analysis Scripts

These scripts have been archived during the refactoring process as they depend on old output files and methodologies that have been corrected.

## Archived Scripts

### camk_correlation_analysis.R
- **Issue**: Depends on `output/CAMK_focused_analysis_results.rds` which has been replaced by `output/CAMK_focused_analysis_results_corrected.rds`
- **Status**: Functionality can be integrated into main analysis scripts when needed
- **Dependencies**: corrplot, pheatmap libraries

### camk_meta_analysis.R  
- **Issue**: Depends on `output/CAMK_focused_DGE_all_datasets.csv` which has been replaced by `output/CAMK_core_DGE_all_datasets.csv`
- **Status**: Functionality can be integrated into main analysis scripts when needed
- **Dependencies**: metafor, meta libraries

## To Reactivate

If these analysis types are needed:

1. Update file paths to use corrected output files
2. Integrate with centralized modules:
   - `source("functions/camk_definitions.R")`
   - `source("functions/common_utils.R")`
3. Use `load_required_libraries()` function instead of individual library calls
4. Use `save_analysis_results()` for standardized output

## Methodological Notes

These scripts were part of the old methodology that had statistical bias from gene pre-filtering. The corrected methodology performs genome-wide DGE first, then filters for CAMK genes, providing proper statistical foundation.