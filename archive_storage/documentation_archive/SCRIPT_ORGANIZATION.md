# CAMK Analysis Scripts Organization

## Current Script Status After Cleanup

### ✅ **PRODUCTION READY** (Methodologically Corrected)

1. **`camk_focused_analysis.R`** - Core CAMK family analysis (CORRECTED)
   - ✅ Genome-wide DGE → CAMK filtering
   - ✅ Includes pathway context genes
   - ✅ Proper statistical foundation

2. **`camk_healthy_vs_disease_analysis.R`** - Disease-specific analysis (CORRECTED)
   - ✅ Genome-wide DGE → CAMK filtering
   - ✅ Focus on healthy vs disease comparisons
   - ✅ Enhanced biological context

3. **`camk2d_independent_analysis.R`** - CAMK2D specialized analysis (VERIFIED CORRECT)
   - ✅ Uses `perform_limma_analysis()` which does genome-wide DGE
   - ✅ Post-DGE CAMK2D extraction
   - ✅ Comprehensive co-expression and functional analysis

### ⚠️ **NEEDS REVIEW/CONSOLIDATION**

4. **`camk_correlation_analysis.R`** - Correlation-focused analysis
   - ❓ Status unknown - needs review
   - 🔄 May contain old methodology

5. **`camk_meta_analysis.R`** - Meta-analysis across datasets
   - ❓ Status unknown - needs review
   - 🔄 May depend on corrected scripts above

6. **`camk_family_interconnection_analysis.R`** - Family interconnection analysis
   - ❓ Status unknown - needs review
   - 🔄 May contain redundant functionality

### 🗂️ **SUPPORT/UTILITY SCRIPTS**

7. **`comprehensive_camk_analysis_report.R`** - Reporting functions
8. **`corrected_camk_analysis_summary.R`** - Summary reporting
9. **`create_camk_mapping.R`** - Gene mapping utilities
10. **`methodology_correction_clinical_interpretation.R`** - Clinical interpretation

### 📊 **GENERAL ANALYSIS SCRIPTS**

11. **`comprehensive_dataset_inventory.R`** - Dataset management
12. **`comprehensive_dataset_investigation.R`** - Dataset exploration
13. **`multi_technology_camk_analysis.R`** - Multi-platform analysis
14. **`enhanced_pipeline_integration.R`** - Pipeline integration

## Consolidation Plan

### Phase 1: Keep Core Production Scripts (3 scripts)
- `camk_focused_analysis.R` (corrected)
- `camk_healthy_vs_disease_analysis.R` (corrected)
- `camk2d_independent_analysis.R` (verified correct)

### Phase 2: Review & Fix/Archive (3 scripts)
- Review `camk_correlation_analysis.R`
- Review `camk_meta_analysis.R`
- Review `camk_family_interconnection_analysis.R`

### Phase 3: Organize Support Scripts
- Move reporting scripts to `functions/reporting.R`
- Move utility scripts to `functions/utilities.R`

## Recommended File Structure

```
analysis_scripts/
├── core/
│   ├── camk_focused_analysis.R           # Main CAMK family analysis
│   ├── camk_healthy_vs_disease_analysis.R # Disease-focused analysis
│   └── camk2d_independent_analysis.R      # CAMK2D specialized analysis
├── specialized/
│   ├── camk_correlation_analysis.R       # Correlation analysis (if needed)
│   ├── camk_meta_analysis.R              # Meta-analysis (if needed)
│   └── camk_family_interconnection_analysis.R # Interconnection (if needed)
└── utilities/
    ├── create_camk_mapping.R              # Gene mapping
    ├── comprehensive_camk_analysis_report.R # Reporting
    └── methodology_correction_clinical_interpretation.R # Clinical
```

---
*Status: In Progress*
*Next: Review specialized scripts for methodology compliance*