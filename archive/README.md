# Archive: Old Methodology Files

This directory contains files that implemented the **old, flawed methodology** where CAMK family genes were filtered **before** differential expression analysis.

## Methodological Issue

The old approach:
1. Load expression matrix
2. **Filter to only CAMK genes** (11 genes)
3. Run DGE analysis on filtered data

This created **statistical bias** because:
- Limited statistical power (11 genes vs full transcriptome)
- Missing biological context (no pathway background)
- Improper normalization (methods work better with larger gene sets)
- No background correction for technical effects

## Corrected Approach

The new methodology:
1. Load full expression matrix (all genes)
2. **Run genome-wide DGE analysis**
3. Apply proper statistical correction with full gene background
4. **Filter CAMK genes from results** (post-DGE filtering)

## Files Archived

- `CAMK2D_Analysis_Documentation_backup_20250813_153343.Rmd` - Contains examples of pre-filtering approach

## Status

❌ **DO NOT USE** these files for analysis - they contain methodological flaws
✅ **Use corrected analysis scripts** in the main directory

---
*Archived on: {Sys.Date()}*
*Reason: Methodological correction to eliminate data skewing*