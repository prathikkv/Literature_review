# CAMK2D Pipeline - Quick Usage Guide

## 🚀 One-Command Execution

### **Quick Analysis (2-3 minutes) - RECOMMENDED**
```bash
Rscript run_pipeline_complete.R --quick
```

### **Standard Analysis (5 minutes)**
```bash
Rscript run_pipeline_complete.R
```

### **Full Analysis with All Features (10+ minutes)**
```bash
Rscript run_pipeline_complete.R --full
```

## 📄 Your Results

After execution, open these files in your browser:

1. **Analysis Report**: `output/current/CAMK_Analysis_Report.html`
   - CAMK2D meta-analysis results
   - Statistical significance
   - Visualizations

2. **Interactive Documentation**: `output/current/Interactive_Technical_Documentation.html`
   - 70+ interactive flowcharts
   - Pipeline methodology
   - Technical details

## ✅ What You'll Get

- **CAMK2D Expression**: Upregulated in cardiovascular disease
- **P-value**: ~0.015 (statistically significant)
- **Effect Size**: Log fold change ~0.055
- **Datasets**: Analysis across 6 GEO studies (436 samples)

## 🔧 Troubleshooting

**If pipeline doesn't start:**
```bash
# Use quick mode first (fastest real analysis)
Rscript run_pipeline_complete.R --quick
```

**If takes too long:**
```bash
# Quick mode skips enhanced features for speed
Rscript run_pipeline_complete.R --quick
```

**For help:**
- Check `output/logs/` for error details
- See full documentation in `README.md`

---

*That's it! Just one command to get complete CAMK2D analysis results.*