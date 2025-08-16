# CAMK2D Analysis Pipeline - Production Ready

## 🎯 Overview

A comprehensive bioinformatics pipeline for CAMK2D cardiovascular meta-analysis, featuring integrated HTML reporting with interactive documentation. This production-ready pipeline analyzes differential gene expression across multiple GEO datasets with a focus on CAMK2D and related genes in cardiovascular disease.

## ✨ Key Features

- **Single Command Execution** - Multiple speed modes from 30 seconds to 10 minutes
- **Dual HTML Output** - Analysis report + Interactive technical documentation  
- **70+ Interactive Flowcharts** - Mermaid.js powered methodology visualization
- **Production Optimized** - Clean folder structure with only essential files
- **Multiple Execution Modes** - Demo, Quick, Standard, and Full analysis options
- **Real-time Progress Tracking** - Visual progress bars with step indicators

## 🚀 Quick Start

### **Fastest Test (30 seconds)**
```bash
Rscript run_pipeline_complete.R --demo
```

### **Recommended Production Run (2-3 minutes)**
```bash
Rscript run_pipeline_complete.R --quick
```

## 📊 Execution Modes

| Mode | Command | Time | Description | Use Case |
|------|---------|------|-------------|----------|
| **Demo** | `--demo` | 30 sec | Sample data, all steps | Testing/Validation |
| **Quick** | `--quick` | 2-3 min | Core analysis only | Production |
| **Standard** | *(default)* | 5 min | Core + essential features | Balanced |
| **Full** | `--full` | 10+ min | All features enabled | Research |

### Mode Details

#### 🎬 Demo Mode
```bash
Rscript run_pipeline_complete.R --demo
```
- Uses subset of data for rapid execution
- Shows complete workflow with progress bars
- Generates both HTML reports
- Perfect for testing pipeline functionality

#### ⚡ Quick Mode (RECOMMENDED)
```bash
Rscript run_pipeline_complete.R --quick
```
- Full analysis on all 6 datasets (436 samples)
- Generates real CAMK2D meta-analysis results
- Skips time-consuming enhanced features
- Optimal for routine analysis

#### 📊 Standard Mode
```bash
Rscript run_pipeline_complete.R
```
- Complete core analysis
- Includes gene family discovery (limited)
- Balanced feature set
- Good for comprehensive reports

#### 🔥 Full Mode
```bash
Rscript run_pipeline_complete.R --full
```
- All features enabled
- Complete gene family discovery
- Dataset discovery and pathway analysis
- Maximum insight generation

## 📁 Pipeline Structure

```
pipeline/
├── run_pipeline_complete.R         # Main runner with progress tracking (RECOMMENDED)
├── run_enhanced_pipeline.R         # Full-featured runner with all modules
├── run_quick_pipeline.R           # Streamlined quick execution
├── generate_interactive_documentation.R  # Documentation generator
├── config.yml                      # Pipeline configuration
├── scripts/
│   ├── pipeline_orchestrator.R     # Core orchestration logic
│   ├── step_01_data_loader.R      # Data loading step
│   ├── step_02_preprocessing.R    # Preprocessing step
│   ├── step_03_dge_analysis.R     # Differential expression
│   ├── step_04_meta_analysis.R    # Meta-analysis
│   ├── step_05_report_generator.R # Report generation
│   └── utilities/                  # Helper functions
├── modules/                        # Enhancement modules
│   ├── auto_download.R            # GEO dataset downloader
│   ├── dataset_discovery.R        # Dataset discovery
│   ├── gene_family_discovery.R    # Gene family analysis
│   ├── pathway_analysis.R         # Pathway enrichment
│   └── validation_framework.R     # Validation utilities
├── templates/
│   └── CAMK_Analysis_Report.Rmd   # Analysis report template
├── cache/                          # Cached GEO datasets
├── output/
│   ├── current/                   # Generated reports location
│   ├── logs/                      # Execution logs
│   └── checkpoints/               # Pipeline checkpoints
└── validation/
    └── baseline_results.yml       # Expected results for validation

## 📄 Output Files

After execution, find your results in `output/current/`:

1. **CAMK_Analysis_Report.html** - Main analysis report containing:
   - CAMK2D differential expression results
   - Meta-analysis statistics
   - Forest plots and visualizations
   - Dataset summaries
   - Cross-linked navigation to documentation

2. **Interactive_Technical_Documentation.html** - Technical documentation with:
   - 70+ interactive Mermaid flowcharts
   - Pipeline methodology
   - Algorithm details
   - Live analysis results summary
   - Cross-linked navigation to analysis report

## 📊 What You'll See During Execution

```
╔══════════════════════════════════════════════════════════════╗
║        CAMK2D PIPELINE - COMPLETE EXECUTION                  ║
║        Mode: QUICK | Estimated Time: 2-3 minutes             ║
╚══════════════════════════════════════════════════════════════╝

[1/8] ⏳ Loading pipeline utilities...
[▓▓▓░░░░░░░░░░░░░░░░░] 13%

[2/8] ✅ Loading cached datasets
[▓▓▓▓▓░░░░░░░░░░░░░░░] 25%

[3/8] ✅ Preprocessing completed
[▓▓▓▓▓▓▓▓░░░░░░░░░░░] 38%

... continues with visual progress ...

═══════════════════════════════════════════════════════════════
📊 PIPELINE EXECUTION COMPLETE
═══════════════════════════════════════════════════════════════
✅ Analysis Report: Generated (245 KB)
   📄 output/current/CAMK_Analysis_Report.html
✅ Interactive Documentation: Generated (148 KB)
   🌐 output/current/Interactive_Technical_Documentation.html

📊 Execution Summary:
   Mode: QUICK
   Total Steps: 8/8
   Time Elapsed: 2.4 minutes
```

## 🔧 Configuration

The pipeline behavior is controlled by `config.yml`:

- **Primary gene**: CAMK2D
- **Gene family**: 11 CAMK genes
- **Datasets**: 6 GEO datasets (GSE57338, GSE41177, GSE79768, GSE115574, GSE31821, GSE14975)
- **Total samples**: 436 (218 disease, 218 control)
- **Diseases**: Heart Failure, Atrial Fibrillation
- **Analysis method**: limma + metafor (fixed-effects)

## 📋 Requirements

- R version 4.0+
- Required packages (auto-checked on startup):
  - tidyverse
  - limma
  - metafor
  - rmarkdown
  - knitr
  - yaml

## 🐛 Troubleshooting

### Pipeline stops at validation
- Check that all required packages are installed
- Verify cache/ directory contains dataset files
- Run `Rscript scripts/utilities/config_validator.R` to validate configuration

### Mermaid diagrams not rendering
- The pipeline automatically handles large diagrams with 200K text limit
- Oversized diagrams show graceful error messages
- All 70+ flowcharts tested and working

### Analysis report not generated
- Use `--quick` mode for faster execution
- Check `output/logs/` for error details
- Ensure RMarkdown is properly installed

### Demo mode for quick testing
If experiencing issues, test with demo mode first:
```bash
Rscript run_pipeline_complete.R --demo
```
This completes in 30 seconds and validates the entire workflow.

## 📊 Expected Results

For CAMK2D in cardiovascular disease:
- **Combined log fold change**: ~0.0552
- **P-value**: ~0.015 (significant)
- **Direction**: Upregulated
- **Consistency**: Across multiple datasets

## 🔬 Scientific Context

This pipeline performs meta-analysis of CAMK2D (Calcium/Calmodulin Dependent Protein Kinase II Delta) expression in cardiovascular disease, combining data from multiple studies to identify consistent expression patterns that may indicate therapeutic targets.

## 📚 Citation

If using this pipeline, please cite:
- Original GEO datasets (GSE IDs listed in reports)
- limma: Ritchie et al., 2015
- metafor: Viechtbauer, 2010

## 💡 Tips for Best Performance

1. **First time**: Run demo mode to verify setup
2. **Regular use**: Use quick mode (2-3 min)
3. **Research**: Use standard or full mode for comprehensive analysis
4. **Debugging**: Check output/logs/ for detailed execution logs
5. **Validation**: Compare with validation/baseline_results.yml

## 🆘 Support

- Check logs in `output/logs/` for execution details
- Review `validation/baseline_results.yml` for expected outputs
- Ensure all R packages are up to date

---

*Pipeline Version: 1.0.0 | Production Ready*