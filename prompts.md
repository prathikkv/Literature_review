# Claude Code Prompts for CAMK2D Multi-Species Analysis Pipeline

## Prompt 1: Project Setup and Comprehensive Data Retrieval

```
Create an R project structure for a comprehensive CAMK2D transcriptomics analysis across multiple species (human, mouse, rat) investigating atrial fibrillation and heart failure. 

Requirements:
1. Create a main RMarkdown document that will serve as the analysis notebook
2. Set up proper directory structure (data/raw/, data/processed/, results/, figures/, reports/, logs/)
3. Install and load all required packages for transcriptomics analysis including:
   - GEOquery for GEO data retrieval
   - ArrayExpress for EBI data retrieval  
   - limma for differential expression
   - DESeq2 for RNA-seq analysis  
   - biomaRt for gene annotation
   - metafor for meta-analysis
   - tidyverse for data manipulation
   - ComplexHeatmap for visualization
   - clusterProfiler for pathway analysis
   - WGCNA for co-expression analysis
   - R.utils for file handling
   - httr for web downloads

4. Create comprehensive data download functions that handle:
   **Human Heart Failure Datasets:**
   - GSE120895: Download expression matrix and sample metadata
   - GSE57338: Retrieve processed data and clinical annotations
   - GSE141910: Large dataset with 366 samples (DCM + controls)
   
   **Human Atrial Fibrillation Datasets:**
   - GSE31821: AF vs SR samples with proper metadata extraction
   - GSE41177: Left atrial samples (16 AF + 3 SR)
   - GSE79768: Left/right atrial samples (7 AF + 6 SR) 
   - GSE115574: Larger AF cohort (14 AF + 15 SR)
   - GSE14975: Balanced AF dataset (5 AF + 5 SR)
   
   **Mouse/Rat Datasets:**
   - E-MTAB-7895: Mouse MI time course (days 0,1,3,7,14,28)
   - GSE132146: Col1a1-GFP+ cardiac fibroblasts
   - GSE155882: TAC-induced heart failure model

5. Implement robust download functions that:
   - Check if data already exists locally before downloading
   - Handle different GEO data formats (SOFT, TXT, CEL files)
   - Download supplementary files when available
   - Extract and parse sample metadata properly
   - Handle missing or corrupted downloads with retry logic
   - Create detailed download logs with timestamps

6. Data validation functions:
   - Verify file integrity after download
   - Check sample sizes match expected numbers
   - Validate gene expression matrix dimensions
   - Ensure metadata consistency
   - Generate download summary reports

7. Include comprehensive error handling:
   - Network timeout handling
   - File permission checks
   - Disk space validation
   - Progress tracking for large downloads
   - Graceful failure with informative error messages

The code should be production-ready with clear comments, follow R best practices, and include example usage for each function.
```

## Prompt 2: Data Preprocessing and Quality Control

```
Build on the previous R project to implement comprehensive data preprocessing and quality control for the downloaded CAMK2D transcriptomics datasets.

Requirements:
1. Create robust data loading functions that handle different GEO formats:
   - Parse GEO SOFT files using GEOquery::getGEO()
   - Load expression matrices from downloaded supplementary files
   - Extract and clean sample metadata (phenotype data)
   - Handle different microarray platforms (GPL570, GPL96, GPL16791, etc.)
   - Process RNA-seq count matrices and TPM/FPKM data
   - Merge technical replicates when present

2. Implement platform-specific preprocessing:
   **Microarray data:**
   - Convert probe IDs to gene symbols using appropriate annotation packages
   - Handle multiple probes per gene (median/max selection)
   - Log2 transformation for non-log data
   - Quantile normalization across samples
   - Background correction using RMA/GCRMA
   
   **RNA-seq data:**
   - Filter low-count genes (CPM threshold)
   - TMM normalization using edgeR
   - Variance stabilizing transformation
   - Handle different count formats (raw, normalized, log-transformed)

3. Data integration and harmonization:
   - Standardize sample naming conventions
   - Map clinical metadata to consistent terminology
   - Merge datasets from same studies when applicable
   - Create unified sample annotation tables
   - Generate study-specific identifiers

4. Comprehensive quality control pipeline:
   - Sample correlation matrices and clustering
   - Principal Component Analysis with variance explained
   - Expression distribution analysis (density plots, boxplots)
   - Outlier detection using PCA distances and correlation
   - Platform effect visualization
   - Missing data assessment and handling

5. Cross-species gene mapping:
   - Use biomaRt to map human-mouse-rat orthologs
   - Handle one-to-many and many-to-one ortholog relationships
   - Create cross-reference tables for downstream analysis
   - Validate mapping success rates

6. Batch effect correction:
   - Detect batch effects using PCA and clustering
   - Apply ComBat correction preserving biological variables
   - Compare before/after batch correction
   - Validate that biological signal is preserved

7. Generate detailed QC reports:
   - Summary statistics for each dataset
   - Sample size and gene coverage tables
   - QC plots with clear pass/fail criteria
   - Outlier sample identification and recommendations
   - Platform comparison and integration success metrics

Include specific code to handle the exact datasets we're downloading, with proper validation that the preprocessing worked correctly for each study.
```

## Prompt 2.5: Dataset Validation and Structure Verification

```
Create a comprehensive dataset validation module to ensure all downloaded transcriptomics data is correct and properly structured before proceeding with analysis.

Requirements:
1. Implement dataset validation functions for each GEO accession:
   **For each dataset, verify:**
   - Correct number of samples matches published studies
   - Expression matrix dimensions are reasonable
   - Sample metadata includes necessary clinical variables
   - Gene identifiers are properly formatted
   - Data ranges are biologically plausible

2. Create specific validation checks for each dataset:
   **GSE120895 (Heart Failure):**
   - Verify ~160+ samples with DCM and control groups
   - Check for heart failure-related clinical variables
   - Validate platform GPL16791 (Illumina HiSeq)
   
   **GSE57338 (Heart Failure):**
   - Confirm presence of HF vs control samples
   - Validate expression data format and range
   
   **GSE41177 (Atrial Fibrillation):**
   - Verify 16 AF samples + 3 SR controls
   - Check left atrial tissue annotation
   - Validate GPL570 platform data
   
   **GSE79768 (Atrial Fibrillation):**
   - Confirm 7 AF + 6 SR samples
   - Verify atrial tissue source annotation
   
   **GSE115574, GSE14975, GSE31821:**
   - Validate sample numbers and AF/SR groups
   - Check expression data completeness

3. Data structure standardization:
   - Create consistent expression matrix format (genes x samples)
   - Standardize sample naming conventions
   - Harmonize clinical variable naming
   - Ensure consistent gene identifier formats
   - Generate unified metadata tables

4. Cross-dataset comparison:
   - Compare gene coverage across datasets
   - Identify shared and unique samples
   - Validate platform compatibility
   - Check for potential sample overlaps between studies

5. Generate validation reports:
   - Dataset summary table with key metrics
   - Sample size verification against published papers
   - Gene coverage and platform comparison
   - Data quality flags and recommendations
   - Missing data assessment

6. Error handling and troubleshooting:
   - Identify common download/parsing issues
   - Provide solutions for data format problems
   - Flag datasets that may need manual curation
   - Generate alternative download strategies if needed

This validation step ensures we have high-quality, properly formatted data before investing time in downstream analysis.
```

## Prompt 3: CAMK2D and CAMK Family Expression Analysis

```
Extend the R analysis pipeline to perform comprehensive expression analysis of CAMK2D and the entire CAMK family across all validated datasets and species.

Requirements:
1. Create gene extraction and analysis functions:
   **CAMK family gene list (use both gene symbols and aliases):**
   - CAMK2D (CaMKII-delta, Camk2d)
   - CAMK2A (CaMKII-alpha, Camk2a) 
   - CAMK2B (CaMKII-beta, Camk2b)
   - CAMK2G (CaMKII-gamma, Camk2g)
   - CAMKK1 (CaMKK-alpha, Camkk1)
   - CAMKK2 (CaMKK-beta, Camkk2)
   - Include all known aliases and Ensembl IDs

2. Robust gene expression extraction:
   - Handle cases where genes are missing from datasets
   - Deal with multiple probes per gene (select best or average)
   - Convert between gene symbols, Ensembl IDs, and probe IDs
   - Cross-validate gene presence across all downloaded datasets
   - Create gene-by-sample expression matrices for CAMK family

3. Differential expression analysis pipeline:
   **For each dataset separately:**
   - AF vs Sinus Rhythm comparisons
   - Heart Failure vs Healthy Control comparisons
   - Use appropriate statistical tests (limma for microarray, DESeq2 for RNA-seq)
   - Calculate fold changes, p-values, and FDR-corrected q-values
   - Generate effect size estimates with confidence intervals

4. Cross-species ortholog analysis:
   - Map human CAMK genes to mouse/rat orthologs using biomaRt
   - Compare expression patterns across species
   - Account for species-specific gene duplications
   - Create cross-species expression comparison tables

5. Meta-analysis implementation:
   - Combine results across studies using random effects models
   - Handle heterogeneity between studies
   - Calculate pooled effect sizes and confidence intervals
   - Generate forest plots for each CAMK family member
   - Assess publication bias and study quality effects

6. Comprehensive visualizations:
   - Volcano plots highlighting CAMK genes for each dataset
   - Box plots showing CAMK2D expression in disease vs control
   - Heatmaps of CAMK family expression across all conditions
   - Forest plots for meta-analysis results
   - Cross-species comparison plots
   - Co-expression networks within CAMK family

7. Statistical robustness checks:
   - Multiple testing correction across all genes and studies
   - Sensitivity analysis excluding individual studies
   - Heterogeneity assessment (I² and Q statistics)
   - Power calculations for effect size detection
   - Bootstrap confidence intervals for meta-analysis

8. Generate detailed results tables:
   - Individual study results for each CAMK gene
   - Meta-analysis summary statistics
   - Cross-species comparison results
   - Sample size and power analysis results

Include specific code to handle the actual gene symbols from our downloaded datasets and validate that CAMK genes are properly detected in each study.
```

## Prompt 3.5: Large-Scale Database Integration (ARCHS4, GTEx, Additional Resources)

```
Integrate large-scale transcriptomics databases to expand the CAMK2D analysis beyond individual GEO studies.

Requirements:
1. ARCHS4 database integration:
   - Install and configure ARCHS4 R package or direct API access
   - Download human and mouse gene expression matrices
   - Filter for cardiac/cardiovascular samples using metadata
   - Extract CAMK family expression across >100,000 samples
   - Identify heart failure and atrial fibrillation samples in ARCHS4

2. GTEx (Genotype-Tissue Expression) integration:
   - Access GTEx heart tissue data through recount3 or direct download
   - Extract normal heart expression patterns for CAMK genes
   - Use as additional healthy control data
   - Compare disease samples to GTEx normals

3. ArrayExpress database access:
   - Implement functions to query and download from ArrayExpress
   - Target E-MTAB-7895 (mouse MI time course) and E-MTAB-6081 (multi-tissue)
   - Handle ArrayExpress-specific data formats and metadata

4. Additional cardiovascular datasets:
   - Identify and download large cardiovascular studies from:
     * Human Protein Atlas cardiac data
     * CZI Heart Cell Atlas
     * Cardiovascular datasets in ENCODE
   - Filter and integrate relevant samples

5. Large-scale data processing:
   - Implement efficient data handling for >50,000 samples
   - Use data.table or similar for memory-efficient processing
   - Batch processing for expression calculations
   - Parallel computing for intensive analyses

6. Quality control for large datasets:
   - Automated outlier detection and removal
   - Platform and batch effect assessment
   - Sample annotation validation and standardization
   - Data completeness assessment

7. Integration with existing analysis:
   - Merge large-scale data with curated GEO studies
   - Weighted analysis accounting for study quality and size
   - Validation of findings across independent large cohorts
   - Power analysis benefits from increased sample size

8. Computational optimization:
   - Memory-efficient data storage and retrieval
   - Parallel processing implementation
   - Progress tracking for long-running analyses
   - Checkpoint saving for large computations

This integration provides orders of magnitude more data for robust CAMK2D analysis and validation.
```

## Prompt 4: Pathway Analysis and Functional Enrichment

```
Develop the pathway analysis and functional enrichment component of the CAMK2D analysis pipeline.

Requirements:
1. Implement comprehensive pathway analysis using:
   - clusterProfiler for GO, KEGG, Reactome enrichment
   - GSEA (Gene Set Enrichment Analysis) for pathway scores
   - Custom pathway analysis for calcium signaling and cardiac function
   - Cross-species pathway conservation analysis

2. Create functions for:
   - Identifying CAMK2D co-expressed genes in each condition
   - Functional enrichment of differentially expressed genes
   - Pathway visualization and interpretation
   - Integration of results across species

3. CAMK2D-specific analysis:
   - Literature mining integration for known CAMK2D targets
   - Prediction of novel phosphorylation substrates
   - Cardiac-specific expression analysis
   - Secretion potential assessment using SignalP predictions

4. Visualization functions:
   - Enrichment plots and dotplots
   - Network diagrams showing CAMK2D interactions
   - Pathway heatmaps across conditions and species
   - GSEA enrichment plots

5. Generate comprehensive reports:
   - Top enriched pathways in AF and HF
   - CAMK2D target gene analysis
   - Cross-species pathway conservation
   - Potential drug targets and biomarkers

Include detailed biological interpretation and clinical relevance of findings.
```

## Prompt 5: Meta-Analysis and Statistical Integration

```
Implement robust meta-analysis methods to combine results across multiple studies and species for the CAMK2D analysis.

Requirements:
1. Create meta-analysis functions using the metafor package:
   - Random effects models for combining effect sizes
   - Heterogeneity assessment and forest plots
   - Subgroup analysis by species, disease type, study design
   - Publication bias assessment using funnel plots and Egger's test

2. Statistical integration methods:
   - Rank-based meta-analysis for combining p-values
   - Effect size standardization across studies
   - Confidence interval calculation and visualization
   - Sensitivity analysis excluding low-quality studies

3. Cross-species integration:
   - Comparative analysis of CAMK2D expression patterns
   - Species-specific vs conserved changes
   - Evolutionary analysis of CAMK2D function
   - Translational relevance scoring

4. Quality assessment:
   - Study quality scoring based on sample size, design, methodology
   - Risk of bias assessment
   - Heterogeneity sources identification
   - Robustness testing

5. Advanced visualizations:
   - Multi-panel forest plots
   - Cross-species comparison plots
   - Meta-analysis summary dashboards
   - Interactive plots using plotly

6. Generate final summary tables and reports:
   - Meta-analysis results table
   - Study characteristics summary
   - Cross-species comparison results
   - Clinical relevance assessment

Ensure all statistical methods are properly documented with citations and biological interpretation.
```

## Prompt 6: Phosphoproteomics and Drug Target Analysis

```
Develop the phosphoproteomics analysis and drug target identification component of the CAMK2D pipeline.

Requirements:
1. Implement phosphorylation target analysis:
   - Literature mining for known CAMK2D substrates
   - Cross-reference with cardiac-specific proteins
   - Analyze substrate expression in AF/HF datasets
   - Predict novel phosphorylation sites using available tools

2. Drug target analysis:
   - Druggability assessment of CAMK2D and substrates
   - Integration with drug databases (DrugBank, ChEMBL)
   - Identification of existing compounds targeting CAMK2D pathway
   - Repurposing opportunities analysis

3. Create analysis functions for:
   - Phosphopeptide sequence extraction and analysis
   - Tryptic digest prediction for mass spectrometry
   - Cardiac-specific expression filtering
   - Secretion potential analysis (signal peptides, membrane proteins)

4. Integration with expression data:
   - Correlation of substrate expression with CAMK2D
   - Co-expression network analysis of CAMK2D-substrate pairs
   - Pathway impact assessment

5. Visualization and reporting:
   - CAMK2D substrate network diagrams
   - Expression correlation plots
   - Drug target priority ranking
   - Phosphorylation site mapping

6. Generate deliverable tables:
   - Comprehensive CAMK2D substrate list with evidence levels
   - Phosphopeptide sequences for mass spectrometry
   - Drug target priority list with druggability scores
   - Biomarker candidate ranking

Include detailed methodology for phosphoproteomics integration and clinical translation potential.
```

## Prompt 7: Final Report Generation and Data Packaging

```
Create the final reporting system and comprehensive data packaging for the CAMK2D analysis with verification of all data sources.

Requirements:
1. Generate a comprehensive HTML report using RMarkdown that includes:
   **Executive Summary:**
   - Key findings for CAMK2D across species and diseases
   - Statistical significance and effect sizes
   - Clinical relevance and translational implications
   - Main conclusions and recommendations

   **Data Sources Verification:**
   - Table of all successfully downloaded datasets with sample sizes
   - Data quality metrics and processing success rates
   - Cross-validation of results across independent studies
   - Coverage analysis for CAMK family genes across datasets

   **Detailed Results Sections:**
   - Individual study results with forest plots
   - Meta-analysis findings with heterogeneity assessment
   - Cross-species comparison and evolutionary insights
   - Pathway analysis and functional enrichment
   - Drug target and biomarker identification

2. Create comprehensive supplementary materials:
   **Excel workbooks with multiple sheets:**
   - Dataset summary with download details and QC metrics
   - CAMK family expression results across all studies
   - Meta-analysis results with confidence intervals
   - Cross-species ortholog mapping and conservation
   - Pathway enrichment results with gene lists
   - Drug target candidates with druggability scores
   - Phosphorylation target analysis results

3. Data verification and validation:
   - Cross-check results against published literature
   - Validate CAMK2D findings with known cardiac biology
   - Confirm statistical significance across multiple studies
   - Verify effect size consistency and biological plausibility

4. Interactive components and visualizations:
   - Interactive volcano plots with gene highlighting
   - Dynamic forest plots with study filtering options
   - Searchable results tables with sorting/filtering
   - Cross-species expression comparison tools
   - Pathway network visualizations

5. Reproducibility package:
   - Complete R code with session info and package versions
   - Step-by-step analysis protocol with timing estimates
   - Data download verification scripts
   - Quality control checkpoints and validation procedures
   - Instructions for replicating analysis

6. Clinical translation section:
   - Biomarker potential assessment for CAMK2D
   - Drug target prioritization with existing compound analysis
   - Clinical trial design recommendations
   - Diagnostic and prognostic applications

7. Publication-ready outputs:
   - High-resolution figures in multiple formats (PDF, PNG, SVG)
   - Publication-quality tables with proper statistical reporting
   - Supplementary data files formatted for journal submission
   - Methods section text with complete statistical procedures

8. Data sharing and archival:
   - Processed datasets formatted for public repositories
   - Metadata files describing all variables and procedures
   - Analysis code repository with version control
   - User guide for interpreting and extending results

9. Quality assurance final checks:
   - Verify all datasets were properly downloaded and processed
   - Cross-validate results across independent analysis methods
   - Check statistical assumptions and model appropriateness
   - Confirm biological interpretation aligns with known literature
   - Validate that all original research objectives were addressed

The final report should be publication-ready and demonstrate clear evidence for CAMK2D's role in cardiovascular disease with actionable insights for drug development and clinical applications.
```

## Additional Data Verification Prompt

```
Before finalizing any analysis, create a comprehensive data verification module that ensures all transcriptomics data was correctly downloaded, processed, and analyzed.

Requirements:
1. **Download Verification:**
   - Confirm all GEO datasets were successfully downloaded
   - Verify file sizes and checksums where available
   - Check that supplementary files contain expected data types
   - Validate that sample sizes match published studies

2. **Data Processing Verification:**
   - Ensure gene symbol mapping was successful for CAMK family
   - Verify that normalization procedures were appropriate for each platform
   - Check that batch correction preserved biological signal
   - Confirm cross-species ortholog mapping accuracy

3. **Analysis Verification:**
   - Cross-validate differential expression results using alternative methods
   - Check statistical assumptions for all analyses
   - Verify meta-analysis calculations manually for key results
   - Confirm pathway enrichment results using multiple databases

4. **Biological Validation:**
   - Compare findings to known CAMK2D biology and literature
   - Verify that effect directions are consistent across studies
   - Check for biologically implausible results
   - Validate cross-species conservation patterns

This verification ensures scientific rigor and reproducibility of all findings.
```

## Usage Instructions for Claude Code

To use these updated prompts with Claude Code:

1. **Start with Prompt 1** - This now includes comprehensive data download with specific error handling
2. **Use Prompt 2** - Enhanced preprocessing with platform-specific handling  
3. **Apply Prompt 2.5** - Critical validation step to ensure data integrity
4. **Continue with Prompt 3** - Now includes specific gene symbol handling
5. **Add Prompt 3.5** - Integrate large-scale databases for more power
6. **Proceed through remaining prompts** - Each building comprehensive analysis capabilities
7. **Finish with verification** - Ensure all data and results are valid

Each prompt now includes specific dataset accession numbers, expected sample sizes, and detailed data handling procedures to ensure successful transcriptomics data acquisition and analysis.