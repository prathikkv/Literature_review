# Claude Code Prompts for Dynamic Pharmaceutical Pipeline

## 🎯 PROMPT 1: Dynamic Dataset Download System
**Goal**: Auto-download new GEO datasets when added to YAML

```
I have an existing R pipeline for microarray analysis. I need to enhance it to automatically download GEO datasets that aren't cached.

CURRENT STRUCTURE:
- config.yml contains dataset IDs
- cache/ directory stores processed datasets  
- Uses GEOquery for data download
- Current datasets: GSE57338, GSE41177, GSE79768, GSE115574

REQUIREMENTS:
1. Create function `auto_download_geo_datasets()` that:
   - Reads dataset list from config.yml
   - Checks if datasets exist in cache/
   - Downloads missing datasets using GEOquery
   - Handles GPL platform detection automatically
   - Validates download integrity
   - Saves in standardized format (.rds)
   - Updates download log

2. Error handling for:
   - Network failures
   - Invalid GSE IDs
   - Platform mismatches
   - Corrupted downloads

3. Integration with existing pipeline without breaking current functionality

Please provide the complete function with proper logging and validation.
```

## 🧬 PROMPT 2: Enhanced YAML Configuration System
**Goal**: Make config.yml fully dynamic for any gene/disease combination

```
Transform my static config.yml into a dynamic configuration system for pharmaceutical research.

CURRENT: Static CAMK2D cardiovascular analysis
NEEDED: Dynamic gene family + disease analysis

CREATE NEW config.yml STRUCTURE:
```yaml
# Dynamic Research Configuration
research_target:
  primary_gene: "CAMK2D"
  gene_family: []  # Auto-populated
  diseases: ["Heart Failure", "Atrial Fibrillation"]
  
datasets:
  manual_list: ["GSE57338", "GSE41177", "GSE79768", "GSE115574"]
  auto_discovery: true
  
analysis_parameters:
  p_value_threshold: 0.05
  log_fc_threshold: 0.1
  min_sample_size: 20
  
output_settings:
  include_pathways: true
  include_drug_targets: true
  include_clinical_trials: true
```

Plus R functions to:
1. Validate configuration
2. Auto-populate gene families
3. Handle parameter updates
4. Maintain backwards compatibility

Ensure this works with my existing pipeline structure.
```

## 🔍 PROMPT 3: Automatic GEO Dataset Discovery
**Goal**: Find relevant datasets automatically based on disease terms

```
I need to build an automated dataset discovery system for my microarray pipeline.

REQUIREMENTS:
1. Search GEO database using disease terms from config.yml
2. Filter datasets by quality criteria
3. Create Excel report with inclusion/exclusion reasons
4. Auto-download qualified datasets

CREATE FUNCTION `discover_geo_datasets()` that:

SEARCH CRITERIA:
- Disease terms: "Heart Failure", "Atrial Fibrillation" 
- Species: Homo sapiens
- Platform: GPL570, GPL96, GPL97
- Study type: Expression profiling by array

QUALITY FILTERS:
- Minimum samples: 20 (configurable)
- Must have case/control groups
- Exclude cell lines and in vitro studies
- Sample size < 1000 (exclude outliers)
- Must have proper group annotations

OUTPUT:
1. Excel file with columns:
   - GSE_ID, Title, Samples, Platform, Disease_Match, 
   - Include_Status, Exclusion_Reason, Quality_Score
2. Auto-download qualified datasets
3. Update config.yml with new dataset list

Include robust error handling and progress reporting.
```

## 🧪 PROMPT 4: GO/KEGG Pathway Analysis Integration
**Goal**: Add comprehensive pathway analysis for pharmaceutical insights

```
Enhance my existing DGE pipeline with comprehensive pathway analysis for drug discovery.

CURRENT: Basic differential gene expression
NEEDED: Full pathway analysis + drug target prediction

ADD TO EXISTING PIPELINE:
1. GO enrichment analysis (Biological Process, Molecular Function, Cellular Component)
2. KEGG pathway analysis
3. Reactome pathway analysis
4. Drug target prediction using DrugBank
5. Druggability assessment

CREATE FUNCTIONS:
- `run_pathway_analysis(dge_results, gene_list)`
- `predict_drug_targets(pathway_results)`
- `assess_druggability(gene_targets)`

OUTPUT REQUIREMENTS:
- Pathway enrichment tables
- Drug target ranking tables
- Druggability scores
- Clinical relevance assessments
- Visualization plots

INTEGRATION:
- Must work with existing meta-analysis results
- Output should feed into dynamic Rmd template
- Results formatted for pharmaceutical teams

Use clusterProfiler, DOSE, and ReactomePA packages.
```

## 🧬 PROMPT 5: Gene Family Discovery System
**Goal**: Auto-discover gene families from primary gene input

```
Build an intelligent gene family discovery system for my pharmaceutical pipeline.

CURRENT: Manual CAMK gene family definition
NEEDED: Automatic family discovery for any input gene

CREATE FUNCTION `discover_gene_family()` that:

INPUT: Primary gene (e.g., "CAMK2D")

DISCOVERY METHODS:
1. NCBI Gene database - search gene family
2. UniProt - protein family classification  
3. KEGG - pathway co-membership
4. String-DB - protein interactions
5. GO terms - functional similarity
6. Literature mining - PubMed co-occurrence

QUALITY SCORING:
- Evidence strength (1-10)
- Functional relatedness
- Expression correlation
- Pathway overlap

OUTPUT:
1. Ranked list of family members with evidence scores
2. Update config.yml with discovered families
3. Gene family report with justifications
4. Family tree visualization

ERROR HANDLING:
- API rate limiting
- Invalid gene symbols
- Network failures
- Empty results

Must integrate seamlessly with existing pipeline and maintain scientific rigor.
```

## 📊 PROMPT 6: Dynamic Rmd Template System
**Goal**: Make report template fully adaptive to any gene/disease combination

```
Transform my static Rmd report into a fully dynamic template that adapts to any analysis.

CURRENT: Fixed CAMK2D cardiovascular report
NEEDED: Universal template for any gene/disease

REQUIREMENTS:
1. Read all parameters from config.yml
2. Adapt content based on:
   - Number of datasets found
   - Gene family size  
   - Pathway analysis results
   - Statistical significance levels

DYNAMIC ELEMENTS:
- Title: Auto-generate from gene + disease
- Introduction: Pull from literature mining
- Methods: Adapt based on datasets used
- Results: Scale based on findings
- Discussion: Generated from pathway analysis
- Conclusions: Based on statistical outcomes

CONDITIONAL CONTENT:
```r
# Example dynamic blocks
if(length(significant_pathways) > 5) {
  # Include detailed pathway section
} else {
  # Brief pathway summary
}

if(drug_targets_found) {
  # Add drug development section
}
```

HYPERLINK INTEGRATION:
- PubMed citations → clickable links
- Gene symbols → NCBI Gene links  
- Pathways → KEGG/Reactome links
- Clinical trials → ClinicalTrials.gov

Create template that maintains scientific quality while being completely adaptive.
```

## 🔗 PROMPT 7: Literature Mining and Hyperlink System
**Goal**: Add intelligent literature mining with dynamic hyperlinks

```
Build a literature mining system that adds context and clickable references to my reports.

REQUIREMENTS:
1. Mine PubMed for gene + disease combinations
2. Extract key research themes
3. Find related clinical trials
4. Generate hyperlinked citations

CREATE FUNCTIONS:
- `mine_pubmed_literature(gene, diseases)`
- `extract_clinical_trials(gene, diseases)`  
- `generate_hyperlinks(text_content)`

PUBMED MINING:
- Search: (gene[Title/Abstract] AND disease[Title/Abstract])
- Extract: Key findings, therapeutic relevance, mechanisms
- Rank: By citation count and recency
- Limit: Top 20 most relevant papers

CLINICAL TRIALS:
- Search ClinicalTrials.gov API
- Filter: Active trials, completed trials
- Extract: Phase, status, intervention, outcomes
- Link: Direct to trial pages

HYPERLINK GENERATION:
- Gene symbols → NCBI Gene database
- Pathways → KEGG/Reactome/GO
- Publications → PubMed
- Clinical trials → ClinicalTrials.gov
- Diseases → Disease ontology

OUTPUT:
- Literature summary with embedded links
- Clinical trials table with links
- Reference list with DOI links
- Context paragraphs for report sections

Must integrate with dynamic Rmd template and handle API rate limits.
```

## 🚀 PROMPT 8: Master Orchestration System
**Goal**: Create single-command execution for entire pipeline

```
Build a master orchestration system that runs the complete pharmaceutical analysis pipeline.

CURRENT: Multiple manual steps
NEEDED: Single command → Complete analysis

CREATE `run_pharma_pipeline.R` that:

EXECUTION FLOW:
1. Validate config.yml parameters
2. Discover gene family (if needed)
3. Search/download datasets (if auto_discovery enabled)
4. Process all datasets
5. Run DGE analysis
6. Perform pathway analysis
7. Mine literature and clinical trials
8. Generate dynamic report
9. Create supplementary files

PARAMETER VALIDATION:
- Valid gene symbols
- Recognizable disease terms
- Reasonable statistical thresholds
- Available computational resources

PROGRESS REPORTING:
- Real-time progress updates
- Estimated completion time
- Error reporting with suggestions
- Success/failure notifications

ERROR RECOVERY:
- Partial failure handling
- Resume from checkpoint
- Alternative dataset suggestions
- Graceful degradation

OUTPUTS:
- HTML report
- Excel data tables
- Quality assessment reports
- Pipeline execution log

COMMAND LINE INTERFACE:
```bash
Rscript run_pharma_pipeline.R --gene CAMK2D --diseases "Heart Failure,Atrial Fibrillation" --config config.yml
```

Must maintain backwards compatibility and scientific accuracy.
```

## ✅ PROMPT 9: Validation Framework
**Goal**: Ensure new pipeline produces identical results to current version

```
Create a validation framework to ensure the dynamic pipeline maintains scientific accuracy.

VALIDATION REQUIREMENTS:
1. Run new pipeline with original CAMK2D parameters
2. Compare all outputs with existing results
3. Validate statistical consistency
4. Check report content accuracy

CREATE VALIDATION FUNCTIONS:
- `validate_dge_results(new_results, original_results)`
- `validate_meta_analysis(new_meta, original_meta)`
- `validate_report_content(new_html, original_html)`
- `validate_statistical_significance(new_stats, original_stats)`

COMPARISON METRICS:
- Log fold changes (tolerance: ±0.001)
- P-values (tolerance: ±0.001)  
- Meta-analysis coefficients
- Sample sizes and dataset counts
- Figure content and captions

TEST CASES:
1. Original CAMK2D + cardiovascular diseases
2. Different gene with same diseases
3. Same gene with different diseases
4. Multiple gene families
5. Various dataset combinations

VALIDATION REPORT:
- Pass/fail for each component
- Detailed discrepancy analysis
- Performance benchmarks
- Recommendations for improvements

OUTPUT:
- Validation summary table
- Detailed comparison plots
- Statistical test results
- Confidence assessment

Must ensure 100% reproducibility of current results before deploying new features.
```

## 🎯 IMPLEMENTATION STRATEGY

### Phase 1 (Week 1): Core Infrastructure
- Prompts 1, 2: Dynamic downloads + YAML system
- Validation: Ensure current pipeline still works

### Phase 2 (Week 2): Discovery Systems  
- Prompts 3, 5: Dataset discovery + gene families
- Integration testing

### Phase 3 (Week 3): Analysis Enhancement
- Prompt 4: Pathway analysis + drug targets
- Scientific validation

### Phase 4 (Week 4): Dynamic Reporting
- Prompts 6, 7: Dynamic templates + literature mining
- Content verification

### Phase 5 (Week 5): Integration & Validation
- Prompts 8, 9: Orchestration + validation
- Production testing

## 🏆 FINAL PRODUCT CAPABILITIES

**INPUT**: Gene symbol + Disease terms + Parameters
**OUTPUT**: Complete pharmaceutical research report

**USER EXPERIENCE**:
```bash
# Single command execution
Rscript run_pharma_pipeline.R --gene TP53 --diseases "Lung Cancer,Breast Cancer"

# Output: Comprehensive analysis in 1-2 hours
```

**REPORT INCLUDES**:
- Automated dataset discovery and quality assessment
- Differential gene expression meta-analysis  
- Pathway analysis and drug target predictions
- Literature review with clickable references
- Clinical trials identification
- Executive summary for pharma teams

This system will revolutionize pharmaceutical target validation and drug discovery research!