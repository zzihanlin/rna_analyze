# RNA/Cytokine Analysis Pipeline

This repository contains a comprehensive R-based analysis pipeline for comparing cytokine expression between treatment groups. The analysis includes statistical testing, visualization, and generation of publication-ready figures and tables.

## Data Description

The analysis works with two main datasets:

### C2 Cytokines (11-19-20-C2.xlsx)
- **IL-15**: Interleukin 15
- **IL-17A/F**: Interleukin 17A/F
- **IL-27p28/IL-30**: Interleukin 27p28/IL-30
- **IL-33**: Interleukin 33
- **IL-9**: Interleukin 9
- **IP-10**: Interferon gamma-induced protein 10
- **MCP-1**: Monocyte chemoattractant protein 1
- **MIP-1α**: Macrophage inflammatory protein 1 alpha
- **MIP-2**: Macrophage inflammatory protein 2

### P2 Cytokines (11-19-20-P2.xlsx)
- **IFN-γ**: Interferon gamma
- **IL-10**: Interleukin 10
- **IL-12p70**: Interleukin 12p70
- **IL-1β**: Interleukin 1 beta
- **IL-2**: Interleukin 2
- **IL-4**: Interleukin 4
- **IL-5**: Interleukin 5
- **IL-6**: Interleukin 6
- **KC/GRO**: Keratinocyte chemoattractant/Growth-regulated oncogene
- **TNF-α**: Tumor necrosis factor alpha

## Sample Groups

The analysis automatically identifies two treatment groups based on sample naming:
- **Group S**: Samples S001-S008 (Treatment group 1)
- **Group U**: Samples U001-U008 (Treatment group 2)

## Installation and Setup

### Prerequisites
- R (version 4.0 or higher)
- RStudio (recommended)

### Install Required Packages

1. **First time setup**: Run the package installation script:
   ```r
   source("install_packages.R")
   ```

2. **Alternative manual installation**:
   ```r
   install.packages(c("tidyverse", "readxl", "pheatmap", "ggplot2", "ggpubr", 
                      "corrplot", "factoextra", "DT", "knitr", "kableExtra", 
                      "reshape2", "viridis", "RColorBrewer"))
   ```

## Running the Analysis

### Quick Start
1. Ensure all data files are in the working directory
2. Run the main analysis script:
   ```r
   source("rna_analysis.R")
   ```

### Step-by-Step Execution
1. **Load and clean data**: The script automatically reads Excel files and cleans the data
2. **Statistical analysis**: Performs t-tests between groups for each cytokine
3. **Visualization**: Generates multiple plot types
4. **Report generation**: Creates comprehensive HTML report and CSV tables

## Output Files

### Generated Directories
- `output/`: Contains summary tables and HTML report
- `plots/`: Contains all visualization files

### Summary Tables
- `c2_cytokine_summary.csv`: Statistical results for C2 cytokines
- `p2_cytokine_summary.csv`: Statistical results for P2 cytokines

### Visualizations
- **Boxplots**: Expression distributions by group
- **Heatmaps**: Expression patterns across samples
- **PCA plots**: Principal component analysis for sample clustering
- **Volcano plots**: Statistical significance vs. fold change
- **Correlation matrices**: Inter-cytokine relationships

### HTML Report
- `rna_analysis_report.html`: Comprehensive analysis report with all results and figures

## Statistical Analysis

### Methods
- **Test**: Two-sample t-test (parametric)
- **Effect size**: Cohen's d
- **Significance threshold**: p < 0.05
- **Multiple testing**: No correction applied (can be modified if needed)

### Output Metrics
- Mean expression by group
- Standard deviation by group
- P-value
- T-statistic
- Effect size (Cohen's d)
- Significance indicator

## Customization Options

### Modify Group Assignments
Edit the `identify_groups()` function in `rna_analysis.R` to change how samples are grouped.

### Adjust Significance Threshold
Change the significance threshold in the `perform_statistical_analysis()` function:
```r
significant = t_test$p.value < 0.01  # Change from 0.05 to 0.01
```

### Modify Plot Aesthetics
Adjust colors, themes, and plot parameters in the visualization functions.

### Add Additional Statistical Tests
Extend the `perform_statistical_analysis()` function to include:
- Non-parametric tests (Wilcoxon rank-sum)
- Multiple testing corrections (FDR, Bonferroni)
- ANOVA for multiple groups
- Mixed-effects models for repeated measures

## Troubleshooting

### Common Issues

1. **Package installation errors**:
   - Ensure R is up to date
   - Try installing packages individually
   - Check internet connection for CRAN access

2. **Data loading errors**:
   - Verify Excel files are in the working directory
   - Check file permissions
   - Ensure Excel files are not corrupted

3. **Memory issues with large datasets**:
   - Increase R memory limit
   - Process data in chunks
   - Use data.table for large datasets

### Data Format Requirements
- Excel files should have cytokine names in row 3
- Sample IDs should be in the first column
- Numeric data should start from row 4
- Missing values should be properly formatted

## Advanced Features

### Gene Annotation Integration
The script loads `genes.txt` which contains gene annotations for different cell populations. This can be used to:
- Categorize cytokines by cell type
- Perform pathway analysis
- Generate cell-type specific visualizations

### Batch Effect Correction
For studies with multiple experimental batches, consider adding:
- ComBat normalization
- SVA (Surrogate Variable Analysis)
- RUVSeq normalization

### Machine Learning Integration
Extend the analysis with:
- Random forest for feature selection
- Support vector machines for classification
- Clustering algorithms (k-means, hierarchical)

## Citation and References

When using this pipeline, please cite:
- R Core Team (2023). R: A language and environment for statistical computing
- Wickham H (2016). ggplot2: Elegant Graphics for Data Analysis
- Kolde R (2019). pheatmap: Pretty Heatmaps

## Support and Contributions

For questions or improvements:
1. Check the troubleshooting section
2. Review R documentation for specific packages
3. Consider contributing improvements to the pipeline

## License

This analysis pipeline is provided as-is for research and educational purposes.

