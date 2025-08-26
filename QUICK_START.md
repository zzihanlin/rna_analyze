# Quick Start Guide - RNA/Cytokine Analysis

## 🚀 Get Started in 3 Steps

### 1. Install R Packages
```r
# In R or RStudio, run:
source("install_packages.R")
```

### 2. Run the Analysis
```bash
# Option A: Use the shell script (recommended)
./run_analysis.sh

# Option B: Run directly in R
source("rna_analysis.R")
```

### 3. View Results
Open `output/rna_analysis_report.html` in your web browser!

## 📊 What You'll Get

- **Statistical Tables**: CSV files with t-test results, p-values, and effect sizes
- **Boxplots**: Expression distributions by treatment group
- **Heatmaps**: Expression patterns across all samples
- **PCA Plots**: Sample clustering visualization
- **Volcano Plots**: Statistical significance vs. fold change
- **Correlation Matrices**: Inter-cytokine relationships
- **HTML Report**: Complete analysis summary with all figures

## 🔬 Data Analyzed

**C2 Cytokines**: IL-15, IL-17A/F, IL-27p28/IL-30, IL-33, IL-9, IP-10, MCP-1, MIP-1α, MIP-2

**P2 Cytokines**: IFN-γ, IL-10, IL-12p70, IL-1β, IL-2, IL-4, IL-5, IL-6, KC/GRO, TNF-α

**Groups**: S001-S008 (Group S) vs U001-U008 (Group U)

## 📁 File Structure
```
rna_analyze-2/
├── 11-19-20-C2.xlsx          # C2 cytokine data
├── 11-19-20-P2.xlsx          # P2 cytokine data  
├── genes.txt                  # Gene annotations
├── rna_analysis.R            # Main analysis script
├── install_packages.R        # Package installer
├── run_analysis.sh           # Automated runner
├── README.md                 # Detailed documentation
└── QUICK_START.md           # This file
```

## ⚡ Need Help?

- Check `README.md` for detailed documentation
- Ensure all data files are in the same directory
- Verify R is installed and accessible
- Check the troubleshooting section in README.md

## 🎯 Expected Output

After running the analysis, you'll have:
- `output/` folder with CSV tables and HTML report
- `plots/` folder with all visualization files
- Console output showing key statistical results

**Happy analyzing! 🧬**

