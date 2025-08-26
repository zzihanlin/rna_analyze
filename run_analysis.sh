#!/bin/bash

# RNA/Cytokine Analysis Pipeline Runner
# This script automates the R analysis process

echo "=========================================="
echo "RNA/Cytokine Analysis Pipeline"
echo "=========================================="

# Check if R is installed
if ! command -v Rscript &> /dev/null; then
    echo "Error: R is not installed or not in PATH"
    echo "Please install R from https://cran.r-project.org/"
    exit 1
fi

# Check if required files exist
if [ ! -f "11-19-20-C2.xlsx" ]; then
    echo "Error: C2 data file (11-19-20-C2.xlsx) not found"
    exit 1
fi

if [ ! -f "11-19-20-P2.xlsx" ]; then
    echo "Error: P2 data file (11-19-20-P2.xlsx) not found"
    exit 1
fi

if [ ! -f "genes.txt" ]; then
    echo "Error: Gene annotations file (genes.txt) not found"
    exit 1
fi

echo "All required data files found ✓"
echo ""

# Create output directories
echo "Creating output directories..."
mkdir -p output plots

# Install R packages if needed
echo "Checking and installing R packages..."
Rscript install_packages.R

# Run the main analysis
echo ""
echo "Running main analysis..."
Rscript rna_analysis.R

echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
echo ""
echo "Output files generated:"
echo "- output/c2_cytokine_summary.csv"
echo "- output/p2_cytokine_summary.csv"
echo "- output/rna_analysis_report.html"
echo "- plots/ (various visualization files)"
echo ""
echo "Open 'output/rna_analysis_report.html' in your browser to view results!"
echo ""

