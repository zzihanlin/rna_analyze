#!/usr/bin/env Rscript

# Install required R packages for RNA/Cytokine analysis
# Run this script first to ensure all dependencies are available

cat("Installing required R packages...\n")

# List of required packages
required_packages <- c(
  "tidyverse",
  "readxl", 
  "pheatmap",
  "ggplot2",
  "ggpubr",
  "corrplot",
  "factoextra",
  "DT",
  "knitr",
  "kableExtra",
  "reshape2",
  "viridis",
  "RColorBrewer"
)

# Function to install packages if not already installed
install_if_missing <- function(package_name) {
  if (!require(package_name, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", package_name, "...\n")
    install.packages(package_name, dependencies = TRUE)
    cat(package_name, "installed successfully!\n")
  } else {
    cat(package_name, "is already installed.\n")
  }
}

# Install all required packages
for (package in required_packages) {
  install_if_missing(package)
}

cat("\nAll required packages have been installed!\n")
cat("You can now run the main analysis script: rna_analysis.R\n")

