#!/usr/bin/env Rscript

# Debug script to investigate the data structure and statistical analysis
library(tidyverse)
library(readxl)

setwd("/Users/suwa/Desktop/rna_analyze-2")

# Load and clean data
cat("Loading and cleaning data...\n")

# Load cytokine data
df <- read_excel("11-19-20-C2.xlsx", skip = 2)

# The first row contains cytokine names, second row contains sample IDs
cytokine_names <- as.character(df[1, -1])  # Exclude first column
sample_ids <- as.character(df[-1, 1])      # First column, excluding header row

# Extract numeric data (starting from row 2, columns 2 onwards)
numeric_data <- df[-1, -1]

# Convert to numeric matrix
numeric_matrix <- as.matrix(numeric_data)
mode(numeric_matrix) <- "numeric"

# Create clean dataframe
clean_df <- data.frame(
  Sample_ID = sample_ids,
  numeric_matrix,
  stringsAsFactors = FALSE
)

# Set column names
colnames(clean_df) <- c("Sample_ID", cytokine_names)

# Remove rows with all NA values
clean_df <- clean_df[rowSums(is.na(clean_df)) < ncol(clean_df), ]

cat("Data loaded successfully!\n")
cat("Data dimensions:", dim(clean_df), "\n")
cat("Sample IDs:", head(clean_df$Sample_ID, 10), "\n")
cat("Cytokine names:", colnames(clean_df)[-1], "\n")

# Check data types
cat("\nData types:\n")
cat("Sample_ID:", class(clean_df$Sample_ID), "\n")
cat("First cytokine column:", class(clean_df[[2]]), "\n")
cat("First few values of first cytokine:", head(clean_df[[2]], 5), "\n")

# Check for NA values
cat("\nNA values in first cytokine column:", sum(is.na(clean_df[[2]])), "\n")
cat("Total rows:", nrow(clean_df), "\n")

# Create groups
groups <- ifelse(grepl("^S", clean_df$Sample_ID), "Group_S", "Group_U")
cat("\nGroup assignments:\n")
cat("Group S count:", sum(groups == "Group_S"), "\n")
cat("Group U count:", sum(groups == "Group_U"), "\n")

# Test statistical analysis on first cytokine
cytokine_col <- 2
cytokine_name <- colnames(clean_df)[cytokine_col]

cat("\nTesting statistical analysis on", cytokine_name, ":\n")

group1_data <- clean_df[groups == "Group_S", cytokine_col]
group2_data <- clean_df[groups == "Group_U", cytokine_col]

cat("Group S data length:", length(group1_data), "\n")
cat("Group U data length:", length(group2_data), "\n")

# Remove NA values
group1_data <- group1_data[!is.na(group1_data)]
group2_data <- group2_data[!is.na(group2_data)]

cat("Group S data length after removing NA:", length(group1_data), "\n")
cat("Group U data length after removing NA:", length(group2_data), "\n")

if (length(group1_data) > 0 && length(group2_data) > 0) {
  cat("Group S mean:", mean(group1_data), "\n")
  cat("Group U mean:", mean(group2_data), "\n")
  
  # Perform t-test
  t_test <- t.test(group1_data, group2_data)
  cat("T-test p-value:", t_test$p.value, "\n")
  
  # Calculate effect size
  pooled_sd <- sqrt(((length(group1_data) - 1) * var(group1_data) + 
                     (length(group2_data) - 1) * var(group2_data)) / 
                    (length(group1_data) + length(group2_data) - 2))
  cohens_d <- (mean(group1_data) - mean(group2_data)) / pooled_sd
  cat("Cohen's d:", cohens_d, "\n")
} else {
  cat("Cannot perform analysis - insufficient data\n")
}

