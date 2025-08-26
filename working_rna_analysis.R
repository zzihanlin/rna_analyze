#!/usr/bin/env Rscript

# Working RNA/Cytokine Analysis Script
# Analysis of treatment effects on cytokine expression
# Author: RNA Analysis Pipeline
# Date: 2024

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggplot2)
  library(corrplot)
  library(RColorBrewer)
  library(reshape2)
})

# Set working directory and create output folder
setwd("/Users/suwa/Desktop/rna_analyze-2")
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("plots")) dir.create("plots")

# Function to clean Excel data - Working approach
clean_excel_data <- function(file_path) {
  # Read the Excel file, skipping the first 2 rows
  df <- read_excel(file_path, skip = 2)
  
  # The first row contains cytokine names, second row contains sample IDs
  cytokine_names <- as.character(df[1, -1])  # Exclude first column
  sample_ids <- as.character(df[-1, 1])      # First column, excluding header row
  
  # Remove the "Row Labels" header if present
  if (sample_ids[1] == "Row Labels") {
    sample_ids <- sample_ids[-1]
    df <- df[-1, ]
  }
  
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
  
  return(clean_df)
}

# Load and clean data
cat("Loading and cleaning data...\n")

# Load cytokine data
c2_data <- clean_excel_data("11-19-20-C2.xlsx")
p2_data <- clean_excel_data("11-19-20-P2.xlsx")

# Load gene annotations
gene_annotations <- read.delim("genes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Display data structure
cat("\nData loaded successfully!\n")
cat("C2 data dimensions:", dim(c2_data), "\n")
cat("P2 data dimensions:", dim(p2_data), "\n")
cat("Gene annotations:", nrow(gene_annotations), "genes\n")

# Function to identify sample groups - Working version
identify_groups <- function(data) {
  # Extract sample IDs from first column
  sample_ids <- data$Sample_ID
  
  # Create group assignments based on sample naming convention
  # S001-S008: Treatment group 1, U001-U026: Treatment group 2
  groups <- character(length(sample_ids))
  
  for (i in 1:length(sample_ids)) {
    sample_id <- sample_ids[i]
    if (substr(sample_id, 1, 1) == "S") {
      groups[i] <- "Group_S"
    } else {
      groups[i] <- "Group_U"
    }
  }
  
  cat("Group assignment check:\n")
  cat("Total samples:", length(sample_ids), "\n")
  cat("Group S count:", sum(groups == "Group_S"), "\n")
  cat("Group U count:", sum(groups == "Group_U"), "\n")
  cat("First few sample IDs:", head(sample_ids, 5), "\n")
  cat("First few groups:", head(groups, 5), "\n")
  
  return(data.frame(
    Sample_ID = sample_ids,
    Group = groups,
    stringsAsFactors = FALSE
  ))
}

# Create sample groups
c2_groups <- identify_groups(c2_data)
p2_groups <- identify_groups(p2_data)

# Function to prepare data for analysis
prepare_analysis_data <- function(data, groups) {
  # Remove sample ID column and keep only numeric data
  numeric_data <- data[, -1, drop = FALSE]
  
  # Create analysis dataframe
  analysis_df <- data.frame(
    Sample_ID = groups$Sample_ID,
    Group = groups$Group,
    numeric_data,
    stringsAsFactors = FALSE
  )
  
  return(analysis_df)
}

# Prepare data for analysis
c2_analysis <- prepare_analysis_data(c2_data, c2_groups)
p2_analysis <- prepare_analysis_data(p2_data, p2_groups)

# Function to perform statistical analysis
perform_statistical_analysis <- function(data, group_col = "Group") {
  results <- list()
  
  # Get numeric columns (excluding sample ID and group)
  numeric_cols <- names(data)[sapply(data, is.numeric)]
  
  cat("Analyzing", length(numeric_cols), "cytokines...\n")
  
  for (col in numeric_cols) {
    # Perform t-test
    group1_data <- data[data[[group_col]] == "Group_S", col]
    group2_data <- data[data[[group_col]] == "Group_U", col]
    
    # Remove NA values
    group1_data <- group1_data[!is.na(group1_data)]
    group2_data <- group2_data[!is.na(group2_data)]
    
    if (length(group1_data) > 0 && length(group2_data) > 0) {
      t_test <- t.test(group1_data, group2_data)
      
      # Calculate effect size (Cohen's d)
      pooled_sd <- sqrt(((length(group1_data) - 1) * var(group1_data) + 
                         (length(group2_data) - 1) * var(group2_data)) / 
                        (length(group1_data) + length(group2_data) - 2))
      cohens_d <- (mean(group1_data) - mean(group2_data)) / pooled_sd
      
      results[[col]] <- list(
        mean_group1 = mean(group1_data, na.rm = TRUE),
        mean_group2 = mean(group2_data, na.rm = TRUE),
        sd_group1 = sd(group1_data, na.rm = TRUE),
        sd_group2 = sd(group2_data, na.rm = TRUE),
        p_value = t_test$p.value,
        t_statistic = t_test$statistic,
        cohens_d = cohens_d,
        significant = t_test$p.value < 0.05
      )
    } else {
      cat("Warning: Insufficient data for", col, "\n")
    }
  }
  
  cat("Statistical analysis completed for", length(results), "cytokines\n")
  return(results)
}

# Perform statistical analysis
cat("\nPerforming statistical analysis...\n")
c2_stats <- perform_statistical_analysis(c2_analysis)
p2_stats <- perform_statistical_analysis(p2_analysis)

# Print results to console
cat("\nC2 Cytokines Analysis Results:\n")
cat("===============================\n")
for (cytokine in names(c2_stats)) {
  stats <- c2_stats[[cytokine]]
  cat(sprintf("%s:\n", cytokine))
  cat(sprintf("  Group S: Mean = %.2f, SD = %.2f\n", stats$mean_group1, stats$sd_group1))
  cat(sprintf("  Group U: Mean = %.2f, SD = %.2f\n", stats$mean_group2, stats$sd_group2))
  cat(sprintf("  P-value = %.4f\n", stats$p_value))
  cat(sprintf("  Effect size (Cohen's d) = %.3f\n", stats$cohens_d))
  cat(sprintf("  Significant: %s\n", ifelse(stats$significant, "YES", "NO")))
  cat("\n")
}

cat("\nP2 Cytokines Analysis Results:\n")
cat("===============================\n")
for (cytokine in names(p2_stats)) {
  stats <- p2_stats[[cytokine]]
  cat(sprintf("%s:\n", cytokine))
  cat(sprintf("  Group S: Mean = %.2f, SD = %.2f\n", stats$mean_group1, stats$sd_group1))
  cat(sprintf("  Group U: Mean = %.2f, SD = %.2f\n", stats$mean_group2, stats$sd_group2))
  cat(sprintf("  P-value = %.4f\n", stats$p_value))
  cat(sprintf("  Effect size (Cohen's d) = %.3f\n", stats$cohens_d))
  cat(sprintf("  Significant: %s\n", ifelse(stats$significant, "YES", "NO")))
  cat("\n")
}

# Function to create summary statistics table
create_summary_table <- function(stats, title) {
  if (length(stats) == 0) return(NULL)
  
  summary_df <- data.frame(
    Cytokine = names(stats),
    Mean_Group_S = sapply(stats, function(x) round(x$mean_group1, 3)),
    Mean_Group_U = sapply(stats, function(x) round(x$mean_group2, 3)),
    SD_Group_S = sapply(stats, function(x) round(x$sd_group1, 3)),
    SD_Group_U = sapply(stats, function(x) round(x$sd_group2, 3)),
    P_Value = sapply(stats, function(x) format.pval(x$p_value, digits = 3)),
    T_Statistic = sapply(stats, function(x) round(x$t_statistic, 3)),
    Cohens_D = sapply(stats, function(x) round(x$cohens_d, 3)),
    Significant = sapply(stats, function(x) ifelse(x$significant, "Yes", "No")),
    stringsAsFactors = FALSE
  )
  
  # Sort by p-value
  summary_df <- summary_df[order(summary_df$P_Value), ]
  
  return(summary_df)
}

# Create summary tables
c2_summary <- create_summary_table(c2_stats, "C2 Cytokines")
p2_summary <- create_summary_table(p2_stats, "P2 Cytokines")

# Save summary tables
write.csv(c2_summary, "output/c2_cytokine_summary.csv", row.names = FALSE)
write.csv(p2_summary, "output/p2_cytokine_summary.csv", row.names = FALSE)

# Function to create boxplots
create_boxplots <- function(data, group_col = "Group", title = "Cytokine Expression by Group") {
  # Melt data for plotting
  plot_data <- melt(data, id.vars = c("Sample_ID", group_col), 
                    variable.name = "Cytokine", value.name = "Expression")
  
  # Remove NA values
  plot_data <- plot_data[!is.na(plot_data$Expression), ]
  
  # Create boxplot
  p <- ggplot(plot_data, aes(x = Cytokine, y = Expression, fill = Group)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 21) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 1) +
    scale_fill_brewer(palette = "Set2", name = "Group") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "bottom") +
    labs(title = title,
         x = "Cytokine",
         y = "Expression Level",
         caption = "Boxplots show median, quartiles, and outliers. Points represent individual samples.")
  
  return(p)
}

# Create boxplots
cat("\nCreating visualizations...\n")

# C2 cytokines boxplots
c2_boxplot <- create_boxplots(c2_analysis, title = "C2 Cytokine Expression by Group")
ggsave("plots/c2_cytokine_boxplots.png", c2_boxplot, width = 12, height = 8, dpi = 300)

# P2 cytokines boxplots
p2_boxplot <- create_boxplots(p2_analysis, title = "P2 Cytokine Expression by Group")
ggsave("plots/p2_cytokine_boxplots.png", p2_boxplot, width = 12, height = 8, dpi = 300)

# Function to create correlation plots
create_correlation_plot <- function(data, title = "Cytokine Correlation Matrix") {
  # Extract numeric data
  numeric_data <- data[, sapply(data, is.numeric)]
  
  # Calculate correlation matrix
  cor_matrix <- cor(numeric_data, use = "pairwise.complete.obs")
  
  # Create correlation plot
  png(paste0("plots/", gsub(" ", "_", tolower(title)), ".png"), 
      width = 10, height = 8, units = "in", res = 300)
  corrplot(cor_matrix, 
           method = "color",
           type = "upper",
           order = "hclust",
           tl.cex = 0.8,
           tl.col = "black",
           addCoef.col = "black",
           number.cex = 0.6,
           col = brewer.pal(n = 8, name = "RdBu"))
  title(title, line = 1)
  dev.off()
}

# Create correlation plots
create_correlation_plot(c2_analysis, "C2 Cytokine Correlation Matrix")
create_correlation_plot(p2_analysis, "P2 Cytokine Correlation Matrix")

# Create simple summary report
cat("\nGenerating summary report...\n")

report_content <- paste0('
RNA/Cytokine Analysis Summary Report
====================================

Analysis Date: ', Sys.Date(), '
Data Files: 11-19-20-C2.xlsx, 11-19-20-P2.xlsx

Sample Groups:
- Group S: Samples S001-S008 (Treatment group 1)
- Group U: Samples U001-U026 (Treatment group 2)

C2 Cytokines Analyzed: ', paste(colnames(c2_data)[-1], collapse = ", "), '

P2 Cytokines Analyzed: ', paste(colnames(p2_data)[-1], collapse = ", "), '

Statistical Analysis:
- Test: Two-sample t-test
- Effect size: Cohen\'s d
- Significance threshold: p < 0.05

Key Findings:
============

C2 Cytokines - Significantly Different (p < 0.05):
', paste(sapply(names(c2_stats), function(cytokine) {
    if (c2_stats[[cytokine]]$significant) {
        paste0("- ", cytokine, ": p = ", format.pval(c2_stats[[cytokine]]$p_value, digits = 3),
               ", Effect size = ", round(c2_stats[[cytokine]]$cohens_d, 3))
    }
}), collapse = "\n"), '

P2 Cytokines - Significantly Different (p < 0.05):
', paste(sapply(names(p2_stats), function(cytokine) {
    if (p2_stats[[cytokine]]$significant) {
        paste0("- ", cytokine, ": p = ", format.pval(p2_stats[[cytokine]]$p_value, digits = 3),
               ", Effect size = ", round(p2_stats[[cytokine]]$cohens_d, 3))
    }
}), collapse = "\n"), '

Files Generated:
- output/c2_cytokine_summary.csv
- output/p2_cytokine_summary.csv
- plots/c2_cytokine_boxplots.png
- plots/p2_cytokine_boxplots.png
- plots/c2_cytokine_correlation_matrix.png
- plots/p2_cytokine_correlation_matrix.png

Analysis complete!
')

# Save summary report
writeLines(report_content, "output/analysis_summary.txt")

# Print final summary
cat("\n==================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("==================================================\n\n")

cat("Files generated:\n")
cat("- output/c2_cytokine_summary.csv\n")
cat("- output/p2_cytokine_summary.csv\n")
cat("- output/analysis_summary.txt\n")
cat("- plots/ (visualization files)\n\n")

cat("Key Results:\n")
cat("C2 Cytokines - Significantly different:\n")
for (cytokine in names(c2_stats)) {
  if (c2_stats[[cytokine]]$significant) {
    cat(sprintf("  %s: p = %.4f, Effect size = %.3f\n", 
                cytokine, c2_stats[[cytokine]]$p_value, c2_stats[[cytokine]]$cohens_d))
  }
}

cat("\nP2 Cytokines - Significantly different:\n")
for (cytokine in names(p2_stats)) {
  if (p2_stats[[cytokine]]$significant) {
    cat(sprintf("  %s: p = %.4f, Effect size = %.3f\n", 
                cytokine, p2_stats[[cytokine]]$p_value, p2_stats[[cytokine]]$cohens_d))
  }
}

cat("\nCheck 'output/analysis_summary.txt' for a complete summary!\n")

