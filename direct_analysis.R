#!/usr/bin/env Rscript

# Direct RNA/Cytokine Analysis Script
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

# Direct approach: read data and process step by step
cat("Loading and processing data...\n")

# Read C2 data
c2_raw <- read_excel("11-19-20-C2.xlsx", skip = 2)
cat("C2 raw data dimensions:", dim(c2_raw), "\n")

# Get cytokine names from first row
c2_cytokines <- as.character(c2_raw[1, -1])
cat("C2 cytokines:", paste(c2_cytokines, collapse = ", "), "\n")

# Get sample IDs from first column, excluding header
c2_samples <- as.character(c2_raw[-1, 1])
# Remove "Row Labels" if present
if (c2_samples[1] == "Row Labels") {
  c2_samples <- c2_samples[-1]
  c2_raw <- c2_raw[-1, ]
}

# Get numeric data
c2_numeric <- c2_raw[-1, -1]
c2_matrix <- as.matrix(c2_numeric)
mode(c2_matrix) <- "numeric"

# Create C2 dataframe
c2_data <- data.frame(
  Sample_ID = c2_samples,
  c2_matrix,
  stringsAsFactors = FALSE
)
colnames(c2_data) <- c("Sample_ID", c2_cytokines)

# Read P2 data
p2_raw <- read_excel("11-19-20-P2.xlsx", skip = 2)
cat("P2 raw data dimensions:", dim(p2_raw), "\n")

# Get cytokine names from first row
p2_cytokines <- as.character(p2_raw[1, -1])
cat("P2 cytokines:", paste(p2_cytokines, collapse = ", "), "\n")

# Get sample IDs from first column, excluding header
p2_samples <- as.character(p2_raw[-1, 1])
# Remove "Row Labels" if present
if (p2_samples[1] == "Row Labels") {
  p2_samples <- p2_samples[-1]
  p2_raw <- p2_raw[-1, ]
}

# Get numeric data
p2_numeric <- p2_raw[-1, -1]
p2_matrix <- as.matrix(p2_numeric)
mode(p2_matrix) <- "numeric"

# Create P2 dataframe
p2_data <- data.frame(
  Sample_ID = p2_samples,
  p2_matrix,
  stringsAsFactors = FALSE
)
colnames(p2_data) <- c("Sample_ID", p2_cytokines)

cat("\nData loaded successfully!\n")
cat("C2 data dimensions:", dim(c2_data), "\n")
cat("P2 data dimensions:", dim(p2_data), "\n")

# Create groups manually using direct string operations
c2_groups <- character(length(c2_samples))
for (i in 1:length(c2_samples)) {
  if (substr(c2_samples[i], 1, 1) == "S") {
    c2_groups[i] <- "Group_S"
  } else {
    c2_groups[i] <- "Group_U"
  }
}

p2_groups <- character(length(p2_samples))
for (i in 1:length(p2_samples)) {
  if (substr(p2_samples[i], 1, 1) == "S") {
    p2_groups[i] <- "Group_S"
  } else {
    p2_groups[i] <- "Group_U"
  }
}

cat("\nGroup assignments:\n")
cat("C2 - Group S:", sum(c2_groups == "Group_S"), ", Group U:", sum(c2_groups == "Group_U"), "\n")
cat("P2 - Group S:", sum(p2_groups == "Group_S"), ", Group U:", sum(p2_groups == "Group_U"), "\n")

# Function to perform statistical analysis
analyze_cytokines <- function(data, groups, cytokine_names, dataset_name) {
  cat("\nAnalyzing", dataset_name, "cytokines...\n")
  
  results <- list()
  
  for (i in 1:length(cytokine_names)) {
    cytokine <- cytokine_names[i]
    values <- data[[cytokine]]
    
    # Split by groups
    group_s_values <- values[groups == "Group_S"]
    group_u_values <- values[groups == "Group_U"]
    
    # Remove NA values
    group_s_values <- group_s_values[!is.na(group_s_values)]
    group_u_values <- group_u_values[!is.na(group_u_values)]
    
    if (length(group_s_values) > 0 && length(group_u_values) > 0) {
      # Perform t-test
      t_test <- t.test(group_s_values, group_u_values)
      
      # Calculate effect size
      pooled_sd <- sqrt(((length(group_s_values) - 1) * var(group_s_values) + 
                         (length(group_u_values) - 1) * var(group_u_values)) / 
                        (length(group_s_values) + length(group_u_values) - 2))
      cohens_d <- (mean(group_s_values) - mean(group_u_values)) / pooled_sd
      
      results[[cytokine]] <- list(
        mean_group_s = mean(group_s_values),
        mean_group_u = mean(group_u_values),
        sd_group_s = sd(group_s_values),
        sd_group_u = sd(group_u_values),
        p_value = t_test$p.value,
        t_statistic = t_test$statistic,
        cohens_d = cohens_d,
        significant = t_test$p.value < 0.05
      )
      
      cat(sprintf("  %s: p = %.4f, Effect size = %.3f, Significant: %s\n", 
                  cytokine, t_test$p.value, cohens_d, ifelse(t_test$p.value < 0.05, "YES", "NO")))
    } else {
      cat("  ", cytokine, ": Insufficient data\n")
    }
  }
  
  return(results)
}

# Perform analysis
c2_stats <- analyze_cytokines(c2_data, c2_groups, c2_cytokines, "C2")
p2_stats <- analyze_cytokines(p2_data, p2_groups, p2_cytokines, "P2")

# Create summary tables
create_summary_table <- function(stats, dataset_name) {
  if (length(stats) == 0) return(NULL)
  
  summary_df <- data.frame(
    Cytokine = names(stats),
    Mean_Group_S = sapply(stats, function(x) round(x$mean_group_s, 3)),
    Mean_Group_U = sapply(stats, function(x) round(x$mean_group_u, 3)),
    SD_Group_S = sapply(stats, function(x) round(x$sd_group_s, 3)),
    SD_Group_U = sapply(stats, function(x) round(x$sd_group_u, 3)),
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

# Create and save summary tables
c2_summary <- create_summary_table(c2_stats, "C2")
p2_summary <- create_summary_table(p2_stats, "P2")

write.csv(c2_summary, "output/c2_cytokine_summary.csv", row.names = FALSE)
write.csv(p2_summary, "output/p2_cytokine_summary.csv", row.names = FALSE)

# Create visualizations
cat("\nCreating visualizations...\n")

# Prepare data for plotting
c2_plot_data <- data.frame(
  Sample_ID = c2_data$Sample_ID,
  Group = c2_groups,
  c2_data[, -1],
  stringsAsFactors = FALSE
)

p2_plot_data <- data.frame(
  Sample_ID = p2_data$Sample_ID,
  Group = p2_groups,
  p2_data[, -1],
  stringsAsFactors = FALSE
)

# Create boxplots
create_boxplots <- function(data, group_col, cytokine_names, title) {
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
         y = "Expression Level")
  
  return(p)
}

# Create and save boxplots
c2_boxplot <- create_boxplots(c2_plot_data, "Group", c2_cytokines, "C2 Cytokine Expression by Group")
ggsave("plots/c2_cytokine_boxplots.png", c2_boxplot, width = 12, height = 8, dpi = 300)

p2_boxplot <- create_boxplots(p2_plot_data, "Group", p2_cytokines, "P2 Cytokine Expression by Group")
ggsave("plots/p2_cytokine_boxplots.png", p2_boxplot, width = 12, height = 8, dpi = 300)

# Create correlation plots
create_correlation_plot <- function(data, cytokine_names, title) {
  # Extract numeric data
  numeric_data <- data[, cytokine_names]
  
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

create_correlation_plot(c2_data, c2_cytokines, "C2 Cytokine Correlation Matrix")
create_correlation_plot(p2_data, p2_cytokines, "P2 Cytokine Correlation Matrix")

# Create summary report
cat("\nGenerating summary report...\n")

report_content <- paste0('
RNA/Cytokine Analysis Summary Report
====================================

Analysis Date: ', Sys.Date(), '
Data Files: 11-19-20-C2.xlsx, 11-19-20-P2.xlsx

Sample Groups:
- Group S: Samples S001-S008 (Treatment group 1)
- Group U: Samples U001-U026 (Treatment group 2)

C2 Cytokines Analyzed: ', paste(c2_cytokines, collapse = ", "), '

P2 Cytokines Analyzed: ', paste(p2_cytokines, collapse = ", "), '

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

