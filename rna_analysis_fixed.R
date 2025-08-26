#!/usr/bin/env Rscript

# RNA/Cytokine Analysis Script - Fixed Version
# Analysis of treatment effects on cytokine expression
# Author: RNA Analysis Pipeline
# Date: 2024

# Load required libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(pheatmap)
  library(ggplot2)
  library(ggpubr)
  library(corrplot)
  library(factoextra)
  library(DT)
  library(knitr)
  library(kableExtra)
  library(reshape2)
  library(viridis)
  library(RColorBrewer)
})

# Set working directory and create output folder
setwd("/Users/suwa/Desktop/rna_analyze-2")
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("plots")) dir.create("plots")

# Function to clean Excel data - Fixed for actual data structure
clean_excel_data <- function(file_path) {
  # Read the Excel file with proper headers
  df <- read_excel(file_path, skip = 2)
  
  # The first row contains cytokine names, second row contains sample IDs
  cytokine_names <- as.character(df[1, -1])  # Exclude first column (Row Labels)
  sample_ids <- as.character(df[-1, 1])      # First column contains sample IDs
  
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

# Function to identify sample groups
identify_groups <- function(data) {
  # Extract sample IDs from first column
  sample_ids <- data$Sample_ID
  
  # Create group assignments based on sample naming convention
  # S001-S008: Treatment group 1, U001-U008: Treatment group 2
  groups <- ifelse(grepl("^S", sample_ids), "Group_S", "Group_U")
  
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
    }
  }
  
  return(results)
}

# Perform statistical analysis
cat("\nPerforming statistical analysis...\n")
c2_stats <- perform_statistical_analysis(c2_analysis)
p2_stats <- perform_statistical_analysis(p2_analysis)

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

# Function to create heatmaps
create_heatmap <- function(data, group_col = "Group", title = "Cytokine Expression Heatmap") {
  # Prepare data for heatmap
  sample_ids <- data$Sample_ID
  groups <- data[[group_col]]
  
  # Extract numeric data
  numeric_data <- data[, sapply(data, is.numeric)]
  
  # Create annotation dataframe
  annotation_df <- data.frame(
    Group = groups,
    row.names = sample_ids
  )
  
  # Create color palette
  colors <- brewer.pal(8, "Set2")
  
  # Create heatmap
  pheatmap(as.matrix(numeric_data),
           annotation_row = annotation_df,
           annotation_colors = list(Group = setNames(colors[1:2], unique(groups))),
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           main = title,
           fontsize = 10,
           fontsize_row = 8,
           fontsize_col = 8,
           filename = paste0("plots/", gsub(" ", "_", tolower(title)), ".png"),
           width = 10,
           height = 8)
}

# Create heatmaps
create_heatmap(c2_analysis, title = "C2 Cytokine Expression Heatmap")
create_heatmap(p2_analysis, title = "P2 Cytokine Expression Heatmap")

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

# Function to create PCA plots
create_pca_plot <- function(data, group_col = "Group", title = "PCA Plot") {
  # Extract numeric data
  numeric_data <- data[, sapply(data, is.numeric)]
  
  # Remove columns with all NA values
  numeric_data <- numeric_data[, colSums(is.na(numeric_data)) < nrow(numeric_data)]
  
  # Perform PCA
  pca_result <- prcomp(numeric_data, scale. = TRUE, na.action = na.omit)
  
  # Create PCA plot
  p <- fviz_pca_ind(pca_result,
                     col.ind = data[[group_col]],
                     palette = "Set2",
                     addEllipses = TRUE,
                     ellipse.level = 0.95,
                     geom = "point",
                     title = title) +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  return(p)
}

# Create PCA plots
c2_pca <- create_pca_plot(c2_analysis, title = "C2 Cytokines PCA Plot")
ggsave("plots/c2_pca_plot.png", c2_pca, width = 10, height = 8, dpi = 300)

p2_pca <- create_pca_plot(p2_analysis, title = "P2 Cytokines PCA Plot")
ggsave("plots/p2_pca_plot.png", p2_pca, width = 10, height = 8, dpi = 300)

# Function to create volcano plots
create_volcano_plot <- function(stats, title = "Volcano Plot") {
  if (length(stats) == 0) return(NULL)
  
  # Create volcano plot data
  volcano_data <- data.frame(
    Cytokine = names(stats),
    Log2FC = sapply(stats, function(x) log2(x$mean_group1 / x$mean_group2)),
    NegLog10P = -log10(sapply(stats, function(x) x$p_value)),
    Significant = sapply(stats, function(x) x$significant),
    stringsAsFactors = FALSE
  )
  
  # Create volcano plot
  p <- ggplot(volcano_data, aes(x = Log2FC, y = NegLog10P, color = Significant)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text(aes(label = ifelse(Significant, Cytokine, "")), 
              hjust = 0, vjust = 0, size = 3) +
    scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "red"),
                       name = "Significant (p < 0.05)") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          legend.position = "bottom") +
    labs(title = title,
         x = "Log2 Fold Change (Group_S / Group_U)",
         y = "-Log10 P-value",
         caption = "Red dashed line: p = 0.05 threshold")
  
  return(p)
}

# Create volcano plots
c2_volcano <- create_volcano_plot(c2_stats, "C2 Cytokines Volcano Plot")
if (!is.null(c2_volcano)) {
  ggsave("plots/c2_volcano_plot.png", c2_volcano, width = 10, height = 8, dpi = 300)
}

p2_volcano <- create_volcano_plot(p2_stats, "P2 Cytokines Volcano Plot")
if (!is.null(p2_volcano)) {
  ggsave("plots/p2_volcano_plot.png", p2_volcano, width = 10, height = 8, dpi = 300)
}

# Create HTML report
cat("\nGenerating HTML report...\n")

html_content <- paste0('
<!DOCTYPE html>
<html>
<head>
    <title>RNA/Cytokine Analysis Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; }
        h1, h2, h3 { color: #2c3e50; }
        .section { margin: 30px 0; padding: 20px; border-left: 4px solid #3498db; }
        .summary-table { width: 100%; border-collapse: collapse; margin: 20px 0; }
        .summary-table th, .summary-table td { border: 1px solid #ddd; padding: 8px; text-align: left; }
        .summary-table th { background-color: #f2f2f2; }
        .plot { text-align: center; margin: 20px 0; }
        .plot img { max-width: 100%; height: auto; }
        .significant { color: #e74c3c; font-weight: bold; }
        .not-significant { color: #7f8c8d; }
    </style>
</head>
<body>
    <h1>RNA/Cytokine Analysis Report</h1>
    
    <div class="section">
        <h2>Analysis Overview</h2>
        <p>This report presents the analysis of cytokine expression data comparing two treatment groups:</p>
        <ul>
            <li><strong>Group S:</strong> Samples S001-S008</li>
            <li><strong>Group U:</strong> Samples U001-U008</li>
        </ul>
        <p>Analysis includes statistical comparisons, visualizations, and identification of differentially expressed cytokines.</p>
    </div>
    
    <div class="section">
        <h2>C2 Cytokines Analysis</h2>
        <p>The following cytokines were analyzed: ', paste(colnames(c2_data)[-1], collapse = ", "), '</p>
        
        <h3>Statistical Summary</h3>
        <table class="summary-table">
            <tr>
                <th>Cytokine</th>
                <th>Mean Group S</th>
                <th>Mean Group U</th>
                <th>P-value</th>
                <th>Effect Size (Cohen\'s d)</th>
                <th>Significant</th>
            </tr>',
            paste(sapply(1:nrow(c2_summary), function(i) {
                row <- c2_summary[i, ]
                sig_class <- ifelse(row$Significant == "Yes", "significant", "not-significant")
                paste0('<tr class="', sig_class, '">
                    <td>', row$Cytokine, '</td>
                    <td>', row$Mean_Group_S, '</td>
                    <td>', row$Mean_Group_U, '</td>
                    <td>', row$P_Value, '</td>
                    <td>', row$Cohens_D, '</td>
                    <td>', row$Significant, '</td>
                </tr>')
            }), collapse = ""),
        '</table>
        
        <h3>Visualizations</h3>
        <div class="plot">
            <h4>Expression Boxplots</h4>
            <img src="plots/c2_cytokine_boxplots.png" alt="C2 Cytokine Boxplots">
        </div>
        
        <div class="plot">
            <h4>Expression Heatmap</h4>
            <img src="plots/c2_cytokine_expression_heatmap.png" alt="C2 Cytokine Heatmap">
        </div>
        
        <div class="plot">
            <h4>PCA Plot</h4>
            <img src="plots/c2_pca_plot.png" alt="C2 PCA Plot">
        </div>
        
        <div class="plot">
            <h4>Volcano Plot</h4>
            <img src="plots/c2_volcano_plot.png" alt="C2 Volcano Plot">
        </div>
    </div>
    
    <div class="section">
        <h2>P2 Cytokines Analysis</h2>
        <p>The following cytokines were analyzed: ', paste(colnames(p2_data)[-1], collapse = ", "), '</p>
        
        <h3>Statistical Summary</h3>
        <table class="summary-table">
            <tr>
                <th>Cytokine</th>
                <th>Mean Group S</th>
                <th>Mean Group U</th>
                <th>P-value</th>
                <th>Effect Size (Cohen\'s d)</th>
                <th>Significant</th>
            </tr>',
            paste(sapply(1:nrow(p2_summary), function(i) {
                row <- p2_summary[i, ]
                sig_class <- ifelse(row$Significant == "Yes", "significant", "not-significant")
                paste0('<tr class="', sig_class, '">
                    <td>', row$Cytokine, '</td>
                    <td>', row$Mean_Group_S, '</td>
                    <td>', row$Mean_Group_U, '</td>
                    <td>', row$P_Value, '</td>
                    <td>', row$Cohens_D, '</td>
                    <td>', row$Significant, '</td>
                </tr>')
            }), collapse = ""),
        '</table>
        
        <h3>Visualizations</h3>
        <div class="plot">
            <h4>Expression Boxplots</h4>
            <img src="plots/p2_cytokine_boxplots.png" alt="P2 Cytokine Boxplots">
        </div>
        
        <div class="plot">
            <h4>Expression Heatmap</h4>
            <img src="plots/p2_cytokine_expression_heatmap.png" alt="P2 Cytokine Heatmap">
        </div>
        
        <div class="plot">
            <h4>PCA Plot</h4>
            <img src="plots/p2_pca_plot.png" alt="P2 PCA Plot">
        </div>
        
        <div class="plot">
            <h4>Volcano Plot</h4>
            <img src="plots/p2_volcano_plot.png" alt="P2 Volcano Plot">
        </div>
    </div>
    
    <div class="section">
        <h2>Key Findings</h2>
        <h3>Significantly Different Cytokines (C2)</h3>
        <ul>',
            paste(sapply(names(c2_stats), function(cytokine) {
                if (c2_stats[[cytokine]]$significant) {
                    paste0('<li><strong>', cytokine, '</strong>: p = ', 
                           format.pval(c2_stats[[cytokine]]$p_value, digits = 3),
                           ', Effect size = ', round(c2_stats[[cytokine]]$cohens_d, 3), '</li>')
                }
            }), collapse = ""),
        '</ul>
        
        <h3>Significantly Different Cytokines (P2)</h3>
        <ul>',
            paste(sapply(names(p2_stats), function(cytokine) {
                if (p2_stats[[cytokine]]$significant) {
                    paste0('<li><strong>', cytokine, '</strong>: p = ', 
                           format.pval(p2_stats[[cytokine]]$p_value, digits = 3),
                           ', Effect size = ', round(p2_stats[[cytokine]]$cohens_d, 3), '</li>')
                }
            }), collapse = ""),
        '</ul>
    </div>
    
    <div class="section">
        <h2>Methods</h2>
        <p><strong>Statistical Analysis:</strong> Two-sample t-tests were performed to compare cytokine expression between groups. Effect sizes were calculated using Cohen\'s d.</p>
        <p><strong>Visualization:</strong> Boxplots show expression distributions, heatmaps display expression patterns, PCA plots show sample clustering, and volcano plots highlight significant differences.</p>
        <p><strong>Significance Threshold:</strong> p < 0.05 was used to identify significantly different cytokines.</p>
    </div>
    
    <div class="section">
        <h2>Files Generated</h2>
        <ul>
            <li>Statistical summary tables (CSV format)</li>
            <li>Expression boxplots (PNG format)</li>
            <li>Expression heatmaps (PNG format)</li>
            <li>PCA plots (PNG format)</li>
            <li>Volcano plots (PNG format)</li>
            <li>Correlation matrices (PNG format)</li>
        </ul>
    </div>
</body>
</html>')

# Save HTML report
writeLines(html_content, "output/rna_analysis_report.html")

# Print summary to console
cat("\n" + paste(rep("=", 50), collapse = "") + "\n")
cat("ANALYSIS COMPLETE!\n")
cat(paste(rep("=", 50), collapse = "") + "\n\n")

cat("Files generated:\n")
cat("- output/c2_cytokine_summary.csv\n")
cat("- output/p2_cytokine_summary.csv\n")
cat("- output/rna_analysis_report.html\n")
cat("- plots/ (various visualization files)\n\n")

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

cat("\nOpen 'output/rna_analysis_report.html' in your browser to view the complete report!\n")
