#!/usr/bin/env python3

# Python-based RNA/Cytokine Analysis Script
# Analysis of treatment effects on cytokine expression
# Author: RNA Analysis Pipeline
# Date: 2024

import pandas as pd
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Set up output directories
if not os.path.exists('output'):
    os.makedirs('output')
if not os.path.exists('plots'):
    os.makedirs('plots')

print("Loading and processing data...")

# Read C2 data
c2_raw = pd.read_excel('11-19-20-C2.xlsx', skiprows=2)
print(f"C2 raw data dimensions: {c2_raw.shape}")

# Get cytokine names from first row
c2_cytokines = c2_raw.iloc[0, 1:].tolist()
print(f"C2 cytokines: {', '.join(c2_cytokines)}")

# Get sample IDs from first column, excluding header
c2_samples = c2_raw.iloc[1:, 0].tolist()
# Remove "Row Labels" if present
if c2_samples[0] == 'Row Labels':
    c2_samples = c2_samples[1:]
    c2_raw = c2_raw.iloc[1:]

# Get numeric data
c2_numeric = c2_raw.iloc[1:, 1:]
c2_matrix = c2_numeric.astype(float)

# Create C2 dataframe
c2_data = pd.DataFrame(c2_matrix.values, index=c2_samples, columns=c2_cytokines)
c2_data.index.name = 'Sample_ID'

# Read P2 data
p2_raw = pd.read_excel('11-19-20-P2.xlsx', skiprows=2)
print(f"P2 raw data dimensions: {p2_raw.shape}")

# Get cytokine names from first row
p2_cytokines = p2_raw.iloc[0, 1:].tolist()
print(f"P2 cytokines: {', '.join(p2_cytokines)}")

# Get sample IDs from first column, excluding header
p2_samples = p2_raw.iloc[1:, 0].tolist()
# Remove "Row Labels" if present
if p2_samples[0] == 'Row Labels':
    p2_samples = p2_samples[1:]
    p2_raw = p2_raw.iloc[1:]

# Get numeric data
p2_numeric = p2_raw.iloc[1:, 1:]
p2_matrix = p2_numeric.astype(float)

# Create P2 dataframe
p2_data = pd.DataFrame(p2_matrix.values, index=p2_samples, columns=p2_cytokines)
p2_data.index.name = 'Sample_ID'

print(f"\nData loaded successfully!")
print(f"C2 data dimensions: {c2_data.shape}")
print(f"P2 data dimensions: {p2_data.shape}")

# Create groups manually using direct string operations
c2_groups = []
for sample_id in c2_samples:
    if str(sample_id).startswith('S'):
        c2_groups.append('Group_S')
    else:
        c2_groups.append('Group_U')

p2_groups = []
for sample_id in p2_samples:
    if str(sample_id).startswith('S'):
        p2_groups.append('Group_S')
    else:
        p2_groups.append('Group_U')

print(f"\nGroup assignments:")
print(f"C2 - Group S: {sum(1 for g in c2_groups if g == 'Group_S')}, Group U: {sum(1 for g in c2_groups if g == 'Group_U')}")
print(f"P2 - Group S: {sum(1 for g in p2_groups if g == 'Group_S')}, Group U: {sum(1 for g in p2_groups if g == 'Group_U')}")

# Function to perform statistical analysis
def analyze_cytokines(data, groups, cytokine_names, dataset_name):
    print(f"\nAnalyzing {dataset_name} cytokines...")
    
    results = {}
    
    for cytokine in cytokine_names:
        values = data[cytokine].dropna()
        
        # Create a mapping between the values index and groups
        value_groups = []
        for sample_id in values.index:
            if str(sample_id).startswith('S'):
                value_groups.append('Group_S')
            else:
                value_groups.append('Group_U')
        
        # Split by groups
        group_s_values = values[np.array(value_groups) == 'Group_S']
        group_u_values = values[np.array(value_groups) == 'Group_U']
        
        if len(group_s_values) > 0 and len(group_u_values) > 0:
            # Perform t-test
            t_stat, p_value = stats.ttest_ind(group_s_values, group_u_values)
            
            # Calculate effect size (Cohen's d)
            pooled_std = np.sqrt(((len(group_s_values) - 1) * group_s_values.var() + 
                                 (len(group_u_values) - 1) * group_u_values.var()) / 
                                (len(group_s_values) + len(group_u_values) - 2))
            cohens_d = (group_s_values.mean() - group_u_values.mean()) / pooled_std
            
            results[cytokine] = {
                'mean_group_s': group_s_values.mean(),
                'mean_group_u': group_u_values.mean(),
                'sd_group_s': group_s_values.std(),
                'sd_group_u': group_u_values.std(),
                'p_value': p_value,
                't_statistic': t_stat,
                'cohens_d': cohens_d,
                'significant': p_value < 0.05
            }
            
            print(f"  {cytokine}: p = {p_value:.4f}, Effect size = {cohens_d:.3f}, Significant: {'YES' if p_value < 0.05 else 'NO'}")
        else:
            print(f"  {cytokine}: Insufficient data")
    
    return results

# Perform analysis
c2_stats = analyze_cytokines(c2_data, c2_groups, c2_cytokines, "C2")
p2_stats = analyze_cytokines(p2_data, p2_groups, p2_cytokines, "P2")

# Create summary tables
def create_summary_table(stats, dataset_name):
    if not stats:
        return None
    
    summary_data = []
    for cytokine, stat in stats.items():
        summary_data.append({
            'Cytokine': cytokine,
            'Mean_Group_S': round(stat['mean_group_s'], 3),
            'Mean_Group_U': round(stat['mean_group_u'], 3),
            'SD_Group_S': round(stat['sd_group_s'], 3),
            'SD_Group_U': round(stat['sd_group_u'], 3),
            'P_Value': f"{stat['p_value']:.3e}",
            'T_Statistic': round(stat['t_statistic'], 3),
            'Cohens_D': round(stat['cohens_d'], 3),
            'Significant': 'Yes' if stat['significant'] else 'No'
        })
    
    summary_df = pd.DataFrame(summary_data)
    # Sort by p-value
    summary_df = summary_df.sort_values('P_Value')
    
    return summary_df

# Create and save summary tables
c2_summary = create_summary_table(c2_stats, "C2")
p2_summary = create_summary_table(p2_stats, "P2")

if c2_summary is not None:
    c2_summary.to_csv('output/c2_cytokine_summary.csv', index=False)
if p2_summary is not None:
    p2_summary.to_csv('output/p2_cytokine_summary.csv', index=False)

# Create visualizations
print("\nCreating visualizations...")

# Prepare data for plotting
c2_plot_data = c2_data.copy()
c2_plot_data['Group'] = c2_groups
c2_plot_data = c2_plot_data.reset_index()

p2_plot_data = p2_data.copy()
p2_plot_data['Group'] = p2_groups
p2_plot_data = p2_plot_data.reset_index()

# Create boxplots
def create_boxplots(data, cytokine_names, title, filename):
    # Melt data for plotting
    plot_data = data.melt(id_vars=['Sample_ID', 'Group'], 
                          value_vars=cytokine_names,
                          var_name='Cytokine', value_name='Expression')
    
    # Remove NA values
    plot_data = plot_data.dropna()
    
    # Create boxplot
    plt.figure(figsize=(12, 8))
    sns.boxplot(data=plot_data, x='Cytokine', y='Expression', hue='Group')
    plt.title(title)
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(f'plots/{filename}.png', dpi=300, bbox_inches='tight')
    plt.close()

# Create and save boxplots
create_boxplots(c2_plot_data, c2_cytokines, 'C2 Cytokine Expression by Group', 'c2_cytokine_boxplots')
create_boxplots(p2_plot_data, p2_cytokines, 'P2 Cytokine Expression by Group', 'p2_cytokine_boxplots')

# Create correlation plots
def create_correlation_plot(data, cytokine_names, title, filename):
    # Extract numeric data
    numeric_data = data[cytokine_names]
    
    # Calculate correlation matrix
    cor_matrix = numeric_data.corr()
    
    # Create correlation plot
    plt.figure(figsize=(10, 8))
    sns.heatmap(cor_matrix, annot=True, cmap='RdBu_r', center=0, 
                square=True, fmt='.2f', cbar_kws={'shrink': 0.8})
    plt.title(title)
    plt.tight_layout()
    plt.savefig(f'plots/{filename}.png', dpi=300, bbox_inches='tight')
    plt.close()

create_correlation_plot(c2_data, c2_cytokines, 'C2 Cytokine Correlation Matrix', 'c2_cytokine_correlation_matrix')
create_correlation_plot(p2_data, p2_cytokines, 'P2 Cytokine Correlation Matrix', 'p2_cytokine_correlation_matrix')

# Create summary report
print("\nGenerating summary report...")

report_content = f"""
RNA/Cytokine Analysis Summary Report
====================================

Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d')}
Data Files: 11-19-20-C2.xlsx, 11-19-20-P2.xlsx

Sample Groups:
- Group S: Samples S001-S008 (Treatment group 1)
- Group U: Samples U001-U026 (Treatment group 2)

C2 Cytokines Analyzed: {', '.join(c2_cytokines)}

P2 Cytokines Analyzed: {', '.join(p2_cytokines)}

Statistical Analysis:
- Test: Two-sample t-test
- Effect size: Cohen's d
- Significance threshold: p < 0.05

Key Findings:
============

C2 Cytokines - Significantly Different (p < 0.05):
"""

if c2_stats:
    for cytokine, stat in c2_stats.items():
        if stat['significant']:
            report_content += f"\n- {cytokine}: p = {stat['p_value']:.3e}, Effect size = {stat['cohens_d']:.3f}"

report_content += "\n\nP2 Cytokines - Significantly Different (p < 0.05):"

if p2_stats:
    for cytokine, stat in p2_stats.items():
        if stat['significant']:
            report_content += f"\n- {cytokine}: p = {stat['p_value']:.3e}, Effect size = {stat['cohens_d']:.3f}"

report_content += """

Files Generated:
- output/c2_cytokine_summary.csv
- output/p2_cytokine_summary.csv
- plots/c2_cytokine_boxplots.png
- plots/p2_cytokine_boxplots.png
- plots/c2_cytokine_correlation_matrix.png
- plots/p2_cytokine_correlation_matrix.png

Analysis complete!
"""

with open('output/analysis_summary.txt', 'w') as f:
    f.write(report_content)

# Print final summary
print("\n" + "="*50)
print("ANALYSIS COMPLETE!")
print("="*50 + "\n")

print("Files generated:")
print("- output/c2_cytokine_summary.csv")
print("- output/p2_cytokine_summary.csv")
print("- output/analysis_summary.txt")
print("- plots/ (visualization files)\n")

print("Key Results:")
print("C2 Cytokines - Significantly different:")
if c2_stats:
    for cytokine, stat in c2_stats.items():
        if stat['significant']:
            print(f"  {cytokine}: p = {stat['p_value']:.4f}, Effect size = {stat['cohens_d']:.3f}")

print("\nP2 Cytokines - Significantly different:")
if p2_stats:
    for cytokine, stat in p2_stats.items():
        if stat['significant']:
            print(f"  {cytokine}: p = {stat['p_value']:.4f}, Effect size = {stat['cohens_d']:.3f}")

print("\nCheck 'output/analysis_summary.txt' for a complete summary!")

