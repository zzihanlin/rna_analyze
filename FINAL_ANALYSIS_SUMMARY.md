# RNA/Cytokine Analysis - Final Results Summary

## Analysis Overview
**Date:** August 26, 2025  
**Data Files:** 11-19-20-C2.xlsx, 11-19-20-P2.xlsx  
**Analysis Method:** Python-based statistical analysis with R backup scripts  

## Sample Groups
- **Group S (Treatment Group 1):** 8 samples (S001-S008)
- **Group U (Treatment Group 2):** 26 samples (U001-U026)

## Statistical Methods
- **Test:** Two-sample t-test (independent samples)
- **Effect Size:** Cohen's d
- **Significance Threshold:** p < 0.05
- **Multiple Testing:** No correction applied

## Key Findings

### C2 Cytokines - Significantly Different (p < 0.05)

| Cytokine | P-Value | Effect Size (Cohen's d) | Group S Mean | Group U Mean | Interpretation |
|----------|---------|------------------------|--------------|--------------|----------------|
| **IL-15** | 0.0119 | **1.145** | 7,481.77 | 110.42 | Large increase in Group S |
| **IL-17A/F** | 0.0099 | **1.171** | - | - | Large increase in Group S |
| **IL-27p28/IL-30** | 0.0097 | **1.173** | - | - | Large increase in Group S |
| **IL-9** | 0.0150 | **1.119** | 709.75 | 5.79 | Large increase in Group S |
| **MCP-1** | 0.0050 | **-1.286** | 95.02 | 341.70 | Large decrease in Group S |

### P2 Cytokines - Significantly Different (p < 0.05)

| Cytokine | P-Value | Effect Size (Cohen's d) | Group S Mean | Group U Mean | Interpretation |
|----------|---------|------------------------|--------------|--------------|----------------|
| **IFN-γ** | 0.0115 | **1.144** | 145.57 | 3.86 | Large increase in Group S |
| **IL-10** | 0.0109 | **1.154** | 505.20 | 5.77 | Large increase in Group S |
| **IL-12p70** | 0.0104 | **1.162** | 5,954.72 | 23.98 | Large increase in Group S |
| **IL-2** | 0.0101 | **1.167** | 463.26 | 2.51 | Large increase in Group S |
| **IL-4** | 0.0099 | **1.170** | - | - | Large increase in Group S |
| **IL-5** | 0.0100 | **1.169** | - | - | Large increase in Group S |
| **TNF-α** | 0.0331 | **0.950** | - | - | Large increase in Group S |

## Effect Size Interpretation
- **Cohen's d > 0.8:** Large effect
- **Cohen's d 0.5-0.8:** Medium effect  
- **Cohen's d < 0.5:** Small effect

## Biological Significance

### C2 Cytokines
- **Pro-inflammatory cytokines (IL-15, IL-17A/F, IL-9):** Significantly elevated in Group S
- **MCP-1:** Significantly reduced in Group S (anti-inflammatory effect)
- **IL-33, IP-10, MIP-1α, MIP-2:** No significant differences

### P2 Cytokines  
- **Th1 cytokines (IFN-γ, IL-12p70, IL-2):** Significantly elevated in Group S
- **Th2 cytokines (IL-4, IL-5, IL-10):** Significantly elevated in Group S
- **Pro-inflammatory (TNF-α):** Significantly elevated in Group S
- **IL-1β, IL-6, KC/GRO:** No significant differences

## Treatment Effect Summary
**Group S shows a strong pro-inflammatory response** with significant increases in:
- Multiple interleukins (IL-2, IL-4, IL-5, IL-9, IL-10, IL-12p70, IL-15, IL-17A/F, IL-27p28/IL-30)
- Interferon-gamma (IFN-γ)
- Tumor necrosis factor-alpha (TNF-α)

**Group S shows reduced MCP-1 levels**, suggesting potential anti-inflammatory modulation.

## Files Generated

### Data Files
- `output/c2_cytokine_summary.csv` - Detailed C2 cytokine statistics
- `output/p2_cytokine_summary.csv` - Detailed P2 cytokine statistics  
- `output/analysis_summary.txt` - Complete analysis report

### Visualizations
- `plots/c2_cytokine_boxplots.png` - C2 cytokine expression by group
- `plots/p2_cytokine_boxplots.png` - P2 cytokine expression by group
- `plots/c2_cytokine_correlation_matrix.png` - C2 cytokine correlations
- `plots/p2_cytokine_correlation_matrix.png` - P2 cytokine correlations

## Technical Notes
- **Data Processing:** Successfully handled Excel file structure with proper group assignment
- **Statistical Analysis:** Used Python with scipy.stats for robust t-test calculations
- **Effect Size:** Cohen's d calculated using pooled standard deviation
- **Data Quality:** All samples included, no missing data issues

## Conclusions
The analysis reveals **significant treatment effects** on cytokine expression profiles:

1. **Group S demonstrates enhanced immune activation** with elevated levels of multiple pro-inflammatory and regulatory cytokines
2. **Effect sizes are consistently large** (Cohen's d > 1.0) for most significant differences
3. **Treatment appears to modulate both Th1 and Th2 responses** in Group S
4. **MCP-1 reduction in Group S** suggests potential anti-inflammatory counter-regulation

This pattern suggests that the treatment in Group S induces a **comprehensive immune response** that may be beneficial for the intended therapeutic outcome.

---
*Analysis completed using Python-based pipeline with statistical validation*
