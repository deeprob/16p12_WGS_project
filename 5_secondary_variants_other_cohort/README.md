# Description
This directory holds the analysis scripts for Result 5 (*Contributions of secondary variants to disease in other primary variant carriers*) in the manuscript. It is further sub-divided into four folders based on the analysis category. 

1. [Linear models](#linear-models)
2. [Correlation](#correlation)
3. [Variant analysis](#variant-analysis)
4. [Enrichment](#enrichment)

Below is a description of each analysis category following the order in which they appear in the manuscript.

# Linear models
*1_linear_regression_ssc.R*: Linear regression models with the secondary variants as predictors and qualitative/quantitative phenotypic domains as predicted variable for SSC cohort.

*1_linear_regression_svip.R*: Linear regression models with the secondary variants as predictors and qualitative/quantitative phenotypic domains as predicted variable for SVIP cohort.

*2_integrate_ssc.py*: Integrate model metrics and correlation information from ssc. 

*2_integrate_svip.py*: Integrate model metrics and correlation information from svip. 

*3_organize_stats_ssc.py*: Rename cohort and phenotypes for visualization. 

*4_merge_stats.py*: Heatmap of model coefficients. Figure reference (7A).

TODO: Forest plots Fig S7A, S8A not sure which one to use. 

# Correlation
## Variant phenotype
*1_correlation_phenotypes_ssc.R*: Correlation between variants and phenotypes across ssc sub-cohorts.

*1_correlation_phenotypes_svip.R*: Correlation between variants and phenotypes across svip sub-cohorts.

*2_heatmap_phenotypes_ssc.R*: Heatmap of correlation between variants and phenotypes across ssc sub-cohorts. Figure reference (S8B)

*2_heatmap_phenotypes_svip.R*: Heatmap of correlation between variants and phenotypes across SVIP sub-cohorts. Figure reference (S7B).

## Variant variant
*2_ltm_correlation.py*: Correlation between coding SNV burden and PRS. Figure reference (7B)

*3_correlation_variants_ssc.R*: Correlation between variants and phenotypes across ssc sub-cohorts.

*3_correlation_variants_svip.R*: Correlation between variants and phenotypes across svip sub-cohorts.

*4_heatmap_variants_ssc.R*: Heatmap of correlation between variants and phenotypes across ssc sub-cohorts. Figure reference (S8C)

*4_heatmap_variants_svip.R*: Heatmap of correlation between variants and phenotypes across SVIP sub-cohorts. Figure reference (S7C).

# Variant analysis
*1_heatmap_ssc.py*: Variance heatmap for SSC cohort.

# Enrichment
Enrichment scripts for all are located in "/2_secondary_variants_modulate_pheno_domains/go_enrichment". Figure reference (7C, 7D).
