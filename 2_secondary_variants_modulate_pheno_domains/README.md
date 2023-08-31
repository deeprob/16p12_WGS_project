This directory holds the analysis scripts for Result 2 (*Classes of secondary variants jointly modulate specific developmental domains*) in the manuscript. It is further sub-divided into three folders based on the analysis category. Below is a description of each analysis category following the order in which they appear in the manuscript.


# Variant traits model
*1_models.Rmd*: Logistic regression  and linear regression models with the secondary variants as predictors and qualitative/quantitative phenotypic domains as predicted variable. 

*2_plot_values.py*: Plot model coefficients and variance values for logistic and linear regression. Figure reference (4A, 4B)

*3_forest_plot_together.py*: Get forest plots for individual predictors. Figure reference (S4A)

# Burden tests
*0_ttest_phenotypes.R*: Burden test to check for significance between phenotypic domains and predictor values. 
 
*1_heatmap_ttest.R*: Heatmap of above. Figure reference (S4B)

# Variant quant trait correlation
*0_correlation_variants_v_phenotypes.R*: Pearson's correlation between variant classes and quant traits.

*1_heatmap_cor_variants_v_phenotypes.R*: Heatmap of above. Figure reference (4C)

*3_scatter_plots.py*: Scatter plot between individual variants and traits. Figure reference (4D)

# Go enrichment
*1_prepare_variants_16p12.py*: Prepare the genes for go enrichment. 

*2_GO_enrichment.py*: Conduct GO enrichment using panther.

*3_download.sh*: Download gmt file for cytoscape figure.

*4_prepare_gmt.py*: Prepare gmt file for cytoscape figure. Figure reference (4E)


# Disease gene set enrichment by phenotype
This script is present in 1_secondary_variants_on_clinical_outcomes -> neuronal_effects -> disease_gene_set_enrichment.
*0_disease_gene_set_enrichment.Rmd*: Enrichment of disease genes by phenotypic domains. Figure reference (4F)
