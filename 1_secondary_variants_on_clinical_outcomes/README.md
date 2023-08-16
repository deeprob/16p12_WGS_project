This folder holds the scripts and data of Result-1 (Differential effects of secondary variants towards clinical outcomes) in the manuscript.

# Description of scripts
The scripts are divided into various subfolders based on the type of analysis and the order in which they appear in the manuscript.

## Phenotypic Difference  


## Power Calculation
0_avengeme.ipynb: Common variant power calculation
 
1_rare_burden_plot.ipynb: Rare variant power calculation


## Variant Burden Difference
0_paired_ttests_max_samples.R: Global differences in variant burden in probands vs parents

1_boxplots.py: Box plots of above

2_heatmap.py: Heatmap of burden tests

3_burden_test.R: Comparing burden with Estonian Biobank

4_boxplots.py: Box plot of above

5_heatmap.R: Heatmap of burden tests

## Mulitple Generation

## Neuronal Effects
0_disease_gene_set_enrichment.Rmd: Neuronal gene set enrichment of secondary variants

1_download_allen_data.sh: Download sc-rnaseq data from allen brain institute

2_prepare_scrna_data.py: script to prepare sc-rna data

3a_fishers_exact.py: run fishers exact to check for enrichment of secondary variants in python

3b_fishers_exact_r_subclass.R: run fishers exact to check for enrichment of neuronal subclass on secondary variants in R

3b_fishers_exact_r_superclass.R: run fishers exact to check for enrichment of neuronal superclass on secondary variants in R

3b_fishers_exact_r.R: run fishers exact to check for enrichment of secondary variants in R

4_all_coding_only.py: run fishers exact to check for enrichment of secondary variants in python for coding variants only

4a_make_barplot.py: Bar plot of odds ration

4b_make_forest_plot_combined.py: Forest plot of above

4b_make_forest_plot.py: Forest plot of above


