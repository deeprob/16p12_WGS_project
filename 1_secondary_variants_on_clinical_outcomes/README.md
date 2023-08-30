This directory holds the analysis scripts for Result 1 (*Differential effects of secondary variants towards clinical outcomes*) in the manuscript. It is further sub-divided into four folders based on the analysis category. Below is a description of each analysis category following the order in which they appear in the manuscript.

# Phenotypic difference  
Phenotype (both categorial and quantitative) comparison between probands and parents. Description of individual scripts are as follows:

*0_histogram_hrs_mat.r*: Smoothed histograms showing distributions of HRS-MAT in proband and parests. Figure reference (2C).

*1_histogram_srs.r*: Smoothed histograms showing distributions of SRS in proband and parests. Figure reference (2C). 

TODO: There are other figures whose script are needed in this paragraph.

# Variant burden difference
Detecting significant differences in rare variant burden. It is further divided into two subsections. Each subsection and their individual scripts are listed below.

## Power calculation
*0_avengeme.ipynb*: Common variant power calculation
 
*1_rare_burden_plot*.ipynb: Rare variant power calculation. Figure reference (S1C).

## Burden comparison
*0_paired_ttests_max_samples.R*: Global differences in variant burden in probands vs parents. 

*1_boxplots.py*: Box plots of above. Figure reference (3A, S2A)

*2_heatmap.py*: Heatmap of burden tests. Figure reference (3A)

*3_burden_test.R*: Comparing burden with Estonian Biobank. 

*4_boxplots.py*: Box plot of above. Figure reference (S2B)

*5_heatmap.R*: Heatmap of burden tests. TODO: Check which heatmap creation script is used.

# Mulitple Generation
Trends of increased phenotypic severity and accumulation of secondary variants over multiple generations of deletion carriers.

## Family history


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


