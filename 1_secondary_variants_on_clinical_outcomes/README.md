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
*1_anova.R*: TODO: find the location
*2_annova_heatmap.R*: TODO find the location
*3_FH_t-test.Rmd*: Significance test to compare rare variant burden in probands with family history and no family history. 
*4_heatmap_ttest.R*: Heatmap of above. Figure reference (3C, S2E)
*11_ttest_family_history.R*: TODO: find the location
*12_heatmap_ttest_family_history.R*: TODO: find the location

# Neuronal Effects
Enrichment of secondary variant genes in various other studies or database. 

## Disease gene set enrichment
*0_disease_gene_set_enrichment.Rmd*: Check enrichment of secondary variants in annotated disease genes for autism, schizophrenia etc.  Figure reference (S3A).


## Brain tissue stage enrichment
*1_expression_timeline_annotations.py*: Get the rare variant genes expressed at different developmental points in brain tissue.

*2_fishers_exact.Rmd*: Calculates fisher's exact OR and p-values for each reltaionship-region-time-type set.

*3_plot_enrichment.py*: Plot enrichment of different variant types by regions.

*4_proband_enrichments_only.py*: Plot enrichment of different variant types by regions in probands only. Figure reference (S3B).

## Neuronal class enrichment
*1_download.sh*: Download sc-rnaseq data from allen brain institute

*2_prepare_scrna_data.py*: binarizes sc-rna data after normalizing

*3a_fishers_exact.py*: run fishers exact to check for enrichment of secondary variants in python

*3b_fishers_exact_r_subclass.R*: run fishers exact to check for enrichment of neuronal subclass on secondary variants in R

3b_fishers_exact_r_superclass.R: run fishers exact to check for enrichment of neuronal superclass on secondary variants in R

*3b_fishers_exact_r.R*: run fishers exact to check for enrichment of secondary variants in R

*4_all_coding_only.py*: run fishers exact to check for enrichment of secondary variants in python for coding variants only

*4a_make_barplot.py: Bar plot of odds ration

*4b_make_forest_plot_combined.py*: Forest plot of above. Figure reference (3D)

*4b_make_forest_plot.py*: Forest plot of above

## LCL DE Enrichment
*0_rnaseq_enrichment.Rmd*: Secondary variant gene enrichment within LCL DE genes. Figure reference (3E)

## Gene network enrichment
*1_download.sh*: Download data for 3q29 network.

*2_get_network_metrics.py*: Get network degree of the loeuf genes. 

*3_figure_degree.py*: Plot the network degree distribution of loeuf genes. 

*4_prepare_variants_estonia.py*: Prepare variants for the estonia cohort.

*5a_simulation_by_phenotype.py*: Simulate the distribution of genes falling within various network degree bins by phenotype. 

*5b_simlulations_by_cohort.py*: Simulate the distribution of genes falling within various network degree bins by cohort.

*6_simulation_ztest.py*: Simulate the distribution of genes falling within various network degree bins by cohort and add z-scores.

*7_simulation_boxplots.py*: Boxplots of the simulation for all cohort. Figure reference (3F, S3C)

