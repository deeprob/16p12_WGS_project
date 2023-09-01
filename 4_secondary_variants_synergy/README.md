# Description
This directory holds the analysis scripts for Result 4 (*Secondary variants synergistically modify clinical phenotypes of deletion carriers*) in the manuscript. It is further sub-divided into four folders based on the analysis category. 

Below is a description of each analysis category following the order in which they appear in the manuscript.

# Other annotated pathogenic variants
*1_identify_pathogenic_SNVs.py*: Flag pathogenic variants annotated in other databases which are present amongst genes with secondary variants.

*2_update_pathogenic_SNVs.py*: Manually review a few. 

*3_upset_plot.py*: Upset plot of the overlap. Figure reference (6A)

# HPO analysis
*1_gene2hpo.py*: Match genes to hpo phenotypes

*2_convert_entrez.Rmd*: Convert entrez id to ensemble id.

*3_format_obo.py*: Parse hpo fomrat to human readable format. 

*4_annotate_MC.txt*: Match hpo terms to current cohort phenos manually.

*5_make_relationships.py*: Add relationships between child and parents 

*6_annotate_variants.py*: Annotate the variants with hpo terms.

*7_annotate_phenotypes.py*: Annotate each person with hpo phenotypes.

*8_make_tableS3.py*: Make a table with sample information about genes, variants and hpo terms. 

*9_find_variability.py*: Get the number of genes in our cohort that are associated with different phenotypes in differnet probands

*10_explained_bygene.py*: Get the genes in our cohort that are already associated with other phenotypes

*11_variant_type_proportions.py*: Get the proportion of phenotypes explained by each class of variants

*12_variant_type_boxplots.py*: Box plot of variant type explaining phenotype proportion. Figure reference (6B)

*13_inheritance_variant_type.py*: Get the proportion of phenotypes explained by each class of variants inherited from each parent

*14_variant_type_inheritance_boxplots.py*: Boxplots showing the proportion of variants explained in probands from different inheritances

*15_proband_phenotype_heatmap.py*: Heatmaps showing proband phenotypes and which can be explained by their genes. Figure reference (S6A, S6B) 

*16_proband_gene_heatmap.py*: Make heatmaps showing probands and the genes they have. Figure reference (6C, S6C)

*17_gene_phenotype_lookup.py*: Quickly lookup the phenotypes explained by a gene for a proband

*18_ClinVar_variants.py*: Get ClinVar variants in probands for manual HPO annotation

*19_ClinVar_only.py*:  Heatmap showing phenotypes associated with genes with ClinVar variants in each proband

*20_get_phenos.py*: Get the phenotypes associated with the ClinVar variant

# RVGDT analysis
*0_preprocess_for_rvgdt.sh*: Exclude rare variants inside structural variant regions for each individual, then convert to the format for RV-GDT

*1_run_rvgdt_each_pheno.sh*: Run RVGDT analysis separately for each phenotypic domain.

*2_analyze.py*: Analyze RVGDT results to find the consensus significant genes across several runs. Figure reference (S6D).

# RareComb analysis
*0_download.sh*: downloads the data to cluster.

*1_prepare_variants.py*: prepare variants in rarecomb format.

*2_create_cases_controls.py*: Create case and control based on presence of phenotype. 

*3_rarecomb.R*: rarecomb script. 

*3_rarecomb.sh*: Run rarecomb on slurm.

*4_filter_significance.py*: Filter combinations by bonferroni.

*5_combine.py*: Combine multiple phenotypes into one file

*6_figure_heatmap.py*: Heatmap of rarecomb results. Figure reference (6D)
