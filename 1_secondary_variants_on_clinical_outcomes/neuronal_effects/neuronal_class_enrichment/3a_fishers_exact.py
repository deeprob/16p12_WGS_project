#!/bin/python3





import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import numpy as np


dropbox = '~/Dropbox'


# load in scrnaseq data
scrna_df = pd.read_csv('allen_institute_tables/cell_type_specific_expression.csv')
scrna_df = scrna_df.set_index('feature')
cell_types = list(scrna_df.columns)

# load in genes with second hits
variant_classes = ['rare_deleterious_snvs', 'dups', 'dels', 'strs_std_2']
variant_df = pd.DataFrame()
for variant_class in variant_classes:
	app = pd.read_csv(dropbox + '/16p12.2 project/Human patients project/WGS paper/19_Variant_GO_enrichment/figure_for_16p12_cohort/variants/{}.csv'.format(variant_class))
	app['variant_class'] = variant_class
	variant_df = variant_df.append(app)



# restrict to probands
master_df = pd.read_excel(dropbox + '/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx')
master_df = master_df.set_index('Sample', drop=False)
master_df = master_df[~master_df['Missense_CADD25'].isna()] # only want samples with WGS 
master_df = master_df[master_df['Relationship'] == 'P']
variant_df = variant_df[variant_df.Sample.isin(master_df.Sample)]


# use the intersection of genes in GENCODE and 
# genes in the scrnaseq table as the gene universe
gencode_df = pd.read_csv(dropbox + '/16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/GENCODE annotations/gencode.v19.parsed.genes.csv')
gencode_df = gencode_df[gencode_df.gene_type == 'protein_coding']
gencode_genes = list(set(gencode_df['gene_name']))
scrna_genes = list(scrna_df.index)
gene_universe = list(set(gencode_genes) & set(scrna_genes))

# restrict scrna_df and variant_df to gene universe
scrna_df = scrna_df.loc[gene_universe]
variant_df = variant_df[variant_df.Gene.isin(gene_universe)]


# for each variant class do a fishers exact test for each cell type
stats = []
for variant_class in variant_classes:
	sub_variant_df = variant_df[variant_df.variant_class == variant_class]
	for cell_type in cell_types:
		genes_with_second_hits = sub_variant_df['Gene'].to_list()
		genes_for_cell_type = list(scrna_df[scrna_df[cell_type] == 1].index)
				
		second_hit_and_pref = list(set(genes_with_second_hits) & set(genes_for_cell_type))
		second_hit_and_not_pref = list(set(genes_with_second_hits) - set(genes_for_cell_type))
		not_second_hit_and_pref = list(set(genes_for_cell_type) - set(genes_with_second_hits))
		not_second_hit_and_not_pref_tmp = list(set(gene_universe) - set(genes_with_second_hits))
		not_second_hit_and_not_pref = list(set(not_second_hit_and_not_pref_tmp) - set(genes_for_cell_type))
		
		contingency_table = [
			[len(second_hit_and_pref), len(second_hit_and_not_pref)],
			[len(not_second_hit_and_pref), len(not_second_hit_and_not_pref)]
		]
		
		
		result = fisher_exact(contingency_table)
		stats.append([variant_class, cell_type, result[0], result[1], len(second_hit_and_pref), len(second_hit_and_not_pref), len(not_second_hit_and_pref), len(not_second_hit_and_not_pref)])
	


stats = pd.DataFrame(stats, columns=['variant_class', 'cell_type', 'oddsratio', 'pvalue', 'second_hit_and_pref', 'second_hit_and_not_pref', 'not_second_hit_and_pref', 'not_second_hit_and_pref'])



# FDR correction for each variant class (BH correction)
stats['FDR'] = np.nan
for variant_class in variant_classes:
	sub_stats = stats[stats.variant_class == variant_class]
	sub_stats['FDR'] = multipletests(sub_stats['pvalue'], method='fdr_bh')[1]
	for i in sub_stats.index:
		stats.at[i, 'FDR'] = sub_stats.at[i, 'FDR']



print(stats.sort_values('FDR'))

stats.to_csv('statistics/fishers_exact.csv', index=False)

















