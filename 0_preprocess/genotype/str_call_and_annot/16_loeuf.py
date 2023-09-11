#!/bin/python3




import pandas as pd




df = pd.read_csv('16p12_cohort.expansions_2SD.annotated.exonic.gene_ids.column_names.csv')



gnomad_df = pd.read_csv('gnomad.v2.1.1.lof_metrics.by_gene.txt', sep='\t')
gnomad_df = gnomad_df.set_index('gene_id')


# remove dot ENSG00000215910.3 > ENSG00000215910
df['Gene_id_'] = df['Gene_id'].apply(lambda s: s.split('.')[0])



df['LOEUF'] = ''
for i, row in df.iterrows():
	gene_id = row['Gene_id_']
	if gene_id not in gnomad_df.index:
		continue
	df.at[i, 'LOEUF'] = gnomad_df.at[gene_id, 'oe_lof_upper']




df.to_csv('16p12_cohort.expansions_2SD.annotated.exonic.gene_ids.column_names.loeuf.csv', index=False)










