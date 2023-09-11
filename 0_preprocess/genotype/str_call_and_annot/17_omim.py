#!/bin/python3


# genemap2.txt is from Dropbox\16p12.2 project\Human patients project\WGS paper\5_SNV pipelines-annotations\OMIM Annotations


import pandas as pd



df = pd.read_csv('16p12_cohort.expansions_2SD.annotated.exonic.gene_ids.column_names.loeuf.csv')



omim_df = pd.read_csv('genemap2.txt', sep='\t', comment='#', header=None)
omim_columns = ['Chromosome', 'Genomic Position Start', 'Genomic Position End', 'Cyto Location', 'Computed Cyto Location', 'MIM Number', 'Gene Symbols', 'Gene Name', 'Approved Gene Symbol', 'Entrez Gene ID', 'Ensembl Gene ID', 'Comments', 'Phenotypes', 'Mouse Gene Symbol/ID']
omim_df.columns = omim_columns



ensembl_genes_in_omim = omim_df['Ensembl Gene ID'].to_list()


df['MIM_number'] = ''
for i, row in df.iterrows():
	gene_ids = row['Gene_id_'].split(';')
	gene_ids_in_omim = [s for s in gene_ids if s in ensembl_genes_in_omim]
	
	if len(gene_ids_in_omim) == 0:
		# none of the gene ids are in omim
		continue
	
	subomim_df = omim_df[omim_df['Ensembl Gene ID'].isin(gene_ids)]
	
	# a gene id may have more than one omim number
	# in which case annotate all seperated by semicolons
	# (600179;607206)
	omim_numbers = subomim_df['MIM Number'].to_list()
	omim_numbers = [str(s) for s in omim_numbers]
	df.at[i, 'MIM_number'] = ';'.join(omim_numbers)




df.to_csv('16p12_cohort.expansions_2SD.annotated.exonic.gene_ids.column_names.loeuf.omim.csv', index=False)










