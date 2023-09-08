import pandas as pd
import numpy as np

# gnomad.v2.1.1.lof_metrics.by_gene.tsv from Dropbox/16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/gnomad_annotations
calls = pd.read_csv('final_calls/calls_by_gene.txt', sep = '\t')

gnomad_df = pd.read_csv('../gnomad.v2.1.1.lof_metrics.by_gene.tsv', sep='\t')

# Get LOEUF score from gnomAD file
# If a gene has 2 gene IDs, use the one with the lower score
print('LOEUF')
def get_loeuf(row):
	# Get IDs
	ids = row['Gene_ID'].split(' ')

	# Restrict to only IDs in gnomAD
	ids2 = [i for i in ids if i in gnomad_df.gene_id.to_list()]

	# No IDs are in gnomAD list
	if len(ids2)==0:
		print(ids)
		return '.'

	loeuf = min(gnomad_df[gnomad_df.gene_id.isin(ids2)]['oe_lof_upper'].to_numpy())
	return loeuf

calls['LOEUF'] = calls.apply(get_loeuf, axis = 1)
print(calls)

# genemap2.txt is from Dropbox\16p12.2 project\Human patients project\WGS paper\5_SNV pipelines-annotations\OMIM Annotations
omim_df = pd.read_csv('../genemap2.txt', sep='\t', comment='#', header=None,
	names = ['Chromosome', 'Genomic Position Start', 'Genomic Position End', 'Cyto Location', 'Computed Cyto Location', 'MIM Number',
		'Gene Symbols', 'Gene Name', 'Approved Gene Symbol', 'Entrez Gene ID', 'Ensembl Gene ID', 'Comments', 'Phenotypes', 'Mouse Gene Symbol/ID'])
ensembl_genes_in_omim = omim_df['Ensembl Gene ID'].to_list()
print('OMIM')
def get_mim(row):
	# Get IDs
	ids = row['Gene_ID'].split(' ')

	# Restrict to only IDs in the MIM set
	ids2 = [i for i in ids if i in ensembl_genes_in_omim]

	# No IDs in MIM
	if len(ids2)==0:
		print(ids)
		return '.'

	mim_nums = omim_df[omim_df["Ensembl Gene ID"].isin(ids)]['MIM Number'].to_list()
	mim = ' '.join([str(i) for i in mim_nums])
	return mim

calls['MIM_number'] = calls.apply(get_mim, axis = 1)
print(calls)

calls.to_csv('final_calls/calls_by_gene_anno.txt', sep = '\t', index = False)
