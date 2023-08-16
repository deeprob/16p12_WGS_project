#!/bin/python3



import pandas as pd
import numpy as np

anno_df = pd.read_csv('~/Dropbox/16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/GENCODE annotations/gencode.v19.parsed.genes.csv')
anno_df = anno_df.drop_duplicates('gene_name')
anno_df = anno_df.set_index('gene_name')


loeuf_df = pd.read_csv('~/Dropbox/16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/gnomad_annotations/gnomad.v2.1.1.lof_metrics.by_gene.tsv', sep='\t')
loeuf_df = loeuf_df.drop_duplicates(['gene'])
loeuf_df = loeuf_df.set_index('gene')


# oe_lof_upper

def get_qunatile(value, lower_quantile, middle_quantile, upper_quantile):
	if value <= lower_quantile:
		return '0-{}'.format(int(lower_quantile))
	if value <= middle_quantile:
		return '{}-{}'.format(int(lower_quantile)+1, int(middle_quantile))
	if value <= upper_quantile:
		return '{}-{}'.format(int(middle_quantile)+1, int(upper_quantile))
	return '{}-Max'.format(int(upper_quantile)+1)


def split_into_4(df):
	df['Degree'] = df['Degree'].astype(float)
	values = df['Degree'].astype(float).to_list()
	
	lower_quantile = np.quantile(values, q=.25)
	middle_quantile = np.quantile(values, q=.5)
	upper_quantile = np.quantile(values, q=.75)
	
	df['Quartile'] = df['Degree'].apply(lambda s: get_qunatile(s, lower_quantile, middle_quantile, upper_quantile))
	print(df['Quartile'].value_counts())
	return df


def get_degree(gene):
	if gene not in counts:
		return 0
	return counts[gene]

def get_loeuf(gene):
	if gene not in loeuf_df.index:
		return ''
	if loeuf_df.at[gene, 'oe_lof_upper'] == '.':
		return ''
	return float(loeuf_df.at[gene, 'oe_lof_upper'])
	

net_df = pd.read_csv('networks/brain.degnorm-ge2.prob-gept02.dat', sep='\t', header=None, names=['gene1','gene2', 'weight'])
net_gene_df = pd.read_csv('networks/brain.genes.protein-coding.txt', sep='\t', header=None, names=['entrez', '1', 'symbol', 'chrom', 'start', 'stop'])


net_gene_df.entrez = net_gene_df.entrez.astype(str)
net_gene_df = net_gene_df.set_index('entrez', drop=False)


net_df.gene1 = net_df.gene1.astype(str)
net_df.gene2 = net_df.gene2.astype(str)
entrez_dict = net_gene_df['symbol'].to_dict()

net_df = net_df[net_df['gene1'].isin(net_gene_df['entrez'].to_list())]
net_df = net_df[net_df['gene2'].isin(net_gene_df['entrez'].to_list())]
net_df['gene1'] = net_df['gene1'].apply(lambda s: entrez_dict[s])
net_df['gene2'] = net_df['gene2'].apply(lambda s: entrez_dict[s])

print(net_df)

counts = pd.Series(net_df['gene1'].to_list() + net_df['gene2'].to_list())
counts = counts.value_counts()
print(counts)


genes = list(counts.index)
genes = list(set(genes))

stats = []
for gene in genes:
	if ';' in gene:
		continue
	degree = get_degree(gene)
	loeuf = get_loeuf(gene)
	app = [gene, degree, loeuf]
	stats.append(app)


stats = pd.DataFrame(stats, columns=['Gene', 'Degree', 'LOEUF'])
stats = split_into_4(stats)
print(stats)
stats.to_csv('statistics/network_degree_loeuf.csv', index=False)





