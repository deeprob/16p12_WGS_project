#!/bin/python3





import pandas as pd
from random import choice
import random


filename = '~/Dropbox/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx'
samp_df = pd.read_excel(filename)
samp_df = samp_df[samp_df.Relationship == 'P']



phenotypes = ['Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial']


variant_classes = ['rare_deleterious_snvs', 'dups', 'dels', 'strs_std_2']
variant_df = pd.DataFrame()
for variant_class in variant_classes:
	app = pd.read_csv('~/Dropbox/16p12.2 project/Human patients project/WGS paper/19_Variant_GO_enrichment/figure_for_16p12_cohort/variants/{}.csv'.format(variant_class))
	app['Variant_class'] = variant_class
	variant_df = variant_df.append(app)

print(variant_df.columns)

variant_class = 'all_variants'




net_df = pd.read_csv('statistics/network_degree_loeuf.csv')
net_df = net_df.set_index('Gene')
gene_universe = list(net_df.index)

def count_degree_bins(genes):
	df = pd.DataFrame(genes, columns=['Gene'])
	df = df[df['Gene'].isin(net_df.index)]
	df['Degree'] = df['Gene'].apply(lambda s: net_df.at[s, 'Degree'])
	df['Quartile'] = df['Gene'].apply(lambda s: net_df.at[s, 'Quartile'])
	
	counts = df['Quartile'].value_counts()
	# count_df.columns = ['Quartile', 'Count']
	
	return counts



def get_random_genes(num_genes):
	return random.sample(gene_universe, num_genes)
	# remaining_genes = gene_universe.copy()
	# genes = []
	# for i in range(num_genes):
	# 	gene = choice(remaining_genes)
	# 	genes.append(gene)
	# 	remaining_genes.remove(gene)
	# return genes
		


summ = []


for phenotype in phenotypes:
	print(phenotype)

	# split at 0
	samples_with_phenotype = samp_df[(~samp_df[phenotype].isna()) & (samp_df[phenotype] > 0)].Sample.to_list()
	
	genes = variant_df[(variant_df.Sample.isin(samples_with_phenotype))]['Gene'].to_list()
	genes = list(set(genes))
	genes = [s for s in genes if ';' not in s]
	genes = [s for s in genes if s in gene_universe]
	print(len(genes))
	
	counts = count_degree_bins(genes)
	app = [variant_class, phenotype, 'True', counts['0-18'], counts['19-95'], counts['96-329'], counts['330-Max']]
	summ.append(app)
	
	for i in range(1000):
		if i % 10 == 0:
			print(i)
		random_genes = get_random_genes(len(genes))
		counts = count_degree_bins(random_genes)
		app = [variant_class, phenotype, i, counts['0-18'], counts['19-95'], counts['96-329'], counts['330-Max']]
		summ.append(app)
			

summ = pd.DataFrame(summ, columns=['Variant_class', 'Phenotype', 'Simulation', '0-18', '19-95', '96-329', '330-Max'])

summ.to_csv('statistics/simulations_by_phenotype.csv', index=False)



















