#!/bin/python3





import pandas as pd
from random import choice
import random




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
		

		


variant_class = 'all_variants'



#======================
# Estonia
#======================

filename = 'D:/Dropbox/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/Estonian_Summary.xlsx'
samp_df = pd.read_excel(filename)
samp_df = samp_df[samp_df.Relationship == 'P']



variant_classes = ['rare_deleterious_snvs', 'dups', 'dels', 'strs_std_2']
variant_df = pd.DataFrame()
for variant_class in variant_classes:
	app = pd.read_csv('variants/{}.csv'.format(variant_class))
	app['Variant_class'] = variant_class
	variant_df = variant_df.append(app)


variant_df = variant_df[variant_df.Sample.isin(samp_df.Sample.to_list())]




summ = []

cohort = 'Estonia'
first_hit = '16p12'

print(cohort, first_hit)

genes = variant_df['Gene'].to_list()
genes = list(set(genes))
genes = [s for s in genes if ';' not in s]
genes = [s for s in genes if s in gene_universe]

counts = count_degree_bins(genes)
app = ['all_variants', cohort, first_hit, 'True', counts['0-18'], counts['19-95'], counts['96-329'], counts['330-Max']]
summ.append(app)

for i in range(1000):
	if i % 10 == 0:
		print(i)
	random_genes = get_random_genes(len(genes))
	counts = count_degree_bins(random_genes)
	app = ['all_variants', cohort, first_hit, i, counts['0-18'], counts['19-95'], counts['96-329'], counts['330-Max']]
	summ.append(app)
			


summ = pd.DataFrame(summ, columns=['Variant_class', 'Cohort', 'First_hit', 'Simulation', '0-18', '19-95', '96-329', '330-Max'])
summ.to_csv('statistics/simulations_{}.csv'.format(cohort), index=False)




#======================
# 16p12
#======================

filename = 'D:/Dropbox/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx'
samp_df = pd.read_excel(filename)
samp_df = samp_df[samp_df.Relationship == 'P']



variant_classes = ['rare_deleterious_snvs', 'dups', 'dels', 'strs_std_2']
variant_df = pd.DataFrame()
for variant_class in variant_classes:
	app = pd.read_csv('D:/Dropbox/16p12.2 project/Human patients project/WGS paper/19_Variant_GO_enrichment/figure_for_all_cohorts_v5/variants/16p12/{}.csv'.format(variant_class))
	app['Variant_class'] = variant_class
	variant_df = variant_df.append(app)


variant_df = variant_df[variant_df.Sample.isin(samp_df.Sample.to_list())]


summ = []

cohort = 'girirajan'
first_hit = '16p12'

print(cohort, first_hit)

genes = variant_df['Gene'].to_list()
genes = list(set(genes))
genes = [s for s in genes if ';' not in s]
genes = [s for s in genes if s in gene_universe]

counts = count_degree_bins(genes)
app = ['all_variants', cohort, first_hit, 'True', counts['0-18'], counts['19-95'], counts['96-329'], counts['330-Max']]
summ.append(app)

for i in range(1000):
	if i % 10 == 0:
		print(i)
	random_genes = get_random_genes(len(genes))
	counts = count_degree_bins(random_genes)
	app = ['all_variants', cohort, first_hit, i, counts['0-18'], counts['19-95'], counts['96-329'], counts['330-Max']]
	summ.append(app)
			


summ = pd.DataFrame(summ, columns=['Variant_class', 'Cohort', 'First_hit', 'Simulation', '0-18', '19-95', '96-329', '330-Max'])
summ.to_csv('statistics/simulations_{}.csv'.format(cohort), index=False)



#======================
# SSC
#======================

first_hits = ['large_rare_duplications','large_rare_deletions','dbd_tier1_snvs']
variant_classes = ['rare_deleterious_snvs', 'dups', 'dels', 'strs_std_2']
variant_class_map = {
	'rare_deleterious_snvs':'rare_deleterious_snvs',
	'dups':'duplications',
	'dels':'deletions',
	'strs_std_2':'strs'
}

cohort = 'SSC'

for first_hit in first_hits:
	filename = f'D:/Dropbox/16p12.2 project/Human patients project/WGS paper/31_SSC analysis/Genotype-phenotype integration/1_variant_preparation/samples/{first_hit}.csv'
	samp_df = pd.read_csv(filename)

	variant_df = pd.DataFrame()
	for variant_class in variant_classes:
		app = pd.read_csv(f'D:/Dropbox/16p12.2 project/Human patients project/WGS paper/31_SSC analysis/Genotype-phenotype integration/1_variant_preparation/variants/{first_hit}/{variant_class_map[variant_class]}.csv')
		app['Variant_class'] = variant_class
		variant_df = variant_df.append(app)


	variant_df = variant_df[variant_df.Sample.isin(samp_df.Sample.to_list())]
	
	summ = []
	samples = samp_df['Sample'].to_list()

	genes = variant_df[variant_df.Sample.isin(samples)]['Gene'].to_list()
	genes = list(set(genes))
	genes = [s for s in genes if ';' not in s]
	genes = [s for s in genes if s in gene_universe]
	
	counts = count_degree_bins(genes)
	app = ['all_variants', cohort, first_hit, 'True', counts['0-18'], counts['19-95'], counts['96-329'], counts['330-Max']]
	summ.append(app)
	
	for i in range(1000):
		if i % 10 == 0:
			print(i)
		random_genes = get_random_genes(len(genes))
		counts = count_degree_bins(random_genes)
		app = ['all_variants', cohort, first_hit, i, counts['0-18'], counts['19-95'], counts['96-329'], counts['330-Max']]
		summ.append(app)
	
	
	summ = pd.DataFrame(summ, columns=['Variant_class', 'Cohort', 'First_hit', 'Simulation', '0-18', '19-95', '96-329', '330-Max'])
	
	print(cohort, first_hit)
	print(summ)
	summ.to_csv(f'statistics/simulations_{cohort}_{first_hit}.csv', index=False)


	




#======================
# SVIP
#======================

filename = 'D:/Dropbox/16p12.2 project/Human patients project/WGS paper/19_Variant_GO_enrichment/figure_for_all_cohorts_v5/cohorts/SVIP_samples.csv'
samp_df = pd.read_csv(filename)
samp_df = samp_df[samp_df.Relationship == 'Proband']



variant_classes = ['rare_deleterious_snvs', 'dups', 'dels']
variant_df = pd.DataFrame()
for variant_class in variant_classes:
	app = pd.read_csv('D:/Dropbox/16p12.2 project/Human patients project/WGS paper/19_Variant_GO_enrichment/figure_for_all_cohorts_v5/variants/SVIP/{}.csv'.format(variant_class))
	app['Variant_class'] = variant_class
	variant_df = variant_df.append(app)
variant_df = variant_df[variant_df.Sample.isin(samp_df.Sample.to_list())]

summ = []

cohort = 'SVIP'
first_hits = ['16p-deletion', '16p-duplication']
for first_hit in first_hits:
	print(cohort, first_hit)
	samples = samp_df[samp_df['Family type'] == first_hit]['Sample'].to_list()

	genes = variant_df[variant_df.Sample.isin(samples)]['Gene'].to_list()
	genes = list(set(genes))
	genes = [s for s in genes if ';' not in s]
	genes = [s for s in genes if s in gene_universe]
	
	counts = count_degree_bins(genes)
	app = ['all_variants', cohort, first_hit, 'True', counts['0-18'], counts['19-95'], counts['96-329'], counts['330-Max']]
	summ.append(app)

	for i in range(1000):
		if i % 10 == 0:
			print(i)
		random_genes = get_random_genes(len(genes))
		counts = count_degree_bins(random_genes)
		app = ['all_variants', cohort, first_hit, i, counts['0-18'], counts['19-95'], counts['96-329'], counts['330-Max']]
		summ.append(app)


summ = pd.DataFrame(summ, columns=['Variant_class', 'Cohort', 'First_hit', 'Simulation', '0-18', '19-95', '96-329', '330-Max'])
summ.to_csv('statistics/simulations_{}.csv'.format(cohort), index=False)


#======================



df1 = pd.read_csv('statistics/simulations_{}.csv'.format('girirajan'))
df2 = pd.read_csv('statistics/simulations_{}.csv'.format('SVIP'))
df3 = pd.read_csv('statistics/simulations_SSC_large_rare_duplications.csv')
df4 = pd.read_csv('statistics/simulations_SSC_large_rare_deletions.csv')
df5 = pd.read_csv('statistics/simulations_SSC_dbd_tier1_snvs.csv')
df6 = pd.read_csv('statistics/simulations_Estonia.csv')


summ = df1.append(df2)
summ = summ.append(df3)
summ = summ.append(df4)
summ = summ.append(df5)
summ = summ.append(df6)

summ.to_csv('statistics/simulations_all_cohorts.csv', index=False)















