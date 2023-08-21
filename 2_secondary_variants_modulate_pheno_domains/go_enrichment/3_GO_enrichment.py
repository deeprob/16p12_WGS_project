#!/bin/python3


import pandas as pd
import numpy as np
from panther_api import *




annotation_dataset = 'ANNOT_TYPE_ID_PANTHER_GO_SLIM_BP'
variant_class = 'all_variants'


#======================
# 16p12
#======================

df = pd.read_csv('cohorts/16p12_samples.csv')
df = df[df.Relationship == 'P']

print(df)



variant_classes = ['rare_deleterious_snvs', 'dups', 'dels', 'strs_std_2']
variant_df = pd.DataFrame()
for variant_class in variant_classes:
	app = pd.read_csv('variants/16p12/{}.csv'.format(variant_class))
	variant_df = variant_df.append(app)


variant_df = variant_df[variant_df.Sample.isin(df.Sample.to_list())]


print(variant_df)


genes = variant_df['Gene'].to_list()
genes = list(set(genes))

print(len(genes))


go_result = panther(genes, annotation_dataset=annotation_dataset)
go_result = go_result[(go_result['fdr'] < 0.05) & (go_result['plus_minus'] == '+')]
output_filename = 'statistics/16p12.csv'
go_result.to_csv(output_filename, index=False)







#======================
# SSC
#======================

	

first_hits = ['large_rare_duplications','large_rare_deletions','dbd_tier1_snvs']

for first_hit in first_hits:
	df = pd.read_csv(f'cohorts/SSC_{first_hit}.csv', header=None)
	df.columns = ['Sample']

	variant_classes = ['rare_deleterious_snvs', 'dups', 'dels', 'strs_std_2']
	variant_df = pd.DataFrame()
	for variant_class in variant_classes:
		app = pd.read_csv(f'variants/SSC/{first_hit}/{variant_class}.csv')
		variant_df = variant_df.append(app)
	variant_df = variant_df[variant_df.Sample.isin(df.Sample.to_list())]


	samples = df['Sample'].to_list()
		
	genes = variant_df[variant_df.Sample.isin(samples)]['Gene'].to_list()
	genes = list(set(genes))
	print(len(samples), len(genes))


	go_result = panther(genes, annotation_dataset=annotation_dataset)
	go_result = go_result[(go_result['fdr'] < 0.05) & (go_result['plus_minus'] == '+')]
	output_filename = 'statistics/SSC_{}.csv'.format(first_hit)
	go_result.to_csv(output_filename, index=False)


#======================
# SVIP
#======================


df = pd.read_csv('cohorts/SVIP_samples.csv')
df = df[df.Relationship == 'Proband']

print(df)


variant_classes = ['rare_deleterious_snvs', 'dups', 'dels']
variant_df = pd.DataFrame()
for variant_class in variant_classes:
	app = pd.read_csv('variants/SVIP/{}.csv'.format(variant_class))
	variant_df = variant_df.append(app)
	print(variant_class)
	print(variant_df)
variant_df = variant_df[variant_df.Sample.isin(df.Sample.to_list())]


print(variant_df)


second_hits = ['16p-deletion', '16p-duplication']

for second_hit in second_hits:
	samples = df[df['Family type'] == second_hit]['Sample'].to_list()
	
	genes = variant_df[variant_df.Sample.isin(samples)]['Gene'].to_list()
	genes = list(set(genes))
	print(len(samples), len(genes))
	




	go_result = panther(genes, annotation_dataset=annotation_dataset)
	go_result = go_result[(go_result['fdr'] < 0.05) & (go_result['plus_minus'] == '+')]
	output_filename = 'statistics/SVIP_{}.csv'.format(second_hit)
	go_result.to_csv(output_filename, index=False)

