#!/bin/python3

import pandas as pd
import sys
import subprocess
import os

subprocess.run('mkdir statistics/5_combine', shell=True)


phenotypes = ['Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial']
combo_lens = [2, 3]



def convert_gene_name(gene):
	gene = gene[len('Input_'):]
	return gene



for combo_len in combo_lens:
	df = pd.DataFrame()
	for phenotype in phenotypes:
		if not os.path.isfile(f'statistics/4_filter/{phenotype}_{combo_len}.csv'):
			continue
		
		app = pd.read_csv(f'statistics/4_filter/{phenotype}_{combo_len}.csv')
		app['Phenotype'] = phenotype
		df = df.append(app)
		
	
	print(df)
	print(df.Phenotype.value_counts())
	
	# reorg columns
	columns = list(df.columns)
	columns = ['Phenotype'] + columns[:-1]
	df = df[columns]
	
	df.to_csv(f'statistics/5_combine/combolen{combo_len}.csv', index=False)









