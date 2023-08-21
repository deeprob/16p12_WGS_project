#!/bin/python3

import pandas as pd
import sys
import os

phenotypes = ['Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial']
combo_lens = [2,3]



def convert_gene_name(gene):
	gene = gene[len('Input_'):]
	return gene


for phenotype in phenotypes:
	for combo_len in combo_lens:
		if not os.path.isfile(f'statistics/3_rarecomb/{phenotype}_{combo_len}.csv'):
			continue
		
		df = pd.read_csv(f'statistics/3_rarecomb/{phenotype}_{combo_len}.csv')
		for combo_num in range(1, combo_len+1):
			df[f'Item_{combo_num}'] = df[f'Item_{combo_num}'].apply(lambda s: convert_gene_name(s))
		

		df = df[df.Case_Adj_Pval_bonf <= 0.05]
		print(df)
		df.to_csv(f'statistics/4_filter/{phenotype}_{combo_len}.csv', index=False)









