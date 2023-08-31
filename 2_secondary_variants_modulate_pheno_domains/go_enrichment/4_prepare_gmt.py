#!/bin/python3


import pandas as pd
import numpy as np




input_filename = '../figure_for_all_cohorts_v2/Human_GOBP_AllPathways_no_GO_iea_November_01_2021_symbol.gmt'
cohorts = ['16p12','SSC_dbd_tier1_snvs','SSC_large_rare_deletions','SSC_large_rare_duplications','SVIP_16p-deletion','SVIP_16p-duplication']



for cohort in cohorts:
	enrichment_file = 'statistics/{}.csv'.format(cohort)
	print(cohort)
	
	df = pd.read_csv(enrichment_file)
	go_terms_to_include = df['term_label'].to_list()
	
	output_filename = 'cytoscape_files/{}.gmt'.format(cohort)

	fin = open(input_filename, 'r')
	fout = open(output_filename, 'w')
	terms_done = []
	
	fout.write('\n')

	for line in fin:
		go_term_label = line.split('\t')[1]
		
		if '%GOBP%' not in line.split('\t')[0]:
			continue
		
		if go_term_label in go_terms_to_include:
			print(line.split('\t')[0], go_term_label)
			fout.write(line)
			terms_done.append(go_term_label)
		
		


	fin.close()
	fout.close()

	terms_not_found = [s for s in go_terms_to_include if s not in terms_done]
	# print(len(terms_done))
	# print(terms_not_found)







