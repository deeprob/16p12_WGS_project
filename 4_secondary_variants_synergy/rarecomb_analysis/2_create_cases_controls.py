#!/bin/python3




import pandas as pd
import subprocess



subprocess.run('mkdir cases_controls', shell=True)




dropbox='~/Dropbox'


samples_df = pd.read_excel(dropbox + '/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx')
samples_df = samples_df[samples_df.Relationship == 'P']
samples_df = samples_df[~samples_df.Missense_CADD25.isna()]






phenotypes = ['Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial']

phenotype2split_value = {
	'Child_ID_DD':2,
	'Child_behav':1,
	'Child_psych':0,
	'Child_nervous_system':0,
	'Child_congenital':1,
	'Child_craniofacial':2,
}



for phenotype in phenotypes:
	subdf = samples_df[~samples_df[phenotype].isna()].copy()
	
	# split at split value
	samples_with_phenotype = subdf[subdf[phenotype] > phenotype2split_value[phenotype]].Sample.to_list()
	samples_control = subdf[subdf[phenotype] <= phenotype2split_value[phenotype]].Sample.to_list()
	
	
	with open(f'cases_controls/{phenotype}_cases.txt', 'w') as f:
		for sample in samples_with_phenotype:
			f.write(sample + '\n')
	
	with open(f'cases_controls/{phenotype}_controls.txt', 'w') as f:
		for sample in samples_control:
			f.write(sample + '\n')
	
	
	
	

















