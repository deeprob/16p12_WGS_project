#!/bin/python3

import pandas as pd
import numpy as np

pd.set_option('display.max_columns', 20)

#---------------------
# Load in data
#---------------------


filename = '~/Dropbox/16p12.2 project/Human patients project/WGS/STR variants/16p12_cohort.expansions_2SD.annotated.protein_coding_and_nearby.xlsx'
df = pd.read_excel(filename)


filename = '~/Dropbox/16p12.2 project/Human patients project/WGS/Structural Variation Analysis/AT_analysis/inheritence_annotation/all_batches.ped'
ped = pd.read_csv(filename, sep='\t', header=None)
ped.columns = ['Family', 'Sample', 'Father', 'Mother', 'Sex', 'Carrier_status']
ped = ped.set_index('Sample')


#---------------------
# Exclude bad samples
#---------------------
exclude_samples = ['SG138', 'SG473', 'SG474', 'SG475',
'SG332', 'SG315','SG138', 'SG224','SG316', 'SG222',
'SG223', 'SG334', 'SG317', 'SG196', 'SG333']

df = df[~df['sample'].isin(exclude_samples)]

#---------------------
# annotate inheritence
# leave a '.' if both parents aren't in pedigree
#---------------------

samples = df['sample'].unique()

df['inheritence'] = '.'

for samp in samples:
	if samp not in ped.index:
		print('{} Sample not in pedigree'.format(samp))
		
	# get mother and father
	# if father or mother == 0 then skip
	father = ped.loc[samp, 'Father']
	mother = ped.loc[samp, 'Mother']
	
	if mother == '0' or father == '0':
		continue
		
	# all STRs for the child, mother, and father
	child_df  = df[df['sample'] == samp]
	mother_df = df[df['sample'] == mother]
	father_df = df[df['sample'] == father]
	
	# create numpy structures for end, start, and svlength
	# (numpy is faster than pandas)
	mother_starts    = mother_df['pos'].to_numpy()
	mother_chroms    = mother_df['chrom'].to_numpy()
	mother_longest   = mother_df['longest_allele'].to_numpy()
	
	father_starts    = father_df['pos'].to_numpy()
	father_chroms    = father_df['chrom'].to_numpy()
	father_longest   = father_df['longest_allele'].to_numpy()
	
	for i, row in child_df.iterrows():
		inh_from_father = False
		inh_from_mother = False
		inheritence = '.'
		
		start    = row['pos']
		chrom    = row['chrom']
		longest  = row['longest_allele']
		
		overlap_mother = (mother_chroms == chrom) & (mother_starts == start) & (mother_longest == longest)
		if overlap_mother.sum() > 0:
			inh_from_mother = True
		
		
		overlap_father = (father_chroms == chrom) & (father_starts == start) & (father_longest == longest)
		if overlap_father.sum() > 0:
			inh_from_father = True
		
		if inh_from_father and inh_from_mother:
			inheritence = 'both'
		elif inh_from_father:
			inheritence = 'father'
		elif inh_from_mother:
			inheritence = 'mother'
		else:
			inheritence = 'de novo'
		
		df.at[i, 'inheritence'] = inheritence
		


# save
filename = '16p12_cohort.expansions_2SD.annotated.protein_coding_and_nearby.inheritence.xlsx'
df.to_excel(filename, index=False)
		
		
