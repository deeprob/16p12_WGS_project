#!/bin/python3

import pandas as pd
import numpy as np

pd.set_option('display.max_columns', 20)

#---------------------
# Load in data
#---------------------

filename = '~/Dropbox/16p12.2 project/Human patients project/WGS/Structural Variation Analysis/SV2 calls final rare coding.txt'
df = pd.read_csv(filename, sep='\t')


filename = 'all_batches.ped'
ped = pd.read_csv(filename, sep='\t', header=None)
ped.columns = ['Family', 'Sample', 'Father', 'Mother', 'Sex', 'Carrier_status']
ped = ped.set_index('Sample')


#---------------------
# Exclude bad samples
#---------------------
exclude_samples = ['SG138', 'SG473', 'SG474', 'SG475']

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
	
	# all structural variants for the child, mother, and father
	child_df  = df[df['sample'] == samp]
	mother_df = df[df['sample'] == mother]
	father_df = df[df['sample'] == father]
	
	
	# create numpy structures for end, start, and svlength
	# (numpy is faster than pandas)
	mother_ends      = mother_df['end'].to_numpy()
	mother_starts    = mother_df['start'].to_numpy()
	mother_svlengths = mother_df['svlength'].to_numpy()
	mother_svtypes   = mother_df['svtype'].to_numpy()
	mother_chroms    = mother_df['chrom'].to_numpy()
	
	father_ends      = father_df['end'].to_numpy()
	father_starts    = father_df['start'].to_numpy()
	father_svlengths = father_df['svlength'].to_numpy()
	father_svtypes   = father_df['svtype'].to_numpy()
	father_chroms    = father_df['chrom'].to_numpy()

	for i, row in child_df.iterrows():
		inh_from_father = False
		inh_from_mother = False
		inheritence = '.'
		
		start    = row['start']
		end      = row['end']
		svlength = row['svlength']
		svtype   = row['svtype']
		chrom    = row['chrom']
		
		min_end    = np.minimum(mother_ends, end)
		max_start  = np.maximum(mother_starts, start)
		max_length = np.maximum(mother_svlengths, svlength)
		
		overlap_mother = (min_end - max_start) > .5 * max_length
		overlap_mother = overlap_mother & (mother_chroms == chrom) & (mother_svtypes == svtype)
		if overlap_mother.sum() > 0:
			inh_from_mother = True
			
		min_end    = np.minimum(father_ends, end)
		max_start  = np.maximum(father_starts, start)
		max_length = np.maximum(father_svlengths, svlength)
		
		overlap_father = (min_end - max_start) > .5 * max_length
		overlap_father = overlap_father & (father_chroms == chrom) & (father_svtypes == svtype)
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
	
	


filename = 'SV2 calls final rare coding.txt'
df.to_csv(filename, sep='\t', index=False)
	

