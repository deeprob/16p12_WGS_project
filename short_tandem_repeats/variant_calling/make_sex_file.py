#!/bin/python3


import pandas as pd

# load in pedigree
ped = pd.read_csv('all_batches.ped', sep='\t', header=None)

# columns in ped file:
#     Family ID
#     Individual ID
#     Paternal ID
#     Maternal ID
#     Sex (1=male; 2=female; other=unknown)
#     Phenotype

def get_sex(s):
# 	if s == 1:
# 		return 'M'
# 	if s == 2:
# 		return 'F'
	return s
	
# make sex column
ped['sex'] = ped[4].apply(get_sex)

# keep only sample and sex columns
ped = ped[[1, 'sex']].copy()

# rename columns
ped.columns = ['sample', 'pedigree']

# save to file
ped.to_csv('sample_sex.txt', sep='\t', index=False)

