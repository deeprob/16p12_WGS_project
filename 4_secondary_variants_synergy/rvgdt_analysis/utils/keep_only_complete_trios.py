#!/bin/python3

# originally prepared by anastasia
# edited by deepro
# edits 
# added in and out dir prefix
# changed file "all_batches.ped" to "Mar2022.fam"
# .fam files have headers unlike .ped files, hence added header=0
# had to change all column calls to index loc calls to make it more generalizable

import pandas as pd

input_file_dir = "/data5/deepro/wgs_16p/rvgdt/data/input_files"
intermediate_file_dir = "/data5/deepro/wgs_16p/rvgdt/data/intermediate_files"


df = pd.read_csv(f'{input_file_dir}/Mar_2022.fam', sep='\t', header=0) # changed header=0


with open(f'{intermediate_file_dir}/tdt_sample_list.txt', 'r') as f:
	samples = f.readlines()

samples = [s.strip() for s in samples]

triodf = df.loc[(df.iloc[:, 2] != '0') & (df.iloc[:, 3] != '0')]


def all_have_wgs(s):
	s = s.to_list()
	if s[0] in samples and s[1] in samples and s[2] in samples:
		return True
	return False



triodf['all_have_wgs'] = triodf.iloc[:, [1,2,3]].apply(all_have_wgs, axis=1)


triodf = triodf[triodf['all_have_wgs']]

samples_trios = triodf.iloc[:, 1].to_list() + triodf.iloc[:, 2].to_list() + triodf.iloc[:, 3].to_list()
samples_trios = list(set(samples_trios))


with open(f'{intermediate_file_dir}/tdt_sample_list_trios.txt', 'w') as f:
	for line in samples_trios:
		f.write(line)
		f.write('\n')
