#!/bin/bash

# create one ped file per trio

import pandas as pd


ped = pd.read_csv('all_batches.ped', sep='\t', header=None)

families = ped[0].to_list()
families = list(set(families))

samples = ped[1].to_list()
samples = list(set(samples))


# for fam in families:
	# fam_ped = ped[ped[0] == fam].copy()

for i, row in ped.iterrows():
	samp = row[1]
	fam  = row[0]

	# skip if not trio
	if row[2] == '0' or row[3] == '0':
		continue

	trio_samples = [samp, row[2], row[3]]


	trio_ped = ped[ped[1].isin(trio_samples)]

	# write out
	trio_ped.to_csv('peds/{}.trio.ped'.format(samp), sep='\t', index=False, header=None)



ls peds/ | cut -f1 -d. > children_of_trios.txt

print('done')


