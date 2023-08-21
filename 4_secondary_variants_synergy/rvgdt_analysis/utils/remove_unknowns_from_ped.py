#!/bin/python3

# originally prepared by anastasia
# edited by deepro
# edits 
# added in and out dir prefix
# changed file "all_batches.ped" to "Mar2022.fam"
# .fam files have headers unlike .ped files, hence added header=0

import pandas as pd

input_file_dir = "/data5/deepro/wgs_16p/rvgdt/data/input_files"
intermediate_file_dir = "/data5/deepro/wgs_16p/rvgdt/data/intermediate_files"

df = pd.read_csv(f'{input_file_dir}/Mar_2022.fam', sep='\t', header=0) # changed to header=0 from header=None
# there was a major error in the line below, value counts was removed to account for that error!
df.iloc[:, 5] = df.iloc[:, 5].apply(lambda s: '0' if s == 'unknown' else s) #.value_counts() # had to change this to index loc to make it more generalizable
df.to_csv(f'{intermediate_file_dir}/all_batches.no_unknowns.ped', sep='\t', header=False, index=False)
