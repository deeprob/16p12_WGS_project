#!/bin/python3

# originally prepared by anastasia
# edited by deepro
# edits 
# added in and out dir prefix


import pandas as pd

input_file_dir = "/data5/deepro/wgs_16p/rvgdt/data/input_files"
intermediate_file_dir = "/data5/deepro/wgs_16p/rvgdt/data/intermediate_files"

filename = f"{input_file_dir}/16p12_All_Participants_v9.xlsx"

df = pd.read_excel(filename)
df = df[df['No_consent_forms'].isna()]
df = df[df.WGS == 'X']

samples  = df.Sample.to_list()

with open(f'{intermediate_file_dir}/tdt_sample_list.txt', 'w') as f:
	for sample in samples:
		f.write(sample)
		f.write('\n')
