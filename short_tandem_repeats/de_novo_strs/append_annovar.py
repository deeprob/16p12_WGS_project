#!/bin/python3


import sys
import pandas as pd

# filename = 'output/filtered/all.filtered4.tab'
filename = sys.argv[1]

# annotations_filename = 'output/filtered/annovar.hg19_multianno.txt'
annotations_filename = sys.argv[2]

outfilename = sys.argv[3]


anno = pd.read_csv(annotations_filename, sep='\t')
anno['var_id'] = anno['Chr'].apply(lambda s: s[3:]) + '_' + anno['Start'].astype(str)
anno = anno.set_index('var_id', drop=False)

# list of columns that will annotate
columns_to_annotate = list(anno.columns)
columns_to_annotate = columns_to_annotate[5:-1]


df = pd.read_csv(filename, sep='\t')

str_ids = df.str_id.to_list()

# append all columns
for col in columns_to_annotate:
	df[col] = anno.loc[str_ids][col].to_list()

df.to_csv(outfilename, sep='\t', index=False)

