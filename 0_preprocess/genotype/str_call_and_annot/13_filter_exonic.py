#!/bin/python3


import pandas as pd



filename = '16p12_cohort.expansions_2SD.annotated.tsv'
df = pd.read_csv(filename, sep='\t')

print(df.columns)

df = df[df['Func.wgEncodeGencodeBasicV19'] == 'exonic']

# remove refgene columns
columns =['chrom', 'pos', 'end', 'sample', 'zscore', 'longest_allele',
       'cohort_mode', 'motif', 'motif_period', 'variant_id', 'Func.wgEncodeGencodeBasicV19',
       'Gene.wgEncodeGencodeBasicV19', 'GeneDetail.wgEncodeGencodeBasicV19',
       'ExonicFunc.wgEncodeGencodeBasicV19',
       'AAChange.wgEncodeGencodeBasicV19']
df = df[columns]


df.to_csv('16p12_cohort.expansions_2SD.annotated.exonic.csv', index=False)












