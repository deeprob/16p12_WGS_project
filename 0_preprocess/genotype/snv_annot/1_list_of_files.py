#!/bin/python3




import pandas as pd
import os


filenames1 = os.listdir('/data3/16p12_WGS/final_calls/annovar/')
filenames1 = [s for s in filenames1 if s.endswith('_input.vcf')]
filenames1 = [f'/data3/16p12_WGS/final_calls/annovar/{s}' for s in filenames1]

filenames2 = os.listdir('/data3/16p12_WGS/final_calls_batch4/annovar/')
filenames2 = [s for s in filenames2 if s.endswith('_input.vcf')]
filenames2 = [f'/data3/16p12_WGS/final_calls_batch4/annovar/{s}' for s in filenames2]



filenames = filenames1 + filenames2

df = pd.DataFrame(index=filenames)
df['Filename'] = df.index.to_series()
df['Sample'] = df['Filename'].apply(lambda s: s.split('/')[-1].split('_')[0])

# remove duplicate 047
df = df[df['Filename'] != '/data3/16p12_WGS/final_calls/annovar/SG047_batch3_input.vcf']



df.to_csv('samples.csv', index=False)
