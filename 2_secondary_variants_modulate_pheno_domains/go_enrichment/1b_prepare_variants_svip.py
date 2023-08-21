#!/bin/python3



import pandas as pd
import numpy as np
import subprocess


dropbox = '~/Dropbox/'
subprocess.run('mkdir -p variants/SVIP', shell=True)


#=========================
# Load mastertable
#=========================


filename = '~/Dropbox/16p12.2 project/Human patients project/WGS paper/32_SVIP analysis/Genotype-phenotype integration/1_variant_preperation/SVIP.csv'
master_df = pd.read_csv(filename)
master_df = master_df.set_index('Sample', drop=False)
master_df = master_df[master_df.Relationship == 'Proband']

#=========================
# Save samples with WGS
#=========================


cols = ['Sample', 'Family', 'Relationship','Sex','Family type']
master_df[cols].to_csv('cohorts/SVIP_samples.csv', index=False)


cols = ['Family', 'Sample', 'Gene', 'Relationship']


#=========================
# WGS SNVs
#=========================


filename = dropbox + '16p12.2 project/Human patients project/WGS paper/32_SVIP analysis/SNVs_indels/SVIP_rare_deleterious_snvs_indels.csv'
df = pd.read_csv(filename)
df = df[df.Sample.isin(master_df.index)]
df = df.set_index('Sample', drop=False)


df['Gene'] = df['Gene_symbol']

df = df[['Sample', 'Gene']]
df.to_csv('variants/SVIP/rare_deleterious_snvs.csv', index=False)



#=========================
# Microarray CNVs
#=========================

# get list of samples with microarray
filename = dropbox + '16p12.2 project/Human patients project/WGS paper/32_SVIP analysis/CNVs/SVIP_CNVs_by_gene_loeuf.csv'
df = pd.read_csv(filename)
df = df.set_index('Sample', drop=False)





# Dups
subdf = df[df['CNV_Type'] == 'Dup']
subdf = subdf[['Sample', 'Gene']]
subdf.to_csv('variants/SVIP/dups.csv', index=False)


# Dels
subdf = df[df['CNV_Type'] == 'Del']
subdf = subdf[['Sample', 'Gene']]
subdf.to_csv('variants/SVIP/dels.csv', index=False)




