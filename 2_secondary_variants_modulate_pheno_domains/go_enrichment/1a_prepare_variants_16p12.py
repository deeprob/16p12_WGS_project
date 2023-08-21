#!/bin/python3



import pandas as pd
import numpy as np
import subprocess



dropbox_location = '~/Dropbox/'

subprocess.run('mkdir -p cohorts', shell=True)
subprocess.run('mkdir -p variants/16p12', shell=True)

#=========================
# Load mastertable
#=========================


filename = '~/Dropbox/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx'
master_df = pd.read_excel(filename)
master_df = master_df.set_index('Sample', drop=False)
print(master_df.columns)
master_df = master_df[~master_df['Rare_Deleterious_SNVs'].isna()] # only want samples with WGS 
samples_with_wgs = list(master_df.Sample.to_list())

#=========================
# Save samples with WGS
#=========================



cols = ['Sample', 'Family', 'Relationship','Sex','Carrier','Child_ID_DD', 'Child_behav', 'Child_psych', 'Child_nervous_system', 'Child_congenital', 'Child_craniofacial', 'De_vrie']
master_df[cols].to_csv('cohorts/16p12_samples.csv', index=False)



#=========================
# DUPs and DELs
#=========================

df = pd.read_csv(dropbox_location+'16p12.2 project/Human patients project/WGS paper/7_Structural variants/sv_calls_combined.txt', sep='\t')
df = df[df.Sample.isin(samples_with_wgs)]


dels_df = df[df.Type == 'DEL']
dups_df = df[df.Type == 'DUP']

dels_df = dels_df[['Sample', 'Gene_Symbol']]
dups_df = dups_df[['Sample', 'Gene_Symbol']]

dels_df.columns = ['Sample', 'Gene']
dups_df.columns = ['Sample', 'Gene']


dups_df.to_csv('variants/16p12/dups.csv', index=False)
dels_df.to_csv('variants/16p12/dels.csv', index=False)




#=========================
# STRs
#=========================


filename = f'{dropbox_location}/16p12.2 project/Human patients project/WGS paper/8_STR variants/Exonic_STR_expansions.csv'


df = pd.read_csv(filename)
df = df[df.Sample.isin(samples_with_wgs)]

df = df[['Sample', 'Gene.wgEncodeGencodeBasicV19']]
df.columns = ['Sample', 'Gene']

# split genes
new_genes_list = []
for i, row in df.iterrows():
	sample = row['Sample']
	genes = row['Gene']
	for gene in genes.split(';'):
		new_genes_list.append([sample, gene])

new_df = pd.DataFrame(new_genes_list, columns=['Sample', 'Gene'])


new_df.to_csv('variants/16p12/strs_std_2.csv', index=False)




#=========================
# Load variant table
#=========================


# New SNV calls
df = pd.read_csv(dropbox_location+'16p12.2 project/Human patients project/WGS paper/9_Rare variant calls/Rare coding SNV calls/Rare_Deleterious_Exonic_Variants.csv')
df = df[df.Sample.isin(samples_with_wgs)]

df = df[['Sample', 'Gene_symbol']]

# split genes
new_genes_list = []
for i, row in df.iterrows():
	sample = row['Sample']
	genes = row['Gene_symbol']
	for gene in genes.split(';'):
		new_genes_list.append([sample, gene])

new_df = pd.DataFrame(new_genes_list, columns=['Sample', 'Gene'])
new_df.to_csv('variants/16p12/rare_deleterious_snvs.csv', index=False)







