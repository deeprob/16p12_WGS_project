#!/bin/python3



import pandas as pd
import numpy as np

dropbox_location = '~/Dropbox/'



part_df = pd.read_csv(dropbox_location+'16p12.2 project/Human patients project/WGS paper/3_Cohort information/Estonian_Participants.csv')
part_df = part_df.sort_values('Sample')
part_df = part_df[part_df['No_consent_forms'].isna()]
part_df = part_df[(part_df.WGS == 'X')]
# part_df = part_df.set_index('Sample', drop=False)
samples_with_wgs = list(set(part_df.Sample.to_list()))


#=========================
# SNVs
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
new_df.to_csv('variants/rare_deleterious_snvs.csv', index=False)

#=========================
# Dups and Dels
#=========================

df = pd.read_csv(dropbox_location+'16p12.2 project/Human patients project/WGS paper/7_Structural variants/sv_calls_combined.txt', sep='\t')
df = df[df.Sample.isin(samples_with_wgs)]


dels_df = df[df.Type == 'DEL']
dups_df = df[df.Type == 'DUP']

dels_df = dels_df[['Sample', 'Gene_Symbol']]
dups_df = dups_df[['Sample', 'Gene_Symbol']]

dels_df.columns = ['Sample', 'Gene']
dups_df.columns = ['Sample', 'Gene']


dups_df.to_csv('variants/dups.csv', index=False)
dels_df.to_csv('variants/dels.csv', index=False)



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


new_df.to_csv('variants/strs_std_2.csv', index=False)


#=========================
# noncoding
#=========================


# enhancer calls
df = pd.read_csv(dropbox_location+'/16p12.2 project/Human patients project/WGS paper/9_Rare variant calls/Non-coding SNVs/Rare_enhancer_variants.txt', sep='\t')
df = df[df['genehancer_doubleelite_target_genesymbol']!='.']
df = df[df.Sample.isin(samples_with_wgs)]

# split genes
new_genes_list = []
for i, row in df.iterrows():
	sample = row['Sample']
	genes = row['genehancer_doubleelite_target_genesymbol']
	for gene in genes.split(','):
		new_genes_list.append([sample, gene])

new_df = pd.DataFrame(new_genes_list, columns=['Sample', 'Gene'])
new_df.to_csv('variants/enhancers.csv', index=False)



# utr5 calls
df = pd.read_csv(dropbox_location+'/16p12.2 project/Human patients project/WGS paper/9_Rare variant calls/Non-coding SNVs/Rare_UTR5_variants.txt', sep='\t')
df = df[['Sample', 'Gene.wgEncodeGencodeBasicV19']]
df.columns = ['Sample', 'Gene']
df = df[df.Sample.isin(samples_with_wgs)]

# split genes
new_genes_list = []
for i, row in df.iterrows():
	sample = row['Sample']
	genes = row['Gene']
	for gene in genes.split(';'):
		new_genes_list.append([sample, gene])

new_df = pd.DataFrame(new_genes_list, columns=['Sample', 'Gene'])
new_df.to_csv('variants/utr5.csv', index=False)


# promoter calls
df = pd.read_csv(dropbox_location+'/16p12.2 project/Human patients project/WGS paper/9_Rare variant calls/Non-coding SNVs/Rare_promoter_variants.txt', sep='\t')
df = df[['Sample', 'Gene.wgEncodeGencodeBasicV19']]
df.columns = ['Sample', 'Gene']
df = df[df.Sample.isin(samples_with_wgs)]

new_genes_list = []
for i, row in df.iterrows():
	sample = row['Sample']
	genes = row['Gene']
	for gene in genes.split(';'):
		new_genes_list.append([sample, gene])

new_df = pd.DataFrame(new_genes_list, columns=['Sample', 'Gene'])
new_df.to_csv('variants/promoters.csv', index=False)





