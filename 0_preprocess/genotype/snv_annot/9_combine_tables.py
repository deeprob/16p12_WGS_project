#!/bin/python3


# combine tables


import pandas as pd


# 16p12_All_Participants_v5.csv is from /Dropbox/16p12.2 project/Human patients project/WGS paper/3_Cohort information/16p12_All_Participants_v5.csv

samples_df = pd.read_csv('16p12_All_Participants_v5.csv')
samples_df = samples_df[samples_df['No_consent_forms'] != 'X']
samples_df = samples_df[samples_df['WGS'] == 'X']


df = pd.DataFrame()

for sample in samples_df.Sample.to_list():
	filename = f'exonic_variants/by_sample/{sample}.txt'
	app = pd.read_csv(filename, header=None, sep='\t')
	df = df.append(app)



columns = ['Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19', 'GeneDetail.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19', 'AAChange.wgEncodeGencodeBasicV19', 'gnomad_exome_AF', 'gnomad_genome_AF', 'CADD_PHRED', 'CADD_RawScore', 'Sample', 'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB']



df.columns = columns



df.to_csv('exonic_variants/rare_exonic_variants.csv', index=False)




