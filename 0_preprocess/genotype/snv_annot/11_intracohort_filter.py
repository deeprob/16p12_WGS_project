#!/bin/python3




import pandas as pd


df = pd.read_csv('exonic_variants/rare_deleterious_exonic_variants.csv')




df['variant_id'] = df['Chrom'] + '_' + df['Pos'].astype(str) + '_' + df['Ref'] + '_' + df['Alt']



# reorder columns
columns = ['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'variant_id', 'Qual', 'Func.refGene', 'Gene.refGene', 'GeneDetail.refGene', 'ExonicFunc.refGene', 'AAChange.refGene', 'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19', 'GeneDetail.wgEncodeGencodeBasicV19', 'ExonicFunc.wgEncodeGencodeBasicV19', 'AAChange.wgEncodeGencodeBasicV19', 'gnomad_exome_AF', 'gnomad_genome_AF', 'CADD_PHRED', 'CADD_RawScore', 'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB', 'Mut_type', 'variant_filter']
df = df[columns]





variant_counts = df['variant_id'].value_counts()
def get_cohort_count(variant_id):
	return variant_counts[variant_id]



df['cohort_count'] = df['variant_id'].apply(get_cohort_count)








# filter intracohort count
df = df[df['cohort_count'] <= 10]

df.to_csv('exonic_variants/rare_deleterious_exonic_variants_cohort_count.csv', index=False)











