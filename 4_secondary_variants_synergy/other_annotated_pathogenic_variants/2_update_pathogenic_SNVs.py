#!/bin/python
import pandas as pd

# Update the pathogenic SNVs by incorporating the manual review
df = pd.read_csv('Intermediate_files/REVIEWED_Rare_Deleterious_Exonic_Variants_Pathogenic_Anno.csv')

df['Pathogenic'] = '.'

df.loc[(df.SFARI_Tier_1_S=='X') | (df.DBD_Tier_1_2=='X') | ((df.CLINVAR_Pathogenic=='X') & (df.ClinVar_Manual_Screen=='X')) | ((df.Dominant_X_OMIM=='X') & (df.OMIM_Manual_Screen=='X')), 'Pathogenic'] = 'X'

print(df.Pathogenic.value_counts())

# Reorder columns
df = df[['Sample', 'Chrom', 'Pos', 'Ref', 'Alt', 'Qual', 'Mut_type',
       'variant_filter', 'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19',
       'GeneDetail.wgEncodeGencodeBasicV19',
       'ExonicFunc.wgEncodeGencodeBasicV19',
       'AAChange.wgEncodeGencodeBasicV19', 'gnomad_exome_AF',
       'gnomad_genome_AF', 'CADD_PHRED', 'CADD_RawScore', 'ClinVar_CLNDN',
       'ClinVar_CLNDISDB', 'ClinVar_CLNREVSTAT', 'ClinVar_CLNSIG',
       'ClinVar_ALLELEID', 'GT', 'DP', 'AD', 'GQ', 'PL', 'RGQ', 'SB',
       'variant_id', 'cohort_count', 'Gene_symbol', 'Gene_id', 'Gene_biotype',
       'Gene_id_', 'LOEUF', 'MIM_number', 'inheritance', 'CLINVAR_Pathogenic', 'ClinVar_Manual_Screen',
       'OMIM_phenotype', 'Dominant_X_OMIM', 'OMIM_Manual_Screen', 'SFARI_gene_score',
       'SFARI_Tier_1_S', 'Geisinger_DBD_Tier', 'DBD_Tier_1_2', 'Pathogenic']]

# Save to file
df.to_csv('Rare_Deleterious_Exonic_Variants_Pathogenic_Anno.csv', index = False)