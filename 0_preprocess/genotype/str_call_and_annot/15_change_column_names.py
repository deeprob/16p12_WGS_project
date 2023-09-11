#!/bin/python3




# gencode.v19.parsed.genes.csv is from Dropbox/16p12.2 project/Human patients project/WGS paper/5_SNV pipelines-annotations/GENCODE annotations



import pandas as pd

df = pd.read_csv('16p12_cohort.expansions_2SD.annotated.exonic.gene_ids.csv')




['chrom', 'pos', 'end', 'sample', 'zscore', 'longest_allele',
       'cohort_mode', 'motif', 'motif_period', 'variant_id',
       'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19',
       'GeneDetail.wgEncodeGencodeBasicV19',
       'ExonicFunc.wgEncodeGencodeBasicV19',
       'AAChange.wgEncodeGencodeBasicV19', 'Gene_symbol', 'Gene_id',
       'Gene_biotype']


# rename columns
columns = ['Chrom', 'Pos', 'End', 'Sample', 'zscore', 'longest_allele',
       'cohort_mode', 'motif', 'motif_period', 'variant_id',
       'Func.wgEncodeGencodeBasicV19', 'Gene.wgEncodeGencodeBasicV19',
       'GeneDetail.wgEncodeGencodeBasicV19',
       'ExonicFunc.wgEncodeGencodeBasicV19',
       'AAChange.wgEncodeGencodeBasicV19', 'Gene_symbol', 'Gene_id',
       'Gene_biotype']
df.columns = columns



df.to_csv('16p12_cohort.expansions_2SD.annotated.exonic.gene_ids.column_names.csv', index=False)



