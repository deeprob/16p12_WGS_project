import pandas as pd

anno = pd.read_csv('output/16p12_cohort.expansions_2SD.hg19_multianno.txt', sep='\t')
anno['variant_id'] = anno['Chr'] + '_' + anno['Start'].astype(str) + '_' + anno['Alt']
anno = anno.set_index('variant_id')
anno = anno.drop_duplicates().copy()

df = pd.read_csv('output/16p12_cohort.expansions_2SD.tsv', sep='\t')
df['variant_id'] = df['chrom'] + '_' + df['pos'].astype(str) + '_' + df['alt_allele']


# filter out intergenic variants
# anno = anno[anno['Func.refGene'] != 'intergenic'].copy()
# df = df[df.variant_id.isin(list(anno.index))].copy()

# change to dict
anno_cols = anno.columns[5:]
anno = anno.to_dict()

cols = ['chrom', 'pos', 'end', 'sample', 'zscore',
       'longest_allele', 'cohort_mode', 'motif_period', 'variant_id']
# drop ref allele and alt allele
df = df[cols].copy()

for col in anno_cols:
	df[col] = '.'

total = df.shape[0]
for i, row in df.iterrows():
	if i % 1000 == 0:
		print('{}/{}'.format(i, total))
	variant_id = row['variant_id']
	for col in anno_cols:
		annotation = anno[col][variant_id]
		df.at[i, col] = annotation


df.to_csv('output/16p12_cohort.expansions_2SD.annotated.tsv', sep='\t', index=False)


