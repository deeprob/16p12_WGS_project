#!/bin/python3




import pandas as pd




dropbox='~/Dropbox'


samples_df = pd.read_excel(dropbox + '/16p12.2 project/Human patients project/WGS paper/11_Variant Integration/16p12_cohort_summary_v17.xlsx')
samples_df = samples_df[samples_df.Relationship == 'P']
samples_df = samples_df[~samples_df.Missense_CADD25.isna()]
samples = samples_df.Sample.to_list()



variants_df = pd.read_csv('tables/rare_exonic_proteincoding_variants.txt', sep='\t')


genes = list(set(variants_df['Gene.refGene'].to_list()))
genes = [s.replace('\\x3b', '_') for s in genes]
genes = sorted(genes)
print(len(genes))

outdf = pd.DataFrame(index=samples)
outdf['Sample_Name'] = outdf.index.to_series()


i = 0
for gene in genes:
	if i % 1000 == 0:
		print(f'{i} {gene}')
	i = i + 1
	
	samples_with_variant_in_gene = variants_df[variants_df['Gene.refGene'] == gene].Sample.to_list()
	
	outdf[f'Input_{gene}'] = outdf['Sample_Name'].isin(samples_with_variant_in_gene)
	outdf[f'Input_{gene}'] = outdf[f'Input_{gene}'].astype(int)



outdf.to_csv('tables/variants.csv', index=False)







