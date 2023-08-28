import pandas as pd






groups = ['16p-deletion', '16p-duplication']
phenotypes = ['Full-scale IQ','ABCL_CBCL external','ABCL_CBCL internal','SRS raw','BSI','BMI Z-score','Head circumference Z-score']

#----------------------------
# variant v. phenotype
#----------------------------

# Group	Test	Direction metric	Other info	Variant class/PRS	Phenotype	P-value	FDR	Direction	Sample size

df = pd.DataFrame()


# Pearson's R
filename = '~/Dropbox/16p12.2 project/Human patients project/WGS paper/32_SVIP analysis/Genotype-phenotype integration/3_correlation/statistics/variant_v_phenotype.csv'
app = pd.read_csv(filename)
app['Test'] = 'Correlation'
app['Direction metric'] = 'Pearson\'s R'
app['Other info'] = ''
df = df.append(app)


for group in groups:
	for phenotype in phenotypes:
		for i in range(1,4):
			print(i)
			model = 'model{}'.format(i)
			filename = '~/Dropbox/16p12.2 project/Human patients project/WGS paper/32_SVIP analysis/Genotype-phenotype integration/4_regression_models/statistics/{}_{}_{}.tsv'.format(group, phenotype, model)
			app = pd.read_csv(filename, sep='\t')
			app.columns = ['Variant class/PRS', 'Direction', '2.5% C.I.', '97.5% C.I.', 'P-value','Test', 'Direction metric', 'Sample size', 'R2']
			app['Group'] = group
			app['Other info'] = model
			app['Phenotype'] = phenotype
			df = df.append(app)








cols = ['Group','Test','Direction metric','Other info','Variant class/PRS','Phenotype','P-value','Direction','Sample size', 'R2']
df = df[cols]

filename = 'SVIP_variant_v_phenotype.csv'
df.to_csv(filename, index=False)

#----------------------------
# variant v. variant
#----------------------------


# Test	Direction metric	Other info	Variant class/PRS1	Variant class/PRS2	Group1	Group2	P-value	FDR	Direction	Sample size



# Pearson's R
filename = '~/Dropbox/16p12.2 project/Human patients project/WGS paper/32_SVIP analysis/Genotype-phenotype integration/3_correlation/statistics/variant_v_variant.csv'
df = pd.read_csv(filename)
df = df.drop('use_in_fdr_calculation', axis=1)
df.columns = ['Variant class/PRS1', 'Variant class/PRS2', 'Group', 'P-value', 'Sample size', 'Direction', 'FDR']


df['Direction metric'] = 'Pearson\'s R'
df['Other info'] = ''
df['Group1'] = 'Probands'
df['Group2'] = 'Probands'
df['Test'] = 'Correlation'




cols = ['Group', 'Test','Direction metric','Other info','Variant class/PRS1','Variant class/PRS2','Group1','Group2','P-value','FDR','Direction','Sample size']
df = df[cols]

filename = 'SVIP_variant_v_variant.csv'
df.to_csv(filename, index=False)



















