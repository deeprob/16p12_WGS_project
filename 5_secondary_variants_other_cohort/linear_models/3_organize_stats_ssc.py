import pandas as pd

# Organize the linear regression statistics for the supplemental table
groups = ['dbd_tier1_snvs', 'large_rare_deletions', 'large_rare_duplications']
phenotypes = ['Full_scale_IQ', 'ABCL_CBCL_external', 'ABCL_CBCL_internal', 'SRS_raw', 'RBS_R', 'DCDQ', 'BMI_zscore']

df = pd.DataFrame()
for group in groups:
	for phenotype in phenotypes:
		for i in range(2,4):
			model = 'model{}'.format(i)
			filename = 'C:/Users/corny/Dropbox/16p12.2 project/Human patients project/WGS paper/31_SSC analysis/Genotype-phenotype integration/4_regression_models/statistics/{}_{}_{}.tsv'.format(group, phenotype, model)
			app = pd.read_csv(filename, sep='\t')
			app.columns = ['Variant class/PRS', 'Direction', '2.5% C.I.', '97.5% C.I.', 'P-value','Test', 'Direction metric', 'Sample size', 'R2']
			app['Group'] = group
			app['Other info'] = model
			app['Phenotype'] = phenotype
			df = df.append(app)
		for v in ['SCZ_PRS', 'Rare_Deleterious_SNVs', 'Genes DEL', 'Genes DUP', 'STRs exonic', 'Rare_Deleterious_SNVs_LOEUF', 'DELs LOEUF0.35', 'DUPs LOEUF0.35', 'STRs exonic LOEUF0.35']:
			filename = 'C:/Users/corny/Dropbox/16p12.2 project/Human patients project/WGS paper/31_SSC analysis/Genotype-phenotype integration/4_regression_models/statistics/{}_{}_{}.tsv'.format(group, phenotype, v)
			app = pd.read_csv(filename, sep='\t')
			app.columns = ['Variant class/PRS', 'Direction', '2.5% C.I.', '97.5% C.I.', 'P-value','Test', 'Direction metric', 'Sample size', 'R2']
			app['Group'] = group
			app['Other info'] = v
			app['Phenotype'] = phenotype
			df = df.append(app)

# Clean up names
pheno_d={'Full_scale_IQ':'Full scale IQ', 'ABCL_CBCL_external':'Externalizing behavior (ABCL/CBCL)', 'ABCL_CBCL_internal':'Internalizing behavior (ABCL/CBCL)', 'SRS_raw':'Social responsiveness (SRS)',
        'RBS_R':'Repetetive behavior (RBS_R)', 'DCDQ':'Coordination disorder (DCDQ)', 'BMI_zscore':'BMI Z-score', 'Full-scale IQ':'Full scale IQ', 'ABCL_CBCL external':'Externalizing behavior (ABCL/CBCL)',
        'ABCL_CBCL internal':'Internalizing behavior (ABCL/CBCL)', 'SRS raw':'Social responsiveness (SRS)', 'BSI':'Autism behavior (BSI)', 'BMI Z-score':'BMI Z-score', 'Head circumference Z-score':'Head circumference Z-score'}
df.Phenotype=df.Phenotype.map(pheno_d)

var_d={'(Intercept)':'Intercept', 'Sex':'Sex', 'SCZ_PRS':'SCZ PRS', 'Rare_Deleterious_SNVs':'All coding SNVs', '`Genes DEL`':'Genes del.', '`Genes DUP`':'Genes dup.', '`STRs exonic`':'STRs', 'Rare_Deleterious_SNVs_LOEUF':'All coding SNVs (LF)',
        '`DELs LOEUF<0.35`':'Genes del. (LF)', '`DUPs LOEUF<0.35`':'Genes dup. (LF)', '`STRs exonic LOEUF<0.35`':'STRs (LF)'}
df['variable']=df['Variant class/PRS'].map(var_d)
df.loc[(df.variable=='SCZ PRS') & (df['Other info']=='model3'), 'variable'] = 'SCZ PRS (model3)'

mod_map={'model2':'Phenotype ~ Sex + All coding SNVs + STRs + Genes del. + Genes dup. + SCZ PRS', 'model3':'Phenotype ~ Sex + All coding SNVs (LF) + STRs (LF) + Genes del. (LF) + Genes dup. (LF) + SCZ PRS',
        'SCZ_PRS':'Phenotype ~ Sex + SCZ PRS', 'Rare_Deleterious_SNVs':'Phenotype ~ Sex + All coding SNVs', 'Genes DEL':'Phenotype ~ Sex + Genes del.', 'Genes DUP':'Phenotype ~ Sex + Genes dup.', 'STRs exonic':'Phenotype ~ Sex + STRs',
        'Rare_Deleterious_SNVs_LOEUF':'Phenotype ~ Sex + All coding SNVs (LF)', 'DELs LOEUF0.35':'Phenotype ~ Sex + Genes del. (LF)', 'DUPs LOEUF0.35':'Phenotype ~ Sex + Genes dup. (LF)', 'STRs exonic LOEUF0.35':'Phenotype ~ Sex + STRs (LF)'}
df['model']=df['Other info'].map(mod_map)

df['Group']=df.Group.map({'dbd_tier1_snvs':'DBD Tier 1 SNVs', 'large_rare_deletions':'Large rare deletions', 'large_rare_duplications':'Large rare duplications'})

df=df[df.variable!='Intercept']
df=df[['Group', 'Phenotype', 'variable', 'Sample size', 'Direction', 'P-value', 'R2', 'model']]

# Save to file
df.to_csv('s6h.csv', index=False)