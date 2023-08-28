#!/bin/python3
import pandas as pd
import numpy as np

# libraries related to plotting
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_style({'font.family':'sans-serif', 'font.sans-serif':'Arial'})
from matplotlib.backends.backend_pdf import PdfPages

df = pd.read_csv('../5_integrate_statistics/SSC_variant_v_phenotype.csv')
df=df[df.Test=='Linear regression']
df=df[~df.Group.str.contains('nejm')]
df=df[df['Other info']!='model1']
df=df[df['Variant class/PRS']!='(Intercept)']
df=df[df['Variant class/PRS']!='Sex']
print(df)

colors=['#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')

phenos={'Full_scale_IQ':'Full scale IQ',
		'ABCL_CBCL_external':'Externalizing behavior',
		'ABCL_CBCL_internal':'Internalizing behavior',
		'SRS_raw':'Social responsiveness',
		'RBS_R':'Repetitive behavior',
		'DCDQ':'Coordination disorder',
		'BMI_zscore':'BMI Z-score'}
		
vars={'model2':'All variants',
			'model3':'All variants (LF)',
			'SCZ_PRS':'SCZ PRS',
			'Rare_Deleterious_SNVs':'All coding SNVs',
			'Genes DEL':'Genes del.',
			'Genes DUP':'Genes dup.',
			'STRs exonic':'STRs',
			'Rare_Deleterious_SNVs_LOEUF':'All coding SNVs (LF)',
			'DELs LOEUF0.35':'Genes del. (LF)',
			'DUPs LOEUF0.35':'Genes dup. (LF)',
			'STRs exonic LOEUF0.35':'STRs (LF)'}

df['variant']=df['Other info'].map(vars)
df['pheno']=df.Phenotype.map(phenos)
df['pheno']=pd.Categorical(df.pheno, ['Full scale IQ', 'Externalizing behavior', 'Internalizing behavior', 'Social responsiveness', 'Repetitive behavior', 'Coordination disorder', 'BMI Z-score'])

mods = [['All variants', 'All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs', 'SCZ PRS'], ['All variants (LF)', 'All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS']]

pdf = PdfPages('variance_heatmap_phenotypes.pdf')
for c in ['dbd_tier1_snvs', 'large_rare_deletions', 'large_rare_duplications']:
	for i in range(2):
		title=c+' '+['Model', 'Model (LF)'][i]
		mod=mods[i]

		subdf = df[(df.Group==c) & (df['variant'].isin(mod))][['variant', 'pheno', 'R2']]
		subdf.drop_duplicates(inplace=True)
		subdf['variant']=pd.Categorical(subdf.variant, mod)
		subdf.sort_values(by=['variant', 'pheno'], inplace=True)
		
		dfplot=subdf.pivot(index='variant', columns='pheno', values='R2')

		sns.set(font_scale=0.7)
		fig = plt.figure(figsize=(2,2))

		g = sns.heatmap(data=dfplot, fmt='', square=True, cmap=cmap, linecolor='#CCCCCC', linewidths=0.75, vmin=0, vmax=0.25)

		g.set_title(title)

		pdf.savefig(fig, bbox_inches='tight')

groups={'dbd_tier1_snvs':'DBD Tier 1 SNVs',
		'large_rare_deletions':'Large rare deletions',
		'large_rare_duplications':'Large rare duplications'}

df['LF']=''
df.loc[df.variant.str.contains('(LF)'), 'LF']=' (LF)'

df['Group']=df.Group.map(groups)
df['label']=df.Group+df.LF

df['label']=pd.Categorical(df.label, ['DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications', 'DBD Tier 1 SNVs (LF)', 'Large rare deletions (LF)', 'Large rare duplications (LF)'])

df=df[(df.variant.str.contains('All variants'))][['pheno', 'label', 'R2']]
df.drop_duplicates(inplace=True)
df.sort_values(by=['label', 'pheno'], inplace=True)

dfplot=df.pivot(index='label', columns='pheno', values='R2')

sns.set(font_scale=0.7)
fig = plt.figure(figsize=(2,2))

g = sns.heatmap(data=dfplot, fmt='', square=True, cmap=cmap, linecolor='#CCCCCC', linewidths=0.75, vmin=0, vmax=0.25)

pdf.savefig(fig, bbox_inches='tight')

pdf.close()





