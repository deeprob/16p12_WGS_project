import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.colors import LinearSegmentedColormap
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
sns.set_style({'font.family':'sans-serif', 'font.sans-serif':'Arial'})

# Plot heatmap of burden comparisons

# Get data and clean up
df=pd.read_csv('statistics/t-tests_16p12_max_samples.tsv', sep='\t')
var_d={'autism_PRS':'Autism PRS', 'educational_attainment_PRS':'Education PRS', 'SCZ_PRS':'SCZ PRS', 'intelligence_PRS':'Intelligence PRS', 'UTR5':"5' UTR", 'promoter':'Promoter', 'enhancer':'Enhancer',
		'Rare_Deleterious_SNVs':'All coding SNVs', 'genes_del':'Genes del.', 'genes_dup':'Genes dup.', 'STRs_exonic':'STRs',
		'Rare_Deleterious_SNVs_LOEUF':'All coding SNVs (LF)', 'dels_loeuf':'Genes del. (LF)', 'dups_loeuf':'Genes dup. (LF)', 'STRs_exonic_LOEUF035':'STRs (LF)',
		'Splice_CADD25_LOEUF035':'Splice (LF)', 'Splice_CADD25':'Splice', 'LOF_LOEUF_035':'LOF (LF)', 'LOF':'LOF', 'Missense_CADD25_LOEUF035':'Missense (LF)', 'Missense_CADD25':'Missense'}
df['variant']=df['Variant Group'].map(var_d)
var_lst=['Missense', 'Missense (LF)', 'LOF', 'LOF (LF)', 'Splice', 'Splice (LF)', 'Genes del.', 'Genes dup.',
			'Genes del. (LF)', 'Genes dup. (LF)', 'STRs', 'STRs (LF)', 'All coding SNVs', 'All coding SNVs (LF)', 'Enhancer', 'Promoter', "5' UTR",
			'Intelligence PRS', 'SCZ PRS', 'Education PRS', 'Autism PRS']
df.variant=pd.Categorical(df.variant, var_lst)

df.sort_values(by='variant', inplace=True)

df['star']=''
df.loc[df['P-value']<=0.05, 'star']='*'
df.loc[df.FDR<=0.05, 'star']='**'

# Plot heatmap
plot_df=pd.pivot(df, columns='group2', index='variant', values='Cohens D')
p_df=pd.pivot(df, columns='group2', index='variant', values='star')

# Plot
colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plot_df, cmap=cmap, vmin=-.28, vmax=.28, square=True, fmt='', linecolor='#CCCCCC', linewidths=0.75, annot=p_df)
plt.savefig('figures/3_heatmap.pdf')
plt.close()
