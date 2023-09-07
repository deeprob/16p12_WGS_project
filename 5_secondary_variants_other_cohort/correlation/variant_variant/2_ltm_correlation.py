import pandas as pd

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rcParams['pdf.fonttype'] = 42

# Merge correlation statistics
ssc=pd.read_csv('../31_SSC analysis/Genotype-phenotype integration/5_integrate_statistics/SSC_variant_v_variant.csv')
ssc=ssc[~ssc.Group.str.contains('nejm')]
ssc['Group']=ssc.Group.map({'dbd_tier1_snvs':'DBD Tier 1 SNVs', 'large_rare_deletions':'Large rare deletions', 'large_rare_duplications':'Large rare duplications'})

svip=pd.read_csv('../32_SVIP analysis/Genotype-phenotype integration/5_integrate_statistics/SVIP_variant_v_variant.csv')
svip['Group']=svip.Group.map({'16p-deletion':'16p11.2 deletion', '16p-duplication':'16p11.2 duplication'})

df=pd.concat([ssc, svip])
df=df[['Group', 'Variant class/PRS1', 'Variant class/PRS2', 'P-value', 'FDR', 'Direction']]
keep_vars=['Rare_Deleterious_SNVs', 'autism_PRS', 'educational_attainment_PRS', 'intelligence_PRS', 'SCZ_PRS']
df=df[(df['Variant class/PRS1']==keep_vars[0]) & (df['Variant class/PRS2'].isin(keep_vars[1:]))]
df.sort_values(by='Variant class/PRS2', inplace=True)

# Annotate significance
df['star']=''
df.loc[df['P-value']<=0.05, 'star']='*'
df.loc[df.FDR<=0.05, 'star']='**'
print(df)

# Reformat to wide for plotting
group_lst=['16p11.2 deletion', '16p11.2 duplication', 'DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications']
prs_lst=['autism_PRS', 'educational_attainment_PRS', 'intelligence_PRS', 'SCZ_PRS']
df.Group=pd.Categorical(df.Group, group_lst)
df['Variant class/PRS2']=pd.Categorical(df['Variant class/PRS2'], prs_lst)
df.sort_values(by=['Group', 'Variant class/PRS2'], inplace=True)

plotdf=df.pivot(columns='Variant class/PRS2', index='Group', values='Direction')
sigdf=df.pivot(columns='Variant class/PRS2', index='Group', values='star')

# Plot
colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plotdf, cmap=cmap, vmin=-0.55, vmax=0.55, square=True, fmt='', linecolor='#CCCCCC', linewidths=0.75, annot=sigdf)
plt.tight_layout()
plt.savefig('ltm_heatmap.pdf')
plt.close()