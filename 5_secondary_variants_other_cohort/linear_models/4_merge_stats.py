import pandas as pd

import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

matplotlib.rcParams['pdf.fonttype'] = 42

# Merge statistics from SSC and SVIP models
ssc=pd.read_csv('../31_SSC analysis/Genotype-phenotype integration/5_integrate_statistics/SSC_variant_v_phenotype.csv')
ssc=ssc[ssc.Test=='Linear regression']
ssc=ssc[ssc['Other info'].isin(['model2', 'model3'])]

ssc=ssc[~ssc.Group.str.contains('nejm')]
ssc['Group']=ssc.Group.map({'dbd_tier1_snvs':'DBD Tier 1 SNVs', 'large_rare_deletions':'Large rare deletions', 'large_rare_duplications':'Large rare duplications'})

svip=pd.read_csv('../32_SVIP analysis/Genotype-phenotype integration/5_integrate_statistics/SVIP_variant_v_phenotype.csv')
svip=svip[svip.Test=='Linear regression']
svip=svip[svip['Other info'].isin(['model2', 'model3'])]

svip['Group']=svip.Group.map({'16p-deletion':'16p11.2 deletion', '16p-duplication':'16p11.2 duplication'})

df=pd.concat([ssc, svip])
df=df[['Group', 'Other info', 'Variant class/PRS', 'Phenotype', 'P-value', 'Direction', 'Sample size', 'R2']]

# Clean up names
pheno_d={'Full_scale_IQ':'Full scale IQ', 'ABCL_CBCL_external':'Externalizing behavior', 'ABCL_CBCL_internal':'Internalizing behavior', 'SRS_raw':'Social responsiveness', 'RBS_R':'Repetetive behavior', 'DCDQ':'Coordination disorder',
        'BMI_zscore':'BMI Z-score', 'Full-scale IQ':'Full scale IQ', 'ABCL_CBCL external':'Externalizing behavior', 'ABCL_CBCL internal':'Internalizing behavior', 'SRS raw':'Social responsiveness', 'BSI':'Autism behaviors', 'BMI Z-score':'BMI Z-score',
        'Head circumference Z-score':'Head circumference Z-score'}
df.Phenotype=df.Phenotype.map(pheno_d)

print(df['Variant class/PRS'].unique())
var_d={'(Intercept)':'Intercept', 'Sex':'Sex', 'SCZ_PRS':'SCZ PRS', 'Rare_Deleterious_SNVs':'All coding SNVs', '`Genes DEL`':'Genes del.', '`Genes DUP`':'Genes dup.', '`STRs exonic`':'STRs', 'Rare_Deleterious_SNVs_LOEUF':'All coding SNVs (LF)',
        '`DELs LOEUF<0.35`':'Genes del. (LF)', '`DUPs LOEUF<0.35`':'Genes dup. (LF)', '`STRs exonic LOEUF<0.35`':'STRs (LF)'}
df['variable']=df['Variant class/PRS'].map(var_d)
df.loc[(df.variable=='SCZ PRS') & (df['Other info']=='model3'), 'variable'] = 'SCZ PRS (model3)'

df['Model']=df['Other info'].map({'model2':'All variants', 'model3':'LF variants'})

df['star']=''
df.loc[df['P-value']<=0.05, 'star']='*'

# Restrict to values we want to plot
df=df[~df.variable.isin(['Sex', 'Intercept'])]

# Organize data
group_lst=['16p11.2 deletion', '16p11.2 duplication', 'DBD Tier 1 SNVs', 'Large rare deletions', 'Large rare duplications']
pheno_lst=['Full scale IQ', 'Social responsiveness', 'Coordination disorder', 'Externalizing behavior', 'Internalizing behavior', 'Repetetive behavior', 'Autism behaviors', 'BMI Z-score', 'Head circumference Z-score']
var_lst=['All coding SNVs', 'Genes del.', 'Genes dup.', 'STRs', 'SCZ PRS', 'All coding SNVs (LF)', 'Genes del. (LF)', 'Genes dup. (LF)', 'STRs (LF)', 'SCZ PRS (model3)']

df.Group=pd.Categorical(df.Group, group_lst)
df.Phenotype=pd.Categorical(df.Phenotype, pheno_lst)
df.variable=pd.Categorical(df.variable, var_lst)

df.sort_values(by=['Phenotype', 'variable', 'Group'], inplace=True)

# Keep only phenotype-var pairs with a significant result
df['p_v']=df.Phenotype.astype(str)+'_'+df.variable.astype(str)

sig=list(df[df['P-value']<=0.05]['p_v'].unique())
df=df[df.p_v.isin(sig)]
print(df)

# Convert from long to wide
plot_df=pd.pivot(df, columns='Group', index=['Phenotype', 'Model', 'variable'], values='Direction')
print(plot_df)

p_df=pd.pivot(df, columns='Group', index=['Phenotype', 'Model', 'variable'], values='star')
p_df.fillna('', inplace=True)
print(p_df)

# Plot
colors=["#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", '#FFFFFF', "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"]
cmap=LinearSegmentedColormap.from_list('BlRd', colors, N=20)
cmap.set_bad('#CCCCCC')
sns.heatmap(data=plot_df, cmap=cmap, vmin=-0.5, vmax=0.5, square=True, fmt='', linecolor='#CCCCCC', linewidths=0.75, annot=p_df)
plt.savefig('summary_heatmap.pdf')
plt.close()
