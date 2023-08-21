#!/bin/python
import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import seaborn.objects as so

matplotlib.rcParams['pdf.fonttype'] = 42

import scipy.stats as stats

# Load in proportion data
df = pd.read_csv('Result_tables/1_variant_type_proportions.csv')
print(df)

# Make boxplots for proportion of phenotypes explained by each variant type
# Get proportion columns
prop_cols = [i for i in df.columns.to_list() if 'proportion' in i]
prop_df = df[['Sample'] + prop_cols].copy()
prop_df.dropna(how = 'any', inplace = True)

# Annotate samples with relationship status
cohort_info = pd.read_csv('../../3_cohort information/16p12_All_Participants_v9.csv')
cohort_info.set_index(cohort_info.Sample, inplace = True)
rel_dict = cohort_info.Relationship.to_dict()
prop_df['Relationship'] = prop_df.Sample.map(rel_dict)

# Get samples with valid consent
samples = cohort_info[(cohort_info.WGS=='X') & (cohort_info.No_consent_forms!='X')]['Sample'].to_list()
prop_df = prop_df[prop_df.Sample.isin(samples)]

# Make boxplots
pdf = PdfPages('Figures/2_variant_type_boxplots.pdf')

# Proband only
prop_pro = prop_df[prop_df.Relationship=='P'].copy()
# Reformat
props=[]
for idx, row in prop_pro.iterrows():
    props+=row[prop_cols].to_list()
pro_df=pd.DataFrame({'relationship':['P']*(len(prop_cols)*prop_pro.shape[0]), 'var_type':prop_cols*prop_pro.shape[0], 'proportion':props, 'sample':list(np.repeat(prop_pro.Sample.to_list(), len(prop_cols)))})
print(pro_df)
cols = sns.color_palette('flare_r', 9)[2:] + sns.color_palette("viridis", 4) + sns.light_palette("seagreen", 3)[1:] + ["white"]
f1 = plt.figure()
sns.boxplot(data=pro_df, x='var_type', y='proportion', fliersize = 0, palette = cols)
sns.stripplot(data=pro_df, x='var_type', y='proportion', marker='$\circ$', color='k', jitter=0.2, size=6)
locs, labs = plt.xticks()
plt.xticks(locs, ['Pathogenic SNVs', 'LOF', 'LOF_LOEUF', 'Missense', 'Missense_LOEUF', 'Splice', 'Splice_LOEUF', 'DEL', 'DEL_LOEUF', 'DUP', 'DUP_LOEUF', 'STR', 'STR_LOEUF', 'All'], rotation = 90)
plt.tight_layout()
pdf.savefig(f1)

# Limit variant types
subpro=pro_df[~pro_df.var_type.str.contains('LOEUF')]
subpro=subpro[~subpro.var_type.str.contains('splice')]
f1 = plt.figure()
sns.boxplot(data=subpro, x='var_type', y='proportion', fliersize = 0, palette = [cols[0], cols[1], cols[3], cols[5], cols[7], cols[9], cols[13]])
sns.stripplot(data=subpro, x='var_type', y='proportion', marker='$\circ$', color='k', jitter=0.2, size=6)
locs, labs = plt.xticks()
plt.xticks(locs, ['Pathogenic SNVs', 'LOF', 'Missense', 'Splice', 'DEL', 'DUP', 'All'], rotation = 90)
plt.tight_layout()
pdf.savefig(f1)

# Save proband proportions
pro_df.to_csv('Result_tables/Proband_proportions.csv', index=False)
# Add mean information using Excel

# Proband and sibling boxplots
# Save stats
stat_lsts = []
names = ['Pathogenic SNVs', 'LOF', 'LOF_LOEUF', 'Missense', 'Missense_LOEUF', 'Splice', 'Splice_LOEUF', 'DEL', 'DEL_LOEUF', 'DUP', 'DUP_LOEUF', 'STR', 'STR_LOEUF', 'All']
kid_df = prop_df[prop_df.Relationship.isin(['P', 'SC', 'SNC'])].copy()
for i, type in enumerate(prop_cols):
    f = plt.figure()
    sns.boxplot(data = kid_df, y = type, x = 'Relationship', fliersize = 0)
    sns.stripplot(data = kid_df, y = type, x = 'Relationship', color = 'k', size = 4, edgecolor = 'white', linewidth = 0.2)
    plt.ylabel('Proportion of phenotypes explained')
    plt.title(names[i])
    pdf.savefig(f)
    
    # Stats
    # Mann Whitney 
    done = []
    for kid1 in ['P', 'SC', 'SNC']:
        for kid2 in ['P', 'SC', 'SNC']:
            if kid1==kid2:
                continue
            kids = [kid1, kid2]
            kids.sort()
            if kids in done:
                continue
            
            t, p = stats.mannwhitneyu(kid_df[kid_df.Relationship==kid1][type].to_numpy(), kid_df[kid_df.Relationship==kid2][type].to_numpy())
            stat_lsts.append([kid1, kid2, names[i], 'two-tailed Mann Whitney U', t, p])
            done.append(kids)
    # Kruskal-Wallis
    t, p = stats.kruskal(kid_df[kid_df.Relationship=='P'][type].to_numpy(), kid_df[kid_df.Relationship=='SC'][type].to_numpy(), kid_df[kid_df.Relationship=='SNC'][type].to_numpy(), nan_policy = 'omit')
    stat_lsts.append(['.', '.', names[i], 'Kruskall-Wallis', t, p])
pdf.close()

stat_df = pd.DataFrame(stat_lsts, columns = ['Group1', 'Group2', 'Variant_type', 'Test', 'Test_statistic', 'p_value'])
stat_df.to_csv('Result_tables/2_variant_type_stats.csv', index = False)

