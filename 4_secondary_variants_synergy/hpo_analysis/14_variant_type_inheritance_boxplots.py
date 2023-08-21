#!/bin/python
import pandas as pd
import numpy as np

import scipy.stats as stats

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

# Open circle marker
pnts = np.linspace(0, np.pi * 2, 24)
circ = np.c_[np.sin(pnts) / 2, -np.cos(pnts) / 2]
vert = np.r_[circ, circ[::-1] * .7]
open_circle = mpl.path.Path(vert)

# Load in proportion data
df = pd.read_csv('Result_tables/3_variant_type_inheritance_proportions.csv')

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
cohort_info = cohort_info[(cohort_info.WGS=='X') & (cohort_info.No_consent_forms!='X')]
samples = cohort_info['Sample'].to_list()
prop_df = prop_df[prop_df.Sample.isin(samples)]

# Filter for only complete trios
# Mother-father trios
def mf_trios(fam):
    members = cohort_info[cohort_info.Family==fam]['Relationship'].to_list()
    if 'P' in members:
        if 'MC' in members or 'MNC' in members:
            if 'FC' in members or 'FNC' in members:
                return True
    return False
fams = pd.Series(cohort_info.Family.unique())
trios_mf = fams[fams.apply(mf_trios)]
mf_trio_pros = cohort_info[(cohort_info.Relationship=='P') & (cohort_info.Family.isin(trios_mf))]['Sample'].to_list()

# Carrier-NonCarrier trios
def cnc_trios(fam):
    members = cohort_info[cohort_info.Family==fam]['Relationship'].to_list()
    if 'P' in members:
        if 'MC' in members and 'FNC' in members:
            return True
        if 'MNC' in members and 'FC' in members:
            return True
    return False
trios_cnc = fams[fams.apply(cnc_trios)]
cnc_trio_pros = cohort_info[(cohort_info.Relationship=='P') & (cohort_info.Family.isin(trios_cnc))]['Sample'].to_list()

# Make boxplots showing the proportion of variants explained in probands from different inheritances
pro_df = prop_df[prop_df.Relationship=='P'].copy()

pdf = PdfPages('Figures/4_variant_type_inheritance_boxplots.pdf')

# Inherited vs. De novo
mf_df = pro_df[pro_df.Sample.isin(mf_trio_pros)]
inh_dn_cols = [i for i in prop_cols if 'inh' in i or 'DN' in i]
use_cols = [i for i in inh_dn_cols if 'STR' not in i]
mf_df = mf_df[['Sample'] + use_cols]
mf_reform = pd.DataFrame(columns = ['Proportion', 'type', 'Inh', 'Sample'])
# Reformat the table for easy plotting
for inh in ['inh', 'DN']:
    for col in [i for i in use_cols if inh in i]:
        name = col.split('_'+inh)[0]
        df = pd.DataFrame(mf_df[[col]+['Sample']])
        df.columns = ['Proportion', 'Sample']
        df['type']=name
        df['Inh'] = inh
        mf_reform = pd.concat([mf_reform, df], axis = 0)
# Sort by sample
mf_reform.sort_values(by = ['Sample', 'type'], inplace = True)

# Make boxplots
type_order = ['pathogenic_SNVs', 'LOF', 'LOF_LOEUF', 'missense', 'missense_LOEUF', 'splice', 'splice_LOEUF', 'DEL', 'DEL_LOEUF', 'DUP', 'DUP_LOEUF', 'All']
f1 = plt.figure()
sns.boxplot(data = mf_reform, x = 'type', y = 'Proportion', hue = 'Inh', fliersize = 0, palette="Set2", order = type_order)
sns.stripplot(data = mf_reform, x = 'type', y = 'Proportion', hue = 'Inh', dodge = True, color = 'k', marker = open_circle, size = 4, linewidth = 0.2, order = type_order)
locs, labs = plt.xticks()
plt.xticks(locs, ['Pathogenic SNVs', 'LOF', 'LOF_LOEUF', 'Missense', 'Missense_LOEUF', 'Splice', 'Splice_LOEUF', 'DEL', 'DEL_LOEUF', 'DUP', 'DUP_LOEUF', 'All'], rotation = 90)
plt.title('Inherited vs. De Novo in Probands')
plt.tight_layout()
pdf.savefig(f1)

# Calculate stats
stat_lst = []
for type in mf_reform.type.unique():
    stat, p = stats.ttest_rel(mf_reform[(mf_reform.type==type) & (mf_reform.Inh=='inh')]['Proportion'].to_numpy(), mf_reform[(mf_reform.type==type) & (mf_reform.Inh=='DN')]['Proportion'].to_numpy(), alternative = 'greater')
    stat_lst.append([type, 'InhvDN', 'one-tailed greater paired t test', stat, p])

# Carrier vs Non-Carrier
cnc_df = pro_df[pro_df.Sample.isin(cnc_trio_pros)]
cnc_cols = [i for i in prop_cols if '_C_' in i or '_NC_' in i]
cnc_df = cnc_df[cnc_cols+['Sample']]
cnc_reform = pd.DataFrame(columns = ['Proportion', 'type', 'Inh', 'Sample'])
# Reformat the table for easy plotting
for inh in ['C', 'NC']:
    for col in [i for i in cnc_cols if '_'+inh+'_' in i]:
        name = col.split('_'+inh+'_')[0]
        df = pd.DataFrame(cnc_df[[col]+['Sample']])
        df.columns = ['Proportion', 'Sample']
        df['type']=name
        df['Inh'] = inh
        cnc_reform = pd.concat([cnc_reform, df], axis = 0)
# Sort by sample
cnc_reform.sort_values(by = ['Sample', 'type'], inplace = True)
cnc_reform.reset_index(inplace=True, drop=True)

# Make boxplots
# type_order2 = ['pathogenic_SNVs', 'LOF', 'LOF_LOEUF', 'missense', 'missense_LOEUF', 'splice', 'splice_LOEUF', 'DEL', 'DEL_LOEUF', 'DUP', 'DUP_LOEUF', 'STR', 'STR_LOEUF', 'All']
# April 18, 2022 - we decided to remove STR from all HPO analysis
type_order2 = ['pathogenic_SNVs', 'LOF', 'LOF_LOEUF', 'missense', 'missense_LOEUF', 'splice', 'splice_LOEUF', 'DEL', 'DEL_LOEUF', 'DUP', 'DUP_LOEUF', 'All']
f2 = plt.figure()
sns.boxplot(data = cnc_reform, x = 'type', y = 'Proportion', hue = 'Inh', fliersize = 0, palette="Set2", order = type_order2)
sns.stripplot(data = cnc_reform, x = 'type', y = 'Proportion', hue = 'Inh', dodge = True, color = 'k', marker = open_circle, size = 4, linewidth = 0.2, order = type_order2)
locs, labs = plt.xticks()
plt.xticks(locs, ['Pathogenic SNVs', 'LOF', 'LOF_LOEUF', 'Missense', 'Missense_LOEUF', 'Splice', 'Splice_LOEUF', 'DEL', 'DEL_LOEUF', 'DUP', 'DUP_LOEUF', 'All'], rotation = 90)
plt.title('Carrier Inherited vs. NonCarrier Inherited in Probands')
plt.tight_layout()
pdf.savefig(f2)

pdf.close()

# Calculate stats
for type in cnc_reform.type.unique():
    stat, p = stats.ttest_rel(cnc_reform[(cnc_reform.type==type) & (cnc_reform.Inh=='C')]['Proportion'].to_numpy(), cnc_reform[(cnc_reform.type==type) & (cnc_reform.Inh=='NC')]['Proportion'].to_numpy())
    stat_lst.append([type, 'CarriervNonCarrier', 'two-tailed paired t test', stat, p])

# Save stats to file
stat_df = pd.DataFrame(stat_lst, columns = ['Variant_Type', 'Comparison', 'Test', 'Test_stat', 'p'])
stat_df.to_csv('Result_tables/4_variant_type_inheritance_stats.csv', index = False)