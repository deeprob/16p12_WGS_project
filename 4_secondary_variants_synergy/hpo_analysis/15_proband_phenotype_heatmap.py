#!/bin/python
import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

# Make heatmaps showing proband phenotypes and which can be explained by their genes

# Load in data
# Get probands
cohort_info = pd.read_csv('../../3_cohort information/16p12_All_Participants_v9.csv')
cohort_info = cohort_info[(cohort_info.WGS=='X') & (cohort_info.No_consent_forms!='X')]
pros = cohort_info[cohort_info.Relationship=='P']['Sample'].to_list()

# Phenotypes
pheno = pd.read_csv('../1_HPO_annotations/Annotated_files/phenotype_hpo.csv')[['Sample', 'HPO_ID']]
pheno = pheno[pheno.Sample.isin(pros)]
# Variants
snvs = pd.read_csv('../1_HPO_annotations/Annotated_files/snvs_hpo_anno.csv')[['Sample', 'Gene_symbol', 'Mut_type', 'Gene_id_', 'GeneHPO_inclusive']]
snvs = snvs[snvs.Sample.isin(pros)]
cnvs = pd.read_csv('../1_HPO_annotations/Annotated_files/cnvs_hpo_anno.csv')[['Sample', 'Gene_Symbol',  'Type', 'Gene_ID', 'GeneHPO_inclusive']]
cnvs = cnvs[cnvs.Sample.isin(pros)]

# Get all phenotypes present in probands
all_pheno_rows = pheno.HPO_ID.to_list()
all_pheno = []
for row in all_pheno_rows:
    all_pheno += [i.strip() for i in row.split(';')]
all_pheno = list(set(all_pheno))
all_pheno.sort()

# Get a list of all probands with variant data
var_pros = [i for i in pros if i in snvs.Sample.to_list() and i in cnvs.Sample.to_list() and i in pheno.Sample.to_list()]

# Make a starting dataframe with phenotypes as rows and probands as columns
pro_pheno = pd.DataFrame(index = all_pheno, columns = var_pros)

# Make heatmaps
pdf = PdfPages('Figures/1_proband_phenotype.pdf')
##############################
# Explained or not explained #
##############################
def explained(col):
    # Get sample ID
    pro = col.name
    
    # Get HPO IDs from genes
    terms = []
    for df in [snvs, cnvs]:
            term_lst = df[df.Sample==pro]['GeneHPO_inclusive'].to_list()
            for i in term_lst:
                terms += [j.strip() for j in i.split(';')]
    terms = list(set(terms))
    
    # Get proband phenotypes
    pheno_lst = pheno[pheno.Sample==pro]['HPO_ID'].to_list()[0]
    phenos = [i.strip() for i in pheno_lst.split(';')]
    
    output = []
    for row in all_pheno:
        # Check if phenotype is present in proband
        if row not in phenos:
            output.append(np.nan)
            continue
        # Check if phenotype is present in their variant HPO terms
        if row in terms:
            output.append(1)
            continue
        output.append(0)
    
    return pd.Series(output, index = col.index)
expl_df = pro_pheno.apply(lambda col: explained(col), axis = 0)

# Replace index with HPO translations
mc2hpo = pd.read_excel('../1_HPO_annotations/Analysis_files/MC_HPO_terms.xlsx')
mc2hpo.dropna(how = 'any', inplace = True)
mc2hpo = mc2hpo[['HPO ID', 'HPO Term', 'Phenotype Category']]
mc2hpo.drop_duplicates(inplace = True, keep = 'first')
mc2hpo.set_index(mc2hpo['HPO ID'], inplace = True)
hpo_dict = mc2hpo['HPO Term'].to_dict()
expl_df.index = expl_df.index.map(hpo_dict)
# Save to file
out_df = expl_df.fillna('.')
out_df.to_csv('Result_tables/1_proband_phenotype_binary.csv')

#Replace NAs with -1 for graphing
expl_df.fillna(-1, inplace = True)

# April 18, 2022
# We decided to group related phenotypes
# Do this using phenotype category information from the Master Checklist-HPO map
domain_defs = dict(zip(mc2hpo['HPO Term'].to_list(), mc2hpo['Phenotype Category'].to_list()))
domains = list(set(mc2hpo['Phenotype Category'].to_list()))
domains.sort()
expl_df['domain'] = expl_df.index.map(domain_defs)

# Function for clustering subgroups
def clustermap_subgroup(df, domains = domains):
    labels = []
    for domain in domains:
        sub_df = df[df.domain==domain].copy()
        if sub_df.shape[0]==0:
            continue
        if sub_df.shape[0]==1:
            labels += sub_df.index.to_list()
            continue
        sub_df.drop('domain', axis = 1, inplace  = True)
        grid = sns.clustermap(sub_df)
        labels += [i.get_text() for i in grid.ax_heatmap.get_ymajorticklabels()]
        plt.close()
    
    return labels
order = clustermap_subgroup(expl_df)
expl_df = expl_df.loc[order]
expl_df.drop('domain', axis = 1, inplace = True)

# Make plot
colors = ['white', 'k', '#9b2226']
cg = sns.heatmap(expl_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True, cbar = False)
cg.set_xticklabels(cg.get_xmajorticklabels(), fontsize = 2)
cg.set_yticklabels(cg.get_ymajorticklabels(), fontsize = 2)
plt.title('Binary: All Probands, All Phenotypes')
pdf.savefig()
plt.close()

# Shorten dataframe
# Limit plot to phenotypes that are present in at least 10 probands
expl_df.replace(-1, np.nan, inplace = True)
expl_df.dropna(thresh=10, axis=0, inplace=True)
expl_df.fillna(-1, inplace = True)
# Make plot
cg2 = sns.heatmap(expl_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True, cbar = False)
cg2.set_xticklabels(cg2.get_xmajorticklabels(), fontsize = 4)
cg2.set_yticklabels(cg2.get_ymajorticklabels(), fontsize = 4)
plt.title('Binary: All Probands, Phenotypes > 10 probands')
pdf.savefig()
plt.close()

# Limit to probands who have at least 3 phenotypes
expl_df.replace(-1, np.nan, inplace = True)
expl_df.dropna(thresh=3, axis=1, inplace=True)
expl_df.fillna(-1, inplace = True)
# Make plot
cg3 = sns.heatmap(expl_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True, cbar = False)
cg3.set_xticklabels(cg3.get_xmajorticklabels(), fontsize = 4)
cg3.set_yticklabels(cg3.get_ymajorticklabels(), fontsize = 4)
plt.title('Binary: Probands > 3 phenotypes, Phenotypes > 10 probands')
pdf.savefig()
plt.close()

# 4/18/2022
# After discussion, we decided to remove ID, borderline (HP:0006889) and Expressive language delay (HP:0002474) from the phenotype list
# 3/27/2022
# We also decided to remove Global Developmental Delay and Intellectual Disability from the list
expl_df = expl_df.loc[~expl_df.index.isin(['Intellectual disability, borderline', 'Expressive language delay', 'Global developmental delay', 'Delayed speech and language development', 'Intellectual disability'])]

# Make plot
cg4 = sns.heatmap(expl_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True, cbar = False)
cg4.set_xticklabels(cg4.get_xmajorticklabels(), fontsize = 4)
cg4.set_yticklabels(cg4.get_ymajorticklabels(), fontsize = 4)
plt.title('Binary: Probands > 3 phenotypes, Phenotypes > 10 probands\nID borderline and Expressive language delay removed')
pdf.savefig()
plt.close()

##########################################
# Number of genes that explain phenotype #
##########################################
def num_explain(col):
    # Get sample ID
    pro = col.name
    
    # Get HPO IDs from genes
    terms = []
    for df in [snvs, cnvs]:
            term_lst = df[df.Sample==pro]['GeneHPO_inclusive'].to_list()
            for i in term_lst:
                terms += [j.strip() for j in i.split(';')]
    terms = list(set(terms))
    
    # Get proband phenotypes
    pheno_lst = pheno[pheno.Sample==pro]['HPO_ID'].to_list()[0]
    phenos = [i.strip() for i in pheno_lst.split(';')]
    
    output = []
    for row in all_pheno:
        # Check if phenotype is present in proband
        if row not in phenos:
            output.append(np.nan)
            continue
        # Check if phenotype is present in their variant HPO terms
        if row not in terms:
            output.append(0)
            continue
        # If it is in their variant HPO terms, find the genes they have that explain that term
        dfs = [snvs, cnvs]
        gene_cols = ['Gene_id_', 'Gene_ID']
        genes = []
        for i in range(2):
            df = dfs[i]
            sub_df = df[(df.Sample==pro) & (df.GeneHPO_inclusive.str.contains(row))]
            if sub_df.shape[0]==0:
                continue
            genes+= sub_df[gene_cols[i]].to_list()
        genes = list(set(genes))
        if len(genes)==0:
            print('No genes that explain phenotype even though phenotype is in variant terms!')
            print(error)
        output.append(len(genes))
    
    return pd.Series(output, index = col.index)
count_df = pro_pheno.apply(lambda col: num_explain(col), axis = 0)
# Save to file
# Replace index with HPO translations
count_df.index = count_df.index.map(hpo_dict)
out_df = count_df.fillna('.')
out_df.to_csv('Result_tables/1_proband_phenotype_count.csv')

#Replace NAs with -1 for graphing
count_df.fillna(-1, inplace = True)

# Order rows
count_df['domain'] = count_df.index.map(domain_defs)
order = clustermap_subgroup(count_df)
count_df = count_df.loc[order]
count_df.drop('domain', axis = 1, inplace = True)

# Make plots
colors = ['white', 'grey'] + sns.color_palette("viridis_r", 20)
cg = sns.heatmap(count_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True)
cg.set_xticklabels(cg.get_xmajorticklabels(), fontsize = 2)
cg.set_yticklabels(cg.get_ymajorticklabels(), fontsize = 2)
plt.title('Count: All Probands, All Phenotypes')
pdf.savefig()
plt.close()

# Shorten dataframe
# Limit plot to phenotypes that are present in at least 10 probands
count_df.replace(-1, np.nan, inplace = True)
count_df.dropna(thresh=10, axis=0, inplace=True)
count_df.fillna(-1, inplace = True)
# Make plot
cg2 = sns.heatmap(count_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True)
cg2.set_xticklabels(cg2.get_xmajorticklabels(), fontsize = 4)
cg2.set_yticklabels(cg2.get_ymajorticklabels(), fontsize = 4)
plt.title('Count: All Probands, Phenotypes > 10 probands')
pdf.savefig()
plt.close()

# Limit to probands who have at least 3 phenotypes
count_df.replace(-1, np.nan, inplace = True)
count_df.dropna(thresh=3, axis=1, inplace=True)
count_df.fillna(-1, inplace = True)
# Make plot
cg3 = sns.heatmap(count_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True)
cg3.set_xticklabels(cg3.get_xmajorticklabels(), fontsize = 4)
cg3.set_yticklabels(cg3.get_ymajorticklabels(), fontsize = 4)
plt.title('Count: Probands > 3 phenotypes, Phenotypes > 10 probands')
pdf.savefig()
plt.close()

# 4/18/2022
# After discussion, we decided to remove ID, borderline (HP:0006889) and Expressive language delay (HP:0002474) from the phenotype list
count_df = count_df.loc[~count_df.index.isin(['Intellectual disability, borderline', 'Expressive language delay', 'Global developmental delay', 'Delayed speech and language development', 'Intellectual disability'])]
# Make plot
cg4 = sns.heatmap(count_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'lightgrey', xticklabels = True, yticklabels = True)
cg4.set_xticklabels(cg4.get_xmajorticklabels(), fontsize = 4)
cg4.set_yticklabels(cg4.get_ymajorticklabels(), fontsize = 4)
plt.title('Count: Probands > 3 phenotypes, Phenotypes > 10 probands\nID borderline and Expressive language delay removed')
pdf.savefig()
plt.close()

# Save the order of phenotypes here so plots used for analysis are consistent
phenotype_order=count_df.index.to_list()

#############################################
# Type of variant that explains a phenotype #
#############################################
def type_explain(col):
    # Get sample ID
    pro = col.name
    
    # Get HPO IDs from genes
    terms = []
    for df in [snvs, cnvs]:
            term_lst = df[df.Sample==pro]['GeneHPO_inclusive'].to_list()
            for i in term_lst:
                terms += [j.strip() for j in i.split(';')]
    terms = list(set(terms))
    
    # Get proband phenotypes
    pheno_lst = pheno[pheno.Sample==pro]['HPO_ID'].to_list()[0]
    phenos = [i.strip() for i in pheno_lst.split(';')]
    
    output = []
    for row in all_pheno:
        # Check if phenotype is present in proband
        if row not in phenos:
            output.append(np.nan)
            continue
        # Check if phenotype is present in their variant HPO terms
        if row not in terms:
            output.append(0)
            continue
        # If it is in their variant HPO terms, find the variant type(s) that explains the term
        dfs = [snvs, cnvs]
        types = []
        for i in range(2):
            df = dfs[i]
            sub_df = df[(df.Sample==pro) & (df.GeneHPO_inclusive.str.contains(row))]
            if sub_df.shape[0]==0:
                # Term is not explained by a variant in that dataframe
                continue
            if i==0:
                for type in ['missense', 'lof', 'splice']:
                    sub_df2 = sub_df[sub_df.Mut_type==type]
                    if sub_df2.shape[0]==0:
                        continue
                    types.append(type)
            elif i==1:
                for type in ['DEL', 'DUP']:
                    sub_df2 = sub_df[sub_df.Type==type]
                    if sub_df2.shape[0]==0:
                        continue
                    types.append(type)
        output.append(';'.join(types))
    
    return pd.Series(output, index = col.index)
type_df = pro_pheno.apply(lambda col: type_explain(col), axis = 0)
type_df.index = type_df.index.map(hpo_dict)
# Save this to a file
out_df = type_df.fillna('.')
out_df.to_csv('Result_tables/1_proband_phenotype_varianttype.csv')

# Make a new dataframe for plotting
# Codes:
# 1. LOF variant
# 2. Missense variant
# 3. Splice variant
# 4. Deletion
# 5. Duplication
# 6. Multiple variants
def type_explain_codes(col):
    # Get sample ID
    pro = col.name
    
    # Get HPO IDs from genes
    terms = []
    for df in [snvs, cnvs]:
            term_lst = df[df.Sample==pro]['GeneHPO_inclusive'].to_list()
            for i in term_lst:
                terms += [j.strip() for j in i.split(';')]
    terms = list(set(terms))
    
    # Get proband phenotypes
    pheno_lst = pheno[pheno.Sample==pro]['HPO_ID'].to_list()[0]
    phenos = [i.strip() for i in pheno_lst.split(';')]
    
    output = []
    for row in all_pheno:
        # Check if phenotype is present in proband
        if row not in phenos:
            output.append(np.nan)
            continue
        # Check if phenotype is present in their variant HPO terms
        if row not in terms:
            output.append(0)
            continue
        # If it is in their variant HPO terms, find the variant type(s) that explains the term
        dfs = [snvs, cnvs]
        types = []
        for i in range(2):
            df = dfs[i]
            sub_df = df[(df.Sample==pro) & (df.GeneHPO_inclusive.str.contains(row))]
            if sub_df.shape[0]==0:
                # Term is not explained by a variant in that dataframe
                continue
            if i==0:
                for type in ['missense', 'lof', 'splice']:
                    sub_df2 = sub_df[sub_df.Mut_type==type]
                    if sub_df2.shape[0]==0:
                        continue
                    types.append(type)
            elif i==1:
                for type in ['DEL', 'DUP']:
                    sub_df2 = sub_df[sub_df.Type==type]
                    if sub_df2.shape[0]==0:
                        continue
                    types.append(type)
        if len(types) > 1:
            val = 6
        elif types==['lof']:
            val = 1
        elif types==['missense']:
            val = 2
        elif types==['splice']:
            val = 3
        elif types==['DEL']:
            val = 4
        elif types==['DUP']:
            val = 5
        output.append(val)
    
    return pd.Series(output, index = col.index)
type_df = pro_pheno.apply(lambda col: type_explain_codes(col), axis = 0)

# Change index
type_df.index = type_df.index.map(hpo_dict)
# Change NA to -2
type_df.fillna(-2, inplace = True)

# Order rows
type_df['domain'] = type_df.index.map(domain_defs)
order = clustermap_subgroup(type_df)
type_df = type_df.loc[order]
type_df.drop('domain', axis = 1, inplace = True)

# Make plots
colors = ['white', 'grey', 'grey'] + sns.color_palette("tab10")[0:5] + ['k']
cg = sns.heatmap(type_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True, cbar = False)
cg.set_xticklabels(cg.get_xmajorticklabels(), fontsize = 2)
cg.set_yticklabels(cg.get_ymajorticklabels(), fontsize = 2)
plt.title('Types: All Probands, All Phenotypes')
pdf.savefig()
plt.close()

# Shorten dataframe
# Limit plot to phenotypes that are present in at least 10 probands
type_df.replace(-2, np.nan, inplace = True)
type_df.dropna(thresh=10, axis=0, inplace=True)
type_df.fillna(-2, inplace = True)
# Make plot
cg2 = sns.heatmap(type_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True, cbar = False)
cg2.set_xticklabels(cg2.get_xmajorticklabels(), fontsize = 4)
cg2.set_yticklabels(cg2.get_ymajorticklabels(), fontsize = 4)
plt.title('Types: All Probands, Phenotypes > 10 probands')
pdf.savefig()
plt.close()

# Limit to probands who have at least 3 phenotypes
type_df.replace(-2, np.nan, inplace = True)
type_df.dropna(thresh=3, axis=1, inplace=True)
type_df.fillna(-2, inplace = True)
# Make plot
cg3 = sns.heatmap(type_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True, cbar = False)
cg3.set_xticklabels(cg3.get_xmajorticklabels(), fontsize = 4)
cg3.set_yticklabels(cg3.get_ymajorticklabels(), fontsize = 4)
plt.title('Types: Probands > 3 phenotypes, Phenotypes > 10 probands')
pdf.savefig()
plt.close()

# 4/18/2022
# After discussion, we decided to remove ID, borderline (HP:0006889) and Expressive language delay (HP:0002474) from the phenotype list
type_df = type_df.loc[phenotype_order]
# Make plot
cg4 = sns.heatmap(type_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'lightgrey', xticklabels = True, yticklabels = True, cbar = False)
cg4.set_xticklabels(cg4.get_xmajorticklabels(), fontsize = 4)
cg4.set_yticklabels(cg4.get_ymajorticklabels(), fontsize = 4)
plt.title('Types: Probands > 3 phenotypes, Phenotypes > 10 probands\nID borderline and Expressive language delay removed')
pdf.savefig()
plt.close()

pdf.close()

