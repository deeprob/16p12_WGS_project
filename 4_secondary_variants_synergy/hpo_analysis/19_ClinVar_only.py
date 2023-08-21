import pandas as pd
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns

matplotlib.rcParams['pdf.fonttype'] = 42

# Make a heatmap showing phenotypes associated with genes with ClinVar variants in each proband

# Load in data
# Get probands
cohort_info = pd.read_csv('../../3_cohort information/16p12_All_Participants_v9.csv')
cohort_info = cohort_info[(cohort_info.WGS=='X') & (cohort_info.No_consent_forms!='X')]
pros = cohort_info[cohort_info.Relationship=='P']['Sample'].to_list()

# Phenotypes
pheno = pd.read_csv('../1_HPO_annotations/Annotated_files/phenotype_hpo.csv')[['Sample', 'HPO_ID']]
pheno = pheno[pheno.Sample.isin(pros)]
# Variants
snvs = pd.read_csv('Analysis_files/4_ClinVar_variants_manualHPO.csv')[['Sample', 'Gene_symbol', 'Mut_type', 'ClinVar_Manual_Screen', 'Gene_id_', 'VarHPO', 'GeneHPO_inclusive']]
snvs = snvs[(snvs.Sample.isin(pros)) & (snvs.ClinVar_Manual_Screen=='X')]

pheno=pheno[pheno.Sample.isin(snvs.Sample.to_list())]

# Get all phenotypes present in probands
all_pheno_rows = pheno.HPO_ID.to_list()
all_pheno = []
for row in all_pheno_rows:
    all_pheno += [i.strip() for i in row.split(';')]
all_pheno = list(set(all_pheno))
all_pheno.sort()

# Make a starting dataframe with phenotypes as rows and probands as columns
pro_pheno = pd.DataFrame(index = all_pheno, columns = list(snvs.Sample.unique()))

# Make heatmaps
pdf = PdfPages('Figures/5_proband_clinvar.pdf')
#######################################################################################################################################
# Associated with varaiant and present, associated with variant and not present, not associated with variant and present, not present #
#######################################################################################################################################
def explained(col):
    # Get sample ID
    pro = col.name
    
    # Get HPO IDs from gene(s)
    var_terms = []
    gene_terms = []
    var_term_lst = snvs[snvs.Sample==pro]['VarHPO'].to_list()
    gene_term_lst=snvs[snvs.Sample==pro]['GeneHPO_inclusive'].to_list()
    for i in range(len(var_term_lst)):
        var_terms += [j.strip() for j in var_term_lst[i].split(';')]
        gene_terms += [j.strip() for j in gene_term_lst[i].split(';')]
    var_terms = list(set(var_terms))
    gene_terms=list(set(gene_terms))
    
    # Get proband phenotypes
    pheno_lst = pheno[pheno.Sample==pro]['HPO_ID'].to_list()[0]
    phenos = [i.strip() for i in pheno_lst.split(';')]
    
    output = []
    for row in all_pheno:
        # Not present
        if row not in phenos:
            output.append(0)
        # Present and not associated with variant terms or gene terms
        elif row not in var_terms and row not in gene_terms:
            output.append(1)
        # Present and associated with variant terms
        elif row in var_terms:
            output.append(2)
        # Present and associated with gene terms
        elif row in gene_terms:
            output.append(3)
    return pd.Series(output, index = col.index)
    
expl_df = pro_pheno.apply(lambda col: explained(col), axis = 0)
print(expl_df)

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
out_df.to_csv('Result_tables/5_clinvar.csv')

# Cluster related phenotypes
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
colors = ['white', 'lightpink', 'black', 'grey']
cg = sns.heatmap(expl_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True, cbar = True)
cg.set_xticklabels(cg.get_xmajorticklabels(), fontsize = 2)
cg.set_yticklabels(cg.get_ymajorticklabels(), fontsize = 2)
plt.title('Probands and ClinVar phenotypes')
pdf.savefig()
plt.close()

# Shorten dataframe
# Limit plot to phenotypes that are present in at least 3 probands
expl_df.replace(0, np.nan, inplace = True)
expl_df.dropna(thresh=3, axis=0, inplace=True)
expl_df.fillna(0, inplace = True)

# Make plot
cg2 = sns.heatmap(expl_df, cmap = colors, square = True, linewidths = 0.3, linecolor = 'k', xticklabels = True, yticklabels = True, cbar = True)
cg2.set_xticklabels(cg2.get_xmajorticklabels(), fontsize = 4)
cg2.set_yticklabels(cg2.get_ymajorticklabels(), fontsize = 4)
plt.title('Probands and ClinVar phenotypes, Phenotypes > 5 probands')
pdf.savefig()
plt.close()

pdf.close()