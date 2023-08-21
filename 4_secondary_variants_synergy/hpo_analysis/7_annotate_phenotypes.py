#!/bin/python
import pandas as pd

dropbox = 'C:/Users/corny/Dropbox/'

# Translate Master Checklist files into HPO terms for each person for easy filtering
mc2hpo = pd.read_excel('Analysis_files/MC_HPO_terms.xlsx')
# Remove any rows with NAs as they are just labels to organize data for readability
mc2hpo.dropna(inplace = True)
print(mc2hpo)

# Child Master Checklist
child_mc = pd.read_excel(dropbox+"16p12.2 project/Human patients project/Phenotypic analysis/de Vries & Psychiatric/WGS/MC_Children/MC_Compiled_working.xlsx", skiprows = 1, index_col = 0, header = None, dtype = str)
# Only keep rows that are in the mc2hpo file and sample ID
keep_rows = [i for i in child_mc.index.to_list() if i in mc2hpo['MC Phenotype'].unique()]
child_mc = child_mc.loc[keep_rows+['Child Code']]
# Reformat
child_mc.columns = child_mc.loc['Child Code'].to_list()
child_mc = child_mc.transpose()

# For each child, annotate the HPO terms for their phenotypes
def add_hpo(row):
    phenos = row[row=='1'].index.to_list()
    hpo = mc2hpo[mc2hpo['MC Phenotype'].isin(phenos)]['HPO ID'].to_list()
    output = list(set(hpo))
    output.sort()
    return ';'.join(output)
child_mc['HPO_ID'] = child_mc.apply(add_hpo, axis = 1)
print(child_mc)

# Add code and phenotypes to new dataframe
df = child_mc[['Child Code', 'HPO_ID']].copy()
df.columns = ['Sample', 'HPO_ID']

# Adult Master Checklist
adult_mc = pd.read_excel(dropbox+"16p12.2 project/Human patients project/Phenotypic analysis/de Vries & Psychiatric/WGS/MC_Adult/AdultMC_FINAL.xlsx", index_col = 0, header = None, sheet_name = 'FINAL final', dtype = str)
adult_mc.fillna('.', inplace = True)
# Only keep rows that are in the mc2hpo file and sample ID
keep_rows = [i for i in adult_mc.index.to_list() if i in mc2hpo['MC Phenotype'].unique()]
adult_mc = adult_mc.loc[keep_rows+['Patient code']]
# Change Sample IDs
# Some samples have 'NO SAMPLE' as an ID - give these all unique numbers
counter = 1
new_ids = []
for i in adult_mc.loc['Patient code'].to_list():
    if i=='NO SAMPLE':
        out = 'NOSAMPLE'+str(counter)
        counter+=1
    else:
        out = i
    new_ids.append(out)
adult_mc.loc['Sample'] = new_ids
print(adult_mc)
# Reformat
adult_mc.columns = adult_mc.loc['Sample'].to_list()
adult_mc = adult_mc.transpose()
# For each person, annotate the HPO terms for their phenotypes
def add_hpo_adult(row):
    phenos = row[(row!='0') & (row!='.')].index.to_list()
    hpo = mc2hpo[mc2hpo['MC Phenotype'].isin(phenos)]['HPO ID'].to_list()
    output = list(set(hpo))
    output.sort()
    return ';'.join(output)
adult_mc['HPO_ID'] = adult_mc.apply(add_hpo_adult, axis = 1)
print(adult_mc)

# Add to dataframe
df = pd.concat([df, adult_mc[['Sample', 'HPO_ID']]], axis = 0)
print(df)

# Fill all empty cells with 0
# If a sample is in this file, it means we have some phenotypic data for them
# Therefore, if they have no HPO terms it is because they have no phenotypes
df.replace('', '0', inplace = True)

# Add in the number of phenotypes for each person
def num_pheno(hpo):
    if hpo=='0':
        return 0
    return len(hpo.split(';'))
df['number_of_phenotypes'] = df.HPO_ID.apply(num_pheno)

# Save to file
df.to_csv('Annotated_files/phenotype_hpo.csv', index = False)