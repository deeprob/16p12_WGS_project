#!/bin/python
import pandas as pd
import numpy as np

# Get the proportion of phenotypes explained by each class of variants inherited from each parent
# Load in phenotype data
pheno = pd.read_csv('../1_HPO_annotations/Annotated_files/phenotype_hpo.csv', dtype = str)
print(pheno)

# Load in variant data
snvs = pd.read_csv('../1_HPO_annotations/Annotated_files/snvs_hpo_anno.csv')

# Add in annotations for the number and IDs of phenotypes explained by the genes
# Important note: as of 1/24/2022, we are only using a strict overlap method to determine explained phenotypes
df = pheno
def snv_explained(row, name, inheritance, patho = False, var_type = False, loeuf = False):
    name = name+'_'+inheritance
    # Phenotype HPO terms
    pheno_hpo = row['HPO_ID'].split(';')
    # UPDATE: 3/28/2023 - Santhosh wanted to remove ID and DD from the HPO calculations
    pheno_hpo=list(set(pheno_hpo)-set(['HP:0001263', 'HP:0001249']))
    
    # If no phenotype terms, return NA
    if pheno_hpo==['0'] or len(pheno_hpo)==0:
        row[name+'_explained_terms'] = np.nan
        row[name+'_number_explained_terms'] = np.nan
        row[name+'_proportion_explained_terms'] = np.nan
        return row
        
    sample = row['Sample']
    # If sample is not in SNV dataframe, return NA (no WGS)
    if sample not in snvs.Sample.to_list():
        row[name+'_explained_terms'] = np.nan
        row[name+'_number_explained_terms'] = np.nan
        row[name+'_proportion_explained_terms'] = np.nan
        return row
    
    # Define variant set
    if patho:
        df = snvs[snvs.Pathogenic=='X']
    else:
        df = snvs
    
    if var_type:
        df = df[df.Mut_type==var_type]
    
    if loeuf:
        df = df[df.LOEUF!='.']
        df = df[df.LOEUF.astype(float) <= 0.35]
    
    if inheritance=='C':
        df = df[df.inheritance.isin(['MC', 'FC'])]
    elif inheritance=='NC':
        df = df[df.inheritance.isin(['MNC', 'FNC'])]
    elif inheritance=='inh':
        df = df[df.inheritance.isin(['MC', 'FC', 'MNC', 'FNC'])]
    elif inheritance=='DN':
        df = df[df.inheritance=='de novo']
    else:
        print('NEED TO DEFINE INHERITANCE')
        print(error)
      
    df = df[df.Sample==sample]
    
    # If there are no variants in a category, return 0
    if df.shape[0]==0:
        row[name+'_explained_terms'] = '.'
        row[name+'_number_explained_terms'] = 0
        row[name+'_proportion_explained_terms'] = 0
        return row
    
    # Get gene HPO terms
    gene_hpo_str = df.GeneHPO_inclusive.to_list()
    gene_hpo = []
    for i in gene_hpo_str:
        gene_hpo += i.split(';')
    gene_hpo = list(set(gene_hpo))
    
    # Get overlap
    overlap = [i for i in gene_hpo if i in pheno_hpo]
    overlap.sort()
    if overlap==[]:
        row[name+'_explained_terms'] = '.'
    else:
        row[name+'_explained_terms'] = ';'.join(overlap)
    row[name+'_number_explained_terms'] = len(overlap)
    row[name+'_proportion_explained_terms'] = len(overlap)/int(row['number_of_phenotypes'])
    
    return row

for inh in ['C', 'NC', 'inh', 'DN']:   
    # Pathogenic SNVs
    df = df.apply(lambda row: snv_explained(row, 'pathogenic_SNVs', inh, patho = True), axis = 1)
    # LOF
    df = df.apply(lambda row: snv_explained(row, 'LOF', inh, var_type = 'lof'), axis = 1)
    # LOF_LOEUF
    df = df.apply(lambda row: snv_explained(row, 'LOF_LOEUF', inh, var_type = 'lof', loeuf = True), axis = 1)
    # Missense
    df = df.apply(lambda row: snv_explained(row, 'missense', inh, var_type = 'missense'), axis = 1)
    # Missense, LOEUF
    df = df.apply(lambda row: snv_explained(row, 'missense_LOEUF', inh, var_type = 'missense', loeuf = True), axis = 1)
    # Splice
    df = df.apply(lambda row: snv_explained(row, 'splice', inh, var_type = 'splice'), axis = 1)
    # Splice, LOEUF
    df = df.apply(lambda row: snv_explained(row, 'splice_LOEUF', inh, var_type = 'splice', loeuf = True), axis = 1)
    print(df.columns)

# CNVs
cnvs = pd.read_csv('../1_HPO_annotations/Annotated_files/cnvs_hpo_anno.csv')
def cnv_explained(row, type, inheritance, loeuf = False):
    # Phenotype HPO terms
    pheno_hpo = row['HPO_ID'].split(';')
    
    if loeuf:
        name = type+'_LOEUF'+'_'+inh
    else:
        name = type+'_'+inh
    
    # If no phenotype terms, return NA
    if pheno_hpo==['0']:
        row[name+'_explained_terms'] = np.nan
        row[name+'_number_explained_terms'] = np.nan
        row[name+'_proportion_explained_terms'] = np.nan
        return row
        
    sample = row['Sample']
    # If sample is not in CNV dataframe, return NA
    if sample not in cnvs.Sample.to_list():
        row[name+'_explained_terms'] = np.nan
        row[name+'_number_explained_terms'] = np.nan
        row[name+'_proportion_explained_terms'] = np.nan
        return row
    
    # Define variant set
    df = cnvs[cnvs.Type==type]
    
    if loeuf:
        df = df[df.LOEUF!='.']
        df = df[df.LOEUF.astype(float) <= 0.35]
    
    if inheritance=='C':
        df = df[df.inheritance.isin(['MC', 'FC'])]
    elif inheritance=='NC':
        df = df[df.inheritance.isin(['MNC', 'FNC'])]
    elif inheritance=='inh':
        df = df[df.inheritance.isin(['MC', 'FC', 'MNC', 'FNC'])]
    elif inheritance=='DN':
        df = df[df.inheritance=='de novo']
    else:
        print('NEED TO DEFINE INHERITANCE')
        print(error)
    
    df = df[df.Sample==sample]
    
    # If there are no variants in a category, return 0
    if df.shape[0]==0:
        row[name+'_explained_terms'] = '.'
        row[name+'_number_explained_terms'] = 0
        row[name+'_proportion_explained_terms'] = 0
        return row
    
    # Get gene HPO terms
    gene_hpo_str = df.GeneHPO_inclusive.to_list()
    gene_hpo = []
    for i in gene_hpo_str:
        gene_hpo += i.split(';')
    gene_hpo = list(set(gene_hpo))
    
    # Get overlap
    overlap = [i for i in gene_hpo if i in pheno_hpo]
    overlap.sort()
    if overlap==[]:
        row[name+'_explained_terms'] = '.'
    else:
        row[name+'_explained_terms'] = ';'.join(overlap)
    row[name+'_number_explained_terms'] = len(overlap)
    row[name+'_proportion_explained_terms'] = len(overlap)/int(row['number_of_phenotypes'])
    
    return row

for inh in ['C', 'NC', 'inh', 'DN']: 
    # DELs
    df = df.apply(lambda row: cnv_explained(row, 'DEL', inh), axis = 1)
    # DEL, LOEUF
    df = df.apply(lambda row: cnv_explained(row, 'DEL', inh, loeuf = True), axis = 1)
    # DUPs
    df = df.apply(lambda row: cnv_explained(row, 'DUP', inh), axis = 1)
    # DUP, LOEUF
    df = df.apply(lambda row: cnv_explained(row, 'DUP', inh, loeuf = True), axis = 1)
    print(df.columns)

# STRs
strs = pd.read_csv('../1_HPO_annotations/Annotated_files/strs_hpo_anno.csv')
def str_explained(row, name, inheritance, loeuf = False):
    name = name+'_'+inheritance
    # Phenotype HPO terms
    pheno_hpo = row['HPO_ID'].split(';')
    
    # If no phenotype terms, return NA
    if pheno_hpo==['0']:
        row[name+'_explained_terms'] = np.nan
        row[name+'_number_explained_terms'] = np.nan
        row[name+'_proportion_explained_terms'] = np.nan
        return row
        
    sample = row['Sample']
    # If sample is not in STR dataframe, return NA
    if sample not in strs.Sample.to_list():
        row[name+'_explained_terms'] = np.nan
        row[name+'_number_explained_terms'] = np.nan
        row[name+'_proportion_explained_terms'] = np.nan
        return row
    
    df = strs
    
    if loeuf:
        df = df[df.LOEUF!='.']
        df = df[df.LOEUF.astype(float) <= 0.35]
    
    if inheritance=='C':
        df = df[df.inheritance.isin(['MC', 'FC'])]
    elif inheritance=='NC':
        df = df[df.inheritance.isin(['MNC', 'FNC'])]
    elif inheritance=='inh':
        df = df[df.inheritance.isin(['MC', 'FC', 'MNC', 'FNC'])]
    else:
        print('NEED TO DEFINE INHERITANCE')
        print(error)
    
    df = df[df.Sample==sample]
    
    # If there are no variants in a category, return 0
    if df.shape[0]==0:
        row[name+'_explained_terms'] = '.'
        row[name+'_number_explained_terms'] = 0
        row[name+'_proportion_explained_terms'] = 0
        return row
    
    # Get gene HPO terms
    gene_hpo_str = df.GeneHPO_inclusive.to_list()
    gene_hpo = []
    for i in gene_hpo_str:
        gene_hpo += i.split(';')
    gene_hpo = list(set(gene_hpo))
    
    # Get overlap
    overlap = [i for i in gene_hpo if i in pheno_hpo]
    overlap.sort()
    if overlap==[]:
        row[name+'_explained_terms'] = '.'
    else:
        row[name+'_explained_terms'] = ';'.join(overlap)
    row[name+'_number_explained_terms'] = len(overlap)
    row[name+'_proportion_explained_terms'] = len(overlap)/int(row['number_of_phenotypes'])
    
    return row

# There are no de novo STRs, so do not annotate for them
for inh in ['C', 'NC', 'inh']:
    # STRs
    df = df.apply(lambda row: str_explained(row, 'STR', inh), axis = 1)
    # STR, LOEUF
    df = df.apply(lambda row: str_explained(row, 'STR_LOEUF', inh, loeuf = True), axis = 1)
print(df)

# Add more columns for "all"
#expl_terms_cols = [i for i in df.columns.to_list() if '_explained_terms' in i and 'number' not in i and 'proportion' not in i]
# April 18, 2022 - we decided to remove STR from all HPO analysis
expl_terms_cols = [i for i in df.columns.to_list() if '_explained_terms' in i and 'number' not in i and 'proportion' not in i and 'STR' not in i]
def add_all(row, inh):
    use_cols = [i for i in expl_terms_cols if inh in i]
    # Phenotype HPO terms
    pheno_hpo = row['HPO_ID'].split(';')
    
    # If no phenotype terms, return NA
    if pheno_hpo==['0']:
        row['All_'+inh+'_explained_terms'] = np.nan
        row['All_'+inh+'_number_explained_terms'] = np.nan
        row['All_'+inh+'_proportion_explained_terms'] = np.nan
        return row

    # If sample is not in snv file, cnv file, or str file, return NA
    sample = row['Sample']
    if sample not in snvs.Sample.to_list() and sample not in cnvs.Sample.to_list() and sample not in strs.Sample.to_list():
        row['All_'+inh+'_explained_terms'] = np.nan
        row['All_'+inh+'_number_explained_terms'] = np.nan
        row['All_'+inh+'_proportion_explained_terms'] = np.nan
        return row

    lst = row[use_cols].to_list()
    total_explained_terms = []
    for i in lst:
        if i==i:
            total_explained_terms += i.split(';')
    total_explained_terms = list(set([i for i in total_explained_terms if i !='.']))
    total_explained_terms.sort()
    
    # Get proportion
    if total_explained_terms==[]:
        row['All_'+inh+'_explained_terms'] = '.'
    else:
        row['All_'+inh+'_explained_terms'] = ';'.join(total_explained_terms)
    row['All_'+inh+'_number_explained_terms'] = len(total_explained_terms)
    row['All_'+inh+'_proportion_explained_terms'] = len(total_explained_terms)/int(row['number_of_phenotypes'])
    
    return row

for inh in ['C', 'NC', 'inh', 'DN']:
    df = df.apply(lambda row: add_all(row, inh), axis = 1)

# Save to file
df.to_csv('Result_tables/3_variant_type_inheritance_proportions.csv', index = False)