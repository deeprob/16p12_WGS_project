#!/bin/python
import pandas as pd
import numpy as np

# Get the proportion of phenotypes explained by each class of variants

# Load in phenotype data
pheno = pd.read_csv('../1_HPO_annotations/Annotated_files/phenotype_hpo.csv', dtype = str)
print(pheno)

# Load in variant data
snvs = pd.read_csv('../1_HPO_annotations/Annotated_files/snvs_hpo_anno.csv')

# Add in annotations for the number and IDs of phenotypes explained by the genes
# Important note: as of 1/24/2022, we are only using  astrict overlap method to determine explained phenotypes

def snv_explained(row, name, patho = False, var_type = False, loeuf = False):
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
    
# Pathogenic SNVs
df = pheno.apply(lambda row: snv_explained(row, 'pathogenic_SNVs', patho = True), axis = 1)
print(df)
# LOF
df = df.apply(lambda row: snv_explained(row, 'LOF', var_type = 'lof'), axis = 1)
print(df)
# LOF_LOEUF
df = df.apply(lambda row: snv_explained(row, 'LOF_LOEUF', var_type = 'lof', loeuf = True), axis = 1)
print(df)
# Missense
df = df.apply(lambda row: snv_explained(row, 'missense', var_type = 'missense'), axis = 1)
print(df)
# Missense, LOEUF
df = df.apply(lambda row: snv_explained(row, 'missense_LOEUF', var_type = 'missense', loeuf = True), axis = 1)
print(df)
# Splice
df = df.apply(lambda row: snv_explained(row, 'splice', var_type = 'splice'), axis = 1)
print(df)
# Splice, LOEUF
df = df.apply(lambda row: snv_explained(row, 'splice_LOEUF', var_type = 'splice', loeuf = True), axis = 1)
print(df)

# CNVs
cnvs = pd.read_csv('../1_HPO_annotations/Annotated_files/cnvs_hpo_anno.csv')
def cnv_explained(row, type, loeuf = False):
    # Phenotype HPO terms
    pheno_hpo = row['HPO_ID'].split(';')
    
    if loeuf:
        name = type+'_LOEUF'
    else:
        name = type
    
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

# DELs
df = df.apply(lambda row: cnv_explained(row, 'DEL'), axis = 1)
print(df)
# DEL, LOEUF
df = df.apply(lambda row: cnv_explained(row, 'DEL', loeuf = True), axis = 1)
print(df)
# DUPs
df = df.apply(lambda row: cnv_explained(row, 'DUP'), axis = 1)
print(df)
# DUP, LOEUF
df = df.apply(lambda row: cnv_explained(row, 'DUP', loeuf = True), axis = 1)
print(df)

# STRs
strs = pd.read_csv('../1_HPO_annotations/Annotated_files/strs_hpo_anno.csv')
def str_explained(row, name, loeuf = False):
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
# STRs
df = df.apply(lambda row: str_explained(row, 'STR'), axis = 1)
print(df)
# STR, LOEUF
df = df.apply(lambda row: str_explained(row, 'STR_LOEUF', loeuf = True), axis = 1)
print(df)

# Add one more column for 'All'
expl_terms_cols = [i for i in df.columns.to_list() if '_explained_terms' in i and 'number' not in i and 'proportion' not in i]
def add_all(row):
    # Phenotype HPO terms
    pheno_hpo = row['HPO_ID'].split(';')
    
    # If no phenotype terms, return NA
    if pheno_hpo==['0']:
        row['All_explained_terms'] = np.nan
        row['All_number_explained_terms'] = np.nan
        row['All_proportion_explained_terms'] = np.nan
        return row

    # If sample is not in snv file, cnv file, or str file, return NA
    sample = row['Sample']
    if sample not in snvs.Sample.to_list() and sample not in cnvs.Sample.to_list() and sample not in strs.Sample.to_list():
        row['All_explained_terms'] = np.nan
        row['All_number_explained_terms'] = np.nan
        row['All_proportion_explained_terms'] = np.nan
        return row

    lst = row[expl_terms_cols].to_list()
    total_explained_terms = []
    for i in lst:
        if i==i:
            total_explained_terms += i.split(';')
    total_explained_terms = list(set([i for i in total_explained_terms if i !='.']))
    total_explained_terms.sort()
    
    # Get proportion
    if total_explained_terms==[]:
        row['All_explained_terms'] = '.'
    else:
        row['All_explained_terms'] = ';'.join(total_explained_terms)
    row['All_number_explained_terms'] = len(total_explained_terms)
    row['All_proportion_explained_terms'] = len(total_explained_terms)/int(row['number_of_phenotypes'])
    
    return row
df = df.apply(add_all, axis = 1)
print(df)

# Save to file
df.to_csv('Result_tables/1_variant_type_proportions.csv', index = False)